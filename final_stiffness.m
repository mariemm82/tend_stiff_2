%%%%%%%%%%%%%%%%%%%%%%%%%%
% final_stiffness
% Marie Moltubakk 16.6.2013
% Read 3+3 trials of MTJ and OTJ scans
% Produce stiffness coeffisients and output arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%



function [fitresult, gof, force_elong_array, final_cutoff_force] = final_stiffness(time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3, time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3, forceintervals, force_cutoff_manual, force_max_trials, at_momentarm) %new2014-04-14

global plot_achilles subject_id plot_conversion % plot_check %plot_norm plot_emg


%% convert input data to cells
MTJ_trials = {time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3};
OTJ_trials = {time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3};
loc_time = 1;
loc_force = 2;
loc_displ = 3;


%% find force onset

% define onset of contraction
%   = time point where the knee extensor torque
%   exceeded baseline + 7.5 Nm (3 = Aagaard P, Simonsen EB, Andersen JL,
%   Magnusson P, and Dyhre-Poulsen P. Increased rate of force development
%   and neural drive of human skeletal muscle following resistance training. J Appl Physiol 93: 1318–1326, 2002.)

% onset as % of max force
threshold_percent = 0.01; %VAR 0.5 percent (0.00-1.00) of highest common force during ramps MMM TODO check BD 102 GM, one trial cut too early with 0.005 --> changed to 0.01
%threshold_fixed = 7.5/at_momentarm; %VAR threshold in Nm, converted to N
threshold_add = threshold_percent * min(force_max_trials);

baseline_frames = 10; %VAR 10 first frames used to define baseline

% preallocate
loc_MTJ_onset(1:3) = NaN;
loc_OTJ_onset(1:3) = NaN;

% locate onset
for i = 1:length(MTJ_trials)
    if length(MTJ_trials{i}) > 1
        threshold = mean(MTJ_trials{i}(1:baseline_frames,loc_force)) + threshold_add;
        loc_MTJ_onset(i) = find(MTJ_trials{i}(:,loc_force) >= threshold, 1, 'first');
    end
end
for i = 1:length(OTJ_trials)
    if length(OTJ_trials{i}) > 1
        threshold = mean(OTJ_trials{i}(1:baseline_frames,loc_force)) + threshold_add;
        loc_OTJ_onset(i) = find(OTJ_trials{i}(:,loc_force) >= threshold, 1, 'first');
    end
end

% plot onset checpoints
if plot_conversion
    plottitle = horzcat('Force onset check, ', subject_id);
    figure('Name',plottitle)
    hold on
    % MTJ
    for i = 1:length(MTJ_trials)
        if length(MTJ_trials{i}) > 1
            plot(MTJ_trials{i}(:,loc_time),MTJ_trials{i}(:,loc_force))
        end
    end
    % OTJ
    for i = 1:length(OTJ_trials)
        if length(OTJ_trials{i}) > 1
            plot(OTJ_trials{i}(:,loc_time),OTJ_trials{i}(:,loc_force))
        end
    end
    % onset x2
    set(gca,'ColorOrderIndex',1)
    for i = 1:length(MTJ_trials)
        if length(MTJ_trials{i}) > 1
            plot([MTJ_trials{i}(loc_MTJ_onset(i),loc_time) MTJ_trials{i}(loc_MTJ_onset(i),loc_time)], [threshold_add*5 -20])
        end
    end
    for i = 1:length(OTJ_trials)
        if length(OTJ_trials{i}) > 1
            plot([OTJ_trials{i}(loc_OTJ_onset(i),loc_time) OTJ_trials{i}(loc_OTJ_onset(i),loc_time)], [threshold_add*5 -20])
        end
    end
    % visual
    xlabel('Time (s)')
    ylabel('Force (N)')
    title(plottitle)
    saveas(gcf, strcat('data_plots_stiff/IND_stiff_onset_', subject_id), 'png')
    axis([0 2 -20 threshold_add+50])
    saveas(gcf, strcat('data_plots_stiff/IND_stiff_onset_', subject_id, '_ZOOM'), 'png')
end

% cut arrays at onset, and offset displacement
for i = 1:length(MTJ_trials)
    if length(MTJ_trials{i}) > 1
        % delete before onset
        MTJ_trials{i}(1:loc_MTJ_onset(i)-1,:) = [];
        % offset displacement and time
        MTJ_trials{i}(:,loc_displ) = MTJ_trials{i}(:,loc_displ) - MTJ_trials{i}(1,loc_displ);
        MTJ_trials{i}(:,loc_time) = MTJ_trials{i}(:,loc_time) - MTJ_trials{i}(1,loc_time);
    end
end
for i = 1:length(OTJ_trials)
    if length(OTJ_trials{i}) > 1
        % delete before onset
        OTJ_trials{i}(1:loc_OTJ_onset(i)-1,:) = [];
        % offset displacement
        OTJ_trials{i}(:,loc_displ) = OTJ_trials{i}(:,loc_displ) - OTJ_trials{i}(1,loc_displ);
        OTJ_trials{i}(:,loc_time) = OTJ_trials{i}(:,loc_time) - OTJ_trials{i}(1,loc_time);
    end
end



%% find force range

% print maximal force from each trial to the screen
cprintf('black', horzcat('Ramp trials, ind max force = '))
cprintf('*blue', horzcat(num2str(round(min(force_max_trials))), ' N'))
cprintf('black', horzcat(', trials = ', num2str(round(force_max_trials(1))), ' - ', num2str(round(force_max_trials(2))), ' - ', num2str(round(force_max_trials(3))), ' - ', num2str(round(force_max_trials(4))), ' - ', num2str(round(force_max_trials(5))), ' - ', num2str(round(force_max_trials(6))), ' N.\n'))

% find lowest common force between all 6 trials
commonforce = min(force_max_trials);

% create array of force values to use for averaging
force_array = (0:forceintervals:commonforce)';



%% calculate tendon elongation

% preallocate
displ_MTJ(1:length(force_array),1:length(MTJ_trials)) = NaN;
displ_OTJ(1:length(force_array),1:length(OTJ_trials)) = NaN;

% keeping 10 data points after commonforce, in order to do proper interpolation around the size of commonforce
keep_datapoints = 10; %VAR

for i = 1:length(MTJ_trials)
    if length(MTJ_trials{i}) > 1
        % detect position of max common force in array
        del_mtj = find(MTJ_trials{i}(:,loc_force) > force_array(end), 1, 'first');
        % delete array contents after commonforce
        MTJ_trials{i}(del_mtj+keep_datapoints:end,:) = [];
        % spline to get displacement values at common force levels
        displ_MTJ(:,i) = interp1(MTJ_trials{i}(:,loc_force),MTJ_trials{i}(:,loc_displ),force_array,'linear'); % ALT: spline(MTJ_trials{i}(:,loc_force),MTJ_trials{i}(:,loc_displ),force_array);
    end
end

for i = 1:length(OTJ_trials)
    if length(OTJ_trials{i}) > 1
        % detect position of max common force in array
        del_OTJ = find(OTJ_trials{i}(:,loc_force) > force_array(end), 1, 'first');
        % delete array contents after commonforce
        OTJ_trials{i}(del_OTJ+keep_datapoints:end,:) = [];
        % spline to get displacement values at common force levels
        displ_OTJ(:,i) = interp1(OTJ_trials{i}(:,loc_force),OTJ_trials{i}(:,loc_displ),force_array,'linear');
    end
end

% average displacement trials, calculate elongation
displ_mtj_mean = nanmean(displ_MTJ,2);
displ_otj_mean = nanmean(displ_OTJ,2);
% manually correct first datapoint (extrapolated by spline, since force is cut above zero N)
displ_mtj_mean(1) = 0;
displ_otj_mean(1) = 0;
tend_elong = displ_mtj_mean - displ_otj_mean;

if plot_achilles
    % MTJ force-disp x 3
    plottitle = horzcat('MTJ stiffness check for ', subject_id);
    figure('Name',plottitle)
    hold on
    % 3 MTJ after cutoff, all datapoints
    for i = 1:length(MTJ_trials)
        if length(MTJ_trials{i}) > 1
            plot(MTJ_trials{i}(:,loc_displ),MTJ_trials{i}(:,loc_force))
        else
            plot([0 0],[0 0])
        end
    end
    % mean MTJ
    plot(displ_mtj_mean, force_array, 'k','LineWidth',2)
    % 3 original MTJ before cutoff
    set(gca,'ColorOrderIndex',1)
    if length(MTJ_trials{1}) > 1
        plot(time_force_displ_mtj1(:,3),time_force_displ_mtj1(:,2),':','LineWidth',0.3)
    else
            plot([0 0],[0 0],':')
    end
    if length(MTJ_trials{2}) > 1
        plot(time_force_displ_mtj2(:,3),time_force_displ_mtj2(:,2),':','LineWidth',0.3)
    else
            plot([0 0],[0 0],':')
    end
    if length(MTJ_trials{3}) > 1
        plot(time_force_displ_mtj3(:,3),time_force_displ_mtj3(:,2),':','LineWidth',0.3)
    else
            plot([0 0],[0 0],':')
    end
    % 3 MTJ - dots at force intervals
    set(gca,'ColorOrderIndex',1)
    for i = 1:length(MTJ_trials)
        if length(MTJ_trials{i}) > 1
            plot(displ_MTJ(:,i), force_array,'.')
        else
            plot([0 0],[0 0])
        end
    end
    % visual
    xlabel('Displacement (mm)'),ylabel('Tendon force (N)'),title(plottitle);
    legend('trial 1 w/onset','trial 2 w/onset','trial 3 w/onset','mean','trial 1 raw','Location','Northwest');
    saveas(gcf, strcat('data_plots_stiff/IND_stiff_MTJ_avg_', subject_id), 'png')
end

if plot_achilles
    % OTJ force-disp x 3
    plottitle = horzcat('OTJ stiffness check for ', subject_id);
    figure('Name',plottitle)
    hold on
    % 3 OTJ after cutoff, all datapoints
    for i = 1:length(OTJ_trials)
        if length(OTJ_trials{i}) > 1
            plot(OTJ_trials{i}(:,loc_displ),OTJ_trials{i}(:,loc_force))
        else
            plot([0 0],[0 0])
        end
    end
    % mean OTJ
    plot(displ_otj_mean, force_array, 'k','LineWidth',2)
    % 3 original OTJ before cutoff
    set(gca,'ColorOrderIndex',1)
    if length(OTJ_trials{1}) > 1
        plot(time_force_displ_otj1(:,3),time_force_displ_otj1(:,2),':','LineWidth',0.3)
    else
            plot([0 0],[0 0],':')
    end
    if length(OTJ_trials{2}) > 1
        plot(time_force_displ_otj2(:,3),time_force_displ_otj2(:,2),':','LineWidth',0.3)
    else
            plot([0 0],[0 0],':')
    end
    if length(OTJ_trials{3}) > 1
        plot(time_force_displ_otj3(:,3),time_force_displ_otj3(:,2),':','LineWidth',0.3) 
    else
            plot([0 0],[0 0],':')
    end
    % 3 OTJ - dots at force intervals
    set(gca,'ColorOrderIndex',1)
    for i = 1:length(OTJ_trials)
        if length(OTJ_trials{i}) > 1
            plot(displ_OTJ(:,i), force_array,'.')
        else
            plot([0 0],[0 0])
        end
    end
    % visual
    xlabel('Displacement (mm)'),ylabel('Tendon force (N)'),title(plottitle);
    legend('trial 1 w/onset','trial 2 w/onset','trial 3 w/onset','mean','trial 1 raw','Location','Northwest');
    saveas(gcf, strcat('data_plots_stiff/IND_stiff_OTJ_avg_', subject_id), 'png')
end




%% cut force (method 4, June 2017: RTD)

% preallocate
time_array = cell(1,6);
force_splined = cell(1,6);
force_dev = cell(1,6);

% common time array:
timeintervals = 0.05;

for i = 1:length(MTJ_trials)
    if length(MTJ_trials{i}) > 1
        time_array{i} = (0:timeintervals:MTJ_trials{i}(end,loc_time))';
        force_splined{i} = spline(MTJ_trials{i}(:,loc_time),MTJ_trials{i}(:,loc_force),time_array{i});
        force_dev{i} = [NaN; diff(force_splined{i})/timeintervals];
    end
end

for i = 1:length(OTJ_trials)
    if length(OTJ_trials{i}) > 1
        time_array{i+length(MTJ_trials)} = (0:timeintervals:OTJ_trials{i}(end,loc_time))';
        force_splined{i+length(MTJ_trials)} = spline(OTJ_trials{i}(:,loc_time),OTJ_trials{i}(:,loc_force),time_array{i+length(MTJ_trials)});
        force_dev{i+length(MTJ_trials)} = [NaN; diff(force_splined{i+length(MTJ_trials)})/timeintervals];
    end
end

if plot_conversion
    plottitle = horzcat('Rate of force development, 6 trials, ', subject_id);
    figure('Name',plottitle)
    hold on
    %yyaxis left
    for i = 1:length(time_array)
        plot(time_array{i},force_splined{i})
    end
    ylabel('Force (N)')
    set(gca,'ColorOrderIndex',1)
    %yyaxis right
    for i = 1:length(time_array)
        plot(time_array{i},force_dev{i},'--')
    end
    plot([0 max(cell2mat(cellfun(@size,time_array,'uni',false))) * timeintervals],[0 0],'k')
    plot([0 max(cell2mat(cellfun(@size,time_array,'uni',false))) * timeintervals],[100 100],'k')
    ylabel('Force (N) / force development (N/s)')
    xlabel('Time (s)')
    axis([-Inf Inf -100 commonforce*1.05])
    title(plottitle)
    saveas(gcf, strcat('data_plots_stiff/IND_RFD_', subject_id), 'png')
end


    
% MMM TODO GOON: RTD calculation for cutoff?
% The best way I think is slope of the force – time curve from onset to +50ms or +100 ms.
% It is important to have a strict criteria for onset. 
% Peak RFD i.e. the steepest point to point slope of the curve is often reported but depends extremely much on filtering and is less reproducible...and valid I think..
%  
% See eg. Enclosed BM et al p. 987 second paragraph. 
% Also enclosed Aagaard et al 2002 with similar description. 
% I think in both these papers it is done as mentioned above + 0-30ms, 0-peak MVC, and peak... 
% also areas under the curve are calculated (‘contractile impulse’). But keep it simple I would suggest...

% MMM goon - define cutoff

% MMM goon - for subjects where one trial is discarded





%% cut force (method 3, June 2014: 90% of max elongation)

% prepare for cutoff at 90% of max elongation
percent_elong = 0.90; %VAR
[max_elong,loc_max_elong] = max(tend_elong);
keep_elong = percent_elong * max_elong;
loc_elong_cut = find(tend_elong>keep_elong,1,'first') - 1; % finds the first value larger than 90% limit - we will use the largest BEFORE 90% limit

% write cutoff to screen
cprintf('black', horzcat('Max elong = ', num2str(round(max_elong,2)), ' mm @ ', num2str(force_array(loc_max_elong)), ' N, ', num2str(percent_elong*100), '%% elong -> cutoff at '))
cprintf('*blue', horzcat(num2str(force_array(loc_elong_cut)), ' N.\n'))

% cut off even earlier if datamaster is set to cut at lower value
if force_array(loc_elong_cut) > str2double(force_cutoff_manual)
    cprintf('red',horzcat('Force cutoff is lowered manually (datamaster), to ', force_cutoff_manual, ' N.\n'));
    loc_elong_cut = find(force_array>str2double(force_cutoff_manual),1,'first') - 1;
end

% save the applied max force level, for output to file
final_cutoff_force = force_array(loc_elong_cut);


%% curve fitting for stiffness
[fitresult, gof] = fit_stiffness(tend_elong, force_array, loc_elong_cut, displ_mtj_mean, displ_otj_mean);

% write stiffness coefficients to screen
cprintf('*blue', horzcat('Stiffness coeffs = ', regexprep(num2str(coeffvalues(fitresult)),' +',' '), '. R2 = ', num2str(gof.rsquare), '\n'))

% create stiffness data array for later use
force_elong_array = [tend_elong(1:loc_elong_cut) force_array(1:loc_elong_cut)];

end