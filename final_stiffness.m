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

% version 1: fixed threshold
%   = time point where the knee extensor torque
%   exceeded baseline + 7.5 Nm (3 = Aagaard P, Simonsen EB, Andersen JL,
%   Magnusson P, and Dyhre-Poulsen P. Increased rate of force development
%   and neural drive of human skeletal muscle following resistance training. J Appl Physiol 93: 1318–1326, 2002.)
%threshold_fixed = 7.5/at_momentarm; %VAR threshold in Nm, converted to N

% version 2: onset = % of max force, above start values
%threshold_percent = 0.005; %VAR percent (0.00-1.00) of highest common force during ramps  --- BD 102 (?) too early cut

% version 3: increase above SD of the first frames
baseline_frames = 6; %VAR X first frames used to define baseline
% remainder is performed inside loop below

% preallocate
threshold_MTJ(1:3) = NaN;
threshold_OTJ(1:3) = NaN;
loc_MTJ_onset(1:3) = NaN;
loc_OTJ_onset(1:3) = NaN;

% locate onset - output variable is the first frame to KEEP
for i = 1:length(MTJ_trials)
    if length(MTJ_trials{i}) > 1
        baseline_mean = mean(MTJ_trials{i}(1:baseline_frames,loc_force));
        baseline_SD = std(MTJ_trials{i}(1:baseline_frames,loc_force));
        threshold_MTJ(i) = mean(MTJ_trials{i}(1:baseline_frames,loc_force)) + 2*baseline_SD;
        loc_MTJ_onset(i) = find(MTJ_trials{i}(:,loc_force) >= threshold_MTJ(i), 1, 'first');
    end
end
for i = 1:length(OTJ_trials)
    if length(OTJ_trials{i}) > 1
        baseline_mean = mean(OTJ_trials{i}(1:baseline_frames,loc_force));
        baseline_SD = std(OTJ_trials{i}(1:baseline_frames,loc_force));
        threshold_OTJ(i) = mean(OTJ_trials{i}(1:baseline_frames,loc_force)) + 2*baseline_SD;
        loc_OTJ_onset(i) = find(OTJ_trials{i}(:,loc_force) >= threshold_OTJ(i), 1, 'first');
    end
end

% plot onset checkpoints
if plot_achilles
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
            plot([MTJ_trials{i}(loc_MTJ_onset(i),loc_time) MTJ_trials{i}(loc_MTJ_onset(i),loc_time)], [max([threshold_MTJ threshold_OTJ])*10 -20])
        end
    end
    for i = 1:length(OTJ_trials)
        if length(OTJ_trials{i}) > 1
            plot([OTJ_trials{i}(loc_OTJ_onset(i),loc_time) OTJ_trials{i}(loc_OTJ_onset(i),loc_time)], [max([threshold_MTJ threshold_OTJ])*10 -20])
        end
    end
    % visual
    xlabel('Time (s)')
    ylabel('Force (N)')
    title(plottitle)
    saveas(gcf, strcat('data_plots_stiff/IND_stiff_onset_', subject_id), 'png')
    axis([0 0.2+MTJ_trials{1}(max([loc_MTJ_onset loc_OTJ_onset]),loc_time) -20 max([threshold_MTJ threshold_OTJ])*10])
    saveas(gcf, strcat('data_plots_stiff/IND_stiff_onset_', subject_id, '_ZOOM'), 'png')
end


%% cut arrays at onset, and offset time+displacement accordingly
for i = 1:length(MTJ_trials)
    if length(MTJ_trials{i}) > 1
        % offset displacement and time
        MTJ_trials{i}(:,loc_displ) = MTJ_trials{i}(:,loc_displ) - mean(MTJ_trials{i}(loc_MTJ_onset(i)-2:loc_MTJ_onset(i),loc_displ)); % average 3 points until and including onset
        MTJ_trials{i}(:,loc_time) = MTJ_trials{i}(:,loc_time) - MTJ_trials{i}(loc_MTJ_onset(i),loc_time);
        % delete before onset
        MTJ_trials{i}(1:loc_MTJ_onset(i)-1,:) = [];
        % warning if first force is larger than ... 20 N?
        if MTJ_trials{i}(1,loc_force) > 20
            cprintf('red', horzcat('WARNING: Cutoff creates force onset > threshold of 20 N: ', num2str(round(MTJ_trials{i}(1,loc_force),1)), '\n'))
        end
    end
end
for i = 1:length(OTJ_trials)
    if length(OTJ_trials{i}) > 1
        % offset displacement and time
        OTJ_trials{i}(:,loc_displ) = OTJ_trials{i}(:,loc_displ) - mean(OTJ_trials{i}(loc_OTJ_onset(i)-2:loc_OTJ_onset(i),loc_displ));
        OTJ_trials{i}(:,loc_time) = OTJ_trials{i}(:,loc_time) - OTJ_trials{i}(loc_OTJ_onset(i),loc_time);
        % delete before onset
        OTJ_trials{i}(1:loc_OTJ_onset(i)-1,:) = [];
        % warning if first force is larger than ... 20 N?
        if OTJ_trials{i}(1,loc_force) > 20
            cprintf('red', horzcat('WARNING: Cutoff creates force onset > threshold of 20 N: ', num2str(round(OTJ_trials{i}(1,loc_force),1)), '\n'))
        end
    end
end




%% find common force level, create common array

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
keep_datapoints = 10; %VAR -- MMM TODO - some adjustment to this?

for i = 1:length(MTJ_trials)
    if length(MTJ_trials{i}) > 1
        % detect position of max common force in array
        loc_del_mtj = find(MTJ_trials{i}(:,loc_force) > force_array(end), 1, 'first');
        % delete array contents after commonforce
        MTJ_trials{i}(loc_del_mtj+keep_datapoints:end,:) = [];
        % spline to get displacement values at common force levels
        displ_MTJ(:,i) = pchip(MTJ_trials{i}(:,loc_force),MTJ_trials{i}(:,loc_displ),force_array); % ALT: interp1(MTJ_trials{i}(:,loc_force),MTJ_trials{i}(:,loc_displ),force_array,'linear'); % ALT: spline(MTJ_trials{i}(:,loc_force),MTJ_trials{i}(:,loc_displ),force_array);
    end
end

for i = 1:length(OTJ_trials)
    if length(OTJ_trials{i}) > 1
        % detect position of max common force in array
        loc_del_OTJ = find(OTJ_trials{i}(:,loc_force) > force_array(end), 1, 'first');
        % delete array contents after commonforce
        OTJ_trials{i}(loc_del_OTJ+keep_datapoints:end,:) = [];
        % spline to get displacement values at common force levels
        displ_OTJ(:,i) = pchip(OTJ_trials{i}(:,loc_force),OTJ_trials{i}(:,loc_displ),force_array);
    end
end

% average displacement trials, calculate elongation
displ_mtj_mean = mean(displ_MTJ,2); % nanmean to average 2 values if 3rd is NaN
displ_otj_mean = mean(displ_OTJ,2);
tend_elong = displ_mtj_mean - displ_otj_mean;
tend_elong(1) = NaN;

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
        plot(time_force_displ_mtj1(:,loc_displ),time_force_displ_mtj1(:,loc_force),':','LineWidth',0.3)
    else
            plot([0 0],[0 0],':')
    end
    if length(MTJ_trials{2}) > 1
        plot(time_force_displ_mtj2(:,loc_displ),time_force_displ_mtj2(:,loc_force),':','LineWidth',0.3)
    else
            plot([0 0],[0 0],':')
    end
    if length(MTJ_trials{3}) > 1
        plot(time_force_displ_mtj3(:,loc_displ),time_force_displ_mtj3(:,loc_force),':','LineWidth',0.3)
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
        plot(time_force_displ_otj1(:,loc_displ),time_force_displ_otj1(:,loc_force),':','LineWidth',0.3)
    else
            plot([0 0],[0 0],':')
    end
    if length(OTJ_trials{2}) > 1
        plot(time_force_displ_otj2(:,loc_displ),time_force_displ_otj2(:,loc_force),':','LineWidth',0.3)
    else
            plot([0 0],[0 0],':')
    end
    if length(OTJ_trials{3}) > 1
        plot(time_force_displ_otj3(:,loc_displ),time_force_displ_otj3(:,loc_force),':','LineWidth',0.3) 
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
RFD_MTJ_raw = cell(1,3);
RFD_OTJ_raw = cell(1,3);
RFD_MTJ_filt = cell(1,3);
RFD_OTJ_filt = cell(1,3);
freq_MTJ(1:3) = NaN;
freq_OTJ(1:3) = NaN;

% create running average filter with x points:
filter_points = 5; %VAR
filter_ravg = 1/filter_points * ones(filter_points,1);
%RFD_window = 0.02; %VAR 0.02s = 20ms. the 20-millisecond sampling window has been shown to be the most reliable (39). https://www.scienceforsport.com/rate-of-force-development-rfd-2/#how-to-calculate-the-rate-of-force-development

% calculate and smooth RFD
for i = 1:length(MTJ_trials)
    if length(MTJ_trials{i}) > 1
%        time_array{i} = (0:timeintervals:MTJ_trials{i}(end,loc_time))';
%        force_splined{i} = pchip(MTJ_trials{i}(:,loc_time),MTJ_trials{i}(:,loc_force),time_array{i});
%        RFD{i} = [NaN; diff(force_splined{i})/timeintervals];
        RFD_MTJ_raw{i} = diff(MTJ_trials{i}(:,loc_force))./diff(MTJ_trials{i}(:,loc_time));
        RFD_MTJ_filt{i} = filtfilt(filter_ravg, 1, RFD_MTJ_raw{i});
        freq_MTJ(i) = 1/mean(diff(MTJ_trials{i}(:,loc_time)));
    end
end

for i = 1:length(OTJ_trials)
    if length(OTJ_trials{i}) > 1
        RFD_OTJ_raw{i} = diff(OTJ_trials{i}(:,loc_force))./diff(OTJ_trials{i}(:,loc_time));
        RFD_OTJ_filt{i} = filtfilt(filter_ravg, 1, RFD_OTJ_raw{i});
        freq_OTJ(i) = 1/mean(diff(OTJ_trials{i}(:,loc_time)));
    end
end

% RFD raw and filt are one frame shorter than MTJ_trials etc, because of
% differentiation.

if plot_achilles
    plottitle = horzcat('Rate of force development, 6 trials, ', subject_id);
    figure('Name',plottitle)
    hold on
    %yyaxis left
    for i = 1:length(MTJ_trials)
        plot(MTJ_trials{i}(:,loc_time),MTJ_trials{i}(:,loc_force))
    end % TODO - if length < 1
    for i = 1:length(OTJ_trials)
        plot(OTJ_trials{i}(:,loc_time),OTJ_trials{i}(:,loc_force))
    end
    ylabel('Force (N)')
    set(gca,'ColorOrderIndex',1)
    %yyaxis right TODO
    for i = 1:length(MTJ_trials)
        plot(MTJ_trials{i}(2:end,loc_time),RFD_MTJ_filt{i},':')
    end
    for i = 1:length(OTJ_trials)
        plot(OTJ_trials{i}(2:end,loc_time),RFD_OTJ_filt{i},':')
    end
    % horizontal lines @ various force levels
    plot([0 4],[0 0],'k')
    plot([0 4],[100 100],'k') %old instead of Inf: max(cell2mat(cellfun(@size,time_array,'uni',false))) * timeintervals
    ylabel('Force (N) / force development (N/s)')
    xlabel('Time (s)')
    axis([-Inf Inf -100 commonforce])
    title(plottitle)
    saveas(gcf, strcat('data_plots_stiff/IND_RFD_', subject_id), 'png')
end


    
% MMM TODO: RTD calculation for cutoff?
% GOON:
% SD smoothed curve, of range from zero to commonforce
% where does smoothed curve change too much




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
% MMM TODO - require going through 0,0?
[fitresult, gof] = fit_stiffness(tend_elong, force_array, loc_elong_cut, displ_mtj_mean, displ_otj_mean);

% write stiffness coefficients to screen
cprintf('*blue', horzcat('Stiffness coeffs = ', regexprep(num2str(coeffvalues(fitresult)),' +',' '), '. R2 = ', num2str(gof.rsquare), '\n'))

% create stiffness data array for later use
force_elong_array = [tend_elong(1:loc_elong_cut) force_array(1:loc_elong_cut)];

end