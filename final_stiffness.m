%%%%%%%%%%%%%%%%%%%%%%%%%%
% final_stiffness
% Marie Moltubakk 16.6.2013
% Read 3+3 trials of MTJ and OTJ scans
% Produce stiffness coeffisients and output arrays
% 
% PLOTS show elongation up to the common force of 6 trials (to force_array_full)
% CALCULATIONS - force_elong_array - uses only 90% of lowest force (to force_array_cut, eventually to manual cutoff point set in datamaster)
% STIFFNESS FIT - uses zero to force_array_cut (eventually manual)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fitresult, gof, force_elong_array] = final_stiffness(time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3, time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3, forceintervals, force_cutoff_manual, force_max_trials) %new2014-04-14
global plot_achilles subject_id % plot_conversion plot_check %plot_norm plot_emg


%% convert input data to cells
MTJ_trials = {time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3};
OTJ_trials = {time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3};
loc_time = 1;
loc_force = 2;
loc_displ = 3;


%% calculate 90% of max force (cutoff method 5, June 2017) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find common FORCE level, create common force array

% find lowest common force between all 6 trials
commonforce = min(force_max_trials);

% define force cutoff level:
force_cutoff = 0.9; %VAR - cut off force at 90% of common max

% create array of force values to use for averaging
force_array_full = (0:forceintervals:commonforce)';
force_array_cut = (0:forceintervals:(commonforce*force_cutoff))';

% print to screen: max force from each trial + final cut force
cprintf('*black', horzcat('Ramps, trials max force: '))
cprintf('black', horzcat(...
    num2str(round(force_max_trials(1))), ' - ', num2str(round(force_max_trials(2))), ' - ', num2str(round(force_max_trials(3))), ' - ', ...
    num2str(round(force_max_trials(4))), ' - ', num2str(round(force_max_trials(5))), ' - ', num2str(round(force_max_trials(6))), ...
    ' N.\n 90%% of common force: '))
cprintf('*blue', horzcat(num2str(force_array_cut(end)), ' N\n'))


%% find ONSET of contraction
% version 3: increase above SD of the first frames
baseline_frames = 10; %VAR X first frames of data is used to define baseline
baseline_SD_multiplier = 2; %VAR onset is when force level increases to X*SD above baseline force
baseline_onset = 2; %VAR  if within baseline determination zone, force goes up by avg of 2N from frame to frame - baseline at the beginning
% remainder is performed inside loop below

% preallocate
threshold_MTJ(1:3) = NaN;
threshold_OTJ(1:3) = NaN;
loc_MTJ_onset(1:3) = NaN;
loc_OTJ_onset(1:3) = NaN;

% locate onset - output variable "loc_MTJ_onset" is the first frame to KEEP
for i = 1:length(MTJ_trials)
    if length(MTJ_trials{i}) > 1
        if mean(diff(MTJ_trials{i}(1:baseline_frames,loc_force))) > baseline_onset
            cprintf('red', horzcat('WARNING: MTJ', num2str(i), ': Force increases from start of file (first ', num2str(baseline_frames), ' frames), by avg 2N/frame. Set first frame as onset = ', num2str(MTJ_trials{i}(1,loc_force)), ' N.\n'))
            loc_MTJ_onset(i) = 1;
            threshold_MTJ(i) = MTJ_trials{i}(loc_MTJ_onset(i),loc_force);
        else
            baseline_SD = std(MTJ_trials{i}(1:baseline_frames,loc_force));
            threshold_MTJ(i) = mean(MTJ_trials{i}(1:baseline_frames,loc_force)) + baseline_SD_multiplier*baseline_SD;
            loc_MTJ_onset(i) = 1 + find(MTJ_trials{i}(2:end,loc_force) >= threshold_MTJ(i), 1, 'first'); % don't check first frame - often contaminated by noise
        end
    end
end
for i = 1:length(OTJ_trials)
    if length(OTJ_trials{i}) > 1
        if mean(diff(OTJ_trials{i}(1:baseline_frames,loc_force))) > baseline_onset
            cprintf('red', horzcat('WARNING: OTJ', num2str(i), ': Force increases from start of file (first ', num2str(baseline_frames), ' frames), by avg 2N/frame. Set first frame as onset =', num2str(OTJ_trials{i}(1,loc_force)), ' N.\n'))
            loc_OTJ_onset(i) = 1;
            threshold_OTJ(i) = OTJ_trials{i}(loc_OTJ_onset(i),loc_force);
        else
            baseline_SD = std(OTJ_trials{i}(1:baseline_frames,loc_force));
            threshold_OTJ(i) = mean(OTJ_trials{i}(1:baseline_frames,loc_force)) + baseline_SD_multiplier*baseline_SD;
            loc_OTJ_onset(i) = 1 + find(OTJ_trials{i}(2:end,loc_force) >= threshold_OTJ(i), 1, 'first'); % don't check first frame - often contaminated by noise
        end
    end
end

% plot onset checkpoints
if plot_achilles
    plottitle = horzcat('Force onset check, 6 trials, ', subject_id);
    figure('Name',plottitle)
    hold on
    % MTJ
    for i = 1:length(MTJ_trials)
        if length(MTJ_trials{i}) > 1
            plot(MTJ_trials{i}(:,loc_time),MTJ_trials{i}(:,loc_force),'.')
        end
    end
    % OTJ
    for i = 1:length(OTJ_trials)
        if length(OTJ_trials{i}) > 1
            plot(OTJ_trials{i}(:,loc_time),OTJ_trials{i}(:,loc_force),'.')
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
    title(plottitle,'Interpreter', 'none')
    saveas(gcf, strcat('data_plots_stiff/IND_force_onset_', subject_id), 'png')
    [latestonset,loc_latestonset] = (max([loc_MTJ_onset loc_OTJ_onset]));
    if loc_latestonset > 3
        axis([0 0.2+OTJ_trials{loc_latestonset-3}(latestonset,loc_time) -20 max([threshold_MTJ threshold_OTJ])*10])
    else
        axis([0 0.2+MTJ_trials{loc_latestonset}(latestonset,loc_time) -20 max([threshold_MTJ threshold_OTJ])*10])
    end
    saveas(gcf, strcat('data_plots_stiff/IND_force_onset_', subject_id, '_ZOOM'), 'png')
end


%% REMOVE frames before onset + offset time+displacement accordingly
for i = 1:length(MTJ_trials)
    if length(MTJ_trials{i}) > 1
        if loc_MTJ_onset(i) == 1
            % do nothing
        elseif loc_MTJ_onset(i) == 2
            MTJ_trials{i}(:,loc_displ) = MTJ_trials{i}(:,loc_displ) - mean(MTJ_trials{i}(1:loc_MTJ_onset(i),loc_displ)); % average 2 points until and including onset
            MTJ_trials{i}(:,loc_time) = MTJ_trials{i}(:,loc_time) - MTJ_trials{i}(loc_MTJ_onset(i),loc_time);
            % delete before onset
            MTJ_trials{i}(1:loc_MTJ_onset(i)-1,:) = [];
        else % USUALLY
            MTJ_trials{i}(:,loc_displ) = MTJ_trials{i}(:,loc_displ) - mean(MTJ_trials{i}(loc_MTJ_onset(i)-2:loc_MTJ_onset(i),loc_displ)); % average 3 points until and including onset
            MTJ_trials{i}(:,loc_time) = MTJ_trials{i}(:,loc_time) - MTJ_trials{i}(loc_MTJ_onset(i),loc_time);
            % delete before onset
            MTJ_trials{i}(1:loc_MTJ_onset(i)-1,:) = [];
        end
    end
end
for i = 1:length(OTJ_trials)
    if length(OTJ_trials{i}) > 1
        if loc_OTJ_onset(i) == 1
            % do nothing
        elseif loc_OTJ_onset(i) == 2
            OTJ_trials{i}(:,loc_displ) = OTJ_trials{i}(:,loc_displ) - mean(OTJ_trials{i}(1:loc_OTJ_onset(i),loc_displ)); % average 2 points until and including onset
            OTJ_trials{i}(:,loc_time) = OTJ_trials{i}(:,loc_time) - OTJ_trials{i}(loc_OTJ_onset(i),loc_time);
            % delete before onset
            OTJ_trials{i}(1:loc_OTJ_onset(i)-1,:) = [];
        else % USUALLY
            OTJ_trials{i}(:,loc_displ) = OTJ_trials{i}(:,loc_displ) - mean(OTJ_trials{i}(loc_OTJ_onset(i)-2:loc_OTJ_onset(i),loc_displ)); % average 3 points until and including onset
            OTJ_trials{i}(:,loc_time) = OTJ_trials{i}(:,loc_time) - OTJ_trials{i}(loc_OTJ_onset(i),loc_time);
            % delete before onset
            OTJ_trials{i}(1:loc_OTJ_onset(i)-1,:) = [];
        end
    end
end


%% calculate tendon ELONGATION

% preallocate
displ_MTJ(1:length(force_array_full),1:length(MTJ_trials)) = NaN;
displ_OTJ(1:length(force_array_full),1:length(OTJ_trials)) = NaN;

% displacement of the MTJ points, at equal force intervals:
for i = 1:length(MTJ_trials)
    if length(MTJ_trials{i}) > 1
        % detect position of max force in array
        [~,loc_del_mtj] = max(MTJ_trials{i}(:,loc_force));
        % delete array contents AFTER (+1) max force
        MTJ_trials{i}(loc_del_mtj+1:end,:) = [];
        % spline to get displacement values at common force levels
        displ_MTJ(:,i) = pchip(MTJ_trials{i}(:,loc_force),MTJ_trials{i}(:,loc_displ),force_array_full);
        % ALTERNATIVE:   interp1(MTJ_trials{i}(:,loc_force),MTJ_trials{i}(:,loc_displ),force_array,'linear');
        % ALTERNATIVE:   spline(MTJ_trials{i}(:,loc_force),MTJ_trials{i}(:,loc_displ),force_array);
        if abs(displ_MTJ(1,i)) > 0.05 % delete first value if spline function creates a point way off (will not be used anyway)
            displ_MTJ(1,i) = NaN;
        end
        if abs(displ_MTJ(end,i)) > 50
            displ_MTJ(end,i) = NaN;
        end
    end
end

% displacement of the OTJ points, at equal force intervals:
for i = 1:length(OTJ_trials)
    if length(OTJ_trials{i}) > 1
        % detect position of max force in array
        [~,loc_del_otj] = max(OTJ_trials{i}(:,loc_force));
        % delete array contents  AFTER (+1) max force
        OTJ_trials{i}(loc_del_otj+1:end,:) = [];
        % spline to get displacement values at common force levels
        displ_OTJ(:,i) = pchip(OTJ_trials{i}(:,loc_force),OTJ_trials{i}(:,loc_displ),force_array_full);
        if abs(displ_OTJ(1,i)) > 0.05
            displ_OTJ(1,i) = NaN;
        end
        if abs(displ_OTJ(end,i)) > 50
            displ_OTJ(end,i) = NaN;
        end
    end
end

% average displacement trials
displ_mtj_mean = nanmean(displ_MTJ,2); % nanmean averages 2 values if 3rd is NaN
displ_otj_mean = nanmean(displ_OTJ,2);
% set elong == 0 at force == 0
displ_mtj_mean(1) = 0;
displ_otj_mean(1) = 0;

% tweak for special trials
% not needed when setting elong == 0 above
%if strcmp(subject_id, 'Control 13 L PRE GM') || strcmp(subject_id,'INT_13_GM_PRE_STR_L')
%    displ_mtj_mean(1) = 0;
%end

% calculate elongation
tend_elong = displ_mtj_mean - displ_otj_mean;
% tend_elong(1) = NaN; % do not include point in stiffness fit

if plot_achilles
    % MTJ force-disp x 3
    plottitle = horzcat('MTJ displacement check, 3 single trials, ', subject_id);
    figure('Name',plottitle)
    hold on
    % 3 MTJ trials, after onset, to trial max force
    for i = 1:length(MTJ_trials)
        if length(MTJ_trials{i}) > 1
            plot(MTJ_trials{i}(:,loc_displ),MTJ_trials{i}(:,loc_force))
        else
            plot([0 0],[0 0])
        end
    end
    % mean MTJ, after onset, to common cutoff force (90%)
    plot(displ_mtj_mean(1:length(force_array_cut)), force_array_cut, 'k','LineWidth',2)
    % 3 MTJ, no onset correction, to trial max force
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
    % 3 MTJ trials, after onset, dots at force intervals, to common cutoff force
    set(gca,'ColorOrderIndex',1)
    for i = 1:length(MTJ_trials)
        if length(MTJ_trials{i}) > 1
            plot(displ_MTJ(1:length(force_array_cut),i), force_array_cut,'.','Markersize',8)
        else
            plot([0 0],[0 0])
        end
    end
    % visual
    xlabel('Displacement (mm)'),ylabel('Tendon force (N)'),title(plottitle,'Interpreter', 'none');
    axis([-1 Inf -100 Inf])
    legend('trial 1 to max','trial 2 to max','trial 3 to max','mean, to cutoff','trial 1 raw','Location','Northwest');
    saveas(gcf, strcat('data_plots_stiff/IND_stiff_avg_MTJ_', subject_id), 'png')
end

if plot_achilles
    % OTJ force-disp x 3
    plottitle = horzcat('OTJ displacement check, 3 single trials, ', subject_id);
    figure('Name',plottitle)
    hold on
    % 3 OTJ trials, after onset, to trial max force
    for i = 1:length(OTJ_trials)
        if length(OTJ_trials{i}) > 1
            plot(OTJ_trials{i}(:,loc_displ),OTJ_trials{i}(:,loc_force))
        else
            plot([0 0],[0 0])
        end
    end
    % mean OTJ, after onset, to common cutoff force (90%)
    plot(displ_otj_mean(1:length(force_array_cut)), force_array_cut, 'k','LineWidth',2)
    % 3 OTJ, no onset correction, to trial max force
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
    % 3 OTJ trials, after onset, dots at force intervals, to common cutoff force
    set(gca,'ColorOrderIndex',1)
    for i = 1:length(OTJ_trials)
        if length(OTJ_trials{i}) > 1
            plot(displ_OTJ(1:length(force_array_cut),i), force_array_cut,'.','Markersize',8)
        else
            plot([0 0],[0 0])
        end
    end
    % visual
    xlabel('Displacement (mm)'),ylabel('Tendon force (N)'),title(plottitle,'Interpreter', 'none');
    axis([-1 Inf -100 Inf])
    legend('trial 1 to max','trial 2 to max','trial 3 to max','mean, to cutoff','trial 1 raw','Location','Northwest');
    saveas(gcf, strcat('data_plots_stiff/IND_stiff_avg_OTJ_', subject_id), 'png')
end


%% choose cutoff method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% keep data up to 90% of 6-trials-common-force  OR  use manually set cutoff (datamster)

if force_array_cut(end) > str2double(force_cutoff_manual)
    cprintf('red',horzcat('Force cutoff is lowered manually (in datamaster), to ', force_cutoff_manual, ' N.\n'));
    loc_cutoff = find(force_array_cut>str2double(force_cutoff_manual),1,'first') - 1;
else % 90% of 6-trials-common-force
    loc_cutoff = length(force_array_cut);
end
% loc_cutoff: KEEP all points up until and including this frame

% create force-elongation array for later use
force_elong_array = [tend_elong(1:loc_cutoff) force_array_full(1:loc_cutoff)];



%% curve fitting for stiffness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% currently stiff fit is NOT forced through zero, but saving plots also through zero
% sending force_array_full, but utilizing for the fit, only up until loc_cutoff
[fitresult, gof] = fit_stiffness(tend_elong, force_array_full, loc_cutoff, displ_mtj_mean, displ_otj_mean);

% check if first coefficient is negative --> red text
coeffs = coeffvalues(fitresult);

% write stiffness coefficients to screen
if coeffs(1) < 0
    cprintf('*red', horzcat('Stiffness coeffs = ', regexprep(num2str(coeffvalues(fitresult)),' +',' '), '. R2 = ', num2str(gof.rsquare), '\n'))
elseif coeffs(1) > 70 %VAR
    cprintf('*red', horzcat('Stiffness coeffs = ', regexprep(num2str(coeffvalues(fitresult)),' +',' '), '. R2 = ', num2str(gof.rsquare), '\n'))
else
    cprintf('*blue', horzcat('Stiffness coeffs = ', regexprep(num2str(coeffvalues(fitresult)),' +',' '), '. R2 = ', num2str(gof.rsquare), '\n'))
end

end