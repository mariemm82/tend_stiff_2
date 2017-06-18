%%%%%%%%%%%%%%%%%%%%%%%%%%
% final_stiffness
% Marie Moltubakk 16.6.2013
% Read 3+3 trials of MTJ and OTJ scans
% Produce stiffness coeffisients and output arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%



function [fitresult, gof, force_elong_array, final_cutoff_force] = final_stiffness(time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3, time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3, forceintervals, force_cutoff_manual, force_max_trials, ~) %new2014-04-14

global plot_achilles subject_id % plot_conversion plot_check %plot_norm plot_emg


%% convert input data to cells
MTJ_trials = {time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3};
OTJ_trials = {time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3};
loc_time = 1;
loc_force = 2;
loc_displ = 3;


%% find common FORCE level, create common force array

% print maximal force from each trial to the screen
cprintf('black', horzcat('Ramp trials, ind max force = '))
cprintf('*blue', horzcat(num2str(round(min(force_max_trials))), ' N'))
cprintf('black', horzcat(', trials = ', num2str(round(force_max_trials(1))), ' - ', num2str(round(force_max_trials(2))), ' - ', num2str(round(force_max_trials(3))), ' - ', num2str(round(force_max_trials(4))), ' - ', num2str(round(force_max_trials(5))), ' - ', num2str(round(force_max_trials(6))), ' N.\n'))

% find lowest common force between all 6 trials
commonforce = min(force_max_trials);

% create array of force values to use for averaging
force_array = (0:forceintervals:commonforce)';


%% find FORCE onset
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
baseline_frames = 10; %VAR X first frames used to define baseline
baseline_SD_multiplier = 2; %VAR onset is when force level increases to 3*SD above baseline force
% remainder is performed inside loop below

% preallocate
threshold_MTJ(1:3) = NaN;
threshold_OTJ(1:3) = NaN;
loc_MTJ_onset(1:3) = NaN;
loc_OTJ_onset(1:3) = NaN;

% locate onset - output variable is the first frame to KEEP
for i = 1:length(MTJ_trials)
    if length(MTJ_trials{i}) > 1
        % in baseline determination zone, if force goes up by an average of 2 N from point to point - baseline is onset
        if mean(diff(MTJ_trials{i}(1:baseline_frames,loc_force))) > 2 %VAR
            cprintf('red', horzcat('WARNING: MTJ', num2str(i), 'Force increasing during baseline, by avg of 2 N/frame. Set baseline = first frame =', num2str(MTJ_trials{i}(1,loc_force)), ' N.\n'))
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
        % in baseline determination zone, if force goes up by an average of 2 N from point to point - baseline is onset
        if mean(diff(OTJ_trials{i}(1:baseline_frames,loc_force))) > 2 %VAR
            cprintf('red', horzcat('WARNING: OTJ', num2str(i), 'Force increasing during baseline, by avg of 2 N/frame. Set baseline = first frame =', num2str(OTJ_trials{i}(1,loc_force)), ' N.\n'))
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
    plottitle = horzcat('Force onset check, ', subject_id);
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
    title(plottitle)
    saveas(gcf, strcat('data_plots_stiff/IND_force_onset_', subject_id), 'png')
    [latestonset,loc_latestonset] = (max([loc_MTJ_onset loc_OTJ_onset]));
    if loc_latestonset > 3
        axis([0 0.2+OTJ_trials{loc_latestonset-3}(latestonset,loc_time) -20 max([threshold_MTJ threshold_OTJ])*10])
    else
        axis([0 0.2+MTJ_trials{loc_latestonset}(latestonset,loc_time) -20 max([threshold_MTJ threshold_OTJ])*10])
    end
    saveas(gcf, strcat('data_plots_stiff/IND_force_onset_', subject_id, '_ZOOM'), 'png')
end


%% CUT CELLS at onset, and offset time+displacement accordingly
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
        % warning if first force is larger than ... 20 N? TODO, delete this warning?
        if MTJ_trials{i}(1,loc_force) > 20 %VAR
            cprintf('red', horzcat('WARNING: High force onset for MTJ', num2str(i), ': ', num2str(round(MTJ_trials{i}(1,loc_force),3)), 'N\n'))
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
        % warning if first force is larger than ... 20 N?
        if OTJ_trials{i}(1,loc_force) > 20 %VAR
            cprintf('red', horzcat('WARNING: High force onset for OTJ', num2str(i), ': ', num2str(round(OTJ_trials{i}(1,loc_force),3)), 'N\n'))
        end
    end
end




%% calculate tendon ELONGATION

% preallocate
displ_MTJ(1:length(force_array),1:length(MTJ_trials)) = NaN;
displ_OTJ(1:length(force_array),1:length(OTJ_trials)) = NaN;

for i = 1:length(MTJ_trials)
    if length(MTJ_trials{i}) > 1
        % detect position of max force in array
        [~,loc_del_mtj] = max(MTJ_trials{i}(:,loc_force));
        % delete array contents AFTER (+1) max force
        MTJ_trials{i}(loc_del_mtj+1:end,:) = [];
        % spline to get displacement values at common force levels
        displ_MTJ(:,i) = pchip(MTJ_trials{i}(:,loc_force),MTJ_trials{i}(:,loc_displ),force_array); % ALT: interp1(MTJ_trials{i}(:,loc_force),MTJ_trials{i}(:,loc_displ),force_array,'linear'); % ALT: spline(MTJ_trials{i}(:,loc_force),MTJ_trials{i}(:,loc_displ),force_array);
        if abs(displ_MTJ(1,i)) > 0.05 % delete first value if spline function creates a point way off (will not be used anyway)
            displ_MTJ(1,i) = NaN;
        end
        if abs(displ_MTJ(end,i)) > 50
            displ_MTJ(end,i) = NaN;
        end
    end
end

for i = 1:length(OTJ_trials)
    if length(OTJ_trials{i}) > 1
        % detect position of max force in array
        [~,loc_del_otj] = max(OTJ_trials{i}(:,loc_force));
        % delete array contents  AFTER (+1) max force
        OTJ_trials{i}(loc_del_otj+1:end,:) = [];
        % spline to get displacement values at common force levels
        displ_OTJ(:,i) = pchip(OTJ_trials{i}(:,loc_force),OTJ_trials{i}(:,loc_displ),force_array);
        if abs(displ_OTJ(1,i)) > 0.05
            displ_OTJ(1,i) = NaN;
        end
        if abs(displ_OTJ(end,i)) > 50
            displ_OTJ(end,i) = NaN;
        end
    end
end

% average displacement trials, calculate elongation
displ_mtj_mean = nanmean(displ_MTJ,2); % nanmean to average 2 values if 3rd is NaN
displ_otj_mean = nanmean(displ_OTJ,2);
tend_elong = displ_mtj_mean - displ_otj_mean;
tend_elong(1) = NaN; % do not include point in stiffness fit

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
    saveas(gcf, strcat('data_plots_stiff/IND_stiff_avg_MTJ_', subject_id), 'png')
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
    saveas(gcf, strcat('data_plots_stiff/IND_stiff_avg_OTJ_', subject_id), 'png')
end


%% calculate rate of force development (RFD) - cutoff method 4, June 2017 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preallocate
RFD_MTJ_raw = cell(1,3);
RFD_OTJ_raw = cell(1,3);
RFD_MTJ_smooth = cell(1,3);
RFD_OTJ_smooth = cell(1,3);
freq_MTJ(1:3) = NaN;
freq_OTJ(1:3) = NaN;

% create running average filter with x points:
filter_points = 5; %VAR
filter_run_avg = 1/filter_points * ones(filter_points,1);

% calculate and smooth RFD
for i = 1:length(MTJ_trials)
    if length(MTJ_trials{i}) > 1
        RFD_MTJ_raw{i} = diff(MTJ_trials{i}(:,loc_force))./diff(MTJ_trials{i}(:,loc_time));
        RFD_MTJ_smooth{i} = filtfilt(filter_run_avg, 1, RFD_MTJ_raw{i});
        freq_MTJ(i) = 1/mean(diff(MTJ_trials{i}(:,loc_time)));
    end
end
for i = 1:length(OTJ_trials)
    if length(OTJ_trials{i}) > 1
        RFD_OTJ_raw{i} = diff(OTJ_trials{i}(:,loc_force))./diff(OTJ_trials{i}(:,loc_time));
        RFD_OTJ_smooth{i} = filtfilt(filter_run_avg, 1, RFD_OTJ_raw{i});
        freq_OTJ(i) = 1/mean(diff(OTJ_trials{i}(:,loc_time)));
    end
end
% RFD raw and filt are one frame shorter than MTJ_trials etc, because of
% differentiation.

% if plot_achilles
% %    % check for matlab version, do not use yyaxis on home computer (2015)
% %    mat_version = version('-release');
%     plottitle = horzcat('Rate of force development, 6 trials, ', subject_id);
%     figure('Name',plottitle)
%     hold on
% %     % left axis
% %     if strcmp(mat_version,'2015b') == 0
% %         yyaxis left
% %     end
%     for i = 1:length(MTJ_trials)
%         if length(MTJ_trials{i}) > 1
%             plot(MTJ_trials{i}(:,loc_time),MTJ_trials{i}(:,loc_force),'-')
%         else
%         end
%     end
%     for i = 1:length(OTJ_trials)
%         if length(OTJ_trials{i}) > 1
%             plot(OTJ_trials{i}(:,loc_time),OTJ_trials{i}(:,loc_force),'-')
%         else
%         end
%     end
%     ylabel('Force (N)')
%     set(gca,'ColorOrderIndex',1)
% %     % right axis if possible
% %     if strcmp(mat_version,'2015b') == 0
% %         yyaxis right
% %         ylabel('Rate of force development (N/s)')
% %     end
%     for i = 1:length(MTJ_trials)
%         if length(MTJ_trials{i}) > 1
%             plot(MTJ_trials{i}(2:end,loc_time),RFD_MTJ_smooth{i},':')
%         else
%         end
%     end
%     for i = 1:length(OTJ_trials)
%         if length(OTJ_trials{i}) > 1
%             plot(OTJ_trials{i}(2:end,loc_time),RFD_OTJ_smooth{i},':')
%         else
%         end
%     end
%     % horizontal lines @ various RFD levels
%     plot([0 4],[200 200],'k-')
%     plot([0 4],[100 100],'k-') %old instead of Inf: max(cell2mat(cellfun(@size,time_array,'uni',false))) * timeintervals
%     % visual
%     xlabel('Time (s)')
%     axis([-Inf Inf -100 commonforce])
%     title(plottitle)
%     saveas(gcf, strcat('data_plots_stiff/IND_RFD_', subject_id), 'png')
% end


%% calculations of cutoff based on RFD

RFD_MTJ_smooth_SD(1:3) = NaN;
loc_startforce_MTJ(1:3) = NaN;
loc_halfforce_MTJ(1:3) = NaN;
RFD_OTJ_smooth_SD(1:3) = NaN;
loc_startforce_OTJ(1:3) = NaN;
loc_halfforce_OTJ(1:3) = NaN;
startforce = 0.1; %VAR 10% of common force
halfforce = 0.5; %VAR 50% of common force

% calculate SD of RFD in stable region of force development
for i = 1:length(RFD_MTJ_smooth)
    if length(RFD_MTJ_smooth{i}) > 1
        loc_startforce_MTJ(i) = find(MTJ_trials{i}(:,loc_force) >= (commonforce*startforce), 1, 'first');
        loc_halfforce_MTJ(i) = find(MTJ_trials{i}(:,loc_force) >= (commonforce*halfforce), 1, 'first');
        RFD_MTJ_smooth_SD(i) = std(RFD_MTJ_smooth{i}(loc_startforce_MTJ(i):loc_halfforce_MTJ(i)));
    end
end
for i = 1:length(RFD_OTJ_smooth)
    if length(RFD_OTJ_smooth{i}) > 1
        loc_startforce_OTJ(i) = find(OTJ_trials{i}(:,loc_force) >= (commonforce*startforce), 1, 'first');
        loc_halfforce_OTJ(i) = find(OTJ_trials{i}(:,loc_force) >= (commonforce*halfforce), 1, 'first');
        RFD_OTJ_smooth_SD(i) = std(RFD_OTJ_smooth{i}(loc_startforce_OTJ(i):loc_halfforce_OTJ(i)));
    end
end

cprintf('black', horzcat('SD for RFD =  ', num2str(round(RFD_MTJ_smooth_SD)), '  N/s for MTJ  -  ', num2str(round(RFD_OTJ_smooth_SD)), '  N/s for OTJ.\n'))

if plot_achilles
    plottitle = horzcat('RFD and SD check, MTJ trials, ', subject_id);
    figure('Name',plottitle)
    hold on
    % force
    for i=1:3
        if length(MTJ_trials{i}) > 1
            plot(MTJ_trials{i}(:,loc_time),MTJ_trials{i}(:,loc_force),':')
        end
    end
    l1 = plot(0, 0,':');
    set(gca,'ColorOrderIndex',1)
    % RFD
    for i=1:3
        if length(MTJ_trials{i}) > 1
            plot(MTJ_trials{i}(2:end,loc_time),RFD_MTJ_smooth{i},'--','Linewidth',2)
        end
    end
    l2 = plot(0,0,'.--');
    % force levels
    [loc_latest_halfforce,loc_latest_trial] = max(loc_halfforce_MTJ);
    l3 = plot([0 MTJ_trials{loc_latest_trial}(loc_latest_halfforce,loc_time)],[(commonforce*startforce) (commonforce*startforce)],'k'); % horizontal line from zero to latest timepoint of halfforce
    plot([0 MTJ_trials{loc_latest_trial}(loc_latest_halfforce,loc_time)],[(commonforce*halfforce) (commonforce*halfforce)],'k');
    % SDs
    set(gca,'ColorOrderIndex',1)
    for i=1:3
        if length(MTJ_trials{i}) > 1
            plot([MTJ_trials{i}(loc_startforce_MTJ(i),loc_time) MTJ_trials{i}(loc_halfforce_MTJ(i),loc_time)],[RFD_MTJ_smooth_SD(i) RFD_MTJ_smooth_SD(i)],'Linewidth',2)
        end
    end
    l4 = plot(0,0,'-','Linewidth',2);
    % visual
    legend([l1 l2 l3 l4], 'Force', 'RFD', 'Region 10-50% force', 'SD of RFD in region', 'Location','Northwest')
    title(plottitle)
    axis([-Inf Inf -100 Inf])
    saveas(gcf, strcat('data_plots_stiff/IND_RFD_SD_', subject_id, ' MTJ'), 'png')
end

if plot_achilles
    plottitle = horzcat('RFD and SD check, OTJ trials, ', subject_id);
    figure('Name',plottitle)
    hold on
    % force
    for i=1:3
        if length(OTJ_trials{i}) > 1
            plot(OTJ_trials{i}(:,loc_time),OTJ_trials{i}(:,loc_force),':')
        end
    end
    l1 = plot(0, 0,':');
    set(gca,'ColorOrderIndex',1)
    % RFD
    for i=1:3
        if length(OTJ_trials{i}) > 1
            plot(OTJ_trials{i}(2:end,loc_time),RFD_OTJ_smooth{i},'--','Linewidth',2)
        end
    end
    l2 = plot(0,0,'.--');
    % force levels
    [loc_latest_halfforce,loc_latest_trial] = max(loc_halfforce_OTJ);
    l3 = plot([0 OTJ_trials{loc_latest_trial}(loc_latest_halfforce,loc_time)],[(commonforce*startforce) (commonforce*startforce)],'k'); % horizontal line from zero to latest timepoint of halfforce
    plot([0 OTJ_trials{loc_latest_trial}(loc_latest_halfforce,loc_time)],[(commonforce*halfforce) (commonforce*halfforce)],'k');
    % SDs
    set(gca,'ColorOrderIndex',1)
    for i=1:3
        if length(OTJ_trials{i}) > 1
            plot([OTJ_trials{i}(loc_startforce_OTJ(i),loc_time) OTJ_trials{i}(loc_halfforce_OTJ(i),loc_time)],[RFD_OTJ_smooth_SD(i) RFD_OTJ_smooth_SD(i)],'Linewidth',2)
        end
    end
    l4 = plot(0,0,'-','Linewidth',2);
    % visual
    legend([l1 l2 l3 l4], 'Force', 'RFD', 'Region 10-50% force', 'SD of RFD in region', 'Location','Northwest')
    title(plottitle)
    saveas(gcf, strcat('data_plots_stiff/IND_RFD_SD_', subject_id, ' OTJ'), 'png')
end

% find cutoff point if RFD drops with X SD ---  MMM TODO

% if RFD_MTJ_smooth_SD(1:3) goes lower than x times RFD_OTJ_smooth_SD(i)
% cut MTJ_trials{i}(:,loc_force)
% OR apply 90% elong

% lowest point of 6 individual trials - or after averaging F-elong?




%% calculate 90% of max elongation - cutoff method 3, June 2014 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find 90% of max elongation
percent_elong = 0.90; %VAR
[max_elong,loc_max_elong] = max(tend_elong);
keep_elong = percent_elong * max_elong;
loc_cutoff_elong = find(tend_elong>keep_elong,1,'first') - 1; % finds the first value larger than 90% limit - we will use the largest BEFORE 90% limit




%% choose cutoff method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% todo: compare 
% loc_cutoff_elong
% vs
% loc_cutoff_RFD
% vs
% loc from manual cutoff
% use lowest incident
% print chosen method + values to screen
loc_cutoff_chosen = loc_cutoff_elong;


% cut off even earlier if datamaster is set to cut at lower value
if force_array(loc_cutoff_chosen) > str2double(force_cutoff_manual)
    cprintf('red',horzcat('Force cutoff is lowered manually (datamaster), to ', force_cutoff_manual, ' N.\n'));
    loc_cutoff_chosen = find(force_array>str2double(force_cutoff_manual),1,'first') - 1;
end

% write cutoff to screen
cprintf('black', horzcat('Max elong = ', num2str(round(max_elong,2)), ' mm @ ', num2str(force_array(loc_max_elong)), ' N, ', num2str(percent_elong*100), '%% elong -> cutoff at '))
cprintf('*blue', horzcat(num2str(force_array(loc_cutoff_chosen)), ' N.\n'))

% save the applied max force level, for output to file
final_cutoff_force = force_array(loc_cutoff_chosen);

% create force-elongation array for later use
force_elong_array = [tend_elong(1:loc_cutoff_chosen) force_array(1:loc_cutoff_chosen)];



%% curve fitting for stiffness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% currently stiff fit not through zero, but saving plots also through zero
[fitresult, gof] = fit_stiffness(tend_elong, force_array, loc_cutoff_chosen, displ_mtj_mean, displ_otj_mean);

% write stiffness coefficients to screen
cprintf('*blue', horzcat('Stiffness coeffs = ', regexprep(num2str(coeffvalues(fitresult)),' +',' '), '. R2 = ', num2str(gof.rsquare), '\n'))


end