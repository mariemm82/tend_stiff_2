%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract_isokinetic
% Marie Moltubakk 2.11.2016
% Read complete, prepared noraxon array
% Produce array with torque and angle
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [torque_max, torque_max_angle, torque_max_velocity, work_max, array_raw_output, array_intervals_output] = extract_isokinetic(noraxondata, side, trial_name, subject_id)
% current output = array containing one best trial. change to array_output_raw = full data series containing 3 trials
    


    %% prepare function settings
    torque_trial_onset_ID_threshold = 0.035; % increase in torque level between two adjacent values, 2x (10 frames apart), for determining that a trial has started
    torque_threshold_PF_initiation = 8; % how low should the avg torque be during the first part of isokinetic phase, in order to discard trial
    
    global column_norm_angle column_norm_torque % column_norm_velocity % column_l_gm column_r_gm column_l_gl column_r_gl column_l_sol column_r_sol column_gonio  column_l_tibant column_r_tibant  column_norm_velocity column_norm_direction column_achilles column_EMG_start column_EMG_end 
    global filepath
    global plot_individual plot_conversion
    global mat_version
    % plot_conversion --- plot figure showing initial IDing of peaks, and removal at each stage
    
    
    %% Read and prepare noraxon arrays
    % Read Noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
    % Produce a new noraxon data array
    freq_resample = 200;
    noraxon_prepped = read_noraxon_active(strcat(filepath, noraxondata), freq_resample, side, horzcat(subject_id, ': ', trial_name));
    
    % create separate arrays for analysis
    array_torque = noraxon_prepped(:,column_norm_torque);
    array_angle = noraxon_prepped(:,column_norm_angle);
    % array_velocity = noraxon_prepped(:,column_norm_velocity);
    % array_output_raw = [array_torque array_angle];
    
    % tweak for BD101 % dirty
    % incorrect offset @ isokinetic trials, due to imprecision during data collection
    if strcmp(subject_id(8:10), '101')
        offset_manual = max(array_angle) - 30;
        array_angle = array_angle - offset_manual;
    end
    
    
        
    %% find peaks
    
    % set minimum distance between peaks, according to velocity
    peakdistance = 80; %- large distance no longer required, since script below is very good at detecting+removing invalid peaks
%     if min(array_velocity) < -85      % trial 90 deg/s + BD 120 deg/s
%         peakdistance = 70;
%     elseif min(array_velocity) < -55  % trial 60 deg/s
%         peakdistance = 100;
%     elseif min(array_velocity) < -40  % trial 45 deg/s
%         peakdistance = 140;
%     elseif min(array_velocity) < -25  % trial 30 deg/s
%         peakdistance = 200; 
%     else
%         peakdistance = 200;
%     end
    
    % using inverted array_angle because this appears to work better for the sometimes short/inexistent stop in max dorsiflexion
    [val_peaks,loc_peaks] = findpeaks(-array_angle,'MinPeakDistance',peakdistance);
    % plot peaks and eliminated peaks
    if plot_conversion
        plottitle = horzcat('Isokinetic phase check, ', subject_id, ', ', trial_name);
        figure('Name',plottitle,'units','normalized','outerposition',[0 0 1 1])
        findpeaks(-array_angle,'MinPeakDistance',peakdistance) % --- will plot figure with peaks IDed
        hold on
        % right axis
        if strcmp(mat_version,'2015b') == 0
            yyaxis right
        end
        plot(array_torque,':')
        ylabel('Isokinetic torque (Nm)')
        if strcmp(mat_version,'2015b') == 0
            yyaxis left
        end
    end
    
    
    
    %% remove "peaks" identified at mid angles (often seen in starting phase)
    delete_peaks = [];
    delete_review = [];
    % define cutoff, what is a peak value (DF/PF) and what is "mid value"
    val_peaks_sort = sort(val_peaks,'descend');
    DF_peak = val_peaks_sort(3) - 4;
    PF_peak =  val_peaks_sort(end-2) + 3;
    % old:
%    DF_peak = max(-array_angle) - 4; % approx +10 deg - 4 deg = +6 deg (DF, inverted)
%    PF_peak = min(-array_angle) + 3; % approx -30 deg + 3 deg = -27 deg (PF, inverted)
    
    for i = 1:length(val_peaks)
        % if not a peak value
        if val_peaks(i) < DF_peak && val_peaks(i) > PF_peak
            % add to list for deletion
            delete_peaks = [delete_peaks, i];
            % if less than 7 degrees from max peak, save angle for potential text report about questionable deletions
            if (abs(val_peaks(i) - max(-array_angle)) < 7) %VAR
                delete_review = [delete_review val_peaks(i)];
            end
        end
    end
    % delete array entries selected above
    if plot_conversion
        plot(loc_peaks(delete_peaks),val_peaks(delete_peaks),'or')
    end
    val_peaks(delete_peaks) = [];
    loc_peaks(delete_peaks) = [];
    delete_peaks = [];
    
    
    
    %% remove multiple peaks at same stop phase
    for i = 2:length(val_peaks)
        peak1 = val_peaks(i-1);
        peak2 = val_peaks(i);
        if abs(peak1-peak2) < 6 && loc_peaks(i)-loc_peaks(i-1) < peakdistance*4 % two consequtive peaks are less than 6 degrees different + close in time --> same peak
            if peak1 < -15
                % delete FIRST for 30/PF positions
                delete_peaks = [delete_peaks, i-1];
            else % 10/DF positions
                % if second one has much higher angle, delete the first one
                if peak2-peak1 > 2
                    delete_peaks = [delete_peaks, i-1];
                else
                    % delete LAST
                    delete_peaks = [delete_peaks, i];
                end
            end
        end
    end
    % delete array entries selected above
    if plot_conversion
        plot(loc_peaks(delete_peaks),val_peaks(delete_peaks),'*c')
    end
    val_peaks(delete_peaks) = [];
    loc_peaks(delete_peaks) = [];
    delete_peaks = [];
    

    
    %% remove DF peak at beginning of DF trial, or PF peak at beginning of PF trial
    if strcmp(trial_name,'isokin DF 30')
        if val_peaks(1) > 0 % value near 10 degrees, e.g. dorsiflexed position (-10 degrees before angles were inverted for peak identification)
            if plot_conversion
                plot(loc_peaks(1),val_peaks(1),'pb')
            end
            val_peaks(1) = [];
            loc_peaks(1) = [];
        end
    else % PF trial
        if val_peaks(1) < -15 % value near -30 degrees, e.g. plantarflexed position (-30 degrees before angles were inverted for peak identification)
            if plot_conversion
                plot(loc_peaks(1),val_peaks(1),'pb')
            end
            val_peaks(1) = [];
            loc_peaks(1) = [];
        end
    end
    
    % write report if less than 6 peaks (3 full sets) remain + there are deletions of peaks near end range
    if length(loc_peaks) < 6 && ~isempty(delete_review)
        cprintf('red', horzcat('WARNING: Less than 3 DF peaks remain, deletions @ ', num2str(round(delete_review,2)), '°, max DF = ', num2str(round(max(-array_angle),2)), '°. Check manually.\n'))
    end
    
    
    
    %% refine identified peaks to the points where movement phases start/stop
    % save as:
    % --- loc_peak_start
    % --- loc_peak_end
    
    loc_peak_start(1:length(loc_peaks)) = zeros;
    loc_peak_end(1:length(loc_peaks)) = zeros;
    
    if loc_peaks(1) < 52 % tweak when first peak is detected VERY soon
        loc_peaks(1) = 52;
    end
    
    for i = 1:length(loc_peaks)
        
        % START of phase
        found = 0;
        j = 0;
        frames_apart = 40;
        while found == 0 && loc_peaks(i)+j+(frames_apart*2)+1 < length(array_angle) % compare two adjacent values, two times (x frames apart), both must change above threshold. E.g. 464-465 & 474-475
            if val_peaks(i) > 0 %+10 = dorsiflexed angle
                if array_angle(loc_peaks(i)+j) - array_angle(loc_peaks(i)+j+1) < -torque_trial_onset_ID_threshold && ...
                        array_angle(loc_peaks(i)+j+frames_apart) - array_angle(loc_peaks(i)+j+frames_apart+1) < -2*torque_trial_onset_ID_threshold && ...
                        min(diff(array_angle(loc_peaks(i)+j:loc_peaks(i)+j+(frames_apart*2)))) > -0.5*torque_trial_onset_ID_threshold % there is no large opposite movement in the next 80 frames
                        % the third requirement is special for DF, not needed for PF trials
                    loc_peak_start(i) = loc_peaks(i)+j-1; % -1 is mathematically correct. Changed to -3 to include the initial force production (phase has started)
                    found = 1;
                else % not found, check next
                    j = j+1;
                end
                
            else % val peak lower than 0 --> -30 = PF angle
                if array_angle(loc_peaks(i)+j) - array_angle(loc_peaks(i)+j+1) > torque_trial_onset_ID_threshold && ...
                        array_angle(loc_peaks(i)+j+frames_apart) - array_angle(loc_peaks(i)+j+frames_apart+1) > 2*torque_trial_onset_ID_threshold
                    loc_peak_start(i) = loc_peaks(i)+j-1; % -1 is mathematically correct. Changed to -3 to include the initial force production (phase has started)
                    found = 1;
                else % not found, check next
                    j = j+1;
                end
            end
        end
        
        if found == 0 % if no startpoint found
            % use end of array
            loc_peak_start(i) = length(array_angle);
        end
        
        % from previously used startpoint: go backwards to find END of previous phase
        found = 0;
        j = 0;
        while found == 0 && loc_peaks(i)-j-(frames_apart+1) > 0 % compare two adjacent values, two times (10 frames apart, both must change above threshold
            if abs(array_angle(loc_peaks(i)-j) - array_angle(loc_peaks(i)-j-1)) > torque_trial_onset_ID_threshold && ...
                    abs(array_angle(loc_peaks(i)-j-frames_apart) - array_angle(loc_peaks(i)-j-(frames_apart+1))) > 2*torque_trial_onset_ID_threshold
                % found a potential endpoint (movement in previous frames)
                % secondary check that 3 values are all ascending (alt. descending), to prevent registering miniscule peak (BD101 120deg/s)
                range_latest = mean(array_angle(loc_peaks(i)-j-1:loc_peaks(i)-j));
                range_mid = mean(array_angle(loc_peaks(i)-j-(frames_apart/2+1):loc_peaks(i)-j-frames_apart/2));
                range_before = mean(array_angle(loc_peaks(i)-j-(frames_apart+1):loc_peaks(i)-j-frames_apart));
                if range_latest > range_mid && range_mid > range_before
                    loc_peak_end(i) = loc_peaks(i)-j+0; % +0 is mathematically correct. Changed to +2 to include the initial force production (phase has started)
                    found = 1;
                elseif range_latest < range_mid && range_mid < range_before
                    loc_peak_end(i) = loc_peaks(i)-j+0; % +0 is mathematically correct. Changed to +2 to include the initial force production (phase has started)
                    found = 1;
                else
                    j = j+1;
                end
            else % not found, check next
                j = j+1;
            end
        end
        
        if found == 0 % if no startpoint found
            % use start of array
            loc_peak_end(i) = 1;
        end
        
    end

    
    
    %% remove double identical peaks (caused by two nearby, similar peaks resulting in same refined peak)
    for i = 1:length(loc_peak_start)-1
        if loc_peak_start(i) == loc_peak_start(i+1)
            delete_peaks = [delete_peaks, i];
        end
    end
	% plot & delete array entries selected above
    if ~isempty(delete_peaks)
        if plot_conversion
            plot(loc_peak_start(delete_peaks),val_peaks(delete_peaks),'xk')
        end
        cprintf('red', horzcat('WARNING: Deleted 2 identical START peaks, POST = ', num2str(loc_peak_end(delete_peaks)), ' vs ', num2str(loc_peak_end(delete_peaks+1)), '.\n'))
        loc_peak_start(delete_peaks) = [];
        loc_peak_end(delete_peaks) = [];
    end
    
    %% at this point, the 3 trials should be correctly identified
    
    % manual correction for one subject where peak ID is not correct:
    if strcmp(subject_id, 'INT_13_CON_PRE_R') == 1 && strcmp(trial_name,'isokin DF 30') % GOON
        loc_peak_start(2) = [];
        loc_peak_end(3) = [];
    end
    
    % plot refined peaks
    if plot_conversion
        for i=1:length(loc_peak_start)
            line([loc_peak_start(i) loc_peak_start(i)], [min(-array_angle) max(-array_angle)],'Color','g','LineWidth',1.5);
        end
        for i=1:2:length(loc_peak_start)-1
            line([loc_peak_start(i) loc_peak_end(i+1)], [min(-array_angle) min(-array_angle)],'Color','g','LineWidth',1.5,'LineStyle','--');
        end
        for i=1:length(loc_peak_end)
            line([loc_peak_end(i) loc_peak_end(i)], [min(-array_angle) max(-array_angle)],'Color','y');
        end
    end
    
    
    
    %% remove PF trials where torque drops low during the first phase (isokinetic movement has started, but subject is not working)
    delete_peaks_start = [];
    delete_peaks_end = [];
    if strcmp(trial_name,'isokin DF 30') == 0
        onset = ceil(peakdistance/3); % # of frames after trial starting point to consider
        for i = 1:2:length(loc_peak_start) % check every second entry = PF phases
            % remove peaks identifying phase which starts at the END of array:
            if loc_peak_start(i) == length(array_torque)
                delete_peaks_start = [delete_peaks_start, i];
            else
                % if torque drops lower than threshold set @ beginning of the script:
                % delete relevant peak entries
                % NB: loc_peak_start(1) belongs to loc_peak_end(2)
                if mean(array_torque(loc_peak_start(i):loc_peak_start(i)+onset)) < torque_threshold_PF_initiation 
                    if i == length(loc_peak_start)
                        % peak is the last: add to list for deletion
                        delete_peaks_start = [delete_peaks_start, i];
                    else
                        % peak is not last: add to list for deletion: current phase (e.g. PF) + next phase (e.g. DF)
                        delete_peaks_start = [delete_peaks_start, i, i+1];
                        if i+1 == length(loc_peak_end)
                            % peak is the last: add to list for deletion
                            delete_peaks_end = [delete_peaks_end, i+1];
                        else
                            % peak is not last: add to list for deletion: current phase (e.g. PF) + next phase (e.g. DF)
                            delete_peaks_end = [delete_peaks_end, i+1, i+2];
                        end
                    end
                end
            end
        end
        % delete array entries selected above
        if plot_conversion && ~isempty(delete_peaks_start)
            plot(loc_peak_start(delete_peaks_start),0,'hexagram','MarkerSize',12,'MarkerEdgeColor','black','MarkerFaceColor','r')
            if ~isempty(delete_peaks_end)
                plot(loc_peak_end(delete_peaks_end),0,'hexagram','MarkerSize',12,'MarkerEdgeColor','black','MarkerFaceColor','y')
            end
            % right axis
            if strcmp(mat_version,'2015b') == 0
                yyaxis right
            end
            line([0 length(array_angle)], [torque_threshold_PF_initiation torque_threshold_PF_initiation],'Color','r')
        end
        loc_peak_start(delete_peaks_start) = [];
        loc_peak_end(delete_peaks_end) = [];
    end    
    
    
    
    %% prepare for output
    
    if plot_conversion
        % PLOT horizontal green lines for trials which remain for analysis
        if strcmp(mat_version,'2015b') == 0
            yyaxis right
        end
        if strcmp(trial_name,'isokin DF 30') == 0 % PF trials
            for i = 1:2:length(loc_peak_end)-1
                plot(loc_peak_start(i):loc_peak_end(i+1), array_torque(loc_peak_start(i):loc_peak_end(i+1)),'r','LineWidth',2)
            end
        else % DF trials
            for i = 1:2:length(loc_peak_end)-1
                plot(loc_peak_start(i):loc_peak_end(i+1), array_torque(loc_peak_start(i):loc_peak_end(i+1)),'r','LineWidth',2)
            end
        end
        if strcmp(mat_version,'2015b') == 0
            yyaxis left
        end
        % axis([loc_peak_start(1)-200 loc_peak_end(end-1)+200 -Inf Inf])
        xlabel('Time (frames)')
        ylabel('Ankle angle (deg) INVERTED')
        title(plottitle,'Interpreter', 'none')
        legend('Angle','Peaks','midval','2cons','init/TQ','location','West')
        
        print(horzcat('data_plots/',plottitle),'-dpng')
    end

    % invert torque for dorsiflexion trial
    if strcmp(trial_name,'isokin DF 30')
        array_torque = -array_torque;
    end
    
        
    
    %% filter each trial
    % variable here is:
    % 20 hz --> 20/(noraxonfreq/2);
    % 0.05 = 37.5/(1500/2) --> 37.5 hz
    [B, A] = butter(8, 0.05, 'low'); %VAR
    
    % collect and filter three trials
    if length(loc_peak_start) >= 5 && length(loc_peak_end) >= 6
        trial1 = [filtfilt(B, A, array_torque(loc_peak_start(1):loc_peak_end(2))) array_angle(loc_peak_start(1):loc_peak_end(2))];
        trial2 = [filtfilt(B, A, array_torque(loc_peak_start(3):loc_peak_end(4))) array_angle(loc_peak_start(3):loc_peak_end(4))];
        trial3 = [filtfilt(B, A, array_torque(loc_peak_start(5):loc_peak_end(6))) array_angle(loc_peak_start(5):loc_peak_end(6))];
        trials = {trial1 trial2 trial3};
    elseif length(loc_peak_start) >= 3 && length(loc_peak_end) >= 4
        trial1 = [filtfilt(B, A, array_torque(loc_peak_start(1):loc_peak_end(2))) array_angle(loc_peak_start(1):loc_peak_end(2))];
        trial2 = [filtfilt(B, A, array_torque(loc_peak_start(3):loc_peak_end(4))) array_angle(loc_peak_start(3):loc_peak_end(4))];
        trials = {trial1 trial2};
    elseif length(loc_peak_start) >= 1 && length(loc_peak_end) >= 2
        trial1 = [filtfilt(B, A, array_torque(loc_peak_start(1):loc_peak_end(2))) array_angle(loc_peak_start(1):loc_peak_end(2))];
        trials = {trial1};
    else % if no valid trials remain...
        cprintf('red*', horzcat('WARNING, ', trial_name, ': No remaining trials! Review peakcheck plot for errors.\n'))
        return
    end

    
    
    %% calculate output variables
    
    % peak torque (best of valid trials)
    
    torque_peaks(1:length(trials),1:2) = zeros;
    for i = 1:length(trials)
        [torque_peaks(i,1),torque_peaks(i,2)] = max(trials{i}(:,1));
    end
    [torque_max,torque_best_i] = max(torque_peaks(:,1));
    torque_max_angle = trials{torque_best_i}(torque_peaks(torque_best_i,2),2);
    torque_max_velocity = str2double(trial_name(11:end)); % or switch to below
    % torque_max_velocity = array_velocity (index); % must add velocity to trials array and repeat as for torque_max_angle, if this is to be used
    
    
    % array of best trial (highest peak torque)
    
    % select the best trial
    array_raw_output = trials{torque_best_i};
    
    % find min and max angles
    angles_min = min(array_raw_output(:,2));
    angles_max = max(array_raw_output(:,2));
    
    % spline for common angles @ 2.5 deg
    angles_common2 = (-7.5:2.5:27.5)';
    array_intervals_output = spline(array_raw_output(:,2), array_raw_output(:,1), angles_common2);

    % if the trial does not have as large range as spline has:
    if angles_min > angles_common2(2)
        % replace interpolated data
        array_intervals_output(1:2) = NaN;
    elseif angles_min > angles_common2(1)
        array_intervals_output(1) = NaN;
    end
    if angles_max < angles_common2(end)
        % replace interpolated data
        array_intervals_output(end) = NaN;
    end
    
     % plot to confirm no smoothing:
%      figure,plot(array_raw_output(:,2),array_raw_output(:,1))
%      hold on
%      plot(angles_common2,array_intervals_output,':ko', 'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0])
    
    
    % work (best work of valid trials)
    
    work(1:length(trials)) = zeros;
    torque_angle_reshaped{length(trials)} = [];
    for i=1:length(trials)
        %reshape torque-angle to 0.1 deg intervals
        work_angle_start = ceil(10*min(trials{i}(:,2)))/10;
        work_angle_stop = floor(10*max(trials{i}(:,2)))/10;
        torque_angle_reshaped{i} = (work_angle_start:0.1:work_angle_stop)'; %VAR
        torque_angle_reshaped{i}(:,2) = spline(trials{i}(:,2), trials{i}(:,1), torque_angle_reshaped{i}(:,1));
        
        % remove entries of negative torque at beginning/end
        torque_pos = torque_angle_reshaped{i}(:,2);
        torque_pos = torque_pos(torque_pos>=0);
        work(i) = mean(torque_pos)*0.5*pi();
    end
    
    work_max = max(work);
    
    
    
    %% checkpoint plot
    
     if plot_individual
        % plotting reshaped torque instead of trials{} which is output to the main method, as a check that data which are basis for WORK calculations are correct
        plottitle = horzcat(horzcat(trial_name, ', ', subject_id));
        figure('Name',plottitle)
        plot(torque_angle_reshaped{1}(:,1),torque_angle_reshaped{1}(:,2))
        hold on
        if length(loc_peak_start) >= 3 && length(loc_peak_end) >= 4
            plot(torque_angle_reshaped{2}(:,1),torque_angle_reshaped{2}(:,2))
        end
        if length(loc_peak_start) >= 5 && length(loc_peak_end) >= 6
            plot(torque_angle_reshaped{3}(:,1),torque_angle_reshaped{3}(:,2))
        end
        plot(angles_common2,array_intervals_output,':ko', 'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'MarkerSize',3)
        ax = gca;
        if strcmp(trial_name,'isokin DF 30')
            set(ax, 'xdir','reverse')
        end
        xlabel('Ankle angle (deg)')
        ylabel('Isokinetic torque (Nm)')
        title(plottitle,'Interpreter', 'none')
        %legend('labels','location','SouthWest')

        print(horzcat('data_plots/',plottitle),'-dpng')
    end
    
 
    
    %% print report
    cprintf('blue', horzcat(trial_name, ': Peak torque = ', num2str(round(torque_max,1)) ,' Nm, angle = ', num2str(round(torque_max_angle,1)) ,'°, work = ', num2str(round(work_max,1)) ,' J.\n'))
    
end