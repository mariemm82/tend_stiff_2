%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract_isokinetic
% Marie Moltubakk 2.11.2016
% Read complete, prepared noraxon array
% Produce array with torque and angle
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [torque_max, torque_max_angle, torque_max_velocity, work_max, array_output] = extract_isokinetic(noraxondata, side, trial_name, subject_id)
% current output = array containing one best trial. change to array_output_raw = full data series containing 3 trials

    global column_norm_angle column_norm_torque column_norm_velocity % column_l_gm column_r_gm column_l_gl column_r_gl column_l_sol column_r_sol column_gonio  column_l_tibant column_r_tibant  column_norm_velocity column_norm_direction column_achilles column_EMG_start column_EMG_end 
    global filepath
    global plot_check plot_individual
    
    
    % Read Noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
    % Produce a new noraxon data array
    freq_resample = 200;
    noraxon_prepped = read_noraxon_active(strcat(filepath, noraxondata), freq_resample, side, horzcat(subject_id, ': ', trial_name));
    
    % array data
    array_torque = noraxon_prepped(:,column_norm_torque);
    array_angle = noraxon_prepped(:,column_norm_angle);
    array_velocity = noraxon_prepped(:,column_norm_velocity);
    % array_output_raw = [array_torque array_angle];

    
    
    %%% identify movement phases
    % separate trials
    % add additional filter per trial (electrical noise removal)
    
    % set minimum distance between peaks, according to velocity
    if min(array_velocity) < -85      % trial 90 deg/s + BD 120 deg/s
        peakdistance = 70;
    elseif min(array_velocity) < -55  % trial 60 deg/s
        peakdistance = 100;
    elseif min(array_velocity) < -40  % trial 45 deg/s
        peakdistance = 140;
    elseif min(array_velocity) < -25  % trial 30 deg/s
        peakdistance = 200;
    else
        peakdistance = 200;
    end
    
    
    
    % find peaks
    
    % using inverted array_angle because this appears to work better for the sometimes short/inexistent stop in max dorsiflexion
    [val_peaks,loc_peaks] = findpeaks(-array_angle,'MinPeakDistance',peakdistance);
    % findpeaks(-array_angle,'MinPeakDistance',peakdistance) % --- will plot figure with peaks IDed
    
    
    
    % remove incorrect peaks

    % remove "peaks" at mid value (e.g. starting phase) in inverted data
    delete_peaks = [];
    delete_review = [];
    % define cutoff for peak value (DF/PF) vs "mid value"
    DF_peak = max(-array_angle) - 4; % approx +10 deg - 4 deg = +6 deg (DF, inverted)
    PF_peak = min(-array_angle) + 3; % approx -30 deg + 4 deg = -26 deg (PF, inverted)
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
    val_peaks(delete_peaks) = [];
    loc_peaks(delete_peaks) = [];
    delete_peaks = [];
    
    % remove multiple peaks at same stop phase
    for i = 2:length(val_peaks)
        peak1 = val_peaks(i-1);
        peak2 = val_peaks(i);
        if abs(peak1-peak2) < 6 % two consequtive peaks are less than 6 degrees different --> same peak
            delete_peaks = [delete_peaks, i-1];
        end
    end
    % delete array entries selected above
    val_peaks(delete_peaks) = [];
    loc_peaks(delete_peaks) = [];
    
    % after invalid peaks have been removed:
    % remove DF peak at beginning of DF trial, or opposite
    if strcmp(trial_name,'isokin DF 30')
        if val_peaks(1) > 0 % value near 10 degrees, e.g. dorsiflexed position (-10 degrees before angles were inverted for peak identification)
            %val_peaks(1) = [];
            loc_peaks(1) = [];
        end
    else % PF trial
        if val_peaks(1) < -15 % value near -30 degrees, e.g. plantarflexed position (-30 degrees before angles were inverted for peak identification)
            %val_peaks(1) = [];
            loc_peaks(1) = [];
        end
    end
    
    % write report if less than 6 peaks (3 full sets) remain + there are deletions of peaks near end range
    if length(loc_peaks) < 6 && ~isempty(delete_review) %VAR
        cprintf('red', horzcat('WARNING: Less than 3 DF peaks remain, deletions @ ', num2str(round(delete_review,2)), '°, max DF = ', num2str(round(max(-array_angle),2)), '°. Check manually.\n'))
    end

    
    % refine peak points to movement phase start/stop
    threshold = 0.035; % VAR
    loc_peak_start(1:length(loc_peaks)) = zeros;
    loc_peak_end(1:length(loc_peaks)) = zeros;
    if loc_peaks(1) < 12
        loc_peaks(1) = 12; % tweak when first peak is detected VERY soon
    end
    for i = 1:length(loc_peaks)
        % START of phase
        found = 0;
        j = 0;
        while found == 0 && loc_peaks(i)+j+11 < length(array_angle) % compare two adjacent values, two times (10 frames apart, both must change above threshold
            if abs(array_angle(loc_peaks(i)+j) - array_angle(loc_peaks(i)+j+1)) > threshold && abs(array_angle(loc_peaks(i)+j+10) - array_angle(loc_peaks(i)+j+11)) > 2*threshold
                loc_peak_start(i) = loc_peaks(i)+j-1; % -1 is mathematically correct. Changed to -3 to include the initial force production (phase has started)
                found = 1;
            else % not found, check next
                j = j+1;
            end
            % if no startpoint found
        end
        if found == 0
            % use end of array
            loc_peak_start(i) = length(array_angle);
        end
        % END of phase
        found = 0;
        j = 0;
        while found == 0 && loc_peaks(i)-j-11 > 0 % compare two adjacent values, two times (10 frames apart, both must change above threshold
            if abs(array_angle(loc_peaks(i)-j) - array_angle(loc_peaks(i)-j-1)) > threshold && abs(array_angle(loc_peaks(i)-j-10) - array_angle(loc_peaks(i)-j-11)) > 2*threshold
                loc_peak_end(i) = loc_peaks(i)-j+0; % +0 is mathematically correct. Changed to +2 to include the initial force production (phase has started)
                found = 1;
            else % not found, check next
                j = j+1;
            end
        end
        % if no startpoint found
        if found == 0
            % use start of array
            loc_peak_end(i) = 1;
        end
    end
    
    % invert data for dorsiflexion trial
    if strcmp(trial_name,'isokin DF 30')
        array_torque = -array_torque;
    end
    
    if plot_individual
        plottitle = horzcat('Isokinetic phase check, ', horzcat(subject_id, ': ', trial_name));
        figure('Name',plottitle)
        hold on
        yyaxis left
        findpeaks(-array_angle,'MinPeakDistance',peakdistance)
        yyaxis right
        plot(array_torque)
        ylabel('Isokinetic torque (Nm)')
        yyaxis left
        for i=1:length(loc_peak_start)
            line([loc_peak_start(i) loc_peak_start(i)], [min(-array_angle) max(-array_angle)],'Color','g');
            line([loc_peak_end(i) loc_peak_end(i)], [min(-array_angle) max(-array_angle)],'Color','y');
        end
        xlabel('Time (frames)')
        ylabel('Ankle angle (deg) INVERTED')
        title(plottitle)
        legend('Angle','Peaks','Torque','location','NorthEast')
    end
    
    % prepare filter
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
    else
        trial1 = [filtfilt(B, A, array_torque(loc_peak_start(1):loc_peak_end(2))) array_angle(loc_peak_start(1):loc_peak_end(2))];
        trials = {trial1};
    end

    
    
    
    %%% calculate output variables
    
    % peak torque
    
    torque_peaks(length(trials),2) = zeros;
    for i = 1:length(trials)
        [torque_peaks(i,1),torque_peaks(i,2)] = max(trials{i}(:,1));
    end
    [torque_max,torque_best_i] = max(torque_peaks(:,1));
    torque_max_angle = trials{torque_best_i}(torque_peaks(torque_best_i,2),2);
    torque_max_velocity = str2double(trial_name(11:end)); % or switch to below
    % torque_max_velocity = array_velocity (index); % must add velocity to trials array and repeat as for torque_max_angle, if this is to be used
    
    array_output = trials{torque_best_i};
    
    
    
    % work
    
    work(length(trials)) = zeros;
    torque_angle_reshaped{length(trials)} = [];
    for i=1:length(trials)
        %reshape
        work_angle_start = ceil(10*min(trials{i}(:,2)))/10;
        work_angle_stop = floor(10*max(trials{i}(:,2)))/10;
        torque_angle_reshaped{i} = (work_angle_start:0.1:work_angle_stop)';
        torque_angle_reshaped{i}(:,2) = spline(trials{i}(:,2), trials{i}(:,1), torque_angle_reshaped{i}(:,1));
        
        % remove negative torques
        torque_pos = torque_angle_reshaped{i}(:,2);
        torque_pos = torque_pos(torque_pos>=0);
        work(i) = mean(torque_pos)*0.5*pi();
    end
    
    work_max = max(work);
    
    
    
    %%% checkpoint plot
    
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
        ax = gca;
        if strcmp(trial_name,'isokin DF 30')
            set(ax, 'xdir','reverse')
        end
        xlabel('Ankle angle (deg)')
        ylabel('Isokinetic torque (Nm)')
        title(plottitle)
        %legend('øeg','location','SouthWest')
        %saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    
    
    %%% print report
    
    cprintf('blue', horzcat(trial_name, ': Peak torque = ', num2str(round(torque_max,2)) ,' Nm, angle = ', num2str(round(torque_max_angle,2)) ,'°, work = ', num2str(round(work_max,2)) ,' J.\n'))
end