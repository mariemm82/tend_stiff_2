%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_angle_constants_active
% Marie Moltubakk 1.12.2014
% Read raw torque data from Norm, filter
% Produce individual conversion factors for ANGLE from ISOKINETIC trial
%
% Note 1: Some variables are defined directly inside this function, see %VAR
%%%%%%%%%%%%%%%%%%%%%%%%%%




function [convert_ind_angle_a, convert_ind_angle_b] = calculate_angle_constants_active(freq_cutoff, inputfile)
    global plot_conversion plot_check subject_id 
    global dm_side
    global column_norm_angle
    global norm_volt_per_degree % norm_volt_per_nm norm_volt_per_velocity norm_mv2nm_a norm_mv2nm_b
        
    % import noraxon data
    noraxondata = importdata(inputfile, '\t', 1);

    % set trigger frame as time zero
    noraxondata.data(:,1) = noraxondata.data(:,1) - noraxondata.data(1,1);
    
    % filter angle
    [B, A] = butter(4, freq_cutoff, 'low');
    norm_angle_filtered = filtfilt(B, A, noraxondata.data(:,column_norm_angle));
    
    % remove 180 degrees offset of input adapter
    angle_offset_180 = norm_volt_per_degree * 180;
    if mean(norm_angle_filtered) > 90 %VAR
        norm_angle_filtered = norm_angle_filtered - angle_offset_180;
        norm_angle_raw = noraxondata.data(:,column_norm_angle) - angle_offset_180;
    else % less than zero
        norm_angle_filtered = norm_angle_filtered + angle_offset_180;
        norm_angle_raw = noraxondata.data(:,column_norm_angle) + angle_offset_180;
    end

   
    
    
    
    %%% Identify movement phases and stop phases
    
    in_phase = 0; % time point is in an area without change in angle
    index = 1;
    phases(2,2) = zeros(); % start and stop oints of time periods without change in angle 
    start = 51; %VAR % point in angle array where we start searching for phases
    tolerance_phase_fluctuation = 1.2; %VAR % fluctuations to allow within a no-change phase, before switching to a movement phase
    
    % find first area where there is movement (discard first phase of zero movement)
    framestep = 30; %VAR
    for i = start:framestep:(length(norm_angle_filtered) - framestep)
        interval = norm_angle_filtered(i:i+framestep-1);
        start_interval = interval(1);
        stop_interval = interval(end);
        
        if abs(start_interval - stop_interval) > tolerance_phase_fluctuation
            start = i;
            break;
        end
    end
    
	% loop through all datapoints, "noframes" at a time
    framestep = 60; %VAR %MMM 2014 120 for passive, 60 for active?
    avgframes = 60/2; %VAR
    for i = start:framestep:(length(norm_angle_filtered) - framestep - avgframes)
        start_interval = mean(norm_angle_filtered(i-avgframes:i+avgframes));
        stop_interval = mean(norm_angle_filtered(i+framestep-avgframes:i+framestep+avgframes));
        
        if abs(start_interval - stop_interval) < tolerance_phase_fluctuation
            if in_phase == 0
                % set startpoint
                in_phase = 1;
                phases(index,1) = i+framestep + 20; %VAR going additional 10 frames into the phase, to avoid slight deviations around the beginning and end
            end
        else
            if in_phase == 1
                % set endpoint
                in_phase = 0;
                phases(index,2) = i-framestep - 20; %VAR
                
                if phases(index,2)-phases(index,1) <= 2*framestep
                   % detected phase is too small, is not real, overwrite with next phase
                   phases(index,:) = 0;
                else
                   index = index + 1;
                end
            end 
        end
    end
    % if the last phase does not have an end point, add END
    len = length(phases(:,2));
    if phases(len,2) == 0
        phases(len,2) = length(norm_angle_filtered)-10;
    end
    
    % Extract average volt values in detected stop phases
    voltvalues(1,1) = zeros();
    for i = 1:length(phases(:,1))
        % averaging RAW angle data
        voltvalues(i) = mean(norm_angle_raw(phases(i,1):phases(i,2)));
    end
    
    
    
    %%% plot phases
    
    if plot_check && plot_conversion
        plottitle = horzcat('Isokin ANGLE zone check, ', subject_id);
        figure('Name', plottitle)
        plot(norm_angle_filtered,'b')
        hold on
        for i = 1:length(phases(:,1))
            plot ([phases(i,1) phases(i,1)], [min(norm_angle_filtered) max(norm_angle_filtered)],'g')
            plot ([phases(i,2) phases(i,2)], [min(norm_angle_filtered) max(norm_angle_filtered)],'r')
        end
        xlabel('Time (frames)'),ylabel('Angle (mV)'),title(plottitle);
        % legend('CPM','Location','Southeast');
    end
    
    
    
    %%% Calculate conversion - section that is common for both methods
    
    % choose only highest and lowest average volt values (REAL stop phases)
    voltvalues_sort = sort(voltvalues);
    len = length(voltvalues_sort);
    tolerance_mV = 0.97; %VAR - how many % mV of the largest/smallest mV phase to accept as another valid phase (phase = max dorsiflexion / max plantarflexion)
    
    % first, third, fifth = dorsiflexed
    if dm_side == 'L' % work with smallest (first) values
        if voltvalues_sort(3) < voltvalues_sort(1) * tolerance_mV
            x2_array = voltvalues_sort(1:3);
        elseif voltvalues_sort(2) < voltvalues_sort(1) * tolerance_mV
            x2_array = voltvalues_sort(1:2);
        else
            x2_array = voltvalues_sort(1);
        end
    else % work with largest values
        if voltvalues_sort(len-2:len) > voltvalues_sort(len) * tolerance_mV
            x2_array = voltvalues_sort(len-2:len);
        elseif voltvalues_sort(len-1:len) > voltvalues_sort(len) * tolerance_mV
            x2_array = voltvalues_sort(len-1:len);
        else
            x2_array = voltvalues_sort(len);
        end
    end
    % second and fourth pair = plantarflexed
    if dm_side == 'L' % work with largest values
        if voltvalues_sort(len-2) > voltvalues_sort(len) * tolerance_mV
            x1_array = voltvalues_sort(len-2:len);
        elseif voltvalues_sort(len-1) > voltvalues_sort(len) * tolerance_mV
            x1_array = voltvalues_sort(len-1:len);
        else
            x1_array = voltvalues_sort(len);
        end
    else % work with smallest values
        if voltvalues_sort(3) > voltvalues_sort(1) * tolerance_mV
            x1_array = voltvalues_sort(1:3);
        elseif voltvalues_sort(2) > voltvalues_sort(1) * tolerance_mV
            x1_array = voltvalues_sort(1:2);
        else
            x1_array = voltvalues_sort(1);
        end
        x1_array = voltvalues_sort(1:3); % MMM TODO delete this entry? why keep 3, given if/else above?
    end
    % TODO MMM shave off excessive zero values? if 8 values are detected, delete those outside x% of mean value...? 
    x1 = mean(x1_array);
    x2 = mean(x2_array);
    
    
    
%     %%% Calculate conversion - old method, y = ax + b
% 
%     % Set corresponding volt values -10 to 30 degrees
%     %   trial contains: incomplete move -> DF, 1st stop = DF, 2nd stop = PF, and so on
%     %   L: high value = -10, low value = +30
%     %   R: high value = +30, low value = -10
%     
%     % Linear regression to determine A and B constants, for later use
%     X = [x1 x2];
%     Y = [y1 y2];
%     coefs = polyfit(X,Y,1);
%     convert_ind_angle_a = coefs(1);
%     convert_ind_angle_b = coefs(2);
%     
%     % Plot figure to confirm conversion
%     if plot_check && plot_conversion
%         plottitle = horzcat('Isokin ANGLE conversion, ', subject_id, ' x1/x2 = ', num2str(length(x1_array)), '/', num2str(length(x2_array)));
%         figure('Name', plottitle)
%         if dm_side == 'L'
%             plot(X,Y,'*',x2-50:10:x1+50,polyval(coefs,x2-50:10:x1+50),'-')
%         else
%             plot(X,Y,'*',x1-50:10:x2+50,polyval(coefs,x1-50:10:x2+50),'-')
%         end
%         % plot vertical lines for min/max values used for calculations
%         hold on
%         for i = 1:length(x1_array)
%             plot ([x1_array(i) x1_array(i)], [y2*1.2 y1*1.2],'g')
%         end
%         for i = 1:length(x2_array)
%             plot ([x2_array(i) x2_array(i)], [y2*1.2 y1*1.2],'g')
%         end
%         xlabel('x (volt)'),ylabel('y (degrees)'),title(plottitle);
%         % legend('isokinetic 45','Location','Southeast');
%     end
%     
%     % find offset value in volt
%     convert_ind_angle_b_volt = convert_ind_angle_b/convert_ind_angle_a;
%     
%     % output angle conversion numbers to screen, as text
%     report = sprintf(horzcat('INDI Angle conversion factors: a = ', num2str(convert_ind_angle_a), ', b = ', num2str(convert_ind_angle_b), '. Offset in millivolt = ', num2str(convert_ind_angle_b_volt), ' mV.' ));
%     disp(report)
    
    
    
    
    
    %%% Calculate conversion - new method, using Norm standard
    
    % plantarflexed
    % x1 - from common section
    if dm_side == 'L'
        y1 = 30;
    else
        y1 = -30;
    end
    z1 = y1 * norm_volt_per_degree;
    q1 = x1 - z1;
    
    % dorsiflexed
    % x2 - from common section
    if dm_side == 'L'
        y2 = -10;
    else
        y2 = 10;
    end
    z2 = y2 * norm_volt_per_degree;
    q2 = x2 - z2;
    
    if dm_side == 'L'
        convert_ind_angle_a = 1/norm_volt_per_degree;
        convert_ind_angle_b = -((q1 + q2)/2) / norm_volt_per_degree;
    else
        convert_ind_angle_a = -1/norm_volt_per_degree;
        convert_ind_angle_b = ((q1 + q2)/2) / norm_volt_per_degree;
    end
    
    % find offset value in volt
    convert_ind_angle_b_volt = convert_ind_angle_b/convert_ind_angle_a;
    
    % output angle conversion numbers to screen, as text
    cprintf('magenta', horzcat('NORM Angle conversion factors: a = ', num2str(convert_ind_angle_a), ', b = ', num2str(convert_ind_angle_b), '. Offset in millivolt = ', num2str(convert_ind_angle_b_volt), ' mV.\n' ));
    
    
    
    
    
    
end

