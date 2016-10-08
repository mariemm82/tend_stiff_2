%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_velocity_constants_active
% Marie Moltubakk 01.12.2014
% Read raw data from Norm, filter
% Produce individual conversion factors for VELOCITY, for ISOKINETIC trial
%
% Note 1: Some variables are defined directly on this method, see %VAR
%%%%%%%%%%%%%%%%%%%%%%%%%%




function [convert_ind_velocity_a, convert_ind_velocity_b] = calculate_velocity_constants_active(freq_cutoff, inputfile)
    global plot_conversion plot_check subject_id
    global dm_side
    global column_norm_velocity
    global norm_volt_per_degree norm_volt_per_nm norm_volt_per_velocity norm_mv2nm_a norm_mv2nm_b
    
    % import noraxon data
    noraxondata = importdata(inputfile, '\t', 1);

    % set trigger frame as time zero
    noraxondata.data(:,1) = noraxondata.data(:,1) - noraxondata.data(1,1);
    
    % filter velocity
    [B, A] = butter(4, freq_cutoff, 'low');
    norm_velocity_filtered = filtfilt(B, A, noraxondata.data(:,column_norm_velocity));
    norm_velocity_raw = noraxondata.data(:,column_norm_velocity);
    
    
    
    %%% Identify movement phases and stop phases
    
    in_phase = 0; % time point is in an area without change in angle
    index = 1;
    phases(2,2) = zeros(); % start and stop points of time periods without change in angle 
    start = 51; %VAR % point in VELOCITY array where we start searching for phases / discard first part with manual movement into dorsiflexion
    tolerance_phase_fluctuation = 1.2; %VAR % fluctuations to allow within a no-change phase, before switching to a movement phase
    
    % loop through all datapoints, "noframes" at a time
    % this is copied from ANGLE calculations. terminology for movement phase vs stop phase is not changed. When angle is changing, velocity is not...
    framestep = 60; %VAR
    avgframes = 60/2; %VAR
    for i = start:framestep:(length(norm_velocity_filtered) - framestep - avgframes)
        start_interval = mean(norm_velocity_filtered(i-avgframes:i+avgframes));
        stop_interval = mean(norm_velocity_filtered(i+framestep-avgframes:i+framestep+avgframes));
        
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
        phases(len,2) = length(norm_velocity_filtered)-10;
    end
    
    % Extract average volt values in stop phases
    voltvalues(1,1) = zeros();
    % delete last entry if non-valid (end of movement not existing)
    if phases(length(phases),2) == 0
        phases = phases(1:length(phases)-1,:);
    end
    for i = 1:length(phases(:,1))
        % averaging RAW data
        voltvalues(i) = mean(norm_velocity_raw(phases(i,1):phases(i,2)));
    end
    
    
    
    %%% plot phases
    
    if plot_check && plot_conversion
        plottitle = horzcat('Active VELOCITY zones, ', subject_id);
        figure('Name', plottitle)
        plot(norm_velocity_filtered,'b')
        hold on
        for i = 1:length(phases(:,1))
            plot ([phases(i,1) phases(i,1)], [min(norm_velocity_filtered) max(norm_velocity_filtered)],'g')
            plot ([phases(i,2) phases(i,2)], [min(norm_velocity_filtered) max(norm_velocity_filtered)],'r')
        end
        text(0.1,100,'Zones of max dorsiflexion velocity are NOT used for calculations.')
        xlabel('Time (frames)'),ylabel('Velocity (mV)'),title(plottitle);
        % legend('CPM','Location','Southeast');
    end
    

    
    
    %%% Calculate conversion - section that is common for both methods
    
    % choose only highest and lowest average volt values (REAL constant velocity phases)
    voltvalues_sort = sort(voltvalues);
    len = length(voltvalues_sort);
    tolerance_mV = 0.97; %VAR - how many % mV of the largest/smallest mV phase to accept as another valid phase (phase = max dorsiflexion / max plantarflexion)
%    tolerance_to_zero = 200; %VAR - mV from extreme value to zero
    
    % movement into plantar flexion
%    y1 = +45;
    if dm_side == 'L' % PF = low (first) values, zero = high values
        if voltvalues_sort(3) < voltvalues_sort(1) * tolerance_mV
            x1_array = voltvalues_sort(1:3);
            x1_no = 3;
        elseif voltvalues_sort(2) < voltvalues_sort(1) * tolerance_mV
            x1_array = voltvalues_sort(1:2);
            x1_no = 2;
        else
            x1_array = voltvalues_sort(1);
            x1_no = 1;
        end
    else % right leg: PF = high (last) values, zero = low values
        if voltvalues_sort(len-2:len) > voltvalues_sort(len) * tolerance_mV
            x1_array = voltvalues_sort(len-2:len);
            x1_no = 3;
        elseif voltvalues_sort(len-1:len) > voltvalues_sort(len) * tolerance_mV
            x1_array = voltvalues_sort(len-1:len);
            x1_no = 2;
        else
            x1_array = voltvalues_sort(len);
            x1_no = 1;
        end
    end
    
    % pauses at zero velocity
%    y2 = 0;
    if dm_side == 'L' % PF = low (first) values, zero = high values
        x2_no = 1;
        i = x1_no+1;
        x2_array = [0 0 0];
        while i <= length(voltvalues_sort) && voltvalues_sort(i) < 80 %VAR
            if voltvalues_sort(i) > -120 && voltvalues_sort(i) < 80 %VAR x2
                x2_array(x2_no) = voltvalues_sort(i);
                x2_no = x2_no + 1;
            end
            i = i + 1;
        end
    else % right leg: PF = high (last) values, zero = low values
        x2_no = 1;
        i = len-x1_no;
        x2_array = [0 0 0];
        while i > 0 && voltvalues_sort(i) > -120  %VAR
            if voltvalues_sort(i) > -120 && voltvalues_sort(i) < 80 %VAR x2
                x2_array(x2_no) = voltvalues_sort(i);
                x2_no = x2_no + 1;
            end
            i = i - 1;
        end
    end
    % TODO MMM shave off excessive zero values? if 8 values are detected, delete those outside x% of mean value...? 
    x1 = mean(x1_array);
    x2 = mean(x2_array);
    
    
    
%     %%% Calculate conversion - old method, y = ax + b
%     
%     % Set corresponding volt values -5 to +5 degrees/sec
%     %   EXTERNALLY: manual move up - passive down - passive up ...
%     %   LEFT trial contains: low volt - zero - high - zero - low - zero - high - zero to end
%     %   LEFT trial: low = movement down = plantar flexion = clockwise = positive velocity
%     %   RIGHT trial contains: high volt - zero - low - zero - high - zero - low - "mess" to end
%     %   RIGHT trial: high = movement down = plantar flexion = COUNTERclockwise = positive velocity
% 
%     % Linear regression to determine A and B constants, for later use
%     X = [x1 x2];
%     Y = [y1 y2];
%     coefs = polyfit(X,Y,1);
%     convert_ind_velocity_a = coefs(1);
%     convert_ind_velocity_b = coefs(2);
%      
%     % Plot figure to confirm conversion
%     if plot_check && plot_conversion
%         plottitle = horzcat('Passive VELOCITY conversion, ', subject_id, ' x1/x2 = ', num2str(length(x1_array)), '/', num2str(length(x2_array)));
%         figure('Name', plottitle)
%         if dm_side == 'L'
%             plot(X,Y,'*',x1-50:10:x2+50,polyval(coefs,x1-50:10:x2+50),'-')
%         else
%             plot(X,Y,'*',x2-50:10:x1+50,polyval(coefs,x2-50:10:x1+50),'-')
%         end
%         % plot vertical lines for min/max values used for calculations
%         hold on
%         for i = 1:length(x1_array)
%             plot ([x1_array(i) x1_array(i)], [y2*1.2 y1*1.2],'g')
%         end
%         for i = 1:length(x2_array)
%             plot ([x2_array(i) x2_array(i)], [y2*1.2 y1*1.2],'g')
%         end
%         xlabel('x (volt)'),ylabel('y (degrees/sec)'),title(plottitle);
%         % legend('isokinetic 45','Location','Southeast');
%     end
%         
%     % find offset value in volt
%     convert_ind_velocity_b_volt = convert_ind_velocity_b/convert_ind_velocity_a;
%     
%     % output angle conversion numbers to screen, as text
%     report = sprintf(horzcat('INDI Velocity conversion factors: a = ', num2str(convert_ind_velocity_a), ', b = ', num2str(convert_ind_velocity_b), '. Offset in millivolt = ', num2str(convert_ind_velocity_b_volt), ' mV.' ));
%     disp(report)
    
    
    
    %%% Calculate conversion - new method, using Norm standard
    
    % movement into plantarflexion
    % x1 - from common section
    if dm_side == 'L'
        y1 = -45;
    else
        y1 = 45;
    end
    z1 = y1 * norm_volt_per_velocity;
    q1 = x1 - z1;
    
    % stop phases at zero velocity
    % x2 - from common section
    if dm_side == 'L'
        y2 = 0;
    else
        y2 = 0;
    end
    z2 = y2 * norm_volt_per_velocity;
    q2 = x2 - z2;
    
    if dm_side == 'L'
        convert_ind_velocity_a = -1/norm_volt_per_velocity;
        convert_ind_velocity_b = ((q1 + q2)/2) / norm_volt_per_velocity;
    else % R side
        convert_ind_velocity_a = 1/norm_volt_per_velocity;
        convert_ind_velocity_b = -((q1 + q2)/2) / norm_volt_per_velocity;
    end
    
    % find offset value in volt
    convert_ind_velocity_b_volt = convert_ind_velocity_b/convert_ind_velocity_a;
    
    % output velocity conversion numbers to screen, as text
    cprintf('magenta', horzcat('NORM Velocity conversion factors: a = ', num2str(convert_ind_velocity_a), ', b = ', num2str(convert_ind_velocity_b), '. Offset in millivolt = ', num2str(convert_ind_velocity_b_volt), ' mV.\n' ));
    
    
    
    
end

