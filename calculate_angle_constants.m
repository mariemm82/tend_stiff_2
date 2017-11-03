%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_angle_constants
% Marie Moltubakk 9.6.2013
% Read raw torque data from Norm, filter
% Produce individual conversion factors for angle from CPM trial
% Note 1: Some variables are defined directly on this method, see %VAR
%%%%%%%%%%%%%%%%%%%%%%%%%%




function [convert_ind_angle_a, convert_ind_angle_b] = calculate_angle_constants(freq_cutoff, inputfile, side)
    global plot_conversion subject_id 
    global column_norm_angle
    global norm_volt_per_degree % norm_volt_per_nm norm_volt_per_velocity norm_mv2nm_a norm_mv2nm_b
        
    
    
    % import noraxon data
    noraxondata = importdata(inputfile, '\t', 1);

    % set trigger frame as time zero
    noraxondata.data(:,1) = noraxondata.data(:,1) - noraxondata.data(1,1);
    
    % filter angle
    [B, A] = butter(4, freq_cutoff, 'low');
    norm_angle_filtered = filtfilt(B, A, noraxondata.data(:,column_norm_angle));
    % norm_angle_raw = noraxondata.data(:,column_norm_angle);
    % remove 180 degrees offset of input adapter
    angle_offset_180 = norm_volt_per_degree * 180;
    if mean(norm_angle_filtered) > 0
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
    tolerance_phase_fluctuation = 0.8; %VAR % fluctuations to allow within a no-change phase, before switching to a movement phase

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
    framestep = 120; %VAR
    avgframes = 100/2; %VAR
    for i = start:framestep:(length(norm_angle_filtered) - framestep - avgframes)
        start_interval = mean(norm_angle_filtered(i-avgframes:i+avgframes));
        stop_interval = mean(norm_angle_filtered(i+framestep-avgframes:i+framestep+avgframes));
        
        if abs(start_interval - stop_interval) < tolerance_phase_fluctuation
            if in_phase == 0
                % set startpoint
                in_phase = 1;
                phases(index,1) = i+framestep + 50; %VAR going additional 10 frames into the phase, to avoid slight deviations around the beginning and end
            end
        else
            if in_phase == 1
                % set endpoint
                in_phase = 0;
                phases(index,2) = i-framestep - 50; %VAR
                
                if phases(index,2)-phases(index,1) <= 2*framestep
                   % detected phase is too small, is not real, overwrite with next phase
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
    
    
    
    
    % Extract average volt values in stop phases
    voltvalues(1,1) = zeros();
    for i = 1:length(phases(:,1))
        % averaging RAW angle data
        voltvalues(i) = mean(norm_angle_raw(phases(i,1):phases(i,2)));
    end
    
    
    
    %%% plot phases
    
    if plot_conversion
        plottitle = horzcat('CPM ANGLE zone check, ', subject_id);
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
    
    
    
%     %%% Calculate conversion - old method, y = ax + b
%     
%     % Set corresponding volt values -5 to 10 degrees
%     %   L leg always high volt first, R leg low volt first
%     %   trial contains: incomplete move -> DF, 1st stop = PF, 2nd stop = DF
%     %   L: high value = -5, low value = +10
%     %   R: high value = +10, low value = -5
%     
%     % first and third pair = plantarflexed
%     x1 = (voltvalues(1)+voltvalues(3))/2;
%     y1 = 10;
%     % second pair = dorsiflexed
%     x2 = voltvalues(2);
%     y2 = -5;
%     
%     % Linear regression to determine A and B constants, for later use
%     X = [x1 x2];
%     Y = [y1 y2];
%     coefs = polyfit(X,Y,1);
%     convert_ind_angle_a = coefs(1);
%     convert_ind_angle_b = coefs(2);
%     
%     % find offset value in volt
%     convert_ind_angle_b_volt = convert_ind_angle_b/convert_ind_angle_a;
%     
%     % output angle conversion numbers to screen, as text
%     report = sprintf(horzcat('INDI Angle conversion factors: a = ', num2str(convert_ind_angle_a), ', b = ', num2str(convert_ind_angle_b), '. Offset in millivolt = ', num2str(convert_ind_angle_b_volt), ' mV.' ));
%     disp(report)
    
    
    
    %%% Calculate conversion - new method, using Norm standard for volt/degree, computing separately only B factore
    
    if isnan(voltvalues(3))
        % special cases where only 2 real values exist:
        % second value = plantarflexed
        x1 = voltvalues(2);
        if side == 'L'
            y1 = 10;
        else
            y1 = -10;
        end
        z1 = y1 * norm_volt_per_degree;
        q1 = x1 - z1;

        % first value = dorsiflexed
        x2 = voltvalues(1);
        if side == 'L'
            y2 = -5;
        else
            y2 = 5;
        end
        z2 = y2 * norm_volt_per_degree;
        q2 = x2 - z2;
        
    else
        % this is the NORMAL procedure, where normally 3-4 real values exist
        % first and third pair = plantarflexed
        x1 = (voltvalues(1)+voltvalues(3))/2;
        if side == 'L'
            y1 = 10;
        else
            y1 = -10;
        end
        z1 = y1 * norm_volt_per_degree;
        q1 = x1 - z1;

        % second pair = dorsiflexed
        x2 = voltvalues(2);
        if side == 'L'
            y2 = -5;
        else
            y2 = 5;
        end
        z2 = y2 * norm_volt_per_degree;
        q2 = x2 - z2;
    end
    if side == 'L'
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

