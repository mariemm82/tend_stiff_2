%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_velocity_constants
% Marie Moltubakk 26.11.2014
% Read raw data from Norm, filter
% Produce individual conversion factors for VELOCITY
%
% Note 1: Some variables are defined directly on this method, see %VAR
%%%%%%%%%%%%%%%%%%%%%%%%%%




function [convert_ind_velocity_a, convert_ind_velocity_b] = calculate_velocity_constants(freq_cutoff, inputfile, side)
    global plot_conversion plot_check plot_achilles subject_id
    global column_norm_velocity
    global norm_volt_per_degree norm_volt_per_nm norm_volt_per_velocity norm_mv2nm_a norm_mv2nm_b
    
    % import noraxon data
    noraxondata = importdata(inputfile, '\t', 1);

    % set trigger frame as time zero
    noraxondata.data(:,1) = noraxondata.data(:,1) - noraxondata.data(1,1);
    
    % filter velocity
    [B, A] = butter(4, freq_cutoff, 'low');
    norm_velocity_filtered = filtfilt(B, A, noraxondata.data(:,column_norm_velocity));

    
    
    %%% Identify movement phases and stop phases
    
    in_phase = 0; % time point is in an area without change in angle
    index = 1;
    phases(2,2) = zeros(); % start and stop points of time periods without change in angle 
    start = 1000; %VAR % point in VELOCITY array where we start searching for phases / discard first part with manual movement into dorsiflexion
    tolerance_phase_fluctuation = 1.2; %VAR % fluctuations to allow within a no-change phase, before switching to a movement phase
    
    % loop through all datapoints, "noframes" at a time
    % this is copied from ANGLE calculations. terminology for movement phase vs stop phase is not changed. When angle is changing, velocity is not...
    framestep = 120; %VAR %MMM 2014 changed from 150 to 120 for F28 pre R sol
    avgframes = 100/2; %VAR
    for i = start:framestep:(length(norm_velocity_filtered) - framestep - avgframes)
        start_interval = mean(norm_velocity_filtered(i-avgframes:i+avgframes));
        stop_interval = mean(norm_velocity_filtered(i+framestep-avgframes:i+framestep+avgframes));
        
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
        phases(len,2) = length(norm_velocity_filtered)-10;
    end
    
    % Extract average volt values in stop phases
    voltvalues(1,1) = zeros();
    % delete last entry if non-valid (end of movement not existing)
    if phases(length(phases),2) == 0
        phases = phases(1:length(phases)-1,:);
    end
    for i = 1:length(phases(:,1))
        % averaging RAW volt data
        voltvalues(i) = mean(noraxondata.data(phases(i,1):phases(i,2),column_norm_velocity));
    end
    
    
    
    %%% plot phases
    
    if plot_check && plot_conversion && plot_achilles
        plottitle = horzcat('Passive VELOCITY zones, ', subject_id);
        figure('Name', plottitle)
        plot(norm_velocity_filtered,'b')
        hold on
        for i = 1:length(phases(:,1))
            plot ([phases(i,1) phases(i,1)], [min(norm_velocity_filtered) max(norm_velocity_filtered)],'g')
            plot ([phases(i,2) phases(i,2)], [min(norm_velocity_filtered) max(norm_velocity_filtered)],'r')
        end
        xlabel('Time (frames)'),ylabel('Velocity (mV)'),title(plottitle);
        % legend('CPM','Location','Southeast');
    end
    
    

    
    %%% Calculate conversion - section that is common for both methods
    
    % TODO MMM: would be better to identify phases by grouping according to volt levels, not assuming that phases are always detected in this order
    % movement into plantar flexion:
    x1 = (voltvalues(1)+voltvalues(5))/2;
    y1 = +5;
    % movement into dorsiflexion
    if length(voltvalues) >= 7
        x2 = (voltvalues(3)+voltvalues(7))/2;
    else
        x2 = voltvalues(3);
    end
    y2 = -5;
    % pauses
    if length(voltvalues) >= 6
        x3 = (voltvalues(2)+voltvalues(4)+voltvalues(6))/3;
    else
        x3 = (voltvalues(2)+voltvalues(4))/2;
    end
    y3 = 0;

    
    
%     %%% Calculate conversion - old method, y = ax + b
%     
%     % Set corresponding volt values -5 to +5 degrees/sec
%     %   EXTERNALLY: manual move up - passive down - passive up ...
%     %   LEFT trial contains: low volt - zero - high - zero - low - zero - high - zero to end
%     %   LEFT trial: low = movement down = plantar flexion = clockwise = positive velocity
%     %   RIGHT trial contains: high volt - zero - low - zero - high - zero - low - "mess" to end
%     %   RIGHT trial: high = movement down = plantar flexion = COUNTERclockwise = positive velocity
% 
%     % prepare array for plotting individual volt values
%     % both L and R start with plantar flexion = positive velocity
%     voltvalues_expected = [+5 0 -5 0 +5 0 -5];
%     
%     % Linear regression to determine A and B constants, for later use
%     X = [x1 x3 x2];
%     Y = [y1 y3 y2];
%     coefs = polyfit(X,Y,1);
%     convert_ind_velocity_a = coefs(1);
%     convert_ind_velocity_b = coefs(2);
%     
%     % Plot figure to confirm conversion
%     if plot_check && plot_conversion
%         plottitle = horzcat('VELOCITY conversion check for ', subject_id);
%         figure('Name', plottitle)
%         if side == 'L'
%             plot(X,Y,'*',x1*1.2:10:x2*0.8,polyval(coefs,x1*1.2:10:x2*0.8),'-')
%             hold on
%             plot(voltvalues(1:7),voltvalues_expected,'*g')
%         else
%             plot(X,Y,'*',x2*1.2:10:x1*0.8,polyval(coefs,x2*1.2:10:x1*0.8),'-')
%             hold on
%             plot(voltvalues(1:7),voltvalues_expected,'*g')
%         end
%         xlabel('x (volt)'),ylabel('y (degrees/sec)'),title(plottitle);
%         % legend('CPM','Location','Southeast');
%     end
%         
%     % find offset value in volt
%     convert_ind_velocity_b_volt = convert_ind_velocity_b/convert_ind_velocity_a;
%     
%     % output velocity conversion numbers to screen, as text
%     report = sprintf(horzcat('INDI Velocity conversion factors: a = ', num2str(convert_ind_velocity_a), ', b = ', num2str(convert_ind_velocity_b), '. Offset in millivolt = ', num2str(convert_ind_velocity_b_volt), ' mV.' ));
%     disp(report)
    
    
    
    
    %%% Calculate conversion - new method, using Norm standard for a constant, calculating individually only the b variable
    
    % plantarflexion
    % x1 - from common section
    if side == 'L'
        y1 = 5;
    else
        y1 = -5;
    end
    z1 = y1 * norm_volt_per_velocity;
    q1 = x1 - z1;
    
    % dorsiflexion
    % x2 - from common section
    if side == 'L'
        y2 = -5;
    else
        y2 = 5;
    end
    z2 = y2 * norm_volt_per_velocity;
    q2 = x2 - z2;
    
    if side == 'L'
        convert_ind_velocity_a = -1/norm_volt_per_velocity;
        convert_ind_velocity_b = ((q1 + q2)/2) / norm_volt_per_velocity;
    else
        convert_ind_velocity_a = 1/norm_volt_per_velocity;
        convert_ind_velocity_b = -((q1 + q2)/2) / norm_volt_per_velocity;
    end
    
    % find offset value in volt
    convert_ind_velocity_b_volt = convert_ind_velocity_b/convert_ind_velocity_a;
    
    % output velocity conversion numbers to screen, as text
    cprintf('magenta', horzcat('NORM Velocity conversion factors: a = ', num2str(convert_ind_velocity_a), ', b = ', num2str(convert_ind_velocity_b), '. Offset in millivolt = ', num2str(convert_ind_velocity_b_volt), ' mV.\n' ));
    
    
    
end

