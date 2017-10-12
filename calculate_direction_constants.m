%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_direction_constants
% Marie Moltubakk 28.11.2014
% Read raw data from Norm, filter
% Produce individual conversion factors for DIRECTION
%
% Note 1: Some variables are defined directly on this method, see %VAR
%%%%%%%%%%%%%%%%%%%%%%%%%%




function [convert_ind_direction_b_volt] = calculate_direction_constants(inputfile)
    global plot_conversion plot_check subject_id
  %  global dm_side
    global column_norm_direction
   % global norm_volt_per_degree norm_volt_per_nm norm_volt_per_velocity norm_mv2nm_a norm_mv2nm_b
    
    % import noraxon data
    noraxondata = importdata(inputfile, '\t', 1);

    % set trigger frame as time zero
    noraxondata.data(:,1) = noraxondata.data(:,1) - noraxondata.data(1,1);
    
     % filter direction - NOT NEEDED
%     [B, A] = butter(4, cutoff, 'low');
%     norm_direction_filtered = filtfilt(B, A, noraxondata.data(:,column_norm_direction));
    norm_direction = noraxondata.data(:,column_norm_direction);
    
    
    
    
    %%% Identify movement phases and stop phases
    % this is really more complex than needed for direction, but is simply copied from velocity
    
    in_phase = 0; % time point is in an area without change in angle
    index = 1;
    phases(2,2) = zeros(); % start and stop points of time periods without change in angle 
    start = 60; %VAR % point in angle array where we start searching for phases
    tolerance = 1.0; %VAR % fluctuations to allow within a no-change phase, before switching to a movement phase
    
    % loop through all datapoints, "noframes" at a time
    % this is copied from ANGLE calculations. terminology for movement phase vs stop phase is not changed. When angle is changing, velocity is not...
    framestep = 120; %VAR %MMM 2014 changed from 150 to 120 for F28 pre R sol
    avgframes = 100/2; %VAR
    for i = start:framestep:(length(norm_direction) - framestep - avgframes)
        start_interval = mean(norm_direction(i-avgframes:i+avgframes));
        stop_interval = mean(norm_direction(i+framestep-avgframes:i+framestep+avgframes));
        
        if abs(start_interval - stop_interval) < tolerance
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
    
    
    
    %%% plot phases
    if plot_check && plot_conversion
        plottitle = horzcat('Norm DIRECTION zone check for ', subject_id);
        figure('Name', plottitle)
        plot(norm_direction,'b')
        hold on
        for i = 1:length(phases(:,1))
            plot ([phases(i,1) phases(i,1)], [min(norm_direction) max(norm_direction)],'g')
            plot ([phases(i,2) phases(i,2)], [min(norm_direction) max(norm_direction)],'r')
        end
        xlabel('Time (frames)'),ylabel('Angle (mV)'),title(plottitle);
    end
    
    
    
    %%% Extract average volt values in stop phases
    
    voltvalues(1,1) = zeros();
    % delete last entry if non-valid (end of movement not existing)
    if phases(length(phases),2) == 0
        phases = phases(1:length(phases)-1,:);
    elseif phases(length(phases),1) >= phases(length(phases),2)
        phases = phases(1:length(phases)-1,:);
    end
    for i = 1:length(phases(:,1))
        % averaging RAW angle data
        voltvalues(i) = mean(noraxondata.data(phases(i,1):phases(i,2),column_norm_direction));
    end

    
    
    %%% Keep only the low volt values (near zero)
    
    voltvalues_low(1,1) = zeros();
    ii = 1;
    for i = 1:length(voltvalues)
        if voltvalues(i) < 300 % VAR
            voltvalues_low(ii) = voltvalues(i);
            ii = ii+1;
        end
    end
    convert_ind_direction_b_volt = -mean(voltvalues_low);
    

    
    %%% output angle conversion numbers to screen, as text
    if plot_conversion
        cprintf('magenta', horzcat('Direction conversion factors: b offset in millivolt = ', num2str(convert_ind_direction_b_volt), ' mV.\n'));
    end
end

