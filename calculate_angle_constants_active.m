%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_angle_constants_active
% Marie Moltubakk 1.12.2014 / 2.11.2016
% 
% Read raw torque data from Norm, filter
% Produce individual conversion factors for ANGLE from ISOMETRIC trials
%
% Note 1: Some variables are defined directly inside this function, see %VAR
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [convert_ind_angle_a, convert_ind_angle_b] = calculate_angle_constants_active(freq_cutoff, inputfile1, inputfile2, side)
    global mute
    global plot_conversion plot_check subject_id 
    global column_norm_angle  % column_norm_velocity
    global norm_volt_per_degree % norm_volt_per_nm norm_volt_per_velocity norm_mv2nm_a norm_mv2nm_b
    
        
    
    %% FILE 1   --- 10 degrees (plantar flexed)
    
    % import noraxon data
    noraxondata = importdata(inputfile1, '\t', 1);

    % set trigger frame as time zero
    noraxondata.data(:,1) = noraxondata.data(:,1) - noraxondata.data(1,1);
    
    % filter angle
    [B, A] = butter(4, freq_cutoff, 'low');
    norm_angle_filtered = filtfilt(B, A, noraxondata.data(:,column_norm_angle));
    
    % remove 180 degrees offset of input adapter
    angle_offset_180 = norm_volt_per_degree * 180;
    if mean(norm_angle_filtered) > 90 %VAR
        norm_angle_filtered = norm_angle_filtered - angle_offset_180;
        %norm_angle_raw = noraxondata.data(:,column_norm_angle) - angle_offset_180;
    else % less than zero
        norm_angle_filtered = norm_angle_filtered + angle_offset_180;
        %norm_angle_raw = noraxondata.data(:,column_norm_angle) + angle_offset_180;
    end

    
    
    %%% Identify movement phases and stop phases
    % if the trial does not contain a movement phase (recording starts @ isometric contraction), this method will correctly use the beginning of the recorded data
    if side == 'L'
        norm_angle_max = min(norm_angle_filtered) + 20; %VAR - very raw method
        % find area of testing angle range
        loc_frame = find(norm_angle_filtered<=norm_angle_max,1,'first');
    else % R
        norm_angle_max = max(norm_angle_filtered) - 20;
        % find area of testing angle range
        loc_frame = find(norm_angle_filtered>=norm_angle_max,1,'first');
    end
    
    % number of frames to step away from min velocity
    framestep = 300;
    % number of frames to USE - very raw method....
    frameuse = 7000;
    if length(norm_angle_filtered) < loc_frame+framestep+frameuse
        frameuse = length(norm_angle_filtered) - loc_frame - framestep;
    end
    
    % average angle data
    x1 = mean(norm_angle_filtered(loc_frame+framestep:loc_frame+framestep+frameuse));
    
    
    
    %%% plot phases
    
    if plot_check && plot_conversion
        plottitle = horzcat('Isometric ANGLE zone check PF10, ', subject_id);
        figure('Name', plottitle)
        plot(norm_angle_filtered,'b')
        hold on
        plot ([loc_frame loc_frame], [min(norm_angle_filtered) max(norm_angle_filtered)],'g')
        plot ([loc_frame+framestep loc_frame+framestep], [min(norm_angle_filtered) max(norm_angle_filtered)],'g')
        plot ([loc_frame+framestep+frameuse loc_frame+framestep+frameuse], [min(norm_angle_filtered) max(norm_angle_filtered)],'r')
        plot ([loc_frame+framestep loc_frame+framestep+frameuse], [x1 x1],'m')
        xlabel('Time (frames)'),ylabel('Angle (mV)'),title(plottitle);
        legend('angle', 'ID angle', 'start', 'stop', 'avg angle', 'Location','Southeast');
    end
    
    
    
    
    %% FILE 2   --- 0 degrees
    
    % import noraxon data
    noraxondata = importdata(inputfile2, '\t', 1);

    % set trigger frame as time zero
    noraxondata.data(:,1) = noraxondata.data(:,1) - noraxondata.data(1,1);
    
    % filter angle
    [B, A] = butter(4, freq_cutoff, 'low');
    norm_angle_filtered = filtfilt(B, A, noraxondata.data(:,column_norm_angle));
    
    % remove 180 degrees offset of input adapter
    angle_offset_180 = norm_volt_per_degree * 180;
    if mean(norm_angle_filtered) > 90 %VAR
        norm_angle_filtered = norm_angle_filtered - angle_offset_180;
        %norm_angle_raw = noraxondata.data(:,column_norm_angle) - angle_offset_180;
    else % less than zero
        norm_angle_filtered = norm_angle_filtered + angle_offset_180;
        %norm_angle_raw = noraxondata.data(:,column_norm_angle) + angle_offset_180;
    end
    
    
    
    %%% Identify movement phases and stop phases
    if side == 'L'
        norm_angle_max = min(norm_angle_filtered) + 20; %VAR - very raw method
        % find area of testing angle range
        loc_frame = find(norm_angle_filtered<=norm_angle_max,1,'first');
    else % R
        norm_angle_max = max(norm_angle_filtered) - 20;
        % find area of testing angle range
        loc_frame = find(norm_angle_filtered>=norm_angle_max,1,'first');
    end
    

    % number of frames to step away from min velocity
    framestep = 300;
    % number of frames to USE % very raw method....
    frameuse = 7000;
    if length(norm_angle_filtered) < loc_frame+framestep+frameuse
        frameuse = length(norm_angle_filtered) - loc_frame - framestep;
    end
    
    % average angle data
    x2 = mean(norm_angle_filtered(loc_frame+framestep:loc_frame+framestep+frameuse));
    
    
    
    %%% plot phases
    
    if plot_check && plot_conversion
        plottitle = horzcat('Isometric ANGLE zone check DF00, ', subject_id);
        figure('Name', plottitle)
        plot(norm_angle_filtered,'b')
        hold on
        plot ([loc_frame loc_frame], [min(norm_angle_filtered) max(norm_angle_filtered)],'g')
        plot ([loc_frame+framestep loc_frame+framestep], [min(norm_angle_filtered) max(norm_angle_filtered)],'g')
        plot ([loc_frame+framestep+frameuse loc_frame+framestep+frameuse], [min(norm_angle_filtered) max(norm_angle_filtered)],'r')
        plot ([loc_frame+framestep loc_frame+framestep+frameuse], [x2 x2],'m')
        xlabel('Time (frames)'),ylabel('Angle (mV)'),title(plottitle);
        legend('angle', 'ID angle', 'start', 'stop', 'avg angle', 'Location','Southeast');
    end
    
    
    
    %% CALCULATE RANGE/SCALE
    % input files =
    % 1 - plantar 10 degrees
    % 2 - dorsi 0 degrees
    
    %%% Calculate conversion - new method, using Norm standard
    
    % plantarflexed
    % x1 - from common section
    if side == 'L'
        y1 = 10;
    else
        y1 = -10;
    end
    z1 = y1 * norm_volt_per_degree;
    q1 = x1 - z1;
    
    % dorsiflexed
    % x2 - from common section
    if side == 'L'
        y2 = 0;
    else
        y2 = 0;
    end
    z2 = y2 * norm_volt_per_degree;
    q2 = x2 - z2;
    
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
    if mute == 0
        if plot_conversion
            cprintf('cyan', horzcat('NORM Angle conversion factors: a = ', num2str(convert_ind_angle_a), ', b = ', num2str(convert_ind_angle_b), '. (((Offset in millivolt = ', num2str(convert_ind_angle_b_volt), ' mV.)))\n' ));
        end
    end
    
    
    
    
end
