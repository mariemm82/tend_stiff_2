%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract_isokinetic
% Marie Moltubakk 2.11.2016
% Read complete, prepared noraxon array
% Produce array with torque and angle
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [torque_max, angle_at_torque_max, velocity_at_torque_max, work, array_output] = extract_isokinetic(noraxondata, side, trial_name)
    
    global column_norm_angle column_norm_torque column_norm_velocity % column_l_gm column_r_gm column_l_gl column_r_gl column_l_sol column_r_sol column_gonio  column_l_tibant column_r_tibant  column_norm_velocity column_norm_direction column_achilles column_EMG_start column_EMG_end 
    global filepath
    
    
    
    % Read Noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
    % Produce a new noraxon data array
    freq_resample = 200;
    noraxon_prepped = read_noraxon_active(strcat(filepath, noraxondata), freq_resample, side, trial_name);
    
    % array data
    array_torque = noraxon_prepped(:,column_norm_torque);
    array_angle = noraxon_prepped(:,column_norm_angle);
    array_velocity = noraxon_prepped(:,column_norm_velocity);
    array_output = [array_torque array_angle];
    
    
    % separate trials, additional filter per trial (electrical noise removal)
    
    % identify movement phases
    [~,loc_peaks] = findpeaks(array_angle,'MinPeakDistance',100);
    
    % TODO MMM: Remove multiple peaks at same level, remove "peaks" at mid
    % value (near angle = zero), etc
    
    if length(loc_peaks) < 7
        % add final part of dorsiflexion phase
        loc_peaks(7) = length(array_angle);
    end
    
    % TMP MMM
    figure,plot(array_angle)
    hold on
    findpeaks(array_angle,'MinPeakDistance',100)
    
    % prepare filter
    [B, A] = butter(8, 0.05, 'low'); %VAR
    if strcmp(trial_name,'isokin DF 30')
        % collect and filter dorsiflexion trials
        df1 = [filtfilt(B, A, array_torque(loc_peaks(1):loc_peaks(2))) array_angle(loc_peaks(1):loc_peaks(2))];
        df2 = [filtfilt(B, A, array_torque(loc_peaks(3):loc_peaks(4))) array_angle(loc_peaks(3):loc_peaks(4))];
        df3 = [filtfilt(B, A, array_torque(loc_peaks(5):loc_peaks(6))) array_angle(loc_peaks(5):loc_peaks(6))];
        
        plottitle = horzcat('Dorsiflexion trials ', trial_name);
        figure('Name',plottitle)
        plot(df1(:,2),df1(:,1))
        hold on
        plot(df2(:,2),df2(:,1))
        plot(df3(:,2),df3(:,1))
        ax = gca;
        set(ax, 'xdir','reverse')
        xlabel('Ankle angle (deg)')
        ylabel('Isokinetic torque (Nm)')
        title(plottitle)
        %legend('øeg','location','SouthWest')
    else
        % collect and filter plantar flexion trials
        df1 = [filtfilt(B, A, array_torque(loc_peaks(2):loc_peaks(3))) array_angle(loc_peaks(2):loc_peaks(3))];
        df2 = [filtfilt(B, A, array_torque(loc_peaks(4):loc_peaks(5))) array_angle(loc_peaks(4):loc_peaks(5))];
        df3 = [filtfilt(B, A, array_torque(loc_peaks(6):loc_peaks(7))) array_angle(loc_peaks(6):loc_peaks(7))];
        
        plottitle = horzcat('Plantar flexion trials ', trial_name);
        figure('Name',plottitle)
        plot(df1(:,2),df1(:,1))
        hold on
        plot(df2(:,2),df2(:,1))
        plot(df3(:,2),df3(:,1))
        ax = gca;
        xlabel('Ankle angle (deg)')
        ylabel('Isokinetic torque (Nm)')
        title(plottitle)
        %legend('øeg','location','SouthWest')
    end
    



    
    
    % max data
    [torque_max,index] = max(array_torque);
    angle_at_torque_max = array_angle(index);
    velocity_at_torque_max = array_velocity (index);
   
    % work data
    work = 0; % MMM TODO
    
    % MMM TODO - deal with dorsiflexion trial
    
    
end