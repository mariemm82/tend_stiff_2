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
    
    % max data
    [torque_max,index] = max(array_torque);
    angle_at_torque_max = array_angle(index);
    velocity_at_torque_max = array_velocity (index);
   
    % work data
    work = 0; % MMM TODO
end