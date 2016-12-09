%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract_isometric
% Marie Moltubakk 2.11.2016
% Read complete, prepared noraxon array
% Produce array with peak torque and angle
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [torque_max,angle_at_torque_max] = extract_isometric(noraxondata, side, trial_name)
    
    global column_norm_angle column_norm_torque % column_l_gm column_r_gm column_l_gl column_r_gl column_l_sol column_r_sol column_gonio  column_l_tibant column_r_tibant  column_norm_velocity column_norm_direction column_achilles column_EMG_start column_EMG_end 
    global filepath
    
    
    
    % Read Noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
    % Produce a new noraxon data array
    freq_resample = 200;
    noraxon_prepped = read_noraxon_active(strcat(filepath, noraxondata), freq_resample, side, trial_name);
    
    [torque_max,index] = max(noraxon_prepped(:,column_norm_torque));
    angle_at_torque_max = noraxon_prepped(index,column_norm_angle);
    
end