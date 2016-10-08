%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script for determining conversion factors Norm/Noraxon, passive&active trials
% 
% passive: using 2 trials: slow CPM +5 to -10 degrees, calc scan and sol scan
% 
% Marie Moltubakk 28.1.2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [convert_norm_angle_a, convert_norm_angle_b, convert_norm_torque_a, convert_norm_torque_b, convert_norm_velocity_a, convert_norm_velocity_b, convert_norm_direction_b] = calculate_norm_constants_complete(subjectno, side, timepoint, line, act_or_pas)
global dm_CPM_calc_NX dm_CPM_sol_NX % dm_MVC_PF dm_MVC_DF  dm_CPM_calc_US dm_CPM_calc_US_frame  dm_CPM_sol_US dm_CPM_sol_US_frame
global filepath

% input:
% subjectno
% side
% timepoint
% act_or_pas






%%% Set constants % PROJECTSPECIFIC
% below variables are copied from other files. those parts that are not needed here are commented out, rather than deleted. 

% sampling frequencies
%global us_zerodispframes noraxonfreq freq_default

% default conversion factors from Norm manuals and Achilles calibration
global norm_volt_per_nm_a norm_volt_per_nm_b % convert_achilles norm_volt_per_degree norm_volt_per_velocity

% conversion factors, calculated per subject (values will be extracted by the script below and inserted to these arrays)
% global convert_norm convert_norm_ind_passive convert_norm_ind_active

% cutoff frequencies for filters
%global emg_bandpass emg_rms_ms mvc_window_ms
global angle_cutoff velocity_cutoff angle_cutoff_active velocity_cutoff_active % torque_cutoff_bandstop torque_cutoff_active

% column placement in Noraxon data
%global column_EMG_start column_EMG_end column_l_gm column_r_gm column_l_gl column_r_gl column_l_sol column_r_sol column_l_tibant column_r_tibant column_gonio column_norm_angle column_norm_torque column_norm_velocity column_norm_direction column_achilles







if strcmpi(act_or_pas,'passive') == 1
    dm_trial = 'PASSIVE';
    % using two CPM trials (from scans on SOL and CALC) for conversion
    
    
    
    %%% Calculate individual Norm offset for DIRECTION
    % trial 1
    [convert_ind_direction_b1_volt] = calculate_direction_constants(horzcat(filepath, dm_CPM_calc_NX{line}));
    % trial 2
    [convert_ind_direction_b2_volt] = calculate_direction_constants(horzcat(filepath, dm_CPM_sol_NX{line}));
    % average
    convert_norm_direction_b = mean([convert_ind_direction_b1_volt convert_ind_direction_b2_volt]);
    
    
    
    %%% Calculate individual Norm conversion factors (y = ax + b) for ANGLE
    % Reads raw data from Norm, filter
    % Produces individual conversion factors for angle, prints report
    
    % trial 1
    [convert_ind_angle_a1, convert_ind_angle_b1] = calculate_angle_constants(angle_cutoff, horzcat(filepath, dm_CPM_calc_NX{line}), side);
    % trial 2
    [convert_ind_angle_a2, convert_ind_angle_b2] = calculate_angle_constants(angle_cutoff, horzcat(filepath, dm_CPM_sol_NX{line}), side);
    % average
    convert_norm_angle_a = mean([convert_ind_angle_a1 convert_ind_angle_a2]);
    convert_norm_angle_b = mean([convert_ind_angle_b1 convert_ind_angle_b2]);
    
    
    
    %%% Calculate individual Norm conversion factors (y = ax + b) for VELOCITY
    % Reads raw data from Norm, filter
    % Produces individual conversion factors for velocity, prints report
    
    % trial 1
    [convert_ind_velocity_a1, convert_ind_velocity_b1] = calculate_velocity_constants(velocity_cutoff, horzcat(filepath, dm_CPM_calc_NX{line}), side);
    % trial 2
    [convert_ind_velocity_a2, convert_ind_velocity_b2] = calculate_velocity_constants(velocity_cutoff, horzcat(filepath, dm_CPM_sol_NX{line}), side);
    % average
    convert_norm_velocity_a = mean([convert_ind_velocity_a1 convert_ind_velocity_a2]);
    convert_norm_velocity_b = mean([convert_ind_velocity_b1 convert_ind_velocity_b2]);
    
    
    
    %%% Calculate individual Norm conversion factors (y = ax + b) for TORQUE
    % Uses constants from direction, angle, velocity
    % Uses default A constant from Norm, computes B constant individually based on default A constant
    % Produces individual conversion factors for torque, and prints report
    
    % torque, A constant
    convert_norm_torque_a = norm_volt_per_nm_a;
    
    % make an array of the three mV values for B constants from direction 7, angle 2, velocity 6
    convert_ind_torque_b_volt = [convert_norm_direction_b convert_norm_angle_b/convert_norm_angle_a convert_norm_velocity_b/convert_norm_velocity_a];
    % create new B constant by combining offset from Noraxon (in mv * A constant) with default offset from Norm system (in mv * A constant)
    convert_norm_torque_b = (mean(convert_ind_torque_b_volt) * norm_volt_per_nm_a) + norm_volt_per_nm_b;
    
    
    
    % output torque conversion numbers to screen, as text
    cprintf('magenta', horzcat('INDI Torque conversion factors: a = ', num2str(convert_norm_torque_a), ', b = ', num2str(convert_norm_torque_b), '. Offset in millivolt = ', num2str(mean(convert_ind_torque_b_volt)), ' mV.\n' ));
    
    
    
    
    
    
elseif strcmpi(act_or_pas,'stiffness') == 1
    dm_trial = 'STIFFNESS';
    % stiffness analyses only need ANGLE data. Using two CPM trials (from scans on SOL and CALC) for conversion.
    
    %%% Calculate individual Norm conversion factors (y = ax + b) for ANGLE
    % Reads raw data from Norm, filter
    % Produces individual conversion factors for angle, prints report
    
    % old line from stiffness, probably obsolete, replaced by lines below as for passive trials.
    %[convert_norm_angle_a, convert_norm_angle_b] = calculate_angle_constants(angle_cutoff, horzcat(filepath, dm_CPM_calc_NX{line}));
    
    % trial 1
    [convert_ind_angle_a1, convert_ind_angle_b1] = calculate_angle_constants(angle_cutoff, horzcat(filepath, dm_CPM_calc_NX{line}), side);
    % trial 2
    [convert_ind_angle_a2, convert_ind_angle_b2] = calculate_angle_constants(angle_cutoff, horzcat(filepath, dm_CPM_sol_NX{line}), side);
    % average
    convert_norm_angle_a = mean([convert_ind_angle_a1 convert_ind_angle_a2]);
    convert_norm_angle_b = mean([convert_ind_angle_b1 convert_ind_angle_b2]);
    
    %%% variables that are not needed for stiffness trials:
    convert_norm_torque_a = 0;
    convert_norm_torque_b = 0;
    convert_norm_velocity_a = 0;
    convert_norm_velocity_b = 0;
    convert_norm_direction_b = 0;
    
    
    
    
else % ACTIVE
    dm_trial = 'ACTIVE';
    %%% note MMM 01.08.2016  ---- TODO MMM
    %%% analyses for active trials not completed
    %%% references to files from which conversion factors will be computed
    %%% should be modified.
    %%% should probably replace:
        % horzcat('data_convert\', dm_isokin45_NX))
        % horzcat(filepath, dm_isokin45_NX{line}), side)
    
    %%% this is the old, manual method:
    %%% retrieve names of files to be used for Norm volt to values conversion
    if strcmpi(subjectno,'1') == 1 && strcmpi(side,'L') == 1 && strcmpi(timepoint,'POST') == 1
        dm_isokin45_NX = '01_POST_isokinetic_5L_trial3.tsv';
    else
        % error!!!
    end
    
    
    
    
    %%% Calculate individual Norm offset for DIRECTION
    convert_norm_direction_b = calculate_direction_constants(horzcat('data_convert\', dm_isokin45_NX));
    
    
    
    
    %%% Calculate individual Norm conversion factors (y = ax + b) for ANGLE
    % Reads raw data from Norm, filter
    % Produces individual conversion factors for angle, prints report
    [convert_norm_angle_a, convert_norm_angle_b] = calculate_angle_constants_active(angle_cutoff_active, horzcat('data_convert\', dm_isokin45_NX));
    
    
    
    %%% Calculate individual Norm conversion factors (y = ax + b) for VELOCITY
    % Reads raw data from Norm, filter
    % Produces individual conversion factors for velocity, prints report
    [convert_norm_velocity_a, convert_norm_velocity_b] = calculate_velocity_constants_active(velocity_cutoff_active, horzcat('data_convert\', dm_isokin45_NX));
    
    
    
    %%% Calculate individual Norm conversion factors (y = ax + b) for TORQUE
    % Uses constants from direction, angle, velocity
    % Uses default A constant from Norm, computes B constant individually based on default A constant
    % Produces individual conversion factors for torque, and prints report
    
    % torque, A constant
    convert_norm_torque_a = norm_volt_per_nm_a;
    
    % make an array of the three mV values for B constants from direction 7, angle 2, velocity 6
    convert_ind_torque_b_volt = [convert_norm_direction_b convert_norm_angle_b/convert_norm_angle_a convert_norm_velocity_b/convert_norm_velocity_a];
    % create new B constant by combining offset from Noraxon (in mv * A constant) with default offset from Norm system (in mv * A constant)
    convert_norm_torque_b = (mean(convert_ind_torque_b_volt) * norm_volt_per_nm_a) + norm_volt_per_nm_b;
    
    % output torque conversion numbers to screen, as text
    cprintf('magenta', horzcat('INDI Torque conversion factors: a = ', num2str(convert_norm_torque_a), ', b = ', num2str(convert_norm_torque_b), '. Offset in millivolt = ', num2str(mean(convert_ind_torque_b_volt)), ' mV.\n' ));
    
    
    
    
end