%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script for determining conversion factors Norm/Noraxon, passive&active trials
% 
% active: using isometric P10 and D00 x2 trials each + isokinetic 45deg/s trial
% passive: using 2 trials: slow CPM +5 to -10 degrees, calc scan and sol scan
% stiffness: using two CPM trials (from scans on SOL and CALC) for conversion.
% 
% Marie Moltubakk 28.1.2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [convert_norm_angle_a, convert_norm_angle_b, convert_norm_torque_a, convert_norm_torque_b, convert_norm_velocity_a, convert_norm_velocity_b, convert_norm_direction_b] = calculate_norm_constants_complete(subjectno, side, timepoint, line, type_of_test)
% input:
% subjectno
% side
% timepoint
% type_of_test (active / passive / stiffness)



%% Set constants % PROJECTSPECIFIC
global dm_isokinP30 dm_isokinP45 dm_CPM_calc_NX dm_CPM_sol_NX
global dm_isomet_P10_2 dm_isomet_D00_2 dm_isomet_P10_1 dm_isomet_D00_1 dm_isomet_D05_1 dm_isomet_D05_2 % dm_isomet_D10_1 dm_isomet_D10_2 dm_isomet_D15_1 dm_isomet_D15_2
global filepath
global plot_conversion

% default conversion factors from Norm manuals and Achilles calibration
global norm_volt_per_nm_a norm_volt_per_nm_b % convert_achilles norm_volt_per_degree norm_volt_per_velocity

% cutoff frequencies for filters
global angle_cutoff velocity_cutoff angle_cutoff_active velocity_cutoff_active % torque_cutoff_bandstop torque_cutoff_active
%global emg_bandpass emg_rms_ms mvc_window_ms



if strcmpi(type_of_test,'passive') == 1
    %% PASSIVE trials
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
    if plot_conversion
        cprintf('magenta', horzcat('INDI Torque conversion factors: a = ', num2str(convert_norm_torque_a), ', b = ', num2str(convert_norm_torque_b), '. Offset in millivolt = ', num2str(mean(convert_ind_torque_b_volt)), ' mV.\n' ));
    end
    
    
    
elseif strcmpi(type_of_test,'stiffness') == 1
    %% STIFFNESS trials
    % stiffness analyses only need ANGLE data. Using two CPM trials (from scans on SOL and CALC) for conversion.
    
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
    
    
    
else
    %% ACTIVE trials
    
    % file/velocity options ISOKINETIC: dm_isokinD30 dm_isokinP30 dm_isokinP45 dm_isokinP60 dm_isokinP90
    if subjectno < 100 % intervention/CON
        file = horzcat(filepath, dm_isokinP45{line});
    else % BD
        file = horzcat(filepath, dm_isokinP30{line});
        % for BD, file named P30 contains 2nd isokinetic trial = 45°/s plantar flexion
    end
    velocity = 45; %VAR
    
    % file options ISOMETRIC: global dm_isomet_P10_1 dm_isomet_P10_2 dm_isomet_D00_1 dm_isomet_D00_2 dm_isomet_D05_1 dm_isomet_D05_2 dm_isomet_D10_1 dm_isomet_D10_2 dm_isomet_D15_1 dm_isomet_D15_2
    file1 = horzcat(filepath, dm_isomet_P10_2{line});
    file1b = horzcat(filepath, dm_isomet_P10_1{line});
    file2 = horzcat(filepath, dm_isomet_D00_2{line});
    file2b = horzcat(filepath, dm_isomet_D00_1{line});
    % hack for INT 22 PRE L, where both trials at P10 are with incorrect angle:
    if subjectno == 22 && strcmpi(side,'L') && strcmpi(timepoint,'PRE') && strcmpi(type_of_test,'active')
        file1 = horzcat(filepath, dm_isomet_D05_2{line});
        file1b = horzcat(filepath, dm_isomet_D05_1{line});
    end     
    
    
    %%% Calculate individual Norm offset for DIRECTION
    convert_norm_direction_b = calculate_direction_constants(file);
    

    
%    %%% Calculate individual Norm conversion factors (y = ax + b) for ANGLE

%    2016-11-FG02: REMOVED use of angle constants from isokinetic file. There
%    is no guarantee that each subject performed the test to maximal
%    dorsiflexion angle. Changed to isometric trials.

    % Reads raw data from Norm, filter
    % Produces individual conversion factors for angle, prints report
    [convert_norm_angle_a1, convert_norm_angle_b1] = calculate_angle_constants_active(angle_cutoff_active, file1, file2, side);
    [convert_norm_angle_a2, convert_norm_angle_b2] = calculate_angle_constants_active(angle_cutoff_active, file1b, file2b, side);
    convert_norm_angle_a = mean([convert_norm_angle_a1 convert_norm_angle_a2]);
    convert_norm_angle_b = mean([convert_norm_angle_b1 convert_norm_angle_b2]);
    if subjectno == 22 && strcmpi(side,'L') && strcmpi(timepoint,'PRE') && strcmpi(type_of_test,'active')
        % angle_a is a constant from Norm manufacturer
        convert_norm_angle_b = convert_norm_angle_b * 0.5; % 0.5 because angle range is 5 not 10 degrees. '-' because angle order is swapped.
    end
    
    
    %%% Calculate individual Norm conversion factors (y = ax + b) for VELOCITY
    % Reads raw data from Norm, filter
    % Produces individual conversion factors for velocity, prints report
    [convert_norm_velocity_a, convert_norm_velocity_b] = calculate_velocity_constants_active(velocity_cutoff_active, file, side, velocity);
    
    
    
    %%% Calculate individual Norm conversion factors (y = ax + b) for TORQUE
    % Uses constants from direction, angle, velocity
    % Uses default A constant from Norm, computes B constant individually based on default A constant
    % Produces individual conversion factors for torque, and prints report
    
    % torque, A constant
    convert_norm_torque_a = norm_volt_per_nm_a;
    
    % make an array of the three mV values for B constants from direction 7, angle 2, velocity 6
    % MMM TODO? include conversion factors from ANGLE in averaging below?
    convert_ind_torque_b_volt = [convert_norm_direction_b convert_norm_velocity_b/convert_norm_velocity_a];
    % create new B constant: combining offset from Noraxon (in mv * A constant) with default offset from Norm system (in mv * A constant)
    convert_norm_torque_b = (mean(convert_ind_torque_b_volt) * norm_volt_per_nm_a) + norm_volt_per_nm_b;
    
    % output torque conversion numbers to screen, as text
    if plot_conversion
        cprintf('magenta', horzcat('INDI Torque conversion factors: a = ', num2str(convert_norm_torque_a), ', b = ', num2str(convert_norm_torque_b), '. Offset in millivolt = ', num2str(mean(convert_ind_torque_b_volt)), ' mV.\n' ));
    end
    
    
    
end