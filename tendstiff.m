%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main file for tendon stiffness
% Marie Moltubakk 17.5.2013
% 
% Note 1:
% The scripts assume a certain structure for input files from Tracker and
% Noraxon. I.e. number of and order of EMG and other channels. If these
% scripts are to be used for other projects, some of the code must be
% modified. These lines are marked % PROJECTSPECIFIC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all



%% PLOTS - determine which plots to display
global plot_achilles plot_norm plot_emg plot_check plot_us plot_conversion subject_id

plot_check = 1; % turn on/off checkpoint plots (leave only stiffness)
plot_achilles = 1; % turn on/off all troubleshoot plots
plot_norm = 0; % show torque before and after initial lowpass filter
plot_emg = 0;  % RMS 3 EMG channels per trial
plot_us = 0;
plot_conversion = 1;


%% Set constants and globals % PROJECTSPECIFIC
% declare for later use:
% variables for NORM conversion factors calculated from actual data
global convert_norm_angle_a convert_norm_angle_b convert_norm_torque_a convert_norm_torque_b convert_norm_velocity_a convert_norm_velocity_b convert_norm_direction_b
    
% sampling frequencies
global us_zerodispframes noraxonfreq freq_default
us_zerodispframes = 1; % No of US frames to average as zero displacement
noraxonfreq = 1500; % sampling frequency of noraxon data
freq_default = 100; % output frequency for noraxon data without US video (MVC, etc)

% default conversion factors from Norm manuals and Achilles calibration
global convert_achilles norm_volt_per_degree norm_volt_per_velocity norm_volt_per_nm_a norm_volt_per_nm_b
convert_achilles = -81.9909;  % conversion factor ACHILLES torque, V -> Nm
norm_volt_per_degree = (2048*((1*16)+0)/1024)*(10000/32768); %   (Gain * (Sampled Position + Offset) / 1024) * (10v/32768) = Volt-value, Sampled Position is in units of 1/16 degree.
norm_volt_per_velocity = (1024*((1*16)+0)/1024)*(10000/32768); %   (Gain * (Sampled Velocity + Offset) / 1024) * (10v/32768) = Volt-value. Sampled Velocity is in units of 1/16 degree.
% norm_volt_per_nm = (1024*((1*(1.355818*32768/500))+0)/1024)*(10000/32768); % 27.116 µV/Nm   Sampled Torque is in units of Foot-Pounds * 32768 / 500
norm_volt_per_nm_a = 0.07490071; % from file "conversion volt-torque 2014-DES NEW DATA.xlsx"
norm_volt_per_nm_b = 0.69283161;

% cutoff frequencies for filters
global emg_bandpass emg_rms_ms mvc_window_ms 
emg_bandpass = [10/(noraxonfreq/2) 500/(noraxonfreq/2)]; % cutoff frequencies for EMG butterworth bandpass filter
emg_rms_ms = 100; % milliseconds RMS window, must be divisible by 4 - ref Basmajian 1985 = 50 ms, ref Aagaard paper Passive tensile stress = 200 ms
mvc_window_ms = 500; % milliseconds window for determining MVC torque and EMG
global angle_cutoff velocity_cutoff torque_cutoff_bandstop torque_cutoff_active angle_cutoff_active velocity_cutoff_active
angle_cutoff = 20/(noraxonfreq/2); % 9/(noraxonfreq/2); % cutoff freq, Norm angle PASSIVE - ref Winter 1990 = 15 hz. Kongsgaard = 8hz
angle_cutoff_active = 20/(noraxonfreq/2); %  20/(noraxonfreq/2);
velocity_cutoff = 20/(noraxonfreq/2); %  12/(noraxonfreq/2);
velocity_cutoff_active = 20/(noraxonfreq/2); %  12/(noraxonfreq/2);
torque_cutoff_bandstop = [0.36/(noraxonfreq/2) 0.42/(noraxonfreq/2)]; % bandstop freq to eliminate noise from Norm engine
torque_cutoff_active = 10/(noraxonfreq/2); % cutoff freq, Norm filtering - ref Winter 1990 = 15 hz. Kongsgaard = 8hz

% Average stiffness across X N
forceintervals = 100; 

% column placement in Noraxon data
global column_EMG_start column_EMG_end column_l_gm column_r_gm column_l_gl column_r_gl column_l_sol column_r_sol column_l_tibant column_r_tibant column_gonio column_norm_angle column_norm_torque column_norm_velocity column_norm_direction column_achilles
column_EMG_start = 2; 
column_EMG_end = 9;
column_l_gm = 2;
column_r_gm = 3;
column_l_gl = 4;
column_r_gl = 5;
column_l_sol = 6;
column_r_sol = 7;
column_l_tibant = 8;
column_r_tibant = 9;
column_gonio = 10;
column_norm_angle = 11;
column_norm_torque = 12;
column_norm_velocity = 13;
column_norm_direction = 14;
column_achilles = 15;


%% Read datamaster file, to connect corresponding data files
% Produces arrays with file names and variables per trial, to be retrieved later

global dm_subjectno dm_timepoint dm_side dm_trial 
global dm_stiff1_NX dm_stiff1_US dm_stiff1_US_frame dm_stiff2_NX dm_stiff2_US dm_stiff2_US_frame dm_stiff3_NX dm_stiff3_US dm_stiff3_US_frame 
global dm_heel1_NX dm_heel1_US dm_heel1_US_frame dm_heel2_NX dm_heel2_US dm_heel2_US_frame dm_heel3_NX dm_heel3_US dm_heel3_US_frame
global dm_MVC_PF dm_MVC_DF dm_CPM_calc_NX dm_CPM_calc_US dm_CPM_calc_US_frame dm_leg_length
global dm_cutforce %new2014-04-14
dm_filename = 'data/datamaster.tsv';
dm_columns = 29; % number of data columns entered per subject % PROJECTSPECIFIC %new2014-04-14 - was 28
linestotal = read_datamaster(dm_filename,dm_columns);


%% preallocate
    %preallocate output arrays
    all_stiff_output_head = {'Subject', 'Time', 'Side', 'Trial', 'Stiff coeff 1', 'Stiff coeff 2', 'Stiff coeff 3', 'Stiff R2', ...
        'Ramp force cutoff (N)', 'Ramp force max (N)', 'Stiffness 80-100 (N/mm)', 'Stiffness 90-100 (N/mm)', 'PF MVC (N)', 'AT moment arm (m)', 'Rotation correction (mm/deg)'}; % PROJECTSPECIFIC
    all_stiff_output = zeros(ceil(linestotal),length(all_stiff_output_head)-4); 
    all_stiff_output_txt = cell(ceil(linestotal),4);


%% Loop through all lines in datamaster file (except header line)
for line = 1:linestotal
    
    
	%% subject/trial identifier
    subject_id = horzcat('subject ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
    cprintf('*black', horzcat('----------------', subject_id, '------------------\n'))
    
    
%     %% calculate NORM conversion factors
%     % retrieve conversion constants for Norm data
%     [convert_norm_angle_a, convert_norm_angle_b, convert_norm_torque_a, convert_norm_torque_b, convert_norm_velocity_a, convert_norm_velocity_b, convert_norm_direction_b] = calculate_norm_constants_complete(dm_subjectno{line}, dm_side{line}, dm_timepoint{line}, line, 'stiffness');
    
    
    %% Calculate conversion factors for tibialis anterior CO-ACTIVATION
    
    % Produce individual conversion factors for angle
    [convert_norm_angle_a, convert_norm_angle_b] = calculate_angle_constants(angle_cutoff, horzcat('data\', dm_CPM_calc_NX{line}), dm_side{line});
    % using only "angle_constants" not "constants_ACTIVE" because this is for NORM data only, and tendstiff uses only passive trials from NORM
    
    % Read co-activation noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
    % Produce a new noraxon data array
    % sending in a length corresponding to 9 seconds (delete anything after 9 sec)
    noraxon_coact = read_noraxon_stiffness(strcat('data/', dm_MVC_DF{line}), freq_default, dm_side{line}, 'MVC dorsi');
    
    % Calculate co-activation constants
    % Read complete, prepared noraxon array + number of frames to average (freq * time)
    % Produce max torque, max EMG constants
    [coact_max_torque,coact_max_EMG] = calculate_coactivation(noraxon_coact, freq_default*(mvc_window_ms/1000), dm_side{line});
    
    
    %% Calculate achilles tendon MOMENT ARM
    
    % Read moment arm trial US data file, determine time stamps, set trigger frame as time = zero
    % Produce US sample frequency, create new US array containing time and displacement
    [usdata_momentarm, usfreq_momentarm] = read_us_file(strcat('data/', dm_CPM_calc_US{line}),str2double(dm_CPM_calc_US_frame{line}), 'CPM calcaneus');
    
    % Read moment arm noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
    % Produce a new noraxon data array
    noraxon_momentarm = read_noraxon_stiffness(strcat('data/', dm_CPM_calc_NX{line}), usfreq_momentarm, dm_side{line}, 'CPM calcaneus');
    
    % Read complete, prepared noraon array + prepared us array
    % Produce AT moment arm constant
    at_momentarm = calculate_momentarm(0,0, dm_leg_length{line});
    
    
    %% Calculate relationship between calc displacement and ANKLE ROTATION
    
    % REUSING us file from moment arm
    % REUSING noraxon file from moment arm
    
    % Read complete, prepared noraxon array + prepared us array
    % Produce ankle rotation constant, mm/degree
    if strcmp(subject_id,'subject 4 R PRE SOL')
        %very special case for this trial...
        at_rotation_const = -0.14;
    else %normally
        at_rotation_const = calculate_rotation_correction(noraxon_momentarm, usdata_momentarm);
    end
    
    
    %% Calculate MVC for plantarflexion

    % Read MVC noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
    % Produce a new noraxon data array
    % sending in a length corresponding to 9 seconds (delete anything after 9 sec)
    noraxon_MVC = read_noraxon_stiffness(strcat('data/', dm_MVC_PF{line}), freq_default, dm_side{line}, 'MVC dorsi');
    
    % Calculate co-activation constants
    % Read complete, prepared noraxon array + number of frames to average (freq * time)
    % Produce max torque, max EMG constants
    [plantflex_max_torque,plantflex_GM_max_EMG] = calculate_MVC(noraxon_MVC, freq_default*(mvc_window_ms/1000), at_momentarm, dm_side{line});
    
    
    %% Calculations for 3x ramp MTJ files
    
    % Allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case in final_stiffness
    trial_force_max(1:6) = 100000;
    if(strcmp(dm_stiff1_NX{line}, 'null'))
        time_force_displ_mtj1 = zeros(0);
    else 
        [time_force_displ_mtj1,trial_force_max(1)] = extract_force_displ_singletrial(dm_stiff1_NX{line}, dm_stiff1_US{line}, dm_stiff1_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, at_rotation_const, dm_side{line}, 'MTJ1');
    end
    if(strcmp(dm_stiff2_NX{line}, 'null'))
        time_force_displ_mtj2 = zeros(0);
    else
        [time_force_displ_mtj2,trial_force_max(2)] = extract_force_displ_singletrial(dm_stiff2_NX{line}, dm_stiff2_US{line}, dm_stiff2_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, at_rotation_const, dm_side{line}, 'MTJ2');
    end
    if(strcmp(dm_stiff3_NX{line}, 'null'))
        time_force_displ_mtj3 = zeros(0);
    else
        [time_force_displ_mtj3,trial_force_max(3)] = extract_force_displ_singletrial(dm_stiff3_NX{line}, dm_stiff3_US{line}, dm_stiff3_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, at_rotation_const, dm_side{line}, 'MTJ3');
    end
    
    
    %% Calculations for 3x ramp calcaneus scans
    if(strcmp(dm_heel1_NX{line}, 'null'))
        time_force_displ_otj1 = zeros(0);
    else
        [time_force_displ_otj1,trial_force_max(4)] = extract_force_displ_singletrial(dm_heel1_NX{line}, dm_heel1_US{line}, dm_heel1_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, 0, dm_side{line}, 'OTJ1');
    end
    if(strcmp(dm_heel2_NX{line}, 'null'))
        time_force_displ_otj2 = zeros(0);
    else
        [time_force_displ_otj2,trial_force_max(5)] = extract_force_displ_singletrial(dm_heel2_NX{line}, dm_heel2_US{line}, dm_heel2_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, 0, dm_side{line}, 'OTJ2');
    end
    if(strcmp(dm_heel3_NX{line}, 'null'))
        time_force_displ_otj3 = zeros(0);
    else
        [time_force_displ_otj3,trial_force_max(6)] = extract_force_displ_singletrial(dm_heel3_NX{line}, dm_heel3_US{line}, dm_heel3_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, 0, dm_side{line}, 'OTJ3');
    end
    
     
    %% Final stiffness
    % Read time-force-displacement data
    % Produce stiffness equation
	[stiff_eq, stiff_gof, force_elong_array, stiff_force_cutoff] = final_stiffness(time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3, time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3, forceintervals, dm_cutforce{line}, trial_force_max);
    
        
    %% calculate stiffness for last 10 and 20% of ind max:
    stiff80 = calculate_stiffness(stiff_eq, stiff_force_cutoff, 0.8, 1.0); % last two variables are percent range, from 0.00 to 1.00
    stiff90 = calculate_stiffness(stiff_eq, stiff_force_cutoff, 0.9, 1.0);
    
    
    %% Strain rate % TODOLATER
    % strain_rate = calculate_strain_rate(time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3, time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3);
    
    
    %% save individual data to common array + force-elong-file

    % add stiffness data points to xls file per subject
    if ispc
        filename_output = strcat('data_output/stiff_forceoutput_', subject_id, '.xls');
        xlswrite(filename_output, {'Elongation (mm)','Force (N)'}, 'Final stiff', 'A1')
        xlswrite(filename_output, force_elong_array, 'Final stiff', 'A2')
    else
        filename_output = strcat('data_output/stiff_forceoutput_', subject_id, '_averaged.csv');
        csvwrite(filename_output, force_elong_array)
    end
    
    % add data to a common array for all subjects    
    all_stiff_output_txt(line,:) = [dm_subjectno(line) dm_timepoint(line) dm_side(line) dm_trial(line)];
    all_stiff_output(line,:) = [coeffvalues(stiff_eq) stiff_gof.rsquare stiff_force_cutoff min(trial_force_max) stiff80 stiff90 plantflex_max_torque at_momentarm at_rotation_const];
    
end
%% LOOP finished


%% calculate stiffness at group common force levels

% select common force
all_stiff_col = 5; % 5 for cutoff force level, 6 for max force level
stiff_common_force_100 = min(all_stiff_output(:,all_stiff_col)); % or y2

% ... for each ... TODO MMM GOON
stiff_coeff1 = all_stiff_output(1,1); % TMP - subj 1
stiff_coeff2 = all_stiff_output(1,2); % TMP

% calculation for a given force range
stiff_common_force_lvl = 0.9; % %VAR
stiff_common_force_percent = stiff_common_force_100 * stiff_common_force_lvl;
stiff_common_x2_pos = (-stiff_coeff2 + (sqrt(stiff_coeff2^2 - (4*stiff_coeff1*-stiff_common_force_100)))) / (2 * stiff_coeff1);
stiff_common_x1_pos = (-stiff_coeff2 + (sqrt(stiff_coeff2^2 - (4*stiff_coeff1*-stiff_common_force_percent)))) / (2 * stiff_coeff1);
stiff_common = (stiff_common_force_100-stiff_common_force_percent) / (stiff_common_x2_pos-stiff_common_x1_pos);


% MMM TODO GOON - add to output array

% TODO - repeat for 70-100%

% repeat for actual force instead of cutoff force


%% other calculations from Melina datasheet?
% MMM TODO GOON


%% save array with individual data to file

% write xls
if ispc
    filename_output = strcat('data_output/all_stiff_output_', datestr(now, 'yyyy-mm-dd HH-MM'));
    xlswrite(filename_output, all_stiff_output_head, 1, 'A1')
    xlswrite(filename_output, all_stiff_output_txt, 1, 'A2')
    xlswrite(filename_output, all_stiff_output, 1, 'E2')
else
    filename_output = strcat('data_output/all_stiff_output_', datestr(now, 'yyyy-mm-dd HH-MM'), '.csv');
    csvwrite(filename_output, all_stiff_output)
end
