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

global plot_achilles plot_norm plot_emg plot_check plot_us subject_id

plot_check = 1; % turn on/off checkpoint plots (leave only fits and stiffness)
plot_achilles = 0; % turn on/off all troubleshoot plots
plot_norm = 0; % show torque before and after initial lowpass filter
plot_emg = 0;  % RMS 3 EMG channels per trial
plot_us = 0;




%%% Set constants % PROJECTSPECIFIC

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
norm_volt_per_nm_a = 0.07490071; % from file "MMM conversion volt-torque 2014-DES NEW DATA.xlsx"
norm_volt_per_nm_b = 0.69283161;

% conversion factors, calculated per subject:
% old data TMP TODO DELETE
% convertion factors Norm, µV -> angle/torque/velocity: 
% convert_norm_old(1) = 0.09918;      %convert_norm_angle_a % no longer needed, because angle factors are calculated per trial
% convert_norm_old(2) = -168.714827;  %convert_norm_angle_b
% convert_norm_old(3) = 3518.257824066; %convert_norm_angle_c
% convert_norm_old(5) = 0.076965;     %convert_norm_torque_a
% convert_norm_old(6) = 7.862109;     %convert_norm_torque_b
% convert_norm_old(7) = 0.20535;      %convert_norm_velocity_a
% convert_norm_old(8) = 10.92228;     %convert_norm_velocity_b
% new:
% convert_norm(1) = convert_ind_angle_a;
% convert_norm(2) = convert_ind_angle_b;
% convert_norm(3) = 0;  % TORQUE A
% convert_norm(4) = (mean(convert_ind_torque_b_volt) * norm_volt_per_nm_a) + norm_volt_per_nm_b; % TORQUE B
% convert_norm(5) = convert_ind_velocity_a;
% convert_norm(6) = convert_ind_velocity_b;
% convert_norm(7) = convert_ind_direction_b_volt;

% cutoff frequencies for filters
global emg_bandpass emg_rms_ms mvc_window_ms 
emg_bandpass = [10/(noraxonfreq/2) 500/(noraxonfreq/2)]; % cutoff frequencies for EMG butterworth bandpass filter
emg_rms_ms = 100; % milliseconds RMS window, must be divisible by 4 - ref Basmajian 1985 = 50 ms, ref Aagaard paper Passive tensile stress = 200 ms
mvc_window_ms = 500; % milliseconds window for determining MVC torque and EMG
global angle_cutoff velocity_cutoff torque_cutoff_bandstop torque_cutoff_active angle_cutoff_active velocity_cutoff_active % TODO MMM add other cutoffs from normconversion?
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




%%% Read datamaster file, to connect corresponding data files
%%% Produces arrays with file names and variables per trial, to be retrieved later

global dm_subjectno dm_timepoint dm_side dm_trial 
global dm_stiff1_NX dm_stiff1_US dm_stiff1_US_frame dm_stiff2_NX dm_stiff2_US dm_stiff2_US_frame dm_stiff3_NX dm_stiff3_US dm_stiff3_US_frame 
global dm_heel1_NX dm_heel1_US dm_heel1_US_frame dm_heel2_NX dm_heel2_US dm_heel2_US_frame dm_heel3_NX dm_heel3_US dm_heel3_US_frame
global dm_MVC_PF dm_MVC_DF dm_CPM_calc_NX dm_CPM_calc_US dm_CPM_calc_US_frame dm_leg_length
global dm_cutforce %new2014-04-14
dm_filename = 'data/datamaster.tsv';
dm_columns = 29; % number of data columns entered per subject % PROJECTSPECIFIC %new2014-04-14 - was 28
linestotal = read_datamaster(dm_filename,dm_columns);




%%% Loop through all lines in datamaster file (except header line)
for line = 1:linestotal;
    
    
    
    
	%%% subject/trial identifier
    subject_id = horzcat('subject ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
    disp(sprintf(horzcat('----------------', subject_id, '------------------')))

    
    
%%% 2015-01-28 - moved into read_noraxon_file
%     %%% Calculate individual Norm conversion factors (a,b) for ANGLE
%     plot_achilles = 0;
%     plot_norm = 0;
%     plot_emg = 0; 
%     plot_us = 0;
%     
%     % Read raw ??? torque data from Norm, filter ???
%     % Produce individual conversion factors for angle
%     [convert_ind_angle_a, convert_ind_angle_b] = calculate_angle_constants(angle_cutoff, horzcat('data\', dm_CPM_calc_NX{line})); % evt add ,dm_side{line});
%     % place individually calculated constant into array with common constants 
%     convert_norm_ind_passive(1) = convert_ind_angle_a;
%     convert_norm_ind_passive(2) = convert_ind_angle_b;
%       function [convert_norm_angle_a, convert_norm_angle_b, convert_norm_torque_a, convert_norm_torque_b, convert_norm_velocity_a, convert_norm_velocity_b, convert_norm_direction_b] = calculate_norm_constants(act_or_pas, FP)
    
    
    
    %%% Calculate conversion factors for tibialis anterior CO-ACTIVATION
    plot_achilles = 0;
    plot_norm = 0; 
    plot_emg = 0;
    plot_us = 0;
    
    
    % Produce individual conversion factors for angle
    global convert_norm_angle_a convert_norm_angle_b
    [convert_norm_angle_a, convert_norm_angle_b] = calculate_angle_constants(angle_cutoff, horzcat('data\', dm_CPM_calc_NX{line}), dm_side{line});
    % place individually calculated constant into array with common constants 

    % Read co-activation noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
    % Produce a new noraxon data array
    % sending in a length corresponding to 9 seconds (delete anything after 9 sec)
    noraxon_coact = read_noraxon_stiffness(strcat('data/', dm_MVC_DF{line}), freq_default, dm_side{line}, 'MVC dorsi');
    
    % Calculate co-activation constants
    % Read complete, prepared noraxon array + number of frames to average (freq * time)
    % Produce max torque, max EMG constants
    [coact_max_torque,coact_max_EMG] = calculate_coactivation(noraxon_coact, freq_default*(mvc_window_ms/1000), dm_side{line});
    
    
    
    %%% Calculate achilles tendon MOMENT ARM
    plot_achilles = 0;
    plot_norm = 0; 
    plot_emg = 0; 
    plot_us = 1;
    
    % Read moment arm trial US data file, determine time stamps, set trigger frame as time = zero
    % Produce US sample frequency, create new US array containing time and displacement
    [usdata_momentarm, usfreq_momentarm] = read_us_file(strcat('data/', dm_CPM_calc_US{line}),str2num(dm_CPM_calc_US_frame{line}), 'CPM calcaneus');
    
    % Read moment arm noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
    % Produce a new noraxon data array
    noraxon_momentarm = read_noraxon_stiffness(strcat('data/', dm_CPM_calc_NX{line}), usfreq_momentarm, dm_side{line}, 'CPM calcaneus');
    
    % Read complete, prepared noraon array + prepared us array
    % Produce AT moment arm constant
    at_momentarm = calculate_momentarm(0,0, dm_leg_length{line}); % TODO MMM now using fixed moment arm. replace 0,0 -> (noraxon_momentarm, usdata_momentarm)
    
    
    
    %%% Calculate relationship between calc displacement and ANKLE ROTATION
    plot_achilles = 0;
    plot_norm = 0; %  turn on for crrection check plot 2
    plot_emg = 0; 
    plot_us = 0;
    
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
    
    
    %%% Calculate MVC for plantarflexion
    plot_us = 0;
    plot_norm = 0;
    
    % Read MVC noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
    % Produce a new noraxon data array
    % sending in a length corresponding to 9 seconds (delete anything after 9 sec)
    noraxon_MVC = read_noraxon_stiffness(strcat('data/', dm_MVC_PF{line}), freq_default, dm_side{line}, 'MVC dorsi');
    
    % Calculate co-activation constants
    % Read complete, prepared noraxon array + number of frames to average (freq * time)
    % Produce max torque, max EMG constants
    [plantflex_max_torque,plantflex_GM_max_EMG] = calculate_MVC(noraxon_MVC, freq_default*(mvc_window_ms/1000), at_momentarm, dm_side{line});

    
    
    %%% Calculations for 3x ramp MTJ files
    plot_achilles = 1;
    plot_norm = 0; 
    plot_emg = 0;
    plot_us = 0;
    plot_check = 1;
    
    % Allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case in final_stiffness
    rampmax(1:6) = 100000;
%    plot_us = 0; % TMP
    if(strcmp(dm_stiff1_NX{line}, 'null'))
        time_force_displ_mtj1 = zeros(0);
    else 
        [time_force_displ_mtj1,rampmax(1)] = extract_force_displ_singletrial(dm_stiff1_NX{line}, dm_stiff1_US{line}, dm_stiff1_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, at_rotation_const, dm_side{line}, 'MTJ1');
    end
%    plot_us = 0; % TMP
    if(strcmp(dm_stiff2_NX{line}, 'null'))
        time_force_displ_mtj2 = zeros(0);
    else
        [time_force_displ_mtj2,rampmax(2)] = extract_force_displ_singletrial(dm_stiff2_NX{line}, dm_stiff2_US{line}, dm_stiff2_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, at_rotation_const, dm_side{line}, 'MTJ2');
    end
%    plot_us = 0; % TMP
    if(strcmp(dm_stiff3_NX{line}, 'null'))
        time_force_displ_mtj3 = zeros(0);
    else
        [time_force_displ_mtj3,rampmax(3)] = extract_force_displ_singletrial(dm_stiff3_NX{line}, dm_stiff3_US{line}, dm_stiff3_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, at_rotation_const, dm_side{line}, 'MTJ3');
    end
    
    
    
    %%% Calculations for 3x ramp calcaneus scans
%    plot_us = 0; % TMP
    if(strcmp(dm_heel1_NX{line}, 'null'))
        time_force_displ_otj1 = zeros(0);
    else
        [time_force_displ_otj1,rampmax(4)] = extract_force_displ_singletrial(dm_heel1_NX{line}, dm_heel1_US{line}, dm_heel1_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, 0, dm_side{line}, 'OTJ1');
    end
%    plot_us = 0; % TMP
    if(strcmp(dm_heel2_NX{line}, 'null'))
        time_force_displ_otj2 = zeros(0);
    else
        [time_force_displ_otj2,rampmax(5)] = extract_force_displ_singletrial(dm_heel2_NX{line}, dm_heel2_US{line}, dm_heel2_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, 0, dm_side{line}, 'OTJ2');
    end
%    plot_us = 0; % TMP
    if(strcmp(dm_heel3_NX{line}, 'null'))
        time_force_displ_otj3 = zeros(0);
    else
        [time_force_displ_otj3,rampmax(6)] = extract_force_displ_singletrial(dm_heel3_NX{line}, dm_heel3_US{line}, dm_heel3_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, 0, dm_side{line}, 'OTJ3');
    end
    
    
        
    %%% Final stiffness
    plot_check = 1;
    plot_achilles = 1;
    plot_norm = 0; 
    plot_emg = 0; 
    plot_us = 0;
    
    % Read time-force-displacement data
    % Produce stiffness equation
	[stiff_eq,stiff_gof,stiff_frames,stiff_usedforce,stiff_maxforce] = final_stiffness(time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3, time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3, forceintervals, dm_cutforce{line}, rampmax);
    
    
    
    %%% calculate stiffness for last 10 and 20% of ind max:
    stiff80 = calculate_stiffness(stiff_eq, stiff_usedforce, 0.8, 1.0); % last two variables are percent range, from 0.00 to 1.00
    stiff90 = calculate_stiffness(stiff_eq, stiff_usedforce, 0.9, 1.0);
    
    
    
    %%% Strain rate % TODOLATER
    % strain_rate = calculate_strain_rate(time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3, time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3);
    
    
        
    %%% Output final data to file

    % add stiffness data points to xls file per subject
    if ispc
        filename_output = strcat('data_output/forceoutput_', subject_id, '.xls');
        xlswrite(filename_output, {'Elongation (mm)','Force (N)'}, 'Final stiff', 'A1')
        xlswrite(filename_output, stiff_frames, 'Final stiff', 'A2')
    else
        filename_output = strcat('data_output/forceoutput_', subject_id, '_averaged.csv');
        csvwrite(filename_output, stiff_frames)
    end
    
    % add data to a common array for all subjects    
    all_output_txt(line,1) = {dm_subjectno{line}};
    all_output_txt(line,2) = {dm_timepoint{line}};
    all_output_txt(line,3) = {dm_side{line}};
    all_output_txt(line,4) = {dm_trial{line}};
    all_output(line,1:3) = coeffvalues(stiff_eq);
    all_output(line,4) = stiff_gof.rsquare;
    all_output(line,5) = stiff_usedforce;
    all_output(line,6) = min(rampmax);
    all_output(line,7) = stiff80;
    all_output(line,8) = stiff90;
    all_output(line,9) = plantflex_max_torque;
    all_output(line,10) = at_momentarm;
    all_output(line,11) = at_rotation_const;
    
end



%%% Output key variables for all subjects to file

% write xls
if ispc
    filename_output = strcat('data_output/all_output_', datestr(now, 'yyyy-mm-dd HH-MM'));
    all_output_head = {'Subject', 'Time', 'Side', 'Trial', 'Stiff coeff 1', 'Stiff coeff 2', 'Stiff coeff 3', 'Stiff R2', ...
        'Stiff cut F (N)', 'Stiff max F (N)', 'Stiffness 80-100 (N/mm)', 'Stiffness 90-100 (N/mm)', 'PF MVC (N)', 'AT moment arm (m)', 'Rotation correction (mm/deg)'}; % PROJECTSPECIFIC
    xlswrite(filename_output, all_output_head, 1, 'A1')
    xlswrite(filename_output, all_output_txt, 1, 'A2')
    xlswrite(filename_output, all_output, 1, 'E2')
else
    filename_output = strcat('data_output/all_output_', datestr(now, 'yyyy-mm-dd HH-MM'), '.csv');
    csvwrite(filename_output, all_output)
end






