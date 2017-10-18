%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEPARATE script for testing determination of conversion factors Norm/Noraxon, passive&active trials, without touching other data
% Marie Moltubakk 27.10.2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all

%%% Plot settings

% select FP to work on here
fp = 10;


global plot_achilles plot_norm plot_emg plot_check plot_us subject_id

plot_check = 0; % turn on/off main checkpoint plots
plot_achilles = 1; % turn on/off troubleshoot plots
plot_norm = 1;
plot_emg = 0;  % RMS 3 EMG channels per trial
plot_us = 1;




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

% conversion factors, calculated per subject (values will be extracted by the script below and inserted to these arrays)
global convert_norm convert_norm_ind_passive convert_norm_ind_active
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
global angle_cutoff velocity_cutoff torque_cutoff_bandstop torque_cutoff_active angle_cutoff_active velocity_cutoff_active
angle_cutoff = 20/(noraxonfreq/2); % cutoff freq, Norm angle PASSIVE - ref Winter 1990 = 15 hz. Kongsgaard = 8hz
angle_cutoff_active = 20/(noraxonfreq/2);
velocity_cutoff = 20/(noraxonfreq/2);
velocity_cutoff_active = 20/(noraxonfreq/2);
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




%%% Read file names of corresponding data files

global dm_subjectno dm_timepoint dm_side dm_trial 

%read list of filenames into variables
filenames

% read subject data from CPM file #1
underscore_indices = strfind(dm_CPM1_NX,'_'); 
dm_subjectno = dm_CPM1_NX(1:2);
dm_timepoint = dm_CPM1_NX(underscore_indices(1)+1:underscore_indices(2)-1);
dm_side = dm_CPM1_NX(underscore_indices(3)+2);




%%%%%%%%%%%%%%%%%%%%%%%%%%%% PASSIVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dm_trial = 'PASSIVE';



%%% subject/trial identifier
subject_id = horzcat('subject ', dm_subjectno, ' ', dm_side, ' ', dm_timepoint, ' ', dm_trial);
disp(sprintf(horzcat('----------------', subject_id, '------------------')))



%%% Calculate individual Norm offset for DIRECTION
% trial 1
[convert_ind_direction_b1_volt] = calculate_direction_constants(horzcat('data_convert\', dm_CPM1_NX));
% trial 2
[convert_ind_direction_b2_volt] = calculate_direction_constants(horzcat('data_convert\', dm_CPM2_NX));
% average
convert_norm_ind_passive(7) = mean([convert_ind_direction_b1_volt convert_ind_direction_b2_volt]);



%%% Calculate individual Norm conversion factors (y = ax + b) for ANGLE
% Reads raw data from Norm, filter
% Produces individual conversion factors for angle, prints report

% trial 1
[convert_ind_angle_a1, convert_ind_angle_b1] = calculate_angle_constants(angle_cutoff, horzcat('data_convert\', dm_CPM1_NX), dm_side); % side = TMP? 
% trial 2
[convert_ind_angle_a2, convert_ind_angle_b2] = calculate_angle_constants(angle_cutoff, horzcat('data_convert\', dm_CPM2_NX), dm_side); % TMP?
% average
convert_norm_ind_passive(1) = mean([convert_ind_angle_a1 convert_ind_angle_a2]);
convert_norm_ind_passive(2) = mean([convert_ind_angle_b1 convert_ind_angle_b2]);



%%% Calculate individual Norm conversion factors (y = ax + b) for VELOCITY
% Reads raw data from Norm, filter
% Produces individual conversion factors for velocity, prints report

% trial 1
[convert_ind_velocity_a1, convert_ind_velocity_b1] = calculate_velocity_constants(velocity_cutoff, horzcat('data_convert\', dm_CPM1_NX));
% trial 2
[convert_ind_velocity_a2, convert_ind_velocity_b2] = calculate_velocity_constants(velocity_cutoff, horzcat('data_convert\', dm_CPM2_NX));
% average
convert_norm_ind_passive(5) = mean([convert_ind_velocity_a1 convert_ind_velocity_a2]);
convert_norm_ind_passive(6) = mean([convert_ind_velocity_b1 convert_ind_velocity_b2]);



%%% Calculate individual Norm conversion factors (y = ax + b) for TORQUE
% Uses constants from direction, angle, velocity
% Uses default A constant from Norm, computes B constant individually based on default A constant
% Produces individual conversion factors for torque, and prints report

% torque, A constant
convert_norm_ind_passive(3) = norm_volt_per_nm_a;

% make an array of the three µV values for B constants from direction 7, angle 2, velocity 6
convert_ind_torque_b_volt = [convert_norm_ind_passive(7) convert_norm_ind_passive(2)/convert_norm_ind_passive(1) convert_norm_ind_passive(6)/convert_norm_ind_passive(5)];
% create new B constant by combining offset from Noraxon (in mv * A constant) with default offset from Norm system (in mv * A constant)
convert_norm_ind_passive(4) = (mean(convert_ind_torque_b_volt) * norm_volt_per_nm_a) + norm_volt_per_nm_b;
% TODO MMM double check that this formula produces the correct conversion calculation

% output torque conversion numbers to screen, as text
report = sprintf(horzcat('INDI Torque conversion factors: a = ', num2str(convert_norm_ind_passive(3)), ', b = ', num2str(convert_norm_ind_passive(4)), '. Offset in millivolt = ', num2str(mean(convert_ind_torque_b_volt)), ' µV.' ));
disp(report)



%%% Test conversions on passive data:
%%% Read stiffness trial Noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
%%% Produce a new noraxon data array
convert_norm = convert_norm_ind_passive;
inputfile = horzcat('data_convert/', dm_ROM1_NX);
noraxon_prepped = read_noraxon_passive(inputfile, freq_default, dm_side, '?trialname?'); %output freq, side, trialname
% code for extracting variables to compare to data directly from Norm
% comp1_maxROM = min(noraxon_prepped(:,column_norm_angle))
% newROM = comp1_maxROM + 2;
% newROM_loc = find(noraxon_prepped(:,column_norm_angle)<newROM,1,'last');
% comp2_newROM_torque = noraxon_prepped(newROM_loc,column_norm_torque)
% zeroROM_loc = newROM_loc - 1 + find(noraxon_prepped(newROM_loc:end,column_norm_angle)<0,1,'last');
% comp3_zeroROM_torque = noraxon_prepped(zeroROM_loc,column_norm_torque)
report = sprintf(horzcat('Passive torque comparison: At ', num2str(maxROM), ' + 2 deg = ', num2str(newROM_torque), ' Nm. At zero deg = ', num2str(zeroROM_torque), ' Nm.' ));
disp(report)


inputfile = horzcat('data_convert/', dm_ROM6_NX);
noraxon_prepped = read_noraxon_passive(inputfile, freq_default, dm_side, '?trialname?'); %output freq, side, trialname
% code for extracting variables to compare to data directly from Norm
% comp4_maxROM = min(noraxon_prepped(:,column_norm_angle))
% newROM = comp4_maxROM + 2;
% newROM_loc = find(noraxon_prepped(:,column_norm_angle)<newROM,1,'last');
% comp5_newROM_torque = noraxon_prepped(newROM_loc,column_norm_torque)
% zeroROM_loc = newROM_loc - 1 + find(noraxon_prepped(newROM_loc:end,column_norm_angle)<0,1,'last');
% comp6_zeroROM_torque = noraxon_prepped(zeroROM_loc,column_norm_torque)
report = sprintf(horzcat('Passive torque comparison: At ', num2str(maxROM), ' + 2 deg = ', num2str(newROM_torque), ' Nm. At zero deg = ', num2str(zeroROM_torque), ' Nm.' ));
disp(report)






%%%%%%%%%%%%%%%%%%%%%%%%%%%% ACTIVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constants to be extracted from dm_isokin45_NX
dm_trial = 'ACTIVE';



%%% subject/trial identifier
subject_id = horzcat('subject ', dm_subjectno, ' ', dm_side, ' ', dm_timepoint, ' ', dm_trial);
disp(sprintf(horzcat('----------------', subject_id, '------------------')))



%%% Calculate individual Norm offset for DIRECTION
[convert_ind_direction_b_volt] = calculate_direction_constants(horzcat('data_convert\', dm_isokin45_NX));
convert_norm_ind_active(7) = convert_ind_direction_b_volt;



%%% Calculate individual Norm conversion factors (y = ax + b) for ANGLE
% Reads raw data from Norm, filter
% Produces individual conversion factors for angle, prints report
[convert_ind_angle_a, convert_ind_angle_b] = calculate_angle_constants_active(angle_cutoff_active, horzcat('data_convert\', dm_isokin45_NX));
convert_norm_ind_active(1) = convert_ind_angle_a;
convert_norm_ind_active(2) = convert_ind_angle_b;



%%% Calculate individual Norm conversion factors (y = ax + b) for VELOCITY
% Reads raw data from Norm, filter
% Produces individual conversion factors for velocity, prints report
[convert_ind_velocity_a, convert_ind_velocity_b] = calculate_velocity_constants_active(velocity_cutoff_active, horzcat('data_convert\', dm_isokin45_NX));
convert_norm_ind_active(5) = convert_ind_velocity_a;
convert_norm_ind_active(6) = convert_ind_velocity_b;



%%% Calculate individual Norm conversion factors (y = ax + b) for TORQUE
% Uses constants from direction, angle, velocity
% Uses default A constant from Norm, computes B constant individually based on default A constant
% Produces individual conversion factors for torque, and prints report

% torque, A constant
convert_norm_ind_active(3) = norm_volt_per_nm_a;

% make an array of the three µV values for B constants from direction 7, angle 2, velocity 6
convert_ind_torque_b_volt = [convert_norm_ind_active(7) convert_norm_ind_active(2)/convert_norm_ind_active(1) convert_norm_ind_active(6)/convert_norm_ind_active(5)];
% create new B constant by combining offset from Noraxon (in mv * A constant) with default offset from Norm system (in mv * A constant)
convert_norm_ind_active(4) = (mean(convert_ind_torque_b_volt) * norm_volt_per_nm_a) + norm_volt_per_nm_b;
% TODO MMM double check that this formula produces the correct conversion calculation

% output torque conversion numbers to screen, as text
report = sprintf(horzcat('INDI Torque conversion factors: a = ', num2str(convert_norm_ind_active(3)), ', b = ', num2str(convert_norm_ind_active(4)), '. Offset in millivolt = ', num2str(mean(convert_ind_torque_b_volt)), ' µV.' ));
disp(report)



%%% Test conversions on active data:
%%% Read stiffness trial Noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
%%% Produce a new noraxon data array
convert_norm = convert_norm_ind_active;
inputfile = horzcat('data_convert/', dm_isokin30_NX);
noraxon_prepped = read_noraxon_file(inputfile, freq_default, dm_side, 'test isokin30'); %output freq, side, trialname
% code for extracting variables to compare to data directly from Norm
%comp7_isokin_PF_max = max(noraxon_prepped(:,column_norm_torque))
%peaktorque_loc = find(noraxon_prepped(:,column_norm_torque)==comp7_isokin_PF_max,1,'last');
%comp71_angle_of_pt = noraxon_prepped(peaktorque_loc,column_norm_angle)

inputfile = horzcat('data_convert/', dm_isokinDF_NX);
noraxon_prepped = read_noraxon_file(inputfile, freq_default, dm_side, 'test isokinDF'); %output freq, side, trialname
% code for extracting variables to compare to data directly from Norm
%comp8_isokin_DF_max = min(noraxon_prepped(:,column_norm_torque))
%peaktorque_loc = find(noraxon_prepped(:,column_norm_torque)==comp8_isokin_DF_max,1,'last');
%comp81_angle_of_pt = noraxon_prepped(peaktorque_loc,column_norm_angle)


%inputfile = horzcat('data_convert/', dm_isokin90_NX);
%noraxon_prepped = read_noraxon_file(inputfile, freq_default, dm_side, 'test isokin90'); %output freq, side, trialname

inputfile = horzcat('data_convert/', dm_isomP10_1_NX);
noraxon_prepped = read_noraxon_file(inputfile, freq_default, dm_side, 'test isomP10-1'); %output freq, side, trialname
% code for extracting variables to compare to data directly from Norm
%comp9_isomP10_1_max = max(noraxon_prepped(:,column_norm_torque))

inputfile = horzcat('data_convert/', dm_isomP10_2_NX);
noraxon_prepped = read_noraxon_file(inputfile, freq_default, dm_side, 'test isomP10-2'); %output freq, side, trialname
% code for extracting variables to compare to data directly from Norm
%comp91_isomP10_2_max = max(noraxon_prepped(:,column_norm_torque))

% add other testing files
%     dm_isomD10_1_NX
%     dm_isomD10_2_NX














