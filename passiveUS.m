%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main file for analysis of passive dorsiflexion with US
% Marie Moltubakk 4.2.2015
% 
% input argument 1 = project selection (1 = BD, 2 = intervent)
% input argument 2 = plot selection (0 = none, 1 = group plots, 2 = ind plots)
% input argument 3 = resume run (0 = start from beginning, 1 = start inloop, 2 = start after loop)
% 
% 2017-12-28: Moving normalization to separate output columns + removing
%    opportunity to analyse BD data from this script. Go back to
%    passiveUS.with_norm.old for reviewing BD data.
% 
% 2018-02-13: Changing "ind_max" angles and forces to be leg independent,
%    i.e. left and right leg may have different individual angles common to
%    PRE and POST.
% 
% The scripts assume a certain structure for input files from Tracker and
% Noraxon. I.e. number of and order of EMG and other channels. If these
% scripts are to be used for other projects, some of the code must be
% modified. These lines are marked % PROJECTSPECIFIC
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [] = passiveUS(~, input_plot, input_resumerun) %  = input project

    %% Startup actions
    close all
    global mute
    mute = 1; % do not output to screen: conversion factors, EMG%, mom arm
    line = 0; % preallocate for resumerun
    AAA_today = date;
    cprintf('*black', horzcat('---------------- PASSIVEUS ', datestr(AAA_today), ' ------------------\n'))
              
    %% Plots: determine which plots to display
    global plot_achilles plot_norm plot_emg plot_check plot_us subject_id plot_licht plot_individual plot_conversion

    plot_check = 0;
    plot_individual = 0;
    plot_norm = 0;
    plot_old = 0; %old plots, data not used for BD study

    if input_plot >= 1 
        plot_check = 1; % LEVEL 1: group plots
    end
    if input_plot >= 2
        plot_individual = 1; % LEVEL 1B: individual plots
    end
    if input_plot >= 3
        plot_norm = 1; % LEVEL 2: main checkpoint plots
    end

    % LEVEL 3:
    plot_conversion = 0; % turn on/off plots for data conversion Norm
    plot_us = 0; % tracked feature vs external marker 
    plot_emg = 0; % RMS 3 EMG channels per trial
    plot_achilles = 0; % turn on/off Achilles machine plots
    plot_licht = 0; % plot averaging of trials from Lichtwark US fascicle tracking


    %% Constants and globals: sampling, conversion, filters, array index % PROJECTSPECIFIC

    % sampling frequencies
    global us_zerodispframes noraxonfreq freq_default
    us_zerodispframes = 1; % No of US frames to average as zero displacement
    noraxonfreq = 1500; % sampling frequency of noraxon data
    freq_default = 100; % output frequency for noraxon data without US video (MVC, etc)
    angle_step_stats_abs = 0.5; % every x degrees
    %angle_step_plots = 0.05; % resampled, averaged data extracted every x degrees for PLOTS
    %angle_step_stats_norm = 10;  % every x percent

    % variables for NORM conversion factors - some from manifacturers, some from calibration trials, some calculated from actual data
    global convert_norm_angle_a convert_norm_angle_b convert_norm_torque_a convert_norm_torque_b convert_norm_velocity_a convert_norm_velocity_b convert_norm_direction_b
    global convert_achilles norm_volt_per_degree norm_volt_per_velocity norm_volt_per_nm_a norm_volt_per_nm_b
    % default conversion factors from Norm manuals and Achilles calibration:
    convert_achilles = -81.9909;  % conversion factor ACHILLES torque, V -> Nm
    norm_volt_per_degree = (2048*((1*16)+0)/1024)*(10000/32768); %   (Gain * (Sampled Position + Offset) / 1024) * (10v/32768) = Volt-value, Sampled Position is in units of 1/16 degree.
    norm_volt_per_velocity = (1024*((1*16)+0)/1024)*(10000/32768); %   (Gain * (Sampled Velocity + Offset) / 1024) * (10v/32768) = Volt-value. Sampled Velocity is in units of 1/16 degree.
    % norm_volt_per_nm = (1024*((1*(1.355818*32768/500))+0)/1024)*(10000/32768); % 27.116 µV/Nm   Sampled Torque is in units of Foot-Pounds * 32768 / 500
    norm_volt_per_nm_a = 0.07490071; % from file "M M M conversion volt-torque 2014-DES NEW DATA.xlsx"
    norm_volt_per_nm_b = 0.69283161;

    % cutoff frequencies for filters
    global emg_bandpass emg_rms_ms mvc_window_ms 
    emg_bandpass = [10/(noraxonfreq/2) 500/(noraxonfreq/2)]; % cutoff frequencies for EMG butterworth bandpass filter
    emg_rms_ms = 100; % milliseconds RMS window, must be divisible by 4 - ref Basmajian 1985 = 50 ms, ref Aagaard paper Passive tensile stress = 200 ms
    mvc_window_ms = 500; % milliseconds window for determining MVC torque and EMG
    global angle_cutoff velocity_cutoff torque_cutoff_bandstop torque_cutoff_active angle_cutoff_active velocity_cutoff_active
    angle_cutoff = 10/(noraxonfreq/2); % cutoff freq, Norm angle PASSIVE - ref Winter 1990 = 15 hz. Kongsgaard = 8hz
    angle_cutoff_active = 20/(noraxonfreq/2);
    velocity_cutoff = 20/(noraxonfreq/2);
    velocity_cutoff_active = 20/(noraxonfreq/2);
    torque_cutoff_bandstop = [0.36/(noraxonfreq/2) 0.42/(noraxonfreq/2)]; % bandstop freq to eliminate noise from Norm engine
    torque_cutoff_active = 10/(noraxonfreq/2); % cutoff freq, Norm filtering - ref Winter 1990 = 15 hz. Kongsgaard = 8hz

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


    %% Constants: axes ++ for plots
    col_lightblue = [0.6 0.8 1];
    col_lightred = [1 0.6 0.8];
    col_orange = [1 0.75 0];
    col_grey = [0.3 0.3 0.3];

    txt_elong = 'Elongation (mm)';
    txt_length = 'Length (mm)';
    txt_strain = 'Strain (% of initial length)';
    txt_displ = 'Displacement (mm)';
    txt_emg = 'EMG (% of MVC)';
    txt_gonio = 'Gonio angle (°)';
    txt_elong_norm = 'Elongation (% of MTU length)';
    txt_length_norm = 'Length (% of MTU length)';
    txt_elong_norm_perc = 'Elongation (% of MTU ELONG)';
    
    angle_max = 36;
    
    axis_EMG_wide = [-1 angle_max 0 30];
    axis_EMG_narrow = [-1 angle_max 0 16];

    axis_displ_GMFAS = [-1 angle_max -2 4];

    axis_el_MTU = [-1 angle_max 20 50];
    axis_str_MTU = [-1 angle_max -0.5 8];
    axis_len_MTU = [-25 angle_max 440 540];

    axis_penn_GMFAS = [-30 angle_max 7 18]; % fasicle
    axis_len_GMFAS = [-30 angle_max 40 80];
    axis_el_GMFAS = [-1 angle_max -10 25];
    axis_str_GMFAS = [-1 angle_max -10 40];

    axis_el_SEE_arch = [-1 angle_max 5 40];
    axis_str_SEE_arch = [-1 angle_max 2 9];
    axis_len_SEE_arch = [-25 angle_max 380 480];
    
    axis_el_percent = [-1 angle_max -5 100];
    
    axis_el_GMmsc_arch = [-1 angle_max 0 30];
    axis_str_GMmsc_arch = [-1 angle_max 0 55];
    axis_len_GMmsc_arch = [-1 angle_max 50 80];

    axis_displ_SOL = [-1 angle_max -3 11];

    axis_el_GMmsc = [-1 angle_max -5 55];
    axis_str_GMmsc = [-1 angle_max -1 22];
    axis_len_GMmsc = [-1 angle_max 250 350];

    axis_el_SOL = [-1 angle_max 10 45];
    axis_str_SOL = [-1 angle_max 3 15];
    axis_len_SOL = [-1 angle_max 230 310];
    
    axis_el_GMapo = [-1 angle_max -10 10];
    axis_str_GMapo = [-1 angle_max -10 10];
    axis_len_GMapo = [-1 angle_max 90 160];

    axis_el_AT = [-1 angle_max -15 20];
    axis_str_AT = [-1 angle_max -20 40];
    axis_len_AT = [-1 angle_max 20 90];

    axis_el_GMtend = [-1 angle_max -18 15];
    axis_str_GMtend = [-1 angle_max -5 10];
    axis_len_GMtend = [-1 angle_max 150 220];

    axis_torque = [-1 angle_max 0 80];

    axis_force = [-1 angle_max 0 1700];
    % axis_PP = [-5 100 0 105];

    axis_ind_elong = [-1 angle_max -8 55];
    axis_ind_strain = [-1 angle_max -6 50];
    axis_ind_len_MTU = [-1 angle_max 0 550];


    %% READ max angles and forces from "create_angles_passive.m"
    % Extracts max angles and forces per trial/subject/common
    % Produces arrays with angles and forces per subject/side/trial/group

    global input_gon_pre_r input_gon_pre_l input_gon_post_r input_gon_post_l input_gon_ind_max input_gon_common_max % gon = goniometer ankle angle
    global input_gon_ind_rmax input_gon_ind_lmax
    global input_for_pre_r input_for_pre_l input_for_post_r input_for_post_l input_for_ind_max input_for_common_max 
    global input_for_ind_rmax input_for_ind_lmax
    dm_filename = 'angles_output.tsv';
    if exist(dm_filename, 'file') == 2
        read_angles_passive(dm_filename);
    else
        cprintf('*red', 'ERROR: Subject ANGLE data file does not exist. Generate by running "create_angles_passive".\n')
        return
    end
    dm_filename = 'forces_output.tsv';
    if exist(dm_filename, 'file') == 2
        read_forces_passive(dm_filename);
    else
        cprintf('*red', 'ERROR: Subject FORCE data file does not exist. Generate by running "create_angles_passive".\n')
        return
    end


    %% READ datamaster file, to connect corresponding data files
    % Produces arrays with file names and variables per trial

    global dm_subjectno dm_timepoint dm_side dm_trial 
    global dm_ROM_sol1_NX dm_ROM_sol1_US dm_ROM_sol1_US_frame dm_ROM_sol2_NX dm_ROM_sol2_US dm_ROM_sol2_US_frame
    global dm_ROM_gmmtj1_NX dm_ROM_gmmtj1_US dm_ROM_gmmtj1_US_frame dm_ROM_gmmtj2_NX dm_ROM_gmmtj2_US dm_ROM_gmmtj2_US_frame
    global dm_ROM_gmfas1_NX dm_ROM_gmfas1_US dm_ROM_gmfas1_US_frame dm_ROM_gmfas2_NX dm_ROM_gmfas2_US dm_ROM_gmfas2_US_frame  dm_ROM_gmfas1_licht dm_ROM_gmfas2_licht
    global dm_MVC_PF %dm_MVC_DF %dm_CPM_calc_NX dm_CPM_calc_US dm_CPM_calc_US_frame dm_CPM_sol_NX dm_CPM_sol_US dm_CPM_sol_US_frame 
    global dm_leg_length dm_at_SOL_length dm_at_GM_length dm_GMmsc_penn dm_GMmsc_faslen dm_ankle_angle_rest
    global dm_at_SOL_zero dm_at_GM_zero
    global at_momentarm 
    global filepath
    dm_filename = 'data/datamaster_passive.tsv';
    
    if input_resumerun == 1
        % resume running of loop
        load all_data_passive_inloop
        input_resumerun = 1; % overwrite loaded variable
        line_start = line+1; % all_data_stiff_inloop was saved at the "end" of a line (trial) - resume with next line
        linestotal = read_datamaster_passive(dm_filename); % allow edits in datamaster (filenames may be edited, line order prior to restart may NOT!)
    elseif input_resumerun == 2
        % do not run loop at all, no need for line_start or linestotal
        load all_data_passive_endloop
        input_resumerun = 2; % overwrite loaded variable
    else
        % == 0, start from beginning
        line_start = 1;
        linestotal = read_datamaster_passive(dm_filename);
    end
    
            
    %% Preallocate output arrays 
    % preallocate only if starting from beginning:
    if input_resumerun == 0
        % common arrays for all subjects:
            all_passive_output_head = {...
'Subject', 'Time', 'Side', 'Trial', ...
'ROM trial (deg)', 'ROM subject (L_R_PRE_POST)', 'ROM common (all subjects)', 'ROM submax - 10deg', 'ROM submax - 67perc ROM', ...
'force at trial ROM (N)', 'force at subject max force', 'force at subject max ROM', 'force at common max ROM', 'force at 0 deg', 'force at 10deg', 'force at 67perc ROM', ...
'torque at trial ROM (N)', 'torque at subject max torque', 'torque at subject max ROM', 'torque at common max ROM', 'torque at 0 deg', 'torque at 10deg', 'torque at 67perc ROM', ...
'angle at trial max force (deg)', 'angle at subject max force', 'angle at common max force', 'angle at ind R max force', 'angle at ind L max force', ...            
'Stiffness at trial ROM (Nmperdeg)', 'Stiffness at subject max ROM', 'Stiffness at common max ROM', 'Stiffness at 15 deg', 'Stiffness at 10deg', 'Stiffness at 67perc ROM',...
'Stiffness index (Nmper_deg_sq)', ... %', 'length,', 'elong,', 'strain
'out_length_AT_trial_max', 'out_length_AT_ind_max', 'out_length_AT_common_max', 'out_length_AT_submax_1', 'out_length_AT_submax_2', 'out_length_AT_max',...
'out_length_GMtend_trial_max', 'out_length_GMtend_ind_max', 'out_length_GMtend_common_max', 'out_length_GMtend_submax_1', 'out_length_GMtend_submax_2', 'out_length_GMtend_max',...
'out_length_GMapo_trial_max', 'out_length_GMapo_ind_max', 'out_length_GMapo_common_max', 'out_length_GMapo_submax_1', 'out_length_GMapo_submax_2', 'out_length_GMapo_max',...
'out_length_msc_GM_trial_max', 'out_length_msc_GM_ind_max', 'out_length_msc_GM_common_max', 'out_length_msc_GM_submax_1', 'out_length_msc_GM_submax_2', 'out_length_msc_GM_max',...
'out_length_msc_SOL_trial_max', 'out_length_msc_SOL_ind_max', 'out_length_msc_SOL_common_max', 'out_length_msc_SOL_submax_1', 'out_length_msc_SOL_submax_2', 'out_length_msc_SOL_max',...
'out_length_leg_trial_max', 'out_length_leg_ind_max', 'out_length_leg_common_max', 'out_length_leg_submax_1', 'out_length_leg_submax_2', 'out_length_leg_max',...
'out_length_SEE_Fuku_trial_max', 'out_length_SEE_Fuku_ind_max', 'out_length_SEE_Fuku_common_max', 'out_length_SEE_Fuku_submax_1', 'out_length_SEE_Fuku_submax_2', 'out_length_SEE_Fuku_max',...
'out_length_msc_GM_Fuku_trial_max', 'out_length_msc_GM_Fuku_ind_max', 'out_length_msc_GM_Fuku_common_max', 'out_length_msc_GM_Fuku_submax_1', 'out_length_msc_GM_Fuku_submax_2', 'out_length_msc_GM_Fuku_max',...
'out_elong_AT_trial_max', 'out_elong_AT_ind_max', 'out_elong_AT_common_max', 'out_elong_AT_submax_1', 'out_elong_AT_submax_2', 'out_elong_AT_max',...
'out_elong_GMtend_trial_max', 'out_elong_GMtend_ind_max', 'out_elong_GMtend_common_max', 'out_elong_GMtend_submax_1', 'out_elong_GMtend_submax_2', 'out_elong_GMtend_max',...
'out_elong_GMapo_trial_max', 'out_elong_GMapo_ind_max', 'out_elong_GMapo_common_max', 'out_elong_GMapo_submax_1', 'out_elong_GMapo_submax_2', 'out_elong_GMapo_max',...
'out_elong_msc_GM_trial_max', 'out_elong_msc_GM_ind_max', 'out_elong_msc_GM_common_max', 'out_elong_msc_GM_submax_1', 'out_elong_msc_GM_submax_2', 'out_elong_msc_GM_max',...
'out_elong_msc_SOL_trial_max', 'out_elong_msc_SOL_ind_max', 'out_elong_msc_SOL_common_max', 'out_elong_msc_SOL_submax_1', 'out_elong_msc_SOL_submax_2', 'out_elong_msc_SOL_max',...
'out_elong_leg_trial_max', 'out_elong_leg_ind_max', 'out_elong_leg_common_max', 'out_elong_leg_submax_1', 'out_elong_leg_submax_2', 'out_elong_leg_max',...
'out_elong_SEE_Fuku_trial_max', 'out_elong_SEE_Fuku_ind_max', 'out_elong_SEE_Fuku_common_max', 'out_elong_SEE_Fuku_submax_1', 'out_elong_SEE_Fuku_submax_2', 'out_elong_SEE_Fuku_max',...
'out_elong_msc_GM_Fuku_trial_max', 'out_elong_msc_GM_Fuku_ind_max', 'out_elong_msc_GM_Fuku_common_max', 'out_elong_msc_GM_Fuku_submax_1', 'out_elong_msc_GM_Fuku_submax_2', 'out_elong_msc_GM_Fuku_max',...
'out_strain_AT_trial_max', 'out_strain_AT_ind_max', 'out_strain_AT_common_max', 'out_strain_AT_submax_1', 'out_strain_AT_submax_2', 'out_strain_AT_max',...
'out_strain_GMtend_trial_max', 'out_strain_GMtend_ind_max', 'out_strain_GMtend_common_max', 'out_strain_GMtend_submax_1', 'out_strain_GMtend_submax_2', 'out_strain_GMtend_max',...
'out_strain_GMapo_trial_max', 'out_strain_GMapo_ind_max', 'out_strain_GMapo_common_max', 'out_strain_GMapo_submax_1', 'out_strain_GMapo_submax_2', 'out_strain_GMapo_max',...
'out_strain_msc_GM_trial_max', 'out_strain_msc_GM_ind_max', 'out_strain_msc_GM_common_max', 'out_strain_msc_GM_submax_1', 'out_strain_msc_GM_submax_2', 'out_strain_msc_GM_max',...
'out_strain_msc_SOL_trial_max', 'out_strain_msc_SOL_ind_max', 'out_strain_msc_SOL_common_max', 'out_strain_msc_SOL_submax_1', 'out_strain_msc_SOL_submax_2', 'out_strain_msc_SOL_max',...
'out_strain_leg_trial_max', 'out_strain_leg_ind_max', 'out_strain_leg_common_max', 'out_strain_leg_submax_1', 'out_strain_leg_submax_2', 'out_strain_leg_max',...
'out_strain_SEE_Fuku_trial_max', 'out_strain_SEE_Fuku_ind_max', 'out_strain_SEE_Fuku_common_max', 'out_strain_SEE_Fuku_submax_1', 'out_strain_SEE_Fuku_submax_2', 'out_strain_SEE_Fuku_max',...
'out_strain_msc_GM_Fuku_trial_max', 'out_strain_msc_GM_Fuku_ind_max', 'out_strain_msc_GM_Fuku_common_max', 'out_strain_msc_GM_Fuku_submax_1', 'out_strain_msc_GM_Fuku_submax_2', 'out_strain_msc_GM_Fuku_max',... %', 'licht', 'GM', 'and', 'SOL
'out_norm_length_leg_trial_max', 'out_norm_length_leg_ind_max', 'out_norm_length_leg_common_max', 'out_norm_length_leg_submax_1', 'out_norm_length_leg_submax_2', 'out_norm_length_leg_max',...
'out_norm_length_SEE_Fuku_trial_max', 'out_norm_length_SEE_Fuku_ind_max', 'out_norm_length_SEE_Fuku_common_max', 'out_norm_length_SEE_Fuku_submax_1', 'out_norm_length_SEE_Fuku_submax_2', 'out_norm_length_SEE_Fuku_max',...
'out_norm_length_msc_GM_Fuku_trial_max', 'out_norm_length_msc_GM_Fuku_ind_max', 'out_norm_length_msc_GM_Fuku_common_max', 'out_norm_length_msc_GM_Fuku_submax_1', 'out_norm_length_msc_GM_Fuku_submax_2', 'out_norm_length_msc_GM_Fuku_max',...
'out_norm_elong_leg_trial_max', 'out_norm_elong_leg_ind_max', 'out_norm_elong_leg_common_max', 'out_norm_elong_leg_submax_1', 'out_norm_elong_leg_submax_2', 'out_norm_elong_leg_max',...
'out_norm_elong_SEE_Fuku_trial_max', 'out_norm_elong_SEE_Fuku_ind_max', 'out_norm_elong_SEE_Fuku_common_max', 'out_norm_elong_SEE_Fuku_submax_1', 'out_norm_elong_SEE_Fuku_submax_2', 'out_norm_elong_SEE_Fuku_max',...
'out_norm_elong_msc_GM_Fuku_trial_max', 'out_norm_elong_msc_GM_Fuku_ind_max', 'out_norm_elong_msc_GM_Fuku_common_max', 'out_norm_elong_msc_GM_Fuku_submax_1', 'out_norm_elong_msc_GM_Fuku_submax_2', 'out_norm_elong_msc_GM_Fuku_max',...
'out_norm_elong_percent_SEE_Fuku_trial_max', 'out_norm_elong_percent_SEE_Fuku_ind_max', 'out_norm_elong_percent_SEE_Fuku_common_max', 'out_norm_elong_percent_SEE_Fuku_submax_1', 'out_norm_elong_percent_SEE_Fuku_submax_2', 'out_norm_elong_percent_SEE_Fuku_max',...
'out_norm_elong_percent_msc_GM_Fuku_trial_max', 'out_norm_elong_percent_msc_GM_Fuku_ind_max', 'out_norm_elong_percent_msc_GM_Fuku_common_max', 'out_norm_elong_percent_msc_GM_Fuku_submax_1', 'out_norm_elong_percent_msc_GM_Fuku_submax_2', 'out_norm_elong_percent_msc_GM_Fuku_max',... %', 'normalized', 'GM', 'fascicle', 'length', 'and', 'elong
'out_licht_pennation_GM_trial_max', 'out_licht_pennation_GM_ind_max', 'out_licht_pennation_GM_common_max', 'out_licht_pennation_GM_submax_1', 'out_licht_pennation_GM_submax_2', 'out_licht_pennation_GM_max', 'out_licht_pennation_GM_zero',...
'out_licht_fas_length_GM_trial_max', 'out_licht_fas_length_GM_ind_max', 'out_licht_fas_length_GM_common_max', 'out_licht_fas_length_GM_submax_1', 'out_licht_fas_length_GM_submax_2', 'out_licht_fas_length_GM_max', 'out_licht_fas_length_GM_zero',...
'out_licht_fas_elong_GM_trial_max', 'out_licht_fas_elong_GM_ind_max', 'out_licht_fas_elong_GM_common_max', 'out_licht_fas_elong_GM_submax_1', 'out_licht_fas_elong_GM_submax_2', 'out_licht_fas_elong_GM_max', 'out_licht_fas_elong_GM_zero',...
'out_licht_fas_strain_GM_trial_max', 'out_licht_fas_strain_GM_ind_max', 'out_licht_fas_strain_GM_common_max', 'out_licht_fas_strain_GM_submax_1', 'out_licht_fas_strain_GM_submax_2', 'out_licht_fas_strain_GM_max', 'out_licht_fas_strain_GM_zero',...
'out_licht_pennation_SOL_trial_max', 'out_licht_pennation_SOL_ind_max', 'out_licht_pennation_SOL_common_max', 'out_licht_pennation_SOL_submax_1', 'out_licht_pennation_SOL_submax_2', 'out_licht_pennation_SOL_max', 'out_licht_pennation_SOL_zero',... %', 'normalized', 'SEE', '&', 'muscle
'out_licht_fas_length_SOL_trial_max', 'out_licht_fas_length_SOL_ind_max', 'out_licht_fas_length_SOL_common_max', 'out_licht_fas_length_SOL_submax_1', 'out_licht_fas_length_SOL_submax_2', 'out_licht_fas_length_SOL_max', 'out_licht_fas_length_SOL_zero',...
'out_norm_licht_fas_length_GM_trial_max', 'out_norm_licht_fas_length_GM_ind_max', 'out_norm_licht_fas_length_GM_common_max', 'out_norm_licht_fas_length_GM_submax_1', 'out_norm_licht_fas_length_GM_submax_2', 'out_norm_licht_fas_length_GM_max', 'out_norm_licht_fas_length_GM_zero',...
'out_norm_licht_fas_elong_GM_trial_max', 'out_norm_licht_fas_elong_GM_ind_max', 'out_norm_licht_fas_elong_GM_common_max', 'out_norm_licht_fas_elong_GM_submax_1', 'out_norm_licht_fas_elong_GM_max', 'out_norm_licht_fas_elong_GM_submax_2', 'out_norm_licht_fas_elong_GM_zero',... %', 'misc
'out_emg_gm_trial_max', 'out_emg_gm_ind_max', 'out_emg_gm_common_max', 'out_emg_gm_submax_1', 'out_emg_gm_submax_2', 'out_emg_gm_max', 'out_emg_gm_zero',...
'out_emg_gl_trial_max', 'out_emg_gl_ind_max', 'out_emg_gl_common_max', 'out_emg_gl_submax_1', 'out_emg_gl_submax_2', 'out_emg_gl_max', 'out_emg_gl_zero',...
'out_emg_sol_trial_max', 'out_emg_sol_ind_max', 'out_emg_sol_common_max', 'out_emg_sol_submax_1', 'out_emg_sol_submax_2', 'out_emg_sol_max', 'out_emg_sol_zero',...
'out_emg_mean_trial_max', 'out_emg_mean_ind_max', 'out_emg_mean_common_max', 'out_emg_mean_submax_1', 'out_emg_mean_submax_2', 'out_emg_mean_max', 'out_emg_mean_zero',...
                }; % PROJECTSPECIFIC
        all_passive_output = zeros(ceil(linestotal),length(all_passive_output_head)-4); 
        all_passive_output_txt = cell(ceil(linestotal),4);

        STR_PRE_count = 0;
        STR_POST_count = 0;
        CON_PRE_count = 0;
        CON_POST_count = 0;

        STR_PRE_ID{ceil(linestotal)} = [];
        STR_POST_ID{ceil(linestotal)} = [];
        CON_PRE_ID{ceil(linestotal)} = [];
        CON_POST_ID{ceil(linestotal)} = [];

        STR_PRE_no(ceil(linestotal)) = zeros;
        STR_PRE_prone(ceil(linestotal),12) = zeros;
        STR_PRE_angle_vars{ceil(linestotal)} = zeros;
        STR_PRE_angle_vars_mean{ceil(linestotal)} = zeros;

        STR_POST_no(ceil(linestotal)) = zeros;
        STR_POST_prone(ceil(linestotal),12) = zeros;
        STR_POST_angle_vars{ceil(linestotal)} = zeros;
        STR_POST_angle_vars_mean{ceil(linestotal)} = zeros;

        CON_PRE_no(ceil(linestotal)) = zeros;
        CON_PRE_prone(ceil(linestotal),12) = zeros;
        CON_PRE_angle_vars{ceil(linestotal)} = zeros;
        CON_PRE_angle_vars_mean{ceil(linestotal)} = zeros;

        CON_POST_no(ceil(linestotal)) = zeros;
        CON_POST_prone(ceil(linestotal),12) = zeros;
        CON_POST_angle_vars{ceil(linestotal)} = zeros;
        CON_POST_angle_vars_mean{ceil(linestotal)} = zeros;
    end
    
    
    %% Prepare for LOOP (if not "resume after loop")
    if input_resumerun < 2
        % 0 = run loop from beginning
        % OR:
        % 1 = continue running loop from "inloop"

              
        %% LOOP through all lines in datamaster file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for line = line_start:linestotal 
            
            
            %% Prepare EMG variables
            clear emg_all;
            emg_all{3,6} = zeros; % 6 EMG channels for 1 subject, to be reused across subjects


            %% Subject/trial identifier
            trial_subjectno = str2double(dm_subjectno{line});
            trial_timepoint = strcmp(dm_timepoint{line},'POST'); % 0 = PRE, 1 = POST
            trial_leg = strcmp(dm_trial{line},'STR'); % 0 = CON, 1 = STR
            trial_calf_length = str2double(dm_leg_length{line}) * 10;  % convert from input in cm to calf length in mm
            prone_GMfas_length = str2double(dm_GMmsc_faslen{line}); % resting (prone) length of GM fascicle 
            prone_GMfas_penn = str2double(dm_GMmsc_penn{line}); % resting (prone) pennation angle of GM 
            prone_GMfas_ankle = -str2double(dm_ankle_angle_rest{line}); % plantarflexed angle - input is positive value, but must be neg for plots

            filepath = 'data\';
            subject_id = horzcat('INT_', dm_subjectno{line}, '_', dm_trial{line}, '_', dm_timepoint{line}, '_', dm_side{line});
            if trial_timepoint == 0 && trial_leg == 1 % PRE, STR
                STR_PRE_count = STR_PRE_count + 1;
                STR_PRE_no(STR_PRE_count) = str2double(dm_subjectno{line});
                STR_PRE_ID{STR_PRE_count} = subject_id;
            elseif trial_timepoint == 1 && trial_leg == 1 % POST, STR
                STR_POST_count = STR_POST_count + 1;
                STR_POST_no(STR_POST_count) = str2double(dm_subjectno{line});
                STR_POST_ID{STR_POST_count} = subject_id;
            elseif trial_timepoint == 0 && trial_leg == 0 % PRE, CON
                CON_PRE_count = CON_PRE_count + 1;
                CON_PRE_no(CON_PRE_count) = str2double(dm_subjectno{line});
                CON_PRE_ID{CON_PRE_count} = subject_id;
                % old: CON_PRE_subject_ID(CON_PRE_count) = trial_subjectno;
            elseif trial_timepoint == 1 && trial_leg == 0 % POST, CON
                CON_POST_count = CON_POST_count + 1;
                CON_POST_no(CON_POST_count) = str2double(dm_subjectno{line});
                CON_POST_ID{CON_POST_count} = subject_id;
            end
            
            cprintf('*black', horzcat('----------------', subject_id, '------------------\n'))
            
            %% extract trial ROM angles (generated by create_angles_passive)
            % retrieve max gonio angle
            if strcmpi(dm_timepoint{line}, 'PRE') == 1
                if strcmpi(dm_side{line}, 'R') == 1
                    out_ROM_trial_max = str2double(input_gon_pre_r{trial_subjectno}) - 0.00000001; % tweak for computer storage of numbers, where 8.3000 is stored as 8.2999999999999999
                else % L
                    out_ROM_trial_max = str2double(input_gon_pre_l{trial_subjectno}) - 0.00000001;
                end
            else % POST
                if strcmpi(dm_side{line}, 'R') == 1
                    out_ROM_trial_max = str2double(input_gon_post_r{trial_subjectno}) - 0.00000001;
                else % L
                    out_ROM_trial_max = str2double(input_gon_post_l{trial_subjectno}) - 0.00000001;
                end
            end

            % error if ROM not existing (incorrect create angles passive?)
            if round(out_ROM_trial_max,1) == 100
                cprintf('*red', 'ERROR: ROM data not calculated for current subject. Run "create_angles_passive".\n')
                return
            end
            

            %% Calculate NORM conversion factors
            % retrieve conversion constants for Norm data
            [convert_norm_angle_a, convert_norm_angle_b, convert_norm_torque_a, convert_norm_torque_b, convert_norm_velocity_a, convert_norm_velocity_b, convert_norm_direction_b] = calculate_norm_constants_complete(dm_subjectno{line}, dm_side{line}, dm_timepoint{line}, line, 'passive');


            %% Calculate MVC from EMG data
            % Read noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample - PLANTAR FLEXION
            % Produce a new noraxon data array
            noraxon_mvc_plantar = read_noraxon_stiffness(strcat(filepath, dm_MVC_PF{line}), freq_default, dm_side{line}, 'MVC plantar');

            % Calculate co-activation constants
            % Read complete, prepared noraxon array + number of frames to average (freq * time)
            % Produce max torque, max EMG constants
            if strcmpi(dm_side{line},'R') == 1
                [~,EMG_max_gm] = calculate_EMG_max(noraxon_mvc_plantar, freq_default*(mvc_window_ms/1000), column_r_gm, 0);
                [~,EMG_max_gl] = calculate_EMG_max(noraxon_mvc_plantar, freq_default*(mvc_window_ms/1000), column_r_gl, 0);
                [~,EMG_max_sol] = calculate_EMG_max(noraxon_mvc_plantar, freq_default*(mvc_window_ms/1000), column_r_sol, 0);
            else % left
                [~,EMG_max_gm] = calculate_EMG_max(noraxon_mvc_plantar, freq_default*(mvc_window_ms/1000), column_l_gm, 0);
                [~,EMG_max_gl] = calculate_EMG_max(noraxon_mvc_plantar, freq_default*(mvc_window_ms/1000), column_l_gl, 0);
                [~,EMG_max_sol] = calculate_EMG_max(noraxon_mvc_plantar, freq_default*(mvc_window_ms/1000), column_l_sol, 0);
            end

            % no calculations for tibialis anterior for passive trials:
            EMG_max_TA = 0;


            %% Calculate ACHILLES TENDON MOMENT ARM
            at_momentarm = calculate_momentarm(0, 0, dm_leg_length{line});


            %% Calculations for 2x SOL trials
            % extract force, torque, gonio, angle, displacement for EACH TRIAL

            if(strcmpi(dm_ROM_sol1_NX{line}, 'null'))
                % allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case in final_stiffness
                SOL_force_1 = zeros;
                SOL_torque_1 = zeros;
                SOL_gonio_1 = zeros;
                SOL_angle_1 = zeros;
                SOL_displacement_1 = zeros;
            else
                [SOL_force_1, SOL_gonio_1, SOL_angle_1, SOL_displacement_1, SOL_emg_gm_1, SOL_emg_gl_1, SOL_emg_sol_1, SOL_time_1, SOL_torque_1] = extract_force_displ_singletrial_passive_EMG(dm_ROM_sol1_NX{line}, dm_ROM_sol1_US{line}, dm_ROM_sol1_US_frame{line}, EMG_max_TA, EMG_max_gm, EMG_max_gl, EMG_max_sol, dm_leg_length{line}, dm_side{line}, line, 'SOL1');
                emg_all{1,1} = [SOL_gonio_1 SOL_emg_gm_1];
                emg_all{2,1} = [SOL_gonio_1 SOL_emg_gl_1];
                emg_all{3,1} = [SOL_gonio_1 SOL_emg_sol_1];
            end
            if(strcmpi(dm_ROM_sol2_NX{line}, 'null'))
                SOL_force_2 = zeros;
                SOL_torque_2 = zeros;
                SOL_gonio_2 = zeros;
                SOL_angle_2 = zeros;
                SOL_displacement_2 = zeros;
            else 
                [SOL_force_2, SOL_gonio_2, SOL_angle_2, SOL_displacement_2, SOL_emg_gm_2, SOL_emg_gl_2, SOL_emg_sol_2, SOL_time_2, SOL_torque_2] = extract_force_displ_singletrial_passive_EMG(dm_ROM_sol2_NX{line}, dm_ROM_sol2_US{line}, dm_ROM_sol2_US_frame{line}, EMG_max_TA, EMG_max_gm, EMG_max_gl, EMG_max_sol, dm_leg_length{line}, dm_side{line}, line, 'SOL2');    
                emg_all{1,2} = [SOL_gonio_2 SOL_emg_gm_2];
                emg_all{2,2} = [SOL_gonio_2 SOL_emg_gl_2];
                emg_all{3,2} = [SOL_gonio_2 SOL_emg_sol_2];
            end

            % average two passive trials
            if(strcmpi(dm_ROM_sol1_NX{line}, 'null')) % trial 1 not existing
                data_SOL = average_passive_trials_EMG(SOL_force_2, SOL_gonio_2, SOL_angle_2, SOL_displacement_2, SOL_emg_gm_2, SOL_emg_gl_2, SOL_emg_sol_2, SOL_time_2, SOL_torque_2);
            elseif(strcmpi(dm_ROM_sol2_NX{line}, 'null'))
                data_SOL = average_passive_trials_EMG(SOL_force_1, SOL_gonio_1, SOL_angle_1, SOL_displacement_1, SOL_emg_gm_1, SOL_emg_gl_1, SOL_emg_sol_1, SOL_time_1, SOL_torque_1);
            else % if 2 trials exist (none of variables are 'null')
                data_SOL = average_passive_trials_EMG(SOL_force_1, SOL_gonio_1, SOL_angle_1, SOL_displacement_1, SOL_emg_gm_1, SOL_emg_gl_1, SOL_emg_sol_1, SOL_time_1, SOL_torque_1, SOL_force_2, SOL_gonio_2, SOL_angle_2, SOL_displacement_2, SOL_emg_gm_2, SOL_emg_gl_2, SOL_emg_sol_2, SOL_time_2, SOL_torque_2);
            end

            % plot summary of 2 trials: Displacement-angle
            if plot_check && plot_individual && plot_old 
                plottitle = horzcat('IND displacement SOL vs angle, ', subject_id);
                figure('Name',plottitle);
                plot(SOL_gonio_1,SOL_displacement_1,'LineWidth',2)
                hold on
                plot(SOL_gonio_2,SOL_displacement_2,'LineWidth',2)
                axis(axis_displ_SOL)
                ylabel(txt_displ)
                xlabel(txt_gonio)
                title(plottitle,'Interpreter', 'none')
                legend('Trial 1','Trial 2','Location','Southeast')
                print(horzcat('data_plots/',plottitle),'-dpng')
            end
         

            %% Calculations for 2x GM MTJ trials

            % extract force, gonio, angle, displacement for EACH TRIAL
            if(strcmpi(dm_ROM_gmmtj1_NX{line}, 'null'))
                % allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case in final_stiffness
                GMMTJ_force_1 = zeros;
                GMMTJ_torque_1 = zeros;
                GMMTJ_angle_1 = zeros;
                GMMTJ_gonio_1 = zeros;
                GMMTJ_displacement_1 = zeros;
            else 
                [GMMTJ_force_1, GMMTJ_gonio_1, GMMTJ_angle_1, GMMTJ_displacement_1, GMMTJ_emg_gm_1, GMMTJ_emg_gl_1, GMMTJ_emg_sol_1, GMMTJ_time_1, GMMTJ_torque_1] = extract_force_displ_singletrial_passive_EMG(dm_ROM_gmmtj1_NX{line}, dm_ROM_gmmtj1_US{line}, dm_ROM_gmmtj1_US_frame{line}, EMG_max_TA, EMG_max_gm, EMG_max_gl, EMG_max_sol, dm_leg_length{line}, dm_side{line}, line, 'GMMTJ1');
                emg_all{1,3} = [GMMTJ_gonio_1 GMMTJ_emg_gm_1];
                emg_all{2,3} = [GMMTJ_gonio_1 GMMTJ_emg_gl_1];
                emg_all{3,3} = [GMMTJ_gonio_1 GMMTJ_emg_sol_1];
            end
            if(strcmpi(dm_ROM_gmmtj2_NX{line}, 'null'))
                GMMTJ_force_2 = zeros;
                GMMTJ_torque_2 = zeros;
                GMMTJ_angle_2 = zeros;
                GMMTJ_gonio_2 = zeros;
                GMMTJ_displacement_2 = zeros;
            else 
                [GMMTJ_force_2, GMMTJ_gonio_2, GMMTJ_angle_2, GMMTJ_displacement_2, GMMTJ_emg_gm_2, GMMTJ_emg_gl_2, GMMTJ_emg_sol_2, GMMTJ_time_2, GMMTJ_torque_2] = extract_force_displ_singletrial_passive_EMG(dm_ROM_gmmtj2_NX{line}, dm_ROM_gmmtj2_US{line}, dm_ROM_gmmtj2_US_frame{line}, EMG_max_TA, EMG_max_gm, EMG_max_gl, EMG_max_sol, dm_leg_length{line}, dm_side{line}, line, 'GMMTJ2');    
                emg_all{1,4} = [GMMTJ_gonio_2 GMMTJ_emg_gm_2];
                emg_all{2,4} = [GMMTJ_gonio_2 GMMTJ_emg_gl_2];
                emg_all{3,4} = [GMMTJ_gonio_2 GMMTJ_emg_sol_2];
            end

            % average two passive trials
            if(strcmpi(dm_ROM_gmmtj1_NX{line}, 'null'))
                data_GMMTJ = average_passive_trials_EMG(GMMTJ_force_2, GMMTJ_gonio_2, GMMTJ_angle_2, GMMTJ_displacement_2, GMMTJ_emg_gm_2, GMMTJ_emg_gl_2, GMMTJ_emg_sol_2, GMMTJ_time_2, GMMTJ_torque_2);
            elseif(strcmpi(dm_ROM_gmmtj2_NX{line}, 'null'))
                data_GMMTJ = average_passive_trials_EMG(GMMTJ_force_1, GMMTJ_gonio_1, GMMTJ_angle_1, GMMTJ_displacement_1, GMMTJ_emg_gm_1, GMMTJ_emg_gl_1, GMMTJ_emg_sol_1, GMMTJ_time_1, GMMTJ_torque_1);
            else % if 2 trials exist (none of variables are 'null')
                data_GMMTJ = average_passive_trials_EMG(GMMTJ_force_1, GMMTJ_gonio_1, GMMTJ_angle_1, GMMTJ_displacement_1, GMMTJ_emg_gm_1, GMMTJ_emg_gl_1, GMMTJ_emg_sol_1, GMMTJ_time_1, GMMTJ_torque_1, GMMTJ_force_2, GMMTJ_gonio_2, GMMTJ_angle_2, GMMTJ_displacement_2, GMMTJ_emg_gm_2, GMMTJ_emg_gl_2, GMMTJ_emg_sol_2, GMMTJ_time_2, GMMTJ_torque_2);
            end

            % plot summary of 2 trials: Displacement-angle
            if plot_check && plot_individual && plot_old 
                plottitle = horzcat('IND displacement GMMTJ vs angle, ', subject_id);
                figure('Name',plottitle);
                plot(GMMTJ_gonio_1,GMMTJ_displacement_1,'LineWidth',2)
                hold on
                plot(GMMTJ_gonio_2,GMMTJ_displacement_2,'LineWidth',2)
                axis(axis_el_GMmsc)
                ylabel(txt_displ)
                xlabel(txt_gonio)
                title(plottitle,'Interpreter', 'none')
                legend('Trial 1','Trial 2','Location','Southeast')
                print(horzcat('data_plots/',plottitle),'-dpng')
            end
        

            %% Calculations for 2x GM fascicle trials - ORIGINAL calculations, as for GMMTJ and SOL

            % extract force, gonio, angle, displacement for EACH TRIAL
            if(strcmpi(dm_ROM_gmfas1_NX{line}, 'null'))
                % allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case in final_stiffness
                GMFAS_force_1 = zeros;
                GMFAS_torque_1 = zeros;
                GMFAS_angle_1 = zeros;
                GMFAS_gonio_1 = zeros;
                GMFAS_displacement_1 = zeros;
            else 
                [GMFAS_force_1, GMFAS_gonio_1, GMFAS_angle_1, GMFAS_displacement_1, GMFAS_emg_gm_1, GMFAS_emg_gl_1, GMFAS_emg_sol_1, GMFAS_time_1, GMFAS_torque_1] = extract_force_displ_singletrial_passive_EMG(dm_ROM_gmfas1_NX{line}, dm_ROM_gmfas1_US{line}, dm_ROM_gmfas1_US_frame{line}, EMG_max_TA, EMG_max_gm, EMG_max_gl, EMG_max_sol, dm_leg_length{line}, dm_side{line}, line, 'GMFAS1');
                emg_all{1,5} = [GMFAS_gonio_1 GMFAS_emg_gm_1];
                emg_all{2,5} = [GMFAS_gonio_1 GMFAS_emg_gl_1];
                emg_all{3,5} = [GMFAS_gonio_1 GMFAS_emg_sol_1];
            end
            if(strcmpi(dm_ROM_gmfas2_NX{line}, 'null'))
                GMFAS_force_2 = zeros;
                GMFAS_torque_2 = zeros;
                GMFAS_angle_2 = zeros;
                GMFAS_gonio_2 = zeros;
                GMFAS_displacement_2 = zeros;
            else 
                [GMFAS_force_2, GMFAS_gonio_2, GMFAS_angle_2, GMFAS_displacement_2, GMFAS_emg_gm_2, GMFAS_emg_gl_2, GMFAS_emg_sol_2, GMFAS_time_2, GMFAS_torque_2] = extract_force_displ_singletrial_passive_EMG(dm_ROM_gmfas2_NX{line}, dm_ROM_gmfas2_US{line}, dm_ROM_gmfas2_US_frame{line}, EMG_max_TA, EMG_max_gm, EMG_max_gl, EMG_max_sol, dm_leg_length{line}, dm_side{line}, line, 'GMFAS2');    
                emg_all{1,6} = [GMFAS_gonio_2 GMFAS_emg_gm_2];
                emg_all{2,6} = [GMFAS_gonio_2 GMFAS_emg_gl_2];
                emg_all{3,6} = [GMFAS_gonio_2 GMFAS_emg_sol_2];
            end

            % average two passive trials
            if(strcmpi(dm_ROM_gmfas1_NX{line}, 'null'))
                data_GMFAS = average_passive_trials_EMG(GMFAS_force_2, GMFAS_gonio_2, GMFAS_angle_2, GMFAS_displacement_2, GMFAS_emg_gm_2, GMFAS_emg_gl_2, GMFAS_emg_sol_2, GMFAS_time_2, GMFAS_torque_2);
            elseif(strcmpi(dm_ROM_gmfas2_NX{line}, 'null'))
                data_GMFAS = average_passive_trials_EMG(GMFAS_force_1, GMFAS_gonio_1, GMFAS_angle_1, GMFAS_displacement_1, GMFAS_emg_gm_1, GMFAS_emg_gl_1, GMFAS_emg_sol_1, GMFAS_time_1, GMFAS_torque_1);
            else % if 2 trials exist (none of variables are 'null')
                data_GMFAS = average_passive_trials_EMG(GMFAS_force_1, GMFAS_gonio_1, GMFAS_angle_1, GMFAS_displacement_1, GMFAS_emg_gm_1, GMFAS_emg_gl_1, GMFAS_emg_sol_1, GMFAS_time_1, GMFAS_torque_1, GMFAS_force_2, GMFAS_gonio_2, GMFAS_angle_2, GMFAS_displacement_2, GMFAS_emg_gm_2, GMFAS_emg_gl_2, GMFAS_emg_sol_2, GMFAS_time_2, GMFAS_torque_2);
            end

            % plot summary of 2 trials: Displacement-angle
            if plot_check && plot_individual && plot_old 
                plottitle = horzcat('IND displacement GMFAS vs angle, ', subject_id);
                figure('Name',plottitle);
                plot(GMFAS_gonio_1,GMFAS_displacement_1,'LineWidth',2)
                hold on
                plot(GMFAS_gonio_2,GMFAS_displacement_2,'LineWidth',2)
                axis(axis_displ_GMFAS)
                ylabel(txt_displ)
                xlabel(txt_gonio)
                title(plottitle,'Interpreter', 'none')
                legend('Trial 1','Trial 2','Location','Southeast')
                print(horzcat('data_plots/',plottitle),'-dpng')
            end
     

            %% Calculations for 2x GM fascicle trials - LICHTWARK analysis

            % read Lichtwark data, trial 1
            GMFAS_licht_SOL_1_exists = 0;
            if(strcmpi(dm_ROM_gmfas1_licht{line}, 'null'))
                % allow for the possibility of discarded trials (null)
            else 
                [GMFAS_licht_data_1] = read_us_licht(strcat(filepath, dm_ROM_gmfas1_licht{line}, '.txt'), str2double(dm_ROM_gmfas1_US_frame{line}), 'GMFAS1_Licht');
                % check if trial has both GM + SOL, or only GM
                if length(GMFAS_licht_data_1(1,:)) == 3 % GM only
                    GMFAS_licht_GM_faslen_1 = GMFAS_licht_data_1(1:length(GMFAS_gonio_1),2);
                    GMFAS_licht_GM_pennation_1 = GMFAS_licht_data_1(1:length(GMFAS_gonio_1),3);
                    GMFAS_licht_SOL_faslen_1 = zeros(0);
                    GMFAS_licht_SOL_pennation_1 = zeros(0);
                else % length = 6, GM + SOL
                    GMFAS_licht_SOL_1_exists = 1;
                    GMFAS_licht_GM_faslen_1 = GMFAS_licht_data_1(1:length(GMFAS_gonio_1),2);
                    GMFAS_licht_GM_pennation_1 = GMFAS_licht_data_1(1:length(GMFAS_gonio_1),3);
                    GMFAS_licht_SOL_faslen_1 = GMFAS_licht_data_1(1:length(GMFAS_gonio_1),5);
                    GMFAS_licht_SOL_pennation_1 = GMFAS_licht_data_1(1:length(GMFAS_gonio_1),6);
                end
            end

            % read Lichtwark data, trial 2
            GMFAS_licht_SOL_2_exists = 0;
            if(strcmpi(dm_ROM_gmfas2_licht{line}, 'null'))
                % allow for the possibility of discarded trials (null)
            else 
                [GMFAS_licht_data_2] = read_us_licht(strcat(filepath, dm_ROM_gmfas2_licht{line}, '.txt'), str2double(dm_ROM_gmfas2_US_frame{line}), 'GMFAS2_Licht');
                % check if trial has both GM + SOL, or only GM
                if length(GMFAS_licht_data_2(1,:)) == 3 % GM only
                    GMFAS_licht_GM_faslen_2 = GMFAS_licht_data_2(1:length(GMFAS_gonio_2),2);
                    GMFAS_licht_GM_pennation_2 = GMFAS_licht_data_2(1:length(GMFAS_gonio_2),3);
                    GMFAS_licht_SOL_faslen_2 = zeros(0);
                    GMFAS_licht_SOL_pennation_2 = zeros(0);
                else % length = 6, GM + SOL
                    GMFAS_licht_SOL_2_exists = 1;
                    GMFAS_licht_GM_faslen_2 = GMFAS_licht_data_2(1:length(GMFAS_gonio_2),2);
                    GMFAS_licht_GM_pennation_2 = GMFAS_licht_data_2(1:length(GMFAS_gonio_2),3);
                    GMFAS_licht_SOL_faslen_2 = GMFAS_licht_data_2(1:length(GMFAS_gonio_2),5);
                    GMFAS_licht_SOL_pennation_2 = GMFAS_licht_data_2(1:length(GMFAS_gonio_2),6);
                end
            end

            % perform averaging of lichtwark data:
            % create two arrays:
            %   data_GMFAS_licht_GM
            %   data_GMFAS_licht_SOL
            % containing:
            %   averaged angle (currently calculated from gonio)
            %   averaged fascicle length
            %   averaged pennation angle
            %   averaged fascicle elongation
            %   averaged fascicle strain
            % OR containing zeros (if nonexistent)

            if(strcmpi(dm_ROM_gmfas1_licht{line}, 'null')==0 && strcmpi(dm_ROM_gmfas2_licht{line}, 'null')==0)

                % if 2 GM trials exist (none of variables are 'null')

                % GM Lichtwark: perform averaging of trial 1 and trial 2
                data_GMFAS_licht_GM = average_passive_trials_licht(GMFAS_gonio_1, GMFAS_angle_1, GMFAS_licht_GM_faslen_1, GMFAS_licht_GM_pennation_1, GMFAS_time_1, GMFAS_gonio_2, GMFAS_angle_2, GMFAS_licht_GM_faslen_2, GMFAS_licht_GM_pennation_2, GMFAS_time_2, out_ROM_trial_max, 'GM fascicles');

                % SOL Lichtwark: check for existence of SOL data
                if GMFAS_licht_SOL_1_exists && GMFAS_licht_SOL_2_exists
                    % average two trials:
                    data_GMFAS_licht_SOL = average_passive_trials_licht(GMFAS_gonio_1, GMFAS_angle_1, GMFAS_licht_SOL_faslen_1, GMFAS_licht_SOL_pennation_1, GMFAS_time_1, GMFAS_gonio_2, GMFAS_angle_2, GMFAS_licht_SOL_faslen_2, GMFAS_licht_SOL_pennation_2, GMFAS_time_2, out_ROM_trial_max,'SOL fascicles');
                elseif GMFAS_licht_SOL_1_exists
                    data_GMFAS_licht_SOL = average_passive_trials_licht(GMFAS_gonio_1, GMFAS_angle_1, GMFAS_licht_SOL_faslen_1, GMFAS_licht_SOL_pennation_1, GMFAS_time_1, out_ROM_trial_max, 'SOL fascicles');
                elseif GMFAS_licht_SOL_2_exists
                    data_GMFAS_licht_SOL = average_passive_trials_licht(GMFAS_gonio_2, GMFAS_angle_2, GMFAS_licht_SOL_faslen_2, GMFAS_licht_SOL_pennation_2, GMFAS_time_2, out_ROM_trial_max, 'SOL fascicles');
                else % no SOL exists
                    data_GMFAS_licht_SOL = zeros(1,3);
                end

            elseif(strcmpi(dm_ROM_gmfas1_licht{line}, 'null')==0) % only trial 1 exists
                % keep GM trial 1
                data_GMFAS_licht_GM = average_passive_trials_licht(GMFAS_gonio_1, GMFAS_angle_1, GMFAS_licht_GM_faslen_1, GMFAS_licht_GM_pennation_1, GMFAS_time_1, out_ROM_trial_max, 'GM fascicles');
                % keep eventual SOL trial 1
                if GMFAS_licht_SOL_1_exists
                    data_GMFAS_licht_SOL = average_passive_trials_licht(GMFAS_gonio_1, GMFAS_angle_1, GMFAS_licht_SOL_faslen_1, GMFAS_licht_SOL_pennation_1, GMFAS_time_1, out_ROM_trial_max, 'SOL fascicles');
                else % no SOL exists
                    data_GMFAS_licht_SOL = zeros(1,3);
                end    

            elseif(strcmpi(dm_ROM_gmfas2_licht{line}, 'null')==0) % only trial 2 exists
                % keep GM trial 2
                data_GMFAS_licht_GM = average_passive_trials_licht(GMFAS_gonio_2, GMFAS_angle_2, GMFAS_licht_GM_faslen_2, GMFAS_licht_GM_pennation_2, GMFAS_time_2, out_ROM_trial_max, 'GM fascicles');
                % keep eventual SOL trial 2
                if GMFAS_licht_SOL_2_exists
                    data_GMFAS_licht_SOL = average_passive_trials_licht(GMFAS_gonio_2, GMFAS_angle_2, GMFAS_licht_SOL_faslen_2, GMFAS_licht_SOL_pennation_2, GMFAS_time_2, out_ROM_trial_max, 'SOL fascicles');
                else % no SOL exists
                    data_GMFAS_licht_SOL = zeros(1,3);
                end    

            else % no trials exist
                data_GMFAS_licht_GM = zeros(1,5);
                data_GMFAS_licht_SOL = zeros(1,3);
            end

            % add GM fascicle elongation and strain to data_GMFAS_licht_GM
            if length(data_GMFAS_licht_GM) > 5
    %             % version 1 - length at zero degrees = zero elong, zero strain
    %             angle_zero = 0 - 0.00000001; % tweak for computer storage of numbers, where 8.3000 is stored as 8.2999999999999999 %VAR
    %             loc_angle_zero = find(data_GMFAS_licht_GM(:,1)>=angle_zero,1,'first');
    %             data_GMFAS_licht_GM(:,4) = data_GMFAS_licht_GM(:,2) - data_GMFAS_licht_GM(loc_angle_zero,2); %elong
    %             data_GMFAS_licht_GM(:,5) = (data_GMFAS_licht_GM(:,2) - data_GMFAS_licht_GM(loc_angle_zero,2)) / data_GMFAS_licht_GM(loc_angle_zero,2) * 100; %strain

                % version 2 - length PRONE is zero elong, zero strain
                data_GMFAS_licht_GM(:,4) = data_GMFAS_licht_GM(:,2) - prone_GMfas_length; %fascicle elong 
                data_GMFAS_licht_GM(:,5) = (data_GMFAS_licht_GM(:,2) - prone_GMfas_length) / prone_GMfas_length * 100; %fascicle strain

                if plot_check && plot_individual
                    plottitle = horzcat('IND Lichtwark GM fascicle vs angle, ', subject_id);
                    figure('Name',plottitle);
                    hold on
                    mat_version = version('-release');
                    % left axis
                    if strcmp(mat_version,'2015b') == 0
                        yyaxis left
                    end
                    plot(data_GMFAS_licht_GM(:,1), data_GMFAS_licht_GM(:,2),'b','LineWidth',2) % fascicle length
                    plot([prone_GMfas_ankle prone_GMfas_ankle],[prone_GMfas_length prone_GMfas_length], 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b') % fas len at rest/prone
                    ylabel(txt_length)
                    if strcmp(mat_version,'2015b') == 0
                        yyaxis right
                    end
                    plot(data_GMFAS_licht_GM(:,1), data_GMFAS_licht_GM(:,3),'r','LineWidth',2) % pennation angle
                    plot([prone_GMfas_ankle prone_GMfas_ankle],[prone_GMfas_penn prone_GMfas_penn], 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r') % penn angle at rest/prone
                    plot(data_GMFAS_licht_GM(:,1), data_GMFAS_licht_GM(:,5),':r','LineWidth',1) % fascicle strain
                    plot(data_GMFAS_licht_GM(:,1), data_GMFAS_licht_GM(:,4),'--r','LineWidth',1) % fascicle elongation
                    ylabel('Elongation (mm) / Strain (%) / Pennation angle (°)')
                    xlabel(txt_gonio)
                    title(plottitle,'Interpreter', 'none')
                    legend('Fasc length','Faslen PRONE', 'Pennation angle', 'Penn ang PRONE', 'Fasc strain', 'Fasc elong', 'Location','West')
                    % not possible to break axis with two y-axes: "Object Copy of Axes with multiple coordinate systems is not supported."
                    % breakxaxis([-15 -5]);
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            % LATER - Fukunaga fascicle/muscle/SEE faslen/elong/strain also for SOL fascicle scans?
            end


            %% Average SOL + GMMTJ + GMFAS trials for force, gonio, angle, EMG
            % data_force_gonio =
            %   force
            %   angle (from gonio)
            %   torque
            %   EMG GM
            %   EMG GL
            %   EMG SOL
            
            % average 3 force arrays into one
            data_force_gonio = average_passive_forces_EMG(data_SOL, data_GMMTJ, data_GMFAS, trial_subjectno);

            % plot all 6 trials separately: Force-angle
            if plot_check && plot_individual
                plottitle = horzcat('IND force vs angle singletrials, ', subject_id);
                figure('Name',plottitle);
                plot(SOL_gonio_1,SOL_force_1,'LineWidth',2, 'Color',[1 0 0])
                hold on
                plot(SOL_gonio_2,SOL_force_2,'LineWidth',2, 'Color',[1 0.6 0])
                plot(GMMTJ_gonio_1,GMMTJ_force_1,'LineWidth',2, 'Color',[1 1 0])
                plot(GMMTJ_gonio_2,GMMTJ_force_2,'LineWidth',2, 'Color',[0 1 0])
                plot(GMFAS_gonio_1,GMFAS_force_1,'LineWidth',2, 'Color',[0 0 1])
                plot(GMFAS_gonio_2,GMFAS_force_2,'LineWidth',2, 'Color',[1 0 1])
                plot(data_force_gonio(:,2), data_force_gonio(:,1), 'LineWidth',2, 'Color','black')
                axis(axis_force)
                ylabel('Force (N)')
                xlabel(txt_gonio)
                title(plottitle,'Interpreter', 'none')
                legend('SOL1','SOL2','GMMTJ1','GMMTJ2','GMFAS1','GMFAS2','mean','Location','Southeast')
                print(horzcat('data_plots/',plottitle),'-dpng')
            end

            % Plot all 6 trials separately: EMG vs angle
            if plot_check && plot_individual
                plottitle = horzcat('IND EMG vs angle, ', subject_id);
                figure('Name',plottitle)
                % tweak to get proper legend:
                hold on
                plot(0,0,'y');
                plot(0,0,'m');
                plot(0,0,'c');
                for i = 1:6 % 6 trials GM
                    if ~isempty(emg_all{1,i})
                        plot(emg_all{1,i}(:,1),emg_all{1,i}(:,2),'y');
                    end
                end
                plot(data_force_gonio(:,2),data_force_gonio(:,4),'y','LineWidth',2); % mean GM

                for i = 1:6 % 6 trials GL
                    if ~isempty(emg_all{1,i})
                        plot(emg_all{2,i}(:,1),emg_all{2,i}(:,2),'m');
                    end
                end
                plot(data_force_gonio(:,2),data_force_gonio(:,5),'m','LineWidth',2); % mean GL

                for i = 1:6 % 6 trials SOL
                    if ~isempty(emg_all{1,i})
                        plot(emg_all{3,i}(:,1),emg_all{3,i}(:,2),'c');
                    end
                end
                plot(data_force_gonio(:,2),data_force_gonio(:,6),'c','LineWidth',2); % mean SOL
                
                plot(data_force_gonio(:,2),data_force_gonio(:,7),'LineWidth',2); % mean all 3 muscles
                
                axis(axis_EMG_wide)
                ylabel('Muscle activation (% of MVC)')
                xlabel(txt_gonio)
                title(plottitle,'Interpreter', 'none')
                legend('Gastr.med.','Gastr.lat.','Soleus','Location','Northwest')
                print(horzcat('data_plots/',plottitle),'-dpng')
            end
            
%             % TMP: separate EMG plots
%             if plot_check && plot_individual
% 
%                 figure('Name','EMG GM')
%                 plot(data_force_gonio(:,2),data_force_gonio(:,4),'y','LineWidth',2); % mean GM
%                 hold on
%                 for i = 1:6 % 6 trials GM
%                     if ~isempty(emg_all{1,i})
%                         plot(emg_all{1,i}(:,1),emg_all{1,i}(:,2));
%                     end
%                 end
% 
%                 figure('Name','EMG GL')
%                 plot(data_force_gonio(:,2),data_force_gonio(:,5),'m','LineWidth',2); % mean GL
%                 hold on
%                 for i = 1:6 % 6 trials GL
%                     if ~isempty(emg_all{1,i})
%                         plot(emg_all{2,i}(:,1),emg_all{2,i}(:,2));
%                     end
%                 end
% 
%                 figure('Name','EMG SOL')
%                 plot(data_force_gonio(:,2),data_force_gonio(:,6),'c','LineWidth',2); % mean SOL
%                 hold on
%                 for i = 1:6 % 6 trials SOL
%                     if ~isempty(emg_all{1,i})
%                         plot(emg_all{3,i}(:,1),emg_all{3,i}(:,2));
%                     end
%                 end
% 
%             end


            %% Check conformation of GONIOMETER to norm angle
            if plot_norm
                plottitle = horzcat('IND Goniometer check 6 trials, ', subject_id);
                figure('Name',plottitle);
                hold on
                plot(SOL_angle_1, SOL_gonio_1, 'LineWidth',2, 'Color',[1 0 0])
                plot(SOL_angle_2, SOL_gonio_2, 'LineWidth',2, 'Color',[1 0.6 0])
                plot(GMMTJ_angle_1, GMMTJ_gonio_1, 'LineWidth',2, 'Color',[1 1 0])
                plot(GMMTJ_angle_2, GMMTJ_gonio_2, 'LineWidth',2, 'Color',[0 1 0])
                plot(GMFAS_angle_1, GMFAS_gonio_1, 'LineWidth',2, 'Color',[0 0 1])
                plot(GMFAS_angle_2, GMFAS_gonio_2, 'LineWidth',2, 'Color',[1 0 1])
                if length(SOL_angle_1) == 1
                    plot(SOL_angle_2, SOL_angle_2, 'LineWidth',2, 'Color',col_grey, 'LineStyle',':') % conformation line
                else
                    plot(SOL_angle_1, SOL_angle_1, 'LineWidth',2, 'Color',col_grey, 'LineStyle',':') % conformation line
                end
                plot([0 0],[0 0], 'Marker','o', 'MarkerSize',10, 'MarkerFaceColor',col_grey)% zero point
                xlabel('Norm angle (°)')
                ylabel(txt_gonio)
                title(plottitle,'Interpreter', 'none')
                legend('SOL1','SOL2','GMMTJ1','GMMTJ2','GMFAS1','GMFAS2','Norm/Norm','Location','Southeast')
                print(horzcat('data_plots/',plottitle),'-dpng')
            end


            %% Extract arrays of MTU LENGTH, ELONGATION, STRAIN

            % MTU_length_array =                            MTU_elong_array =      MTU_strain_array =
            %1       angle
            %2       calc to SOL ins = free tendon                 
            %3       calc go GM ins = GM tendon
            %4       calc to knee = entire calf
            %5       DISPLACEMENT from GMFAS tracking
            %6       SOL to GM = GM apo
            %7       GM to knee = GM msc
            %8       SOL length (H&H) minus free AT = SOL msc
            %9       SEE Fukunaga
            %10      GM muscle Fukunaga (longitudinal aspect of one fascicle)

            % column choices MTU ELONGATION / LENGTH / STRAIN
            col_angle = 1;
            col_AT = 2;
            col_GMtend = 3;
            col_leg = 4;
            col_GMFAS = 5;
            col_GMapo = 6;
            col_GMmsc = 7;
            col_SOLmsc = 8;
            col_GMmsc_Fukunaga = 9;
            col_SEE_Fukunaga = 10;
            col_penn_ang = 11;
            col_GMfaslen = 12;

            % calculate MTU lengths
            [MTU_length_array, MTU_elong_array, MTU_strain_array, MTU_percentelong_array, MTU_prone_vars] = calculate_mtu_length(data_SOL(:,2:3), data_GMMTJ(:,2:3), data_GMFAS(:,2:3), data_GMFAS_licht_GM, dm_at_SOL_length{line}, dm_at_GM_length{line}, trial_calf_length, prone_GMfas_penn, prone_GMfas_length, out_ROM_trial_max, prone_GMfas_ankle, dm_at_SOL_zero{line}, dm_at_GM_zero{line});
             % NB: prone_vars contains pennation angle and fascicle length
             %     in the columns normally used for GMFAS displacement and apo


            %% NORMALIZE Lichtwark/Fukunaga length/elongation to prone leg length
            % MTU NORMALIZED:
            %   1 = length MTU %len
            %   2 = length GMmsc %len
            %   3 = length SEE %len
            %   4 = elong MTU %len
            %   5 = elong GMmsc %len
            %   6 = elong SEE %len
            %   7 = elong GMmsc %elong
            %   8 = elong SEE %elong
            % MTU normalized licht:
            %   2 = length GMfas %len
            %   3 = elong GMfas %len
            
            % resting MTU length is found in variable:
            MTU_length_rest = MTU_prone_vars(col_leg);
                
            % Normalization of MTU, GMmsc, SEE
            % to MTU length (prone rest) and MTU elongation (from prone rest)

            % GMmsc/SEE/MTU LENGTH in % of prone MTU LENGTH
            MTU_normalized(:,1) = MTU_length_array(:,col_leg)/MTU_length_rest*100;
            MTU_normalized(:,2) = MTU_length_array(:,col_GMmsc_Fukunaga)/MTU_length_rest*100;
            MTU_normalized(:,3) = MTU_length_array(:,col_SEE_Fukunaga)/MTU_length_rest*100;
            % MTU elongation in % of prone MTU LENGTH
            MTU_normalized(:,4) = MTU_elong_array(:,col_leg)/MTU_length_rest*100;
            % GMmsc/SEE elongation in % of prone MTU LENGTH
            MTU_normalized(:,5) = MTU_elong_array(:,col_GMmsc_Fukunaga)/MTU_length_rest*100;
            MTU_normalized(:,6) = MTU_elong_array(:,col_SEE_Fukunaga)/MTU_length_rest*100;
            % GMmsc/SEE elongation in % of MTU ELONGATION
            % version 1: % of elongation from resting length
            %MTU_normalized(:,7) = MTU_elong_array(:,col_GMmsc_Fukunaga)./MTU_elong_array(:,col_leg)*100;
            %MTU_normalized(:,8) = MTU_elong_array(:,col_SEE_Fukunaga)./MTU_elong_array(:,col_leg)*100;
            % version 2: % of elongation from zero length
            MTU_normalized(:,7) = MTU_percentelong_array(:,1); % GMmsc
            MTU_normalized(:,8) = MTU_percentelong_array(:,2); % SEE

            % Normalization of fascicle length/elongation 
            % to prone MTU length
            
            MTU_normalized_licht(:,1) = data_GMFAS_licht_GM(:,1);
            MTU_normalized_licht(:,2) = data_GMFAS_licht_GM(:,2)/MTU_length_rest*100; %    2 averaged fasicle length
            MTU_normalized_licht(:,3) = data_GMFAS_licht_GM(:,4)/MTU_length_rest*100; %    4 averaged fasicle elongation
                

            %% Plot arrays of MTU LENGTH, ELONGATION, STRAIN
            if plot_individual && plot_old
                
                % raw lengths (mm)
                plottitle = horzcat('IND MTU length vs angle v1, ', subject_id);
                figure('Name',plottitle);
                hold on
                plot(MTU_length_array(:,col_angle),MTU_length_array(:,col_leg),'LineWidth',2,'Color','blue') % AT + GM apo + GM msc = full GM MTU / calf length
                plot(MTU_length_array(:,col_angle),MTU_length_array(:,col_SEE_Fukunaga),'LineWidth',2,'LineStyle','-','Color','red') % SEE from anthropometry Lichtwark/Fukunaga
                plot(MTU_length_array(:,col_angle),MTU_length_array(:,col_SOLmsc)+MTU_length_array(:,col_AT),'LineWidth',0.8,'Color','black') % AT + SOL msc = SOL MTU
                plot(MTU_length_array(:,col_angle),MTU_length_array(:,col_GMtend),'LineWidth',0.8,'Color','red') % AT + GM apo = GM tendon (linear)
                plot(MTU_length_array(:,col_angle),MTU_length_array(:,col_AT),'LineWidth',0.8,'Color','green') % AT
                % plot vertical lines to show color coding of subsequent figs
                plot([0,0], [MTU_length_array(1,col_leg), MTU_length_array(1,col_GMtend)],'Color','cyan')
                plot([1,1], [MTU_length_array(1,col_leg), MTU_length_array(1,col_SEE_Fukunaga)],'Color','cyan','LineWidth',2)
                plot([1,1], [MTU_length_array(1,col_SOLmsc)+MTU_length_array(1,col_AT), MTU_length_array(1,col_GMtend)],'Color','black')
                plot([0,0], [MTU_length_array(1,col_GMtend), MTU_length_array(1,col_AT)],'Color',col_orange)
                plot([0,0], [MTU_length_array(1,col_AT), 0],'Color','green')

                axis(axis_ind_len_MTU)
                ylabel(txt_length)
                xlabel(txt_gonio)
                title(plottitle,'Interpreter', 'none')
                legend('Full GM MTU', 'SEE (archi)', 'SOL MTU', 'GM tendon (linear)', 'AT free', 'GM msc. (linear)', 'GM msc. (from archi)', 'SOL msc.', 'GM apo.', 'Location','Southeast');
                print(horzcat('data_plots/',plottitle),'-dpng')

                % raw elongation (mm) - linear
                plottitle = horzcat('IND MTU elongation (linear) vs angle, ', subject_id);
                figure('Name',plottitle);
                hold on
                plot(MTU_elong_array(:,col_angle),MTU_elong_array(:,col_leg),'LineWidth',1.5,'Color','blue') % MTU = AT + GM apo + msc
                plot(MTU_elong_array(:,col_angle),MTU_elong_array(:,col_AT),'LineWidth',0.8,'Color','green') % AT
                plot(MTU_elong_array(:,col_angle),MTU_elong_array(:,col_GMtend),'LineWidth',0.8,'Color','red') % GM tendon = AT + GM apo
                plot(MTU_elong_array(:,col_angle),MTU_elong_array(:,col_GMmsc),'LineWidth',0.8,'Color','cyan') % GM msc
                plot(MTU_elong_array(:,col_angle),MTU_elong_array(:,col_SOLmsc),'LineWidth',0.8,'Color','black') % SOL msc
                plot(MTU_elong_array(:,col_angle),MTU_elong_array(:,col_GMFAS),'LineWidth',0.8,'Color','yellow') % GM FAS displacement
                plot(MTU_elong_array(:,col_angle),MTU_elong_array(:,col_GMapo),'LineWidth',0.8,'Color',col_orange) % GM apo
                axis(axis_ind_elong)
                ylabel(txt_elong)
                xlabel(txt_gonio)
                title(plottitle,'Interpreter', 'none')
                legend('GM MTU', 'AT free', 'GM tendon (linear)', 'GM msc. (linear)', 'SOL msc.', 'GM fasc.DISPL.', 'GM apo.', 'Location','Southeast');
                print(horzcat('data_plots/',plottitle),'-dpng')
                xlim([-0.5 inf]) 
                ylim([-inf inf]) 
                print(horzcat('data_plots/',plottitle, ' ZOOM'),'-dpng')
                
                % strain (percent of initial length) - linear
                plottitle = horzcat('IND MTU strain vs angle v1, ', subject_id);
                figure('Name',plottitle);
                hold on
                plot(MTU_strain_array(:,col_angle),MTU_strain_array(:,col_leg),'LineWidth',1.5,'Color','blue') % MTU = AT + GM apo + msc
                plot(MTU_strain_array(:,col_angle),MTU_strain_array(:,col_AT),'LineWidth',0.8,'Color','green') % AT
                plot(MTU_strain_array(:,col_angle),MTU_strain_array(:,col_GMtend),'LineWidth',0.8,'Color','red') % GM tendon = AT + GM apo
                plot(MTU_strain_array(:,col_angle),MTU_strain_array(:,col_GMmsc),'LineWidth',0.8,'Color','cyan') % GM msc
                plot(MTU_strain_array(:,col_angle),MTU_strain_array(:,col_SOLmsc),'LineWidth',0.8,'Color','black') % SOL msc
                plot(MTU_strain_array(:,col_angle),MTU_strain_array(:,col_GMapo),'LineWidth',0.8,'Color',col_orange) % GM apo
                if max(MTU_strain_array(:,col_AT)) > axis_ind_strain(4)
                    text(2,axis_ind_strain(4)-1, horzcat('Max AT str = ', num2str(round(max(MTU_strain_array(:,col_AT)),1))),'Color','green') % TEXT: max AT strain
                end
                axis(axis_ind_strain)
                ylabel(txt_strain)
                xlabel(txt_gonio)
                title(plottitle,'Interpreter', 'none')
                legend('GM MTU', 'AT free', 'GM tendon (linear)', 'GM msc. (linear)', 'SOL msc.', 'GM apo.', 'Location','Southeast');
                print(horzcat('data_plots/',plottitle),'-dpng')
                xlim([-0.5 inf]) 
                ylim([-inf inf]) 
                print(horzcat('data_plots/', plottitle,' ZOOM'),'-dpng')
            end
            
            if plot_individual
                % raw elongation (mm) - Lichtwark/Fukunaga
                plottitle = horzcat('IND MTU elongation (Lichtwark) vs angle, ', subject_id);
                figure('Name',plottitle);
                hold on
                plot(MTU_elong_array(:,col_angle),MTU_elong_array(:,col_leg),'LineWidth',1.5,'Color','blue') % MTU = AT + GM apo + msc
                plot(MTU_elong_array(:,col_angle),MTU_elong_array(:,col_SEE_Fukunaga),'LineWidth',1.5,'Color','red','LineStyle','-') % SEE from anthropometry Lichtwark/Fukunaga
                plot(MTU_elong_array(:,col_angle),MTU_elong_array(:,col_GMmsc_Fukunaga),'LineWidth',1.5,'Color','cyan','LineStyle','-') % GM msc from anthropometry Lichtwark/Fukunaga
                plot(data_GMFAS_licht_GM(:,1), data_GMFAS_licht_GM(:,4),'LineWidth',1.5,'Color','yellow') % fascicle elongation
                axis(axis_ind_elong)
                ylabel(txt_elong)
                xlabel(txt_gonio)
                title(plottitle,'Interpreter', 'none')
                legend('GM MTU', 'SEE (archi)', 'GM msc. (from archi)', 'GM fascicle', 'Location','Southeast');
                print(horzcat('data_plots/',plottitle),'-dpng')

                % strain (percent of initial length) - Lichtwark
                plottitle = horzcat('IND MTU strain v2 vs angle, ', subject_id);
                figure('Name',plottitle);
                hold on
                plot(MTU_strain_array(:,col_angle),MTU_strain_array(:,col_leg),'LineWidth',1.5,'Color','blue') % MTU = AT + GM apo + msc
                plot(MTU_strain_array(:,col_angle),MTU_strain_array(:,col_SEE_Fukunaga),'LineWidth',1.5,'Color','red','LineStyle','-') % SEE from anthropometry Lichtwark/Fukunaga
                plot(MTU_strain_array(:,col_angle),MTU_strain_array(:,col_GMmsc_Fukunaga),'LineWidth',1.5,'Color','cyan','LineStyle','-') % GM msc from anthropometry Lichtwark/Fukunaga
                plot(data_GMFAS_licht_GM(:,1), data_GMFAS_licht_GM(:,5),'LineWidth',1.5,'Color','yellow') % fascicle strain
                axis(axis_ind_strain)
                ylabel(txt_strain)
                xlabel(txt_gonio)
                title(plottitle,'Interpreter', 'none')
                legend('GM MTU', 'SEE (archi)', 'GM msc. (from archi)', 'GM fascicle', 'Location','Southeast');
                print(horzcat('data_plots/',plottitle),'-dpng')
            end


            %% Extract FORCE, ANGLE, LENGTH/ELONGATION/STRAIN, EMG @ various joint angles
            %   force, angle, EMG: from all 3 scan locations/6 trials, averaged
            %   displacement of SOL, GMMTJ, GMFAS: from 2 trials per scan location, averaged
            %   length, elongation, strain of AT, GM apo, etc: from 2 trials per scan location, averaged

            %   joint angles = 
            %       out_ROM_trial_max = trial max (different PRE, POST, L, R)
            %       out_ROM_ind_max = subject ind max (lowest PRE, POST, L, R)
            %       out_ROM_common_max = common max (lowest of all subjects)
            %       out_ROM_submax_1 OLD = additional predetermined angle: 1/3 of trial max ROM
            %       out_ROM_submax_2 = additional predetermined angle: 2/3 of trial max ROM
            %       out_ROM_submax_1 = 10 degrees for those subjects which have this ROM available 

            % column choices force, torque, displacement
            col_force = 1;      % in data_force_gonio / data_SOL etc
            col_angle_DFG = 2;  % in data_force_gonio
            col_torque = 3;     % in data_force_gonio

            % print error if max angles do not exist --> create_angles_passive needs to be run
            if str2double(input_gon_ind_max{trial_subjectno}) == 100
                cprintf('*red', 'ERROR: Max ROM values are not calculated for subject. Run create_angles_passive first.\n')
            end

            % find relevant goniometer angles 
            if strcmpi(dm_timepoint{line}, 'PRE') == 1
                if strcmpi(dm_side{line}, 'R') == 1
                    out_ROM_trial_max = str2double(input_gon_pre_r{trial_subjectno}) - 0.00000001; % tweak for computer storage of numbers, where 8.3000 is stored as 8.2999999999999999
                    out_ROM_ind_max = str2double(input_gon_ind_rmax{trial_subjectno}) - 0.00000001;
                else % L
                    out_ROM_trial_max = str2double(input_gon_pre_l{trial_subjectno}) - 0.00000001;
                    out_ROM_ind_max = str2double(input_gon_ind_lmax{trial_subjectno}) - 0.00000001;
                end
            else % POST
                if strcmpi(dm_side{line}, 'R') == 1
                    out_ROM_trial_max = str2double(input_gon_post_r{trial_subjectno}) - 0.00000001;
                    out_ROM_ind_max = str2double(input_gon_ind_rmax{trial_subjectno}) - 0.00000001;
                else % L
                    out_ROM_trial_max = str2double(input_gon_post_l{trial_subjectno}) - 0.00000001;
                    out_ROM_ind_max = str2double(input_gon_ind_lmax{trial_subjectno}) - 0.00000001;
                end
            end
            out_ROM_common_max = str2double(input_gon_common_max{trial_subjectno}) - 0.00000001;
            out_ROM_submax_1 = 10; %VAR
            out_ROM_submax_2 =  out_ROM_trial_max * 2/3; %VAR

            % forces (using data_force_gonio = averaged data from 3 scan locations / 6 trials)
            loc_frame = find(data_force_gonio(:,col_angle_DFG)>=out_ROM_trial_max,1,'first'); % when averaging angles across trials, intervals of 0.05 degrees are used. This means that any required angle will exist in all data series
            out_F_trial_max_ROM = data_force_gonio(loc_frame,col_force); % force at highest angle in array
            out_T_trial_max_ROM = data_force_gonio(loc_frame,col_torque);
            out_F_trial_max_F = max(data_force_gonio(:,col_force)); % highest force in array
            out_T_trial_max_F = max(data_force_gonio(:,col_torque));
            loc_frame = find(data_force_gonio(:,col_angle_DFG)>=out_ROM_ind_max,1,'first'); 
            out_F_ind_max = data_force_gonio(loc_frame,col_force);
            out_T_ind_max = data_force_gonio(loc_frame,col_torque);
            loc_frame = find(data_force_gonio(:,col_angle_DFG)>=out_ROM_common_max,1,'first'); 
            out_F_common_max = data_force_gonio(loc_frame,col_force);
            out_T_common_max = data_force_gonio(loc_frame,col_torque);
            loc_frame = find(data_force_gonio(:,col_angle_DFG)>=0,1,'first'); % zero angle
            out_F_zero = data_force_gonio(loc_frame,col_force);
            out_T_zero = data_force_gonio(loc_frame,col_torque);
            loc_frame = find(data_force_gonio(:,col_angle_DFG)>=out_ROM_submax_1,1,'first'); 
            if isempty(loc_frame) == 0
                out_F_submax_1 = data_force_gonio(loc_frame,col_force);
                out_T_submax_1 = data_force_gonio(loc_frame,col_torque);
            else
                out_F_submax_1 = NaN;
                out_T_submax_1 = NaN;
            end
            loc_frame = find(data_force_gonio(:,col_angle_DFG)>=out_ROM_submax_2,1,'first'); 
            out_F_submax_2 = data_force_gonio(loc_frame,col_force);
            out_T_submax_2 = data_force_gonio(loc_frame,col_torque);

%             % displacements SOL
%             loc_frame = find(data_SOL(:,col_angle_DFG)>=out_ROM_trial_max,1,'first'); 
%             out_displ_SOL_trial_max = data_SOL(loc_frame,col_displ);
%             loc_frame = find(data_SOL(:,col_angle_DFG)>=out_ROM_ind_max,1,'first'); 
%             out_displ_SOL_ind_max = data_SOL(loc_frame,col_displ); 
%             loc_frame = find(data_SOL(:,col_angle_DFG)>=out_ROM_common_max,1,'first'); 
%             out_displ_SOL_common_max = data_SOL(loc_frame,col_displ); 
%             loc_frame = find(data_SOL(:,col_angle_DFG)>=out_ROM_submax_1,1,'first'); 
%             if isempty(loc_frame) == 0
%                 out_displ_SOL_submax_1 = data_SOL(loc_frame,col_displ); 
%             else
%                 out_displ_SOL_submax_1 = NaN;
%             end
%             loc_frame = find(data_SOL(:,col_angle_DFG)>=out_ROM_submax_2,1,'first'); 
%             out_displ_SOL_submax_2 = data_SOL(loc_frame,col_displ); 
% 
%             % displacements GMMTJ
%             loc_frame = find(data_GMMTJ(:,col_angle_DFG)>=out_ROM_trial_max,1,'first'); 
%             out_displ_GMMTJ_trial_max = data_GMMTJ(loc_frame,col_displ); 
%             loc_frame = find(data_GMMTJ(:,col_angle_DFG)>=out_ROM_ind_max,1,'first'); 
%             out_displ_GMMTJ_ind_max = data_GMMTJ(loc_frame,col_displ); 
%             loc_frame = find(data_GMMTJ(:,col_angle_DFG)>=out_ROM_common_max,1,'first'); 
%             out_displ_GMMTJ_common_max = data_GMMTJ(loc_frame,col_displ); 
%             loc_frame = find(data_GMMTJ(:,col_angle_DFG)>=out_ROM_submax_1,1,'first'); 
%             if isempty(loc_frame) == 0
%                 out_displ_GMMTJ_submax_1 = data_GMMTJ(loc_frame,col_displ); 
%             else
%                 out_displ_GMMTJ_submax_1 = NaN;
%             end
%             loc_frame = find(data_GMMTJ(:,col_angle_DFG)>=out_ROM_submax_2,1,'first'); 
%             out_displ_GMMTJ_submax_2 = data_GMMTJ(loc_frame,col_displ); 
% 
%             % displacements GMFAS
%             loc_frame = find(data_GMFAS(:,col_angle_DFG)>=out_ROM_trial_max,1,'first'); 
%             out_displ_GMFAS_trial_max = data_GMFAS(loc_frame,col_displ); 
%             loc_frame = find(data_GMFAS(:,col_angle_DFG)>=out_ROM_ind_max,1,'first'); 
%             out_displ_GMFAS_ind_max = data_GMFAS(loc_frame,col_displ); 
%             loc_frame = find(data_GMFAS(:,col_angle_DFG)>=out_ROM_common_max,1,'first'); 
%             out_displ_GMFAS_common_max = data_GMFAS(loc_frame,col_displ); 
%             loc_frame = find(data_GMFAS(:,col_angle_DFG)>=out_ROM_submax_1,1,'first'); 
%             if isempty(loc_frame) == 0
%                 out_displ_GMFAS_submax_1 = data_GMFAS(loc_frame,col_displ); 
%             else
%                 out_displ_GMFAS_submax_1 = NaN;
%             end
%             loc_frame = find(data_GMFAS(:,col_angle_DFG)>=out_ROM_submax_2,1,'first'); 
%             out_displ_GMFAS_submax_2 = data_GMFAS(loc_frame,col_displ); 
            
            % elongations, strains @ trial max
            loc_frame = find(MTU_elong_array(:,col_angle)>=out_ROM_trial_max,1,'first'); 
            
            out_elong_AT_trial_max = MTU_elong_array(loc_frame,col_AT);
            out_elong_GMtend_trial_max = MTU_elong_array(loc_frame,col_GMtend);
            out_elong_leg_trial_max = MTU_elong_array(loc_frame,col_leg);
            out_elong_GMapo_trial_max = MTU_elong_array(loc_frame,col_GMapo);
            out_elong_msc_GM_trial_max = MTU_elong_array(loc_frame,col_GMmsc);
            out_elong_msc_SOL_trial_max = MTU_elong_array(loc_frame,col_SOLmsc);
            out_strain_AT_trial_max = MTU_strain_array(loc_frame,col_AT);
            out_strain_GMtend_trial_max = MTU_strain_array(loc_frame,col_GMtend);
            out_strain_leg_trial_max = MTU_strain_array(loc_frame,col_leg);
            out_strain_GMapo_trial_max = MTU_strain_array(loc_frame,col_GMapo);
            out_strain_msc_GM_trial_max = MTU_strain_array(loc_frame,col_GMmsc);
            out_strain_msc_SOL_trial_max = MTU_strain_array(loc_frame,col_SOLmsc);
            out_length_AT_trial_max = MTU_length_array(loc_frame,col_AT);
            out_length_GMtend_trial_max = MTU_length_array(loc_frame,col_GMtend);
            out_length_leg_trial_max = MTU_length_array(loc_frame,col_leg);
            out_length_GMapo_trial_max = MTU_length_array(loc_frame,col_GMapo);
            out_length_msc_GM_trial_max = MTU_length_array(loc_frame,col_GMmsc);
            out_length_msc_SOL_trial_max = MTU_length_array(loc_frame,col_SOLmsc);
            
            out_elong_SEE_Fuku_trial_max = MTU_elong_array(loc_frame,col_SEE_Fukunaga);
            out_elong_msc_GM_Fuku_trial_max = MTU_elong_array(loc_frame,col_GMmsc_Fukunaga);
            out_strain_SEE_Fuku_trial_max = MTU_strain_array(loc_frame,col_SEE_Fukunaga);
            out_strain_msc_GM_Fuku_trial_max = MTU_strain_array(loc_frame,col_GMmsc_Fukunaga);
            out_length_SEE_Fuku_trial_max = MTU_length_array(loc_frame,col_SEE_Fukunaga);
            out_length_msc_GM_Fuku_trial_max = MTU_length_array(loc_frame,col_GMmsc_Fukunaga);
            
                % normalized versions of above:
            out_norm_length_leg_trial_max = MTU_normalized(loc_frame,1);
            out_norm_length_msc_GM_Fuku_trial_max = MTU_normalized(loc_frame,2);
            out_norm_length_SEE_Fuku_trial_max = MTU_normalized(loc_frame,3);
            out_norm_elong_leg_trial_max = MTU_normalized(loc_frame,4);
            out_norm_elong_msc_GM_Fuku_trial_max = MTU_normalized(loc_frame,5);
            out_norm_elong_SEE_Fuku_trial_max = MTU_normalized(loc_frame,6);
            out_norm_elong_percent_msc_GM_Fuku_trial_max = MTU_normalized(loc_frame,7);
            out_norm_elong_percent_SEE_Fuku_trial_max = MTU_normalized(loc_frame,8);
            
            % elongations, strains @ ind max
            loc_frame = find(MTU_elong_array(:,col_angle)>=out_ROM_ind_max,1,'first'); 
            out_elong_AT_ind_max = MTU_elong_array(loc_frame,col_AT); 
            out_elong_GMtend_ind_max = MTU_elong_array(loc_frame,col_GMtend);
            out_elong_leg_ind_max = MTU_elong_array(loc_frame,col_leg);
            out_elong_GMapo_ind_max = MTU_elong_array(loc_frame,col_GMapo); 
            out_elong_msc_GM_ind_max = MTU_elong_array(loc_frame,col_GMmsc); 
            out_elong_msc_SOL_ind_max = MTU_elong_array(loc_frame,col_SOLmsc); 
            out_strain_AT_ind_max = MTU_strain_array(loc_frame,col_AT); 
            out_strain_GMtend_ind_max = MTU_strain_array(loc_frame,col_GMtend);
            out_strain_leg_ind_max = MTU_strain_array(loc_frame,col_leg);
            out_strain_GMapo_ind_max = MTU_strain_array(loc_frame,col_GMapo); 
            out_strain_msc_GM_ind_max = MTU_strain_array(loc_frame,col_GMmsc); 
            out_strain_msc_SOL_ind_max = MTU_strain_array(loc_frame,col_SOLmsc); 
            out_length_AT_ind_max = MTU_length_array(loc_frame,col_AT); 
            out_length_GMtend_ind_max = MTU_length_array(loc_frame,col_GMtend);
            out_length_leg_ind_max = MTU_length_array(loc_frame,col_leg);
            out_length_GMapo_ind_max = MTU_length_array(loc_frame,col_GMapo); 
            out_length_msc_GM_ind_max = MTU_length_array(loc_frame,col_GMmsc); 
            out_length_msc_SOL_ind_max = MTU_length_array(loc_frame,col_SOLmsc); 
            out_elong_SEE_Fuku_ind_max = MTU_elong_array(loc_frame,col_SEE_Fukunaga);
            out_elong_msc_GM_Fuku_ind_max = MTU_elong_array(loc_frame,col_GMmsc_Fukunaga);
            out_strain_SEE_Fuku_ind_max = MTU_strain_array(loc_frame,col_SEE_Fukunaga);
            out_strain_msc_GM_Fuku_ind_max = MTU_strain_array(loc_frame,col_GMmsc_Fukunaga);
            out_length_SEE_Fuku_ind_max = MTU_length_array(loc_frame,col_SEE_Fukunaga);
            out_length_msc_GM_Fuku_ind_max = MTU_length_array(loc_frame,col_GMmsc_Fukunaga);
            
                % normalized versions of above:
            out_norm_length_leg_ind_max = MTU_normalized(loc_frame,1);
            out_norm_length_msc_GM_Fuku_ind_max = MTU_normalized(loc_frame,2);
            out_norm_length_SEE_Fuku_ind_max = MTU_normalized(loc_frame,3);
            out_norm_elong_leg_ind_max = MTU_normalized(loc_frame,4);
            out_norm_elong_msc_GM_Fuku_ind_max = MTU_normalized(loc_frame,5);
            out_norm_elong_SEE_Fuku_ind_max = MTU_normalized(loc_frame,6);
            out_norm_elong_percent_msc_GM_Fuku_ind_max = MTU_normalized(loc_frame,7);
            out_norm_elong_percent_SEE_Fuku_ind_max = MTU_normalized(loc_frame,8);

            % elongations, strains @ common max
            loc_frame = find(MTU_elong_array(:,col_angle)>=out_ROM_common_max,1,'first'); 
            out_elong_AT_common_max = MTU_elong_array(loc_frame,col_AT); 
            out_elong_GMtend_common_max = MTU_elong_array(loc_frame,col_GMtend);
            out_elong_leg_common_max = MTU_elong_array(loc_frame,col_leg);
            out_elong_GMapo_common_max = MTU_elong_array(loc_frame,col_GMapo); 
            out_elong_msc_GM_common_max = MTU_elong_array(loc_frame,col_GMmsc); 
            out_elong_msc_SOL_common_max = MTU_elong_array(loc_frame,col_SOLmsc); 
            out_strain_AT_common_max = MTU_strain_array(loc_frame,col_AT); 
            out_strain_GMtend_common_max = MTU_strain_array(loc_frame,col_GMtend);
            out_strain_leg_common_max = MTU_strain_array(loc_frame,col_leg);
            out_strain_GMapo_common_max = MTU_strain_array(loc_frame,col_GMapo); 
            out_strain_msc_GM_common_max = MTU_strain_array(loc_frame,col_GMmsc); 
            out_strain_msc_SOL_common_max = MTU_strain_array(loc_frame,col_SOLmsc); 
            out_length_AT_common_max = MTU_length_array(loc_frame,col_AT); 
            out_length_GMtend_common_max = MTU_length_array(loc_frame,col_GMtend);
            out_length_leg_common_max = MTU_length_array(loc_frame,col_leg);
            out_length_GMapo_common_max = MTU_length_array(loc_frame,col_GMapo); 
            out_length_msc_GM_common_max = MTU_length_array(loc_frame,col_GMmsc); 
            out_length_msc_SOL_common_max = MTU_length_array(loc_frame,col_SOLmsc); 
            out_elong_SEE_Fuku_common_max = MTU_elong_array(loc_frame,col_SEE_Fukunaga);
            out_elong_msc_GM_Fuku_common_max = MTU_elong_array(loc_frame,col_GMmsc_Fukunaga);
            out_strain_SEE_Fuku_common_max = MTU_strain_array(loc_frame,col_SEE_Fukunaga);
            out_strain_msc_GM_Fuku_common_max = MTU_strain_array(loc_frame,col_GMmsc_Fukunaga);
            out_length_SEE_Fuku_common_max = MTU_length_array(loc_frame,col_SEE_Fukunaga);
            out_length_msc_GM_Fuku_common_max = MTU_length_array(loc_frame,col_GMmsc_Fukunaga);
            
                % normalized versions of above:
            out_norm_length_leg_common_max = MTU_normalized(loc_frame,1);
            out_norm_length_msc_GM_Fuku_common_max = MTU_normalized(loc_frame,2);
            out_norm_length_SEE_Fuku_common_max = MTU_normalized(loc_frame,3);
            out_norm_elong_leg_common_max = MTU_normalized(loc_frame,4);
            out_norm_elong_msc_GM_Fuku_common_max = MTU_normalized(loc_frame,5);
            out_norm_elong_SEE_Fuku_common_max = MTU_normalized(loc_frame,6);
            out_norm_elong_percent_msc_GM_Fuku_common_max = MTU_normalized(loc_frame,7);
            out_norm_elong_percent_SEE_Fuku_common_max = MTU_normalized(loc_frame,8);

            % elongations, strains @ submax1
            loc_frame = find(MTU_elong_array(:,col_angle)>=out_ROM_submax_1,1,'first'); 
            if isempty(loc_frame)
                out_elong_AT_submax_1 = NaN;
                out_elong_GMtend_submax_1 = NaN;
                out_elong_leg_submax_1 = NaN;
                out_elong_GMapo_submax_1 = NaN;
                out_elong_msc_GM_submax_1 = NaN;
                out_elong_msc_SOL_submax_1 = NaN;
                out_strain_AT_submax_1 = NaN;
                out_strain_GMtend_submax_1 = NaN;
                out_strain_leg_submax_1 = NaN;
                out_strain_GMapo_submax_1 = NaN;
                out_strain_msc_GM_submax_1 = NaN;
                out_strain_msc_SOL_submax_1 = NaN;
                out_length_AT_submax_1 = NaN;
                out_length_GMtend_submax_1 = NaN;
                out_length_leg_submax_1 = NaN;
                out_length_GMapo_submax_1 = NaN;
                out_length_msc_GM_submax_1 = NaN;
                out_length_msc_SOL_submax_1 = NaN;
                out_elong_SEE_Fuku_submax_1 = NaN;
                out_elong_msc_GM_Fuku_submax_1 = NaN;
                out_strain_SEE_Fuku_submax_1 = NaN;
                out_strain_msc_GM_Fuku_submax_1 = NaN;
                out_length_SEE_Fuku_submax_1 = NaN;
                out_length_msc_GM_Fuku_submax_1 = NaN;
                out_norm_length_leg_submax_1 = NaN;
                out_norm_length_msc_GM_Fuku_submax_1 = NaN;
                out_norm_length_SEE_Fuku_submax_1 = NaN;
                out_norm_elong_leg_submax_1 = NaN;
                out_norm_elong_msc_GM_Fuku_submax_1 = NaN;
                out_norm_elong_SEE_Fuku_submax_1 = NaN;
                out_norm_elong_percent_msc_GM_Fuku_submax_1 = NaN;
                out_norm_elong_percent_SEE_Fuku_submax_1 = NaN;

            else
                out_elong_AT_submax_1 = MTU_elong_array(loc_frame,col_AT); 
                out_elong_GMtend_submax_1 = MTU_elong_array(loc_frame,col_GMtend);
                out_elong_leg_submax_1 = MTU_elong_array(loc_frame,col_leg);
                out_elong_GMapo_submax_1 = MTU_elong_array(loc_frame,col_GMapo); 
                out_elong_msc_GM_submax_1 = MTU_elong_array(loc_frame,col_GMmsc); 
                out_elong_msc_SOL_submax_1 = MTU_elong_array(loc_frame,col_SOLmsc); 
                out_strain_AT_submax_1 = MTU_strain_array(loc_frame,col_AT); 
                out_strain_GMtend_submax_1 = MTU_strain_array(loc_frame,col_GMtend);
                out_strain_leg_submax_1 = MTU_strain_array(loc_frame,col_leg);
                out_strain_GMapo_submax_1 = MTU_strain_array(loc_frame,col_GMapo); 
                out_strain_msc_GM_submax_1 = MTU_strain_array(loc_frame,col_GMmsc); 
                out_strain_msc_SOL_submax_1 = MTU_strain_array(loc_frame,col_SOLmsc); 
                out_length_AT_submax_1 = MTU_length_array(loc_frame,col_AT); 
                out_length_GMtend_submax_1 = MTU_length_array(loc_frame,col_GMtend);
                out_length_leg_submax_1 = MTU_length_array(loc_frame,col_leg);
                out_length_GMapo_submax_1 = MTU_length_array(loc_frame,col_GMapo); 
                out_length_msc_GM_submax_1 = MTU_length_array(loc_frame,col_GMmsc); 
                out_length_msc_SOL_submax_1 = MTU_length_array(loc_frame,col_SOLmsc); 
                out_elong_SEE_Fuku_submax_1 = MTU_elong_array(loc_frame,col_SEE_Fukunaga);
                out_elong_msc_GM_Fuku_submax_1 = MTU_elong_array(loc_frame,col_GMmsc_Fukunaga);
                out_strain_SEE_Fuku_submax_1 = MTU_strain_array(loc_frame,col_SEE_Fukunaga);
                out_strain_msc_GM_Fuku_submax_1 = MTU_strain_array(loc_frame,col_GMmsc_Fukunaga);
                out_length_SEE_Fuku_submax_1 = MTU_length_array(loc_frame,col_SEE_Fukunaga);
                out_length_msc_GM_Fuku_submax_1 = MTU_length_array(loc_frame,col_GMmsc_Fukunaga);
                out_norm_length_leg_submax_1 = MTU_normalized(loc_frame,1);
                out_norm_length_msc_GM_Fuku_submax_1 = MTU_normalized(loc_frame,2);
                out_norm_length_SEE_Fuku_submax_1 = MTU_normalized(loc_frame,3);
                out_norm_elong_leg_submax_1 = MTU_normalized(loc_frame,4);
                out_norm_elong_msc_GM_Fuku_submax_1 = MTU_normalized(loc_frame,5);
                out_norm_elong_SEE_Fuku_submax_1 = MTU_normalized(loc_frame,6);
                out_norm_elong_percent_msc_GM_Fuku_submax_1 = MTU_normalized(loc_frame,7);
                out_norm_elong_percent_SEE_Fuku_submax_1 = MTU_normalized(loc_frame,8);
            end
                
            % elongations, strains @ submax2
            loc_frame = find(MTU_elong_array(:,col_angle)>=out_ROM_submax_2,1,'first');
            out_elong_AT_submax_2 = MTU_elong_array(loc_frame,col_AT); 
            out_elong_GMtend_submax_2 = MTU_elong_array(loc_frame,col_GMtend);
            out_elong_leg_submax_2 = MTU_elong_array(loc_frame,col_leg);
            out_elong_GMapo_submax_2 = MTU_elong_array(loc_frame,col_GMapo); 
            out_elong_msc_GM_submax_2 = MTU_elong_array(loc_frame,col_GMmsc); 
            out_elong_msc_SOL_submax_2 = MTU_elong_array(loc_frame,col_SOLmsc); 
            out_strain_AT_submax_2 = MTU_strain_array(loc_frame,col_AT); 
            out_strain_GMtend_submax_2 = MTU_strain_array(loc_frame,col_GMtend);
            out_strain_leg_submax_2 = MTU_strain_array(loc_frame,col_leg);
            out_strain_GMapo_submax_2 = MTU_strain_array(loc_frame,col_GMapo); 
            out_strain_msc_GM_submax_2 = MTU_strain_array(loc_frame,col_GMmsc); 
            out_strain_msc_SOL_submax_2 = MTU_strain_array(loc_frame,col_SOLmsc); 
            out_length_AT_submax_2 = MTU_length_array(loc_frame,col_AT); 
            out_length_GMtend_submax_2 = MTU_length_array(loc_frame,col_GMtend);
            out_length_leg_submax_2 = MTU_length_array(loc_frame,col_leg);
            out_length_GMapo_submax_2 = MTU_length_array(loc_frame,col_GMapo); 
            out_length_msc_GM_submax_2 = MTU_length_array(loc_frame,col_GMmsc); 
            out_length_msc_SOL_submax_2 = MTU_length_array(loc_frame,col_SOLmsc); 
            out_elong_SEE_Fuku_submax_2 = MTU_elong_array(loc_frame,col_SEE_Fukunaga);
            out_elong_msc_GM_Fuku_submax_2 = MTU_elong_array(loc_frame,col_GMmsc_Fukunaga);
            out_strain_SEE_Fuku_submax_2 = MTU_strain_array(loc_frame,col_SEE_Fukunaga);
            out_strain_msc_GM_Fuku_submax_2 = MTU_strain_array(loc_frame,col_GMmsc_Fukunaga);
            out_length_SEE_Fuku_submax_2 = MTU_length_array(loc_frame,col_SEE_Fukunaga);
            out_length_msc_GM_Fuku_submax_2 = MTU_length_array(loc_frame,col_GMmsc_Fukunaga);

                % normalized versions of above:
            out_norm_length_leg_submax_2 = MTU_normalized(loc_frame,1);
            out_norm_length_msc_GM_Fuku_submax_2 = MTU_normalized(loc_frame,2);
            out_norm_length_SEE_Fuku_submax_2 = MTU_normalized(loc_frame,3);
            out_norm_elong_leg_submax_2 = MTU_normalized(loc_frame,4);
            out_norm_elong_msc_GM_Fuku_submax_2 = MTU_normalized(loc_frame,5);
            out_norm_elong_SEE_Fuku_submax_2 = MTU_normalized(loc_frame,6);
            out_norm_elong_percent_msc_GM_Fuku_submax_2 = MTU_normalized(loc_frame,7);
            out_norm_elong_percent_SEE_Fuku_submax_2 = MTU_normalized(loc_frame,8);


            % EMG (from all 3 scan locations / 6 trials, averaged)

            emg_step = 9; %VAR - number of EMG values BEFORE relevant angle, to include in average. 9 values = a span of 0.5 degrees.
            % cannot average values AROUND relevant angle, because for a few (least flexible) trials, we want the data AT the last available angle.

            % identify locations of the relevant goniometer angles, in the array with all 2+2+2 trials averaged:
            loc_frame_trial_max = find(data_force_gonio(:,col_angle_DFG)>=out_ROM_trial_max,1,'first'); % when averaging angles across trials, intervals of 0.05 degrees are used. This means that any required angle will exist in all data series
            loc_frame_ind_max = find(data_force_gonio(:,col_angle_DFG)>=out_ROM_ind_max,1,'first'); 
            loc_frame_common_max = find(data_force_gonio(:,col_angle_DFG)>=out_ROM_common_max,1,'first'); 
            loc_frame_submax_1 = find(data_force_gonio(:,col_angle_DFG)>=out_ROM_submax_1,1,'first'); 
            loc_frame_submax_2 = find(data_force_gonio(:,col_angle_DFG)>=out_ROM_submax_2,1,'first'); 
            loc_frame_zero = 1; % find(data_force_gonio(:,col_angle)>=0,1,'first'); 

            % stop running if frame not found
            if isempty(loc_frame_ind_max)
                % cprintf('error', horzcat('ERROR: Computed trial max ROM (', num2str(out_ROM_trial_max), ') does not exist in the current data series (max = ', num2str(max(data_force_gonio(:,col_angle))), '). Check "create_angles_passive" vs current data.\n' ));
                error(horzcat('ERROR: Computed trial max ROM (', num2str(out_ROM_trial_max), ') does not exist in the current data series (max = ', num2str(max(data_force_gonio(:,col_angle_DFG))), '). Check "create_angles_passive" vs current data.\n' ))
            end

            % EMG gm
            out_emg_gm_trial_max = mean(data_force_gonio(loc_frame_trial_max-emg_step:loc_frame_trial_max,4)); % EMG gm = column 4
            out_emg_gm_ind_max = mean(data_force_gonio(loc_frame_ind_max-emg_step:loc_frame_ind_max,4));
            out_emg_gm_common_max = mean(data_force_gonio(loc_frame_common_max-emg_step:loc_frame_common_max,4));
            out_emg_gm_submax_2 = mean(data_force_gonio(loc_frame_submax_2-emg_step:loc_frame_submax_2,4));
            out_emg_gm_zero = data_force_gonio(1,4); % mean(data_force_gonio(loc_frame_zero:loc_frame_zero,4));
            out_emg_gm_max = max(data_force_gonio(:,4));

            % EMG gl
            out_emg_gl_trial_max = mean(data_force_gonio(loc_frame_trial_max-emg_step:loc_frame_trial_max,5)); % EMG gl = column 5
            out_emg_gl_ind_max = mean(data_force_gonio(loc_frame_ind_max-emg_step:loc_frame_ind_max,5)); 
            out_emg_gl_common_max = mean(data_force_gonio(loc_frame_common_max-emg_step:loc_frame_common_max,5));
            out_emg_gl_submax_2 = mean(data_force_gonio(loc_frame_submax_2-emg_step:loc_frame_submax_2,5));
            out_emg_gl_zero = data_force_gonio(1,5); % mean(data_force_gonio(loc_frame_zero:loc_frame_zero,5));
            out_emg_gl_max = max(data_force_gonio(:,5));

            % EMG sol
            out_emg_sol_trial_max = mean(data_force_gonio(loc_frame_trial_max-emg_step:loc_frame_trial_max,6)); % EMG sol = column 6
            out_emg_sol_ind_max = mean(data_force_gonio(loc_frame_ind_max-emg_step:loc_frame_ind_max,6));
            out_emg_sol_common_max = mean(data_force_gonio(loc_frame_common_max-emg_step:loc_frame_common_max,6));
            out_emg_sol_submax_2 = mean(data_force_gonio(loc_frame_submax_2-emg_step:loc_frame_submax_2,6));
            out_emg_sol_zero = data_force_gonio(1,6); % mean(data_force_gonio(emg_step:loc_frame_zero,6));
            out_emg_sol_max = max(data_force_gonio(:,6));
            
            % EMG mean 3 msc
            out_emg_mean_trial_max = mean(data_force_gonio(loc_frame_trial_max-emg_step:loc_frame_trial_max,7)); % EMG mean = column 7
            out_emg_mean_ind_max = mean(data_force_gonio(loc_frame_ind_max-emg_step:loc_frame_ind_max,7));
            out_emg_mean_common_max = mean(data_force_gonio(loc_frame_common_max-emg_step:loc_frame_common_max,7));
            out_emg_mean_submax_2 = mean(data_force_gonio(loc_frame_submax_2-emg_step:loc_frame_submax_2,7));
            out_emg_mean_zero = data_force_gonio(1,7); % mean(data_force_gonio(emg_step:loc_frame_zero,7));
            out_emg_mean_max = max(data_force_gonio(:,7));
            
            if isempty(loc_frame_submax_1)
                out_emg_gm_submax_1 = NaN;
                out_emg_gl_submax_1 = NaN;
                out_emg_sol_submax_1 = NaN;
                out_emg_mean_submax_1 = NaN;
            else
                out_emg_gm_submax_1 = mean(data_force_gonio(loc_frame_submax_1-emg_step:loc_frame_submax_1,4));
                out_emg_gl_submax_1 = mean(data_force_gonio(loc_frame_submax_1-emg_step:loc_frame_submax_1,5));
                out_emg_sol_submax_1 = mean(data_force_gonio(loc_frame_submax_1-emg_step:loc_frame_submax_1,6));
                out_emg_mean_submax_1 = mean(data_force_gonio(loc_frame_submax_1-emg_step:loc_frame_submax_1,7));
            end


            %% Extract FASCICLE LENGTH, PENNATION ANGLE, GM MSC ELONG @ various joint angles
            %     data_GMFAS_licht_GM
            %     data_GMFAS_licht_SOL
            % containing:
            %   averaged angle (currently calculated from gonio)
            %   averaged fasicle length
            %   averaged pennation angle
            %   averaged fasicle elong
            %   averaged fasicle strain
            % OR containing zeros (if nonexistent)
            %    data_GMFAS_elong_Lichtwark
            % containing:
            %  angle
            %  GM muscle elongation (Lichtwark/Fukunaga)

            % set columns:
            col_licht_angle = 1;   % in data_GMFAS_licht_x
            col_licht_faslen = 2;
            col_licht_penn = 3;
            col_licht_fas_elong = 4;
            col_licht_fas_strain = 5;

            % preallocate / set zero values:

            % --- GM pennation and fascicle length
            out_licht_fas_length_GM_trial_max = NaN;
            out_licht_pennation_GM_trial_max = NaN;
            out_licht_fas_length_GM_common_max = NaN;
            out_licht_pennation_GM_common_max = NaN;
            out_licht_fas_length_GM_ind_max = NaN;
            out_licht_pennation_GM_ind_max = NaN;
            out_licht_fas_length_GM_zero = NaN;
            out_licht_pennation_GM_zero = NaN;
            out_licht_fas_length_GM_submax_1 = NaN;
            out_licht_pennation_GM_submax_1 = NaN;
            out_licht_fas_length_GM_submax_2 = NaN;
            out_licht_pennation_GM_submax_2 = NaN;
            
            % --- GM elongation and strain
            out_licht_fas_elong_GM_trial_max = NaN;
            out_licht_fas_strain_GM_trial_max = NaN;
            out_licht_fas_elong_GM_common_max = NaN;
            out_licht_fas_strain_GM_common_max = NaN;
            out_licht_fas_elong_GM_ind_max = NaN;
            out_licht_fas_strain_GM_ind_max = NaN;
            out_licht_fas_elong_GM_zero = NaN;
            out_licht_fas_strain_GM_zero = NaN;
            out_licht_fas_elong_GM_submax_1 = NaN;
            out_licht_fas_strain_GM_submax_1 = NaN;
            out_licht_fas_elong_GM_submax_2 = NaN;
            out_licht_fas_strain_GM_submax_2 = NaN;

            % --- GM length and elong normalized
            out_norm_licht_fas_elong_GM_trial_max = NaN;
            out_norm_licht_fas_length_GM_trial_max = NaN;
            out_norm_licht_fas_elong_GM_common_max = NaN;
            out_norm_licht_fas_length_GM_common_max = NaN;
            out_norm_licht_fas_elong_GM_ind_max = NaN;
            out_norm_licht_fas_length_GM_ind_max = NaN;
            out_norm_licht_fas_elong_GM_zero = NaN;
            out_norm_licht_fas_length_GM_zero = NaN;
            out_norm_licht_fas_elong_GM_submax_1 = NaN;
            out_norm_licht_fas_length_GM_submax_1 = NaN;
            out_norm_licht_fas_elong_GM_submax_2 = NaN;
            out_norm_licht_fas_length_GM_submax_2 = NaN;
            
            % --- SOL pennation and fascicle length
            out_licht_fas_length_SOL_trial_max = NaN;
            out_licht_pennation_SOL_trial_max = NaN;
            out_licht_fas_length_SOL_common_max = NaN;
            out_licht_pennation_SOL_common_max = NaN;
            out_licht_fas_length_SOL_ind_max = NaN;
            out_licht_pennation_SOL_ind_max = NaN;
            out_licht_fas_length_SOL_zero = NaN;
            out_licht_pennation_SOL_zero = NaN;
            out_licht_fas_length_SOL_submax_1 = NaN;
            out_licht_pennation_SOL_submax_1 = NaN;
            out_licht_fas_length_SOL_submax_2 = NaN;
            out_licht_pennation_SOL_submax_2 = NaN;


            if data_GMFAS_licht_GM == 0
                % no licht data existing
            else
                % at trial max angle:
                loc_frame = find(data_GMFAS(:,col_angle_DFG)>=out_ROM_trial_max,1,'first'); 
                out_licht_fas_length_GM_trial_max = data_GMFAS_licht_GM(loc_frame,col_licht_faslen);
                out_licht_pennation_GM_trial_max = data_GMFAS_licht_GM(loc_frame,col_licht_penn);
                out_licht_fas_elong_GM_trial_max = data_GMFAS_licht_GM(loc_frame,col_licht_fas_elong);
                out_licht_fas_strain_GM_trial_max = data_GMFAS_licht_GM(loc_frame,col_licht_fas_strain);
                out_norm_licht_fas_elong_GM_trial_max = MTU_normalized_licht(loc_frame,3);
                out_norm_licht_fas_length_GM_trial_max = MTU_normalized_licht(loc_frame,2);
                if (length((data_GMFAS_licht_SOL)) == 3) == 0 && loc_frame <= length(data_GMFAS_licht_SOL)
                    % SOL exists
                    out_licht_fas_length_SOL_trial_max = data_GMFAS_licht_SOL(loc_frame,col_licht_faslen); 
                    out_licht_pennation_SOL_trial_max = data_GMFAS_licht_SOL(loc_frame,col_licht_penn); 
                end

                % at individual max angle (across sides/timepoints):
                loc_frame = find(data_GMFAS(:,col_angle_DFG)>=out_ROM_ind_max,1,'first'); 
                out_licht_fas_length_GM_ind_max = data_GMFAS_licht_GM(loc_frame,col_licht_faslen); 
                out_licht_pennation_GM_ind_max = data_GMFAS_licht_GM(loc_frame,col_licht_penn); 
                out_licht_fas_elong_GM_ind_max = data_GMFAS_licht_GM(loc_frame,col_licht_fas_elong);
                out_licht_fas_strain_GM_ind_max = data_GMFAS_licht_GM(loc_frame,col_licht_fas_strain);
                out_norm_licht_fas_elong_GM_ind_max = MTU_normalized_licht(loc_frame,3);
                out_norm_licht_fas_length_GM_ind_max = MTU_normalized_licht(loc_frame,2);
                if (length((data_GMFAS_licht_SOL)) == 3) == 0 && loc_frame <= length(data_GMFAS_licht_SOL)
                    % SOL exists
                    out_licht_fas_length_SOL_ind_max = data_GMFAS_licht_SOL(loc_frame,col_licht_faslen); 
                    out_licht_pennation_SOL_ind_max = data_GMFAS_licht_SOL(loc_frame,col_licht_penn); 
                end

                % at subject common max angle:
                loc_frame = find(data_GMFAS(:,col_angle_DFG)>=out_ROM_common_max,1,'first'); 
                out_licht_fas_length_GM_common_max = data_GMFAS_licht_GM(loc_frame,col_licht_faslen); 
                out_licht_pennation_GM_common_max = data_GMFAS_licht_GM(loc_frame,col_licht_penn); 
                out_licht_fas_elong_GM_common_max = data_GMFAS_licht_GM(loc_frame,col_licht_fas_elong);
                out_licht_fas_strain_GM_common_max = data_GMFAS_licht_GM(loc_frame,col_licht_fas_strain);
                out_norm_licht_fas_elong_GM_common_max = MTU_normalized_licht(loc_frame,3);
                out_norm_licht_fas_length_GM_common_max = MTU_normalized_licht(loc_frame,2);
                if (length((data_GMFAS_licht_SOL)) == 3) == 0 && loc_frame <= length(data_GMFAS_licht_SOL)
                    % SOL exists
                    out_licht_fas_length_SOL_common_max = data_GMFAS_licht_SOL(loc_frame,col_licht_faslen); 
                    out_licht_pennation_SOL_common_max = data_GMFAS_licht_SOL(loc_frame,col_licht_penn); 
                end

                % at zero angle:
                loc_frame = find(data_GMFAS(:,2)>=0,1,'first'); 
                out_licht_fas_length_GM_zero = data_GMFAS_licht_GM(loc_frame,col_licht_faslen); 
                out_licht_pennation_GM_zero = data_GMFAS_licht_GM(loc_frame,col_licht_penn); 
                out_licht_fas_elong_GM_zero = data_GMFAS_licht_GM(loc_frame,col_licht_fas_elong);
                out_licht_fas_strain_GM_zero = data_GMFAS_licht_GM(loc_frame,col_licht_fas_strain);
                out_norm_licht_fas_elong_GM_zero = MTU_normalized_licht(loc_frame,3);
                out_norm_licht_fas_length_GM_zero = MTU_normalized_licht(loc_frame,2);
                if (length((data_GMFAS_licht_SOL)) == 3) == 0
                    % SOL exists
                    out_licht_fas_length_SOL_zero = data_GMFAS_licht_SOL(loc_frame,col_licht_faslen); 
                    out_licht_pennation_SOL_zero = data_GMFAS_licht_SOL(loc_frame,col_licht_penn); 
                end

                % at submax_1 angle:
                loc_frame = find(data_GMFAS(:,col_angle_DFG)>=out_ROM_submax_1,1,'first'); 
                % does angle exist in data?
                if isempty(loc_frame)
                    %
                elseif loc_frame <= length(data_GMFAS_licht_GM)
                    % GM data exist:
                    out_licht_fas_length_GM_submax_1 = data_GMFAS_licht_GM(loc_frame,col_licht_faslen); 
                    out_licht_pennation_GM_submax_1 = data_GMFAS_licht_GM(loc_frame,col_licht_penn); 
                    out_licht_fas_elong_GM_submax_1 = data_GMFAS_licht_GM(loc_frame,col_licht_fas_elong);
                    out_licht_fas_strain_GM_submax_1 = data_GMFAS_licht_GM(loc_frame,col_licht_fas_strain);
                    out_norm_licht_fas_elong_GM_submax_1 = MTU_normalized_licht(loc_frame,3);
                    out_norm_licht_fas_length_GM_submax_1 = MTU_normalized_licht(loc_frame,2);
                    if (length((data_GMFAS_licht_SOL)) == 3) == 0 && loc_frame <= length(data_GMFAS_licht_SOL)
                        % SOL exists
                        out_licht_fas_length_SOL_submax_1 = data_GMFAS_licht_SOL(loc_frame,col_licht_faslen); 
                        out_licht_pennation_SOL_submax_1 = data_GMFAS_licht_SOL(loc_frame,col_licht_penn); 
                    end
                end

                % at submax_2 angle:
                loc_frame = find(data_GMFAS(:,col_angle_DFG)>=out_ROM_submax_2,1,'first'); 
                out_licht_fas_length_GM_submax_2 = data_GMFAS_licht_GM(loc_frame,col_licht_faslen); 
                out_licht_pennation_GM_submax_2 = data_GMFAS_licht_GM(loc_frame,col_licht_penn); 
                out_licht_fas_elong_GM_submax_2 = data_GMFAS_licht_GM(loc_frame,col_licht_fas_elong);
                out_licht_fas_strain_GM_submax_2 = data_GMFAS_licht_GM(loc_frame,col_licht_fas_strain);
                out_norm_licht_fas_elong_GM_submax_2 = MTU_normalized_licht(loc_frame,3);
                out_norm_licht_fas_length_GM_submax_2 = MTU_normalized_licht(loc_frame,2);
                if (length((data_GMFAS_licht_SOL)) == 3) == 0 && loc_frame <= length(data_GMFAS_licht_SOL)
                    % SOL exists
                    out_licht_fas_length_SOL_submax_2 = data_GMFAS_licht_SOL(loc_frame,col_licht_faslen); 
                    out_licht_pennation_SOL_submax_2 = data_GMFAS_licht_SOL(loc_frame,col_licht_penn); 
                end
            end


            %% Extract MAXIMAL ELONGATIONS (not angle specific) 
            % find zero angle:
            loc_zero = find(data_GMFAS(:,2)>=0,1,'first'); 
            
            out_elong_AT_max = max(MTU_elong_array(loc_zero:end,col_AT)); 
            out_elong_GMtend_max = max(MTU_elong_array(loc_zero:end,col_GMtend));
            out_elong_leg_max = max(MTU_elong_array(loc_zero:end,col_leg));
            out_elong_GMapo_max = max(MTU_elong_array(loc_zero:end,col_GMapo)); 
            out_elong_msc_GM_max = max(MTU_elong_array(loc_zero:end,col_GMmsc)); 
            out_elong_msc_SOL_max = max(MTU_elong_array(loc_zero:end,col_SOLmsc)); 
            out_strain_AT_max = max(MTU_strain_array(loc_zero:end,col_AT)); 
            out_strain_GMtend_max = max(MTU_strain_array(loc_zero:end,col_GMtend));
            out_strain_leg_max = max(MTU_strain_array(loc_zero:end,col_leg));
            out_strain_GMapo_max = max(MTU_strain_array(loc_zero:end,col_GMapo)); 
            out_strain_msc_GM_max = max(MTU_strain_array(loc_zero:end,col_GMmsc)); 
            out_strain_msc_SOL_max = max(MTU_strain_array(loc_zero:end,col_SOLmsc)); 
            out_length_AT_max = max(MTU_length_array(loc_zero:end,col_AT)); 
            out_length_GMtend_max = max(MTU_length_array(loc_zero:end,col_GMtend));
            out_length_leg_max = max(MTU_length_array(loc_zero:end,col_leg));
            out_length_GMapo_max = max(MTU_length_array(loc_zero:end,col_GMapo)); 
            out_length_msc_GM_max = max(MTU_length_array(loc_zero:end,col_GMmsc)); 
            out_length_msc_SOL_max = max(MTU_length_array(loc_zero:end,col_SOLmsc)); 
            
            out_elong_SEE_Fuku_max = max(MTU_elong_array(loc_zero:end,col_SEE_Fukunaga));
            out_elong_msc_GM_Fuku_max = max(MTU_elong_array(loc_zero:end,col_GMmsc_Fukunaga));
            out_strain_SEE_Fuku_max = max(MTU_strain_array(loc_zero:end,col_SEE_Fukunaga));
            out_strain_msc_GM_Fuku_max = max(MTU_strain_array(loc_zero:end,col_GMmsc_Fukunaga));
            out_length_SEE_Fuku_max = max(MTU_length_array(loc_zero:end,col_SEE_Fukunaga));
            out_length_msc_GM_Fuku_max = max(MTU_length_array(loc_zero:end,col_GMmsc_Fukunaga));
            if (length((data_GMFAS_licht_GM)) == 5) == 0
                out_licht_fas_length_GM_max = max(data_GMFAS_licht_GM(loc_zero:end,col_licht_faslen));
                out_licht_pennation_GM_max = max(data_GMFAS_licht_GM(loc_zero:end,col_licht_penn));
                out_licht_fas_elong_GM_max = max(data_GMFAS_licht_GM(loc_zero:end,col_licht_fas_elong));
                out_licht_fas_strain_GM_max = max(data_GMFAS_licht_GM(loc_zero:end,col_licht_fas_strain));
            else
                out_licht_fas_length_GM_max = NaN;
                out_licht_pennation_GM_max = NaN;
                out_licht_fas_elong_GM_max = NaN;
                out_licht_fas_strain_GM_max = NaN;
            end
            if (length((data_GMFAS_licht_SOL)) == 3) == 0
                % SOL exists
                out_licht_fas_length_SOL_max = max(data_GMFAS_licht_SOL(loc_zero:end,col_licht_faslen));
                out_licht_pennation_SOL_max = max(data_GMFAS_licht_SOL(loc_zero:end,col_licht_penn));
            else
                out_licht_fas_length_SOL_max = NaN;
                out_licht_pennation_SOL_max = NaN;
            end
            
            out_norm_length_leg_max = max(MTU_normalized(loc_zero:end,1));
            out_norm_length_msc_GM_Fuku_max = max(MTU_normalized(loc_zero:end,2));
            out_norm_length_SEE_Fuku_max = max(MTU_normalized(loc_zero:end,3));
            out_norm_elong_leg_max = max(MTU_normalized(loc_zero:end,4));
            out_norm_elong_msc_GM_Fuku_max = max(MTU_normalized(loc_zero:end,5));
            out_norm_elong_SEE_Fuku_max = max(MTU_normalized(loc_zero:end,6));
            out_norm_elong_percent_msc_GM_Fuku_max = max(MTU_normalized(loc_zero:end,7));
            out_norm_elong_percent_SEE_Fuku_max = max(MTU_normalized(loc_zero:end,8));
            if (length((MTU_normalized_licht)) == 3) == 0
                out_norm_licht_fas_elong_GM_max = max(MTU_normalized_licht(loc_zero:end,3));
                out_norm_licht_fas_length_GM_max = max(MTU_normalized_licht(loc_zero:end,2));
            else
                out_norm_licht_fas_elong_GM_max = NaN;
                out_norm_licht_fas_length_GM_max = NaN;
            end


            %% Calculate PASSIVE STIFFNESS and STIFFNESS INDEX (Nordez 2006)
            % gonio angle = data_force_gonio(:,col_angle)
            % torque = data_force_gonio(:,col_torque)

            %%% PASSIVE STIFFNESS:

            % passive stiffness = delta torque / delta angle, at various angles
            % fit 4th order polynomial to averaged torque-angle curve, using data from zero angle to TRIAL max ROM
            fit_ind_max = polyfit(data_force_gonio(loc_frame_zero:loc_frame_trial_max,col_angle_DFG), data_force_gonio(loc_frame_zero:loc_frame_trial_max,col_torque), 4);

            % extract passive stiffness (derivate of 4th order poly) at:
            out_pstiff_trial_max = (4 * fit_ind_max(1) * out_ROM_trial_max^3) + (3 * fit_ind_max(2) * out_ROM_trial_max^2) + (2 * fit_ind_max(3) * out_ROM_trial_max) + fit_ind_max(4);
            out_pstiff_ind_max = (4 * fit_ind_max(1) * out_ROM_ind_max^3) + (3 * fit_ind_max(2) * out_ROM_ind_max^2) + (2 * fit_ind_max(3) * out_ROM_ind_max) + fit_ind_max(4);
            out_pstiff_common_max = (4 * fit_ind_max(1) * out_ROM_common_max^3) + (3 * fit_ind_max(2) * out_ROM_common_max^2) + (2 * fit_ind_max(3) * out_ROM_common_max) + fit_ind_max(4);
            % if out_ROM_trial_max > out_ROM_submax_1 
            %    out_pstiff_submax_1 = (4 * fit_ind_max(1) * out_ROM_submax_1^3) + (3 * fit_ind_max(2) * out_ROM_submax_1^2) + (2 * fit_ind_max(3) * out_ROM_submax_1) + fit_ind_max(4);
            %else
            %    out_pstiff_submax_1 = NaN;
            %end
            out_pstiff_submax_1 = (4 * fit_ind_max(1) * out_ROM_submax_1^3) + (3 * fit_ind_max(2) * out_ROM_submax_1^2) + (2 * fit_ind_max(3) * out_ROM_submax_1) + fit_ind_max(4);
            out_pstiff_submax_2 = (4 * fit_ind_max(1) * out_ROM_submax_2^3) + (3 * fit_ind_max(2) * out_ROM_submax_2^2) + (2 * fit_ind_max(3) * out_ROM_submax_2) + fit_ind_max(4);
            out_pstiff_angle = 15; %VAR
            if out_ROM_trial_max > out_pstiff_angle
                out_pstiff_15 = (4 * fit_ind_max(1) * out_pstiff_angle^3) + (3 * fit_ind_max(2) * out_pstiff_angle^2) + (2 * fit_ind_max(3) * out_pstiff_angle) + fit_ind_max(4);
            else
                out_pstiff_15 = NaN;
            end
            
            %%% STIFFNESS INDEX:

            % fit 2nd order polynomial to averaged torque-angle curve, using data from zero angle to trial max ROM:
            fit_ind_max = polyfit(data_force_gonio(loc_frame_zero:loc_frame_trial_max,col_angle_DFG), data_force_gonio(loc_frame_zero:loc_frame_trial_max,col_torque), 2);

            % extract stiffness index as 2 * a
            out_pstiff_index = 2 * fit_ind_max(1);




            %% Extract ANGLES @ specific FORCE levels
            %   force = 6 trials with max force -> 1 lowest max force
            %   force levels = 
            %       ind max force
            %       common max force 
            %       ind R max force 
            %       ind L max force
            %   extracting gonio angle in degrees

            % data = input_for_pre_r 
            %        input_for_pre_l 
            %        input_for_post_r 
            %        input_for_post_l 
            %        input_for_ind_max 
            %        input_for_common_max 
            %        input_for_ind_rmax 
            %        input_for_ind_lmax

            % print error if forces do not exist --> create_angles_passive needs to be run
            if str2double(input_for_ind_max{trial_subjectno}) == 10000
                cprintf('*red', 'ERROR: Max ROM values are not calculated for current subject. Run create_angles_passive first.\n')
            end

            % load the predetermined forces to be searched for
            if strcmpi(dm_timepoint{line}, 'PRE') == 1
                if strcmpi(dm_side{line}, 'R') == 1
                    loc_F_trial_max = str2double(input_for_pre_r{trial_subjectno}) - 0.0001; % tweak for computer storage of numbers, where 8.3000 is stored as 8.2999999999999999
                    loc_F_ind_max = str2double(input_for_ind_rmax{trial_subjectno}) - 0.0001;
                else % L
                    loc_F_trial_max = str2double(input_for_pre_l{trial_subjectno}) - 0.0001;
                    loc_F_ind_max = str2double(input_for_ind_lmax{trial_subjectno}) - 0.0001;
                end
            else % POST
                if strcmpi(dm_side{line}, 'R') == 1
                    loc_F_trial_max = str2double(input_for_post_r{trial_subjectno}) - 0.0001;
                    loc_F_ind_max = str2double(input_for_ind_rmax{trial_subjectno}) - 0.0001;
                else % L
                    loc_F_trial_max = str2double(input_for_post_l{trial_subjectno}) - 0.0001;
                    loc_F_ind_max = str2double(input_for_ind_lmax{trial_subjectno}) - 0.0001;
                end
            end
            loc_F_common_max = str2double(input_for_common_max{trial_subjectno}) - 0.0001;
            loc_F_ind_rmax = str2double(input_for_ind_rmax{trial_subjectno}) - 0.0001; % will not be used, left in script for simplicity
            loc_F_ind_lmax = str2double(input_for_ind_lmax{trial_subjectno}) - 0.0001; % will not be used, left in script for simplicity

            % find goniometer angles (from all 3 scan locations / 6 trials, averaged)
            loc_frame = find(data_force_gonio(:,col_force)>=loc_F_trial_max,1,'first'); % when averaging angles across trials, intervals of 0.05 degrees are used. This means that any required angle will exist in all data series
            if isempty(loc_frame) % catch error of FORCE not found:
                loc_frame = find(data_force_gonio(:,col_force)>=max(data_force_gonio(:,col_force)),1,'first'); 
                cprintf('*red', 'ERROR: Recorded ind max FORCE not found in data. Run create_angles_passive?\n')
            end
            out_angle_trial_max = data_force_gonio(loc_frame,col_angle_DFG); 
            loc_frame = find(data_force_gonio(:,col_force)>=loc_F_ind_max,1,'first'); 
            if isempty(loc_frame) % catch error of FORCE not found:
                loc_frame = find(data_force_gonio(:,col_force)>=max(data_force_gonio(:,col_force)),1,'first'); 
                cprintf('*red', 'ERROR: Recorded ind max FORCE not found in data. Run create_angles_passive?\n')
            end
            out_angle_ind_max = data_force_gonio(loc_frame,col_angle_DFG);
            loc_frame = find(data_force_gonio(:,col_force)>=loc_F_common_max,1,'first');
            out_angle_common_max = data_force_gonio(loc_frame,col_angle_DFG);

            % the following will not be used - left in script for simplicity
            % for current trial (L or R), find the angle at the max force for
            % the other leg, if the other leg has a lower max force as lowest
            % PRE/POST
            if str2double(input_for_ind_rmax(trial_subjectno)) > 9000
                % data do not exist for the right side (array preloaded with "empty" values of 10000)
                out_angle_ind_rmax = NaN;
            else
                loc_frame = find(data_force_gonio(:,col_force)>=loc_F_ind_rmax,1,'first'); 
                if isempty(loc_frame) % force level does not exist
                    out_angle_ind_rmax = NaN;
                else
                    out_angle_ind_rmax = data_force_gonio(loc_frame,col_angle_DFG);
                end
            end

            if str2double(input_for_ind_lmax(trial_subjectno)) > 9000
                % data do not exist for the left side (array preloaded with "empty" values of 10000)
                out_angle_ind_lmax = NaN;
            else
                loc_frame = find(data_force_gonio(:,col_force)>=loc_F_ind_lmax,1,'first'); 
                if isempty(loc_frame) % force level does not exist
                    out_angle_ind_lmax = NaN;
                else
                    out_angle_ind_lmax = data_force_gonio(loc_frame,col_angle_DFG);
                end
            end


            %% Prepare ARRAYS angle_vars - group plots & stats

            % data_force_gonio
            % from M-file average_passive_forces_EMG
            % average of 6 trials
            % output array "data_force_gonio" contains: 
            %   average_force_gonio
            %   average_angle_array
            %   average_emg_gm_gonio
            %   average_emg_gl_gonio
            %	average_emg_sol_gonio

            % data_SOL / data_GMMTJ / data_GMFAS
            % from M-file average_passive_trials_EMG
            % average of 2 trials
            % output array "data_SOL" etc contains:
            %   average_force_gonio
            %   average_angle_array
            %   average_displ_gonio
            %   average_emg_gm_gonio
            %   average_emg_gl_gonio
            %	average_emg_sol_gonio

            % data_GMFAS_licht_GM
            % average of 2 trials
            % contains:
            %   averaged angle (currently calculated from gonio)
            %   averaged fasicle length
            %   averaged pennation angle
            %   fascicle elongation
            %   fascicle strain
            %   --- OR containing 5x zeros (if nonexistent)


            % extract angle range common to all trials for current subject
            angle_start = 0 - 0.00000001; % tweak for computer storage of numbers, where 8.3000 is stored as 8.2999999999999999 %VAR
            angle_stop = out_ROM_trial_max;

            % identify locations of start/stop angles in above mentioned arrays
            loc_angle_start = find(data_force_gonio(:,col_angle_DFG)>=angle_start,1,'first');
            loc_angle_stop = find(data_force_gonio(:,col_angle_DFG)>=angle_stop,1,'first');
            % check if licht data exist
            if length(data_GMFAS_licht_GM) > 5
                loc_angle_licht_start = find(data_GMFAS_licht_GM(:,col_licht_angle)>=angle_start,1,'first');
                loc_angle_licht_stop = find(data_GMFAS_licht_GM(:,col_licht_angle)>=angle_stop,1,'first');
            else % change tables to NaN
                loc_angle_licht_start = 1;
                loc_angle_licht_stop = loc_angle_stop-loc_angle_start+1;
                MTU_normalized_licht(loc_angle_licht_start:loc_angle_licht_stop,:) = NaN;
                data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,:) = NaN;
            end

            % contents of below angle_vars arrays /// angle_vars contain:
            col_AV_angle = 1;
            col_AV_F = 2;
            col_AV_T = 3;

            col_AV_EMG_gm = 4;
            col_AV_EMG_gl = 5;
            col_AV_EMG_sol = 6;
            
            col_AV_len_AT = 7;
            col_AV_len_GMtend = 8;
            col_AV_len_leg = 9;
            % 15 - no len GMFAS
            col_AV_len_GMapo = 11;
            col_AV_len_msc_GM = 12;
            col_AV_len_msc_SOL = 13;
            col_AV_len_msc_GM_fuku = 14;
            col_AV_len_SEE_fuku = 15;
            
            col_AV_elong_AT = 16;
            col_AV_elong_GMtend = 17;
            col_AV_elong_leg = 18;
            col_AV_elong_GMFAS = 19;
            col_AV_elong_GMapo = 20;
            col_AV_elong_msc_GM = 21;
            col_AV_elong_msc_SOL = 22;
            col_AV_elong_msc_GM_fuku = 23;
            col_AV_elong_SEE_fuku = 24;
            
            col_AV_strain_AT = 25;
            col_AV_strain_GMtend = 26;
            col_AV_strain_leg = 27;
            % 22 - no strain GMFAS
            col_AV_strain_GMapo = 29;
            col_AV_strain_msc_GM = 30;
            col_AV_strain_msc_SOL = 31;
            col_AV_strain_msc_GM_fuku = 32;
            col_AV_strain_SEE_fuku = 33;
            
            col_AV_len_GMfas_licht = 34;
            col_AV_elong_GMfas_licht = 35;
            col_AV_strain_GMfas_licht = 36;
            col_AV_pennation_GMfas_licht = 37;
            
            col_AV_norm_len_leg = 38;
            col_AV_norm_len_msc_GM_fuku = 39;
            col_AV_norm_len_SEE_fuku = 40;
            col_AV_norm_elong_leg = 41;
            col_AV_norm_elong_msc_GM_fuku = 42;
            col_AV_norm_elong_SEE_fuku = 43;
            col_AV_norm_percent_elong_msc_GM_fuku = 44;
            col_AV_norm_percent_elong_SEE_fuku = 45;
            
            col_AV_norm_len_GMfas_licht = 46; % MTU_normalized_licht
            col_AV_norm_elong_GMfas_licht = 47; % MTU_normalized_licht
            
            %col_AV_angle2 = 48; %(repeated - without normalization)
            col_AV_EMG_mean = 49;

            
            if trial_timepoint == 0 && trial_leg == 1 % PRE, STR
                %% STR PRE
                % all data in ONE cell, common angles, RAW data:
                STR_PRE_angle_vars{STR_PRE_count} = [ ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_angle_DFG) ...  1
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_force) ...  2
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_torque) ... % 18old
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_AV_EMG_gm) ...          3old
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_AV_EMG_gl)...           4old
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_AV_EMG_sol) ...          5old
                    ...
                    MTU_length_array(:,2) ...                                       12
                    MTU_length_array(:,3) ...                                       13
                    MTU_length_array(:,4) ...                                       14
                    MTU_length_array(:,5) ...                                       15
                    MTU_length_array(:,6) ...                                       16
                    MTU_length_array(:,7) ...                                       17
                    MTU_length_array(:,8) ...                                       20
                    MTU_length_array(:,9) ...                                       36
                    MTU_length_array(:,10) ...                                      37
                    ...
                    MTU_elong_array(:,2) ...                                        6
                    MTU_elong_array(:,3) ...                                        7
                    MTU_elong_array(:,4) ...                                        8
                    MTU_elong_array(:,5) ...                                        9
                    MTU_elong_array(:,6) ...                                        10
                    MTU_elong_array(:,7) ...                                        11
                    MTU_elong_array(:,8) ...                                        19
                    MTU_elong_array(:,9) ...                                        28
                    MTU_elong_array(:,10) ...                                       29
                    ...
                    MTU_strain_array(:,2) ...                                       21
                    MTU_strain_array(:,3) ...                                       22
                    MTU_strain_array(:,4) ...                                       23
                    MTU_strain_array(:,5) ...                                       24
                    MTU_strain_array(:,6) ...                                       25
                    MTU_strain_array(:,7) ...                                       26
                    MTU_strain_array(:,8) ...                                       27
                    MTU_strain_array(:,9) ...                                       30
                    MTU_strain_array(:,10) ...                                      31
                    ...
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_licht_faslen) ...   32
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_licht_fas_elong) ...   34
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_licht_fas_strain) ...   35
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_licht_penn) ...  33
                    ...
                    MTU_normalized(:,1) ...%leg length
                    MTU_normalized(:,2) ...%col_AV_norm_len_msc_GM_fuku = 100;
                    MTU_normalized(:,3) ...%col_AV_norm_len_SEE_fuku = 100;
                    MTU_normalized(:,4) ...%leg elong
                    MTU_normalized(:,5) ...%col_AV_norm_elong_msc_GM_fuku = 100;
                    MTU_normalized(:,6) ...%col_AV_norm_elong_SEE_fuku = 100;
                    MTU_normalized(:,7) ...%col_AV_norm_percent_elong_msc_GM_fuku = 100;
                    MTU_normalized(:,8) ...%col_AV_norm_percent_elong_SEE_fuku = 100;
                    ...
                    MTU_normalized_licht(loc_angle_licht_start:loc_angle_licht_stop,2) ... % norm fas len
                    MTU_normalized_licht(loc_angle_licht_start:loc_angle_licht_stop,3) ... % norm fas elong
                    ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_angle_DFG) ...   
                    ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,7) ... % mean EMG 3 msc
                    ];


            elseif trial_timepoint == 1 && trial_leg == 1 % POST, STR
                %% STR POST
                % all data in ONE cell, up to each subject's max angle, RAW data:
                STR_POST_angle_vars{STR_POST_count} = [ ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_angle_DFG) ...  1
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_force) ...  2
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_torque) ... % 18
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_AV_EMG_gm) ...          3
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_AV_EMG_gl)...           4
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_AV_EMG_sol) ...          5
                    ...
                    MTU_length_array(:,2) ...                                       12
                    MTU_length_array(:,3) ...                                       13
                    MTU_length_array(:,4) ...                                       14
                    MTU_length_array(:,5) ...                                       15
                    MTU_length_array(:,6) ...                                       16
                    MTU_length_array(:,7) ...                                       17
                    MTU_length_array(:,8) ...                                       20
                    MTU_length_array(:,9) ...                                       36
                    MTU_length_array(:,10) ...                                      37
                    ...
                    MTU_elong_array(:,2) ...                                        6
                    MTU_elong_array(:,3) ...                                        7
                    MTU_elong_array(:,4) ...                                        8
                    MTU_elong_array(:,5) ...                                        9
                    MTU_elong_array(:,6) ...                                        10
                    MTU_elong_array(:,7) ...                                        11
                    MTU_elong_array(:,8) ...                                        19
                    MTU_elong_array(:,9) ...                                        28
                    MTU_elong_array(:,10) ...                                       29
                    ...
                    MTU_strain_array(:,2) ...                                       21
                    MTU_strain_array(:,3) ...                                       22
                    MTU_strain_array(:,4) ...                                       23
                    MTU_strain_array(:,5) ...                                       24
                    MTU_strain_array(:,6) ...                                       25
                    MTU_strain_array(:,7) ...                                       26
                    MTU_strain_array(:,8) ...                                       27
                    MTU_strain_array(:,9) ...                                       30
                    MTU_strain_array(:,10) ...                                      31
                    ...
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_licht_faslen) ...   32
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_licht_fas_elong) ...   34
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_licht_fas_strain) ...   35
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_licht_penn) ...  33
                    ...
                    MTU_normalized(:,1) ...%leg length
                    MTU_normalized(:,2) ...%col_AV_norm_len_msc_GM_fuku = 100;
                    MTU_normalized(:,3) ...%col_AV_norm_len_SEE_fuku = 100;
                    MTU_normalized(:,4) ...%leg elong
                    MTU_normalized(:,5) ...%col_AV_norm_elong_msc_GM_fuku = 100;
                    MTU_normalized(:,6) ...%col_AV_norm_elong_SEE_fuku = 100;
                    MTU_normalized(:,7) ...%col_AV_norm_percent_elong_msc_GM_fuku = 100;
                    MTU_normalized(:,8) ...%col_AV_norm_percent_elong_SEE_fuku = 100;
                    ...
                    MTU_normalized_licht(loc_angle_licht_start:loc_angle_licht_stop,2) ... 
                    MTU_normalized_licht(loc_angle_licht_start:loc_angle_licht_stop,3) ... 
                    ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_angle_DFG) ...   
                    ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,7) ... % mean EMG 3 msc
                    ];


            elseif trial_timepoint == 0 && trial_leg == 0 % PRE, CON
                %% CON PRE
                % all data in ONE cell, up to each subject's max angle, RAW data:
                CON_PRE_angle_vars{CON_PRE_count} = [ ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_angle_DFG) ...  1
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_force) ...  2
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_torque) ... % 18
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_AV_EMG_gm) ...          3
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_AV_EMG_gl)...           4
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_AV_EMG_sol) ...          5
                    ...
                    MTU_length_array(:,2) ...                                       12
                    MTU_length_array(:,3) ...                                       13
                    MTU_length_array(:,4) ...                                       14
                    MTU_length_array(:,5) ...                                       15
                    MTU_length_array(:,6) ...                                       16
                    MTU_length_array(:,7) ...                                       17
                    MTU_length_array(:,8) ...                                       20
                    MTU_length_array(:,9) ...                                       36
                    MTU_length_array(:,10) ...                                      37
                    ...
                    MTU_elong_array(:,2) ...                                        6
                    MTU_elong_array(:,3) ...                                        7
                    MTU_elong_array(:,4) ...                                        8
                    MTU_elong_array(:,5) ...                                        9
                    MTU_elong_array(:,6) ...                                        10
                    MTU_elong_array(:,7) ...                                        11
                    MTU_elong_array(:,8) ...                                        19
                    MTU_elong_array(:,9) ...                                        28
                    MTU_elong_array(:,10) ...                                       29
                    ...
                    MTU_strain_array(:,2) ...                                       21
                    MTU_strain_array(:,3) ...                                       22
                    MTU_strain_array(:,4) ...                                       23
                    MTU_strain_array(:,5) ...                                       24
                    MTU_strain_array(:,6) ...                                       25
                    MTU_strain_array(:,7) ...                                       26
                    MTU_strain_array(:,8) ...                                       27
                    MTU_strain_array(:,9) ...                                       30
                    MTU_strain_array(:,10) ...                                      31
                    ...
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_licht_faslen) ...   32
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_licht_fas_elong) ...   34
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_licht_fas_strain) ...   35
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_licht_penn) ...  33
                    ...
                    MTU_normalized(:,1) ...%leg length
                    MTU_normalized(:,2) ...%col_AV_norm_len_msc_GM_fuku = 100;
                    MTU_normalized(:,3) ...%col_AV_norm_len_SEE_fuku = 100;
                    MTU_normalized(:,4) ...%leg elong
                    MTU_normalized(:,5) ...%col_AV_norm_elong_msc_GM_fuku = 100;
                    MTU_normalized(:,6) ...%col_AV_norm_elong_SEE_fuku = 100;
                    MTU_normalized(:,7) ...%col_AV_norm_percent_elong_msc_GM_fuku = 100;
                    MTU_normalized(:,8) ...%col_AV_norm_percent_elong_SEE_fuku = 100;
                    ...
                    MTU_normalized_licht(loc_angle_licht_start:loc_angle_licht_stop,2) ... 
                    MTU_normalized_licht(loc_angle_licht_start:loc_angle_licht_stop,3) ... 
                    ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_angle_DFG) ...   
                    ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,7) ... % mean EMG 3 msc
                    ];

            elseif trial_timepoint == 1 && trial_leg == 0 % POST, CON
                %% CON POST
                % all data in ONE cell, up to each subject's max angle, RAW data:
                CON_POST_angle_vars{CON_POST_count} = [ ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_angle_DFG) ...  1
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_force) ...  2
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_torque) ... % 18
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_AV_EMG_gm) ...          3
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_AV_EMG_gl)...           4
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_AV_EMG_sol) ...          5
                    ...
                    MTU_length_array(:,2) ...                                       12
                    MTU_length_array(:,3) ...                                       13
                    MTU_length_array(:,4) ...                                       14
                    MTU_length_array(:,5) ...                                       15
                    MTU_length_array(:,6) ...                                       16
                    MTU_length_array(:,7) ...                                       17
                    MTU_length_array(:,8) ...                                       20
                    MTU_length_array(:,9) ...                                       36
                    MTU_length_array(:,10) ...                                      37
                    ...
                    MTU_elong_array(:,2) ...                                        6
                    MTU_elong_array(:,3) ...                                        7
                    MTU_elong_array(:,4) ...                                        8
                    MTU_elong_array(:,5) ...                                        9
                    MTU_elong_array(:,6) ...                                        10
                    MTU_elong_array(:,7) ...                                        11
                    MTU_elong_array(:,8) ...                                        19
                    MTU_elong_array(:,9) ...                                        28
                    MTU_elong_array(:,10) ...                                       29
                    ...
                    MTU_strain_array(:,2) ...                                       21
                    MTU_strain_array(:,3) ...                                       22
                    MTU_strain_array(:,4) ...                                       23
                    MTU_strain_array(:,5) ...                                       24
                    MTU_strain_array(:,6) ...                                       25
                    MTU_strain_array(:,7) ...                                       26
                    MTU_strain_array(:,8) ...                                       27
                    MTU_strain_array(:,9) ...                                       30
                    MTU_strain_array(:,10) ...                                      31
                    ...
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_licht_faslen) ...   32
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_licht_fas_elong) ...   34
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_licht_fas_strain) ...   35
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_licht_penn) ...  33
                    ...
                    MTU_normalized(:,1) ...%leg length
                    MTU_normalized(:,2) ...%col_AV_norm_len_msc_GM_fuku = 100;
                    MTU_normalized(:,3) ...%col_AV_norm_len_SEE_fuku = 100;
                    MTU_normalized(:,4) ...%leg elong
                    MTU_normalized(:,5) ...%col_AV_norm_elong_msc_GM_fuku = 100;
                    MTU_normalized(:,6) ...%col_AV_norm_elong_SEE_fuku = 100;
                    MTU_normalized(:,7) ...%col_AV_norm_percent_elong_msc_GM_fuku = 100;
                    MTU_normalized(:,8) ...%col_AV_norm_percent_elong_SEE_fuku = 100;
                    ...
                    MTU_normalized_licht(loc_angle_licht_start:loc_angle_licht_stop,2) ... 
                    MTU_normalized_licht(loc_angle_licht_start:loc_angle_licht_stop,3) ... 
                    ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_angle_DFG) ... 
                    ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,7) ... % mean EMG 3 msc
                    ];

            end


            %% Prepare ARRAYS prone - group plots & stats
            if trial_timepoint == 0 && trial_leg == 1 % PRE, STR
                STR_PRE_prone(STR_PRE_count,:) = MTU_prone_vars;
            elseif trial_timepoint == 1 && trial_leg == 1 % POST, STR
                STR_POST_prone(STR_POST_count,:) = MTU_prone_vars;
            elseif trial_timepoint == 0 && trial_leg == 0 % PRE, CON
                CON_PRE_prone(CON_PRE_count,:) = MTU_prone_vars;
            elseif trial_timepoint == 1 && trial_leg == 0 % POST, CON
                CON_POST_prone(CON_POST_count,:) = MTU_prone_vars;
            end


            %% Prepare data table - individual trial data to file

            % txt trial ID
            all_passive_output_txt(line,:) = [dm_subjectno(line) dm_timepoint(line) dm_side(line) dm_trial(line)];

            % add data to a common array for all subjects    
            all_passive_output(line,:) = [...
                out_ROM_trial_max out_ROM_ind_max out_ROM_common_max out_ROM_submax_1 out_ROM_submax_2...
                out_F_trial_max_ROM out_F_trial_max_F out_F_ind_max out_F_common_max out_F_zero out_F_submax_1 out_F_submax_2...
                out_T_trial_max_ROM out_T_trial_max_F out_T_ind_max out_T_common_max out_T_zero out_T_submax_1 out_T_submax_2...
                out_angle_trial_max out_angle_ind_max out_angle_common_max out_angle_ind_rmax out_angle_ind_lmax...
                out_pstiff_trial_max out_pstiff_ind_max out_pstiff_common_max out_pstiff_15 out_pstiff_submax_1 out_pstiff_submax_2 out_pstiff_index...
                ...                % length, elong, strain
                out_length_AT_trial_max out_length_AT_ind_max out_length_AT_common_max out_length_AT_submax_1 out_length_AT_submax_2 out_length_AT_max...
                out_length_GMtend_trial_max out_length_GMtend_ind_max out_length_GMtend_common_max out_length_GMtend_submax_1 out_length_GMtend_submax_2 out_length_GMtend_max...
                out_length_GMapo_trial_max out_length_GMapo_ind_max out_length_GMapo_common_max out_length_GMapo_submax_1 out_length_GMapo_submax_2 out_length_GMapo_max...
                out_length_msc_GM_trial_max out_length_msc_GM_ind_max out_length_msc_GM_common_max out_length_msc_GM_submax_1 out_length_msc_GM_submax_2 out_length_msc_GM_max...
                out_length_msc_SOL_trial_max out_length_msc_SOL_ind_max out_length_msc_SOL_common_max out_length_msc_SOL_submax_1 out_length_msc_SOL_submax_2 out_length_msc_SOL_max...
                out_length_leg_trial_max out_length_leg_ind_max out_length_leg_common_max out_length_leg_submax_1 out_length_leg_submax_2 out_length_leg_max...
                out_length_SEE_Fuku_trial_max out_length_SEE_Fuku_ind_max out_length_SEE_Fuku_common_max out_length_SEE_Fuku_submax_1 out_length_SEE_Fuku_submax_2 out_length_SEE_Fuku_max...
                out_length_msc_GM_Fuku_trial_max out_length_msc_GM_Fuku_ind_max out_length_msc_GM_Fuku_common_max out_length_msc_GM_Fuku_submax_1 out_length_msc_GM_Fuku_submax_2 out_length_msc_GM_Fuku_max...
                ...
                out_elong_AT_trial_max out_elong_AT_ind_max out_elong_AT_common_max out_elong_AT_submax_1 out_elong_AT_submax_2 out_elong_AT_max...
                out_elong_GMtend_trial_max out_elong_GMtend_ind_max out_elong_GMtend_common_max out_elong_GMtend_submax_1 out_elong_GMtend_submax_2 out_elong_GMtend_max...
                out_elong_GMapo_trial_max out_elong_GMapo_ind_max out_elong_GMapo_common_max out_elong_GMapo_submax_1 out_elong_GMapo_submax_2 out_elong_GMapo_max...
                out_elong_msc_GM_trial_max out_elong_msc_GM_ind_max out_elong_msc_GM_common_max out_elong_msc_GM_submax_1 out_elong_msc_GM_submax_2 out_elong_msc_GM_max...
                out_elong_msc_SOL_trial_max out_elong_msc_SOL_ind_max out_elong_msc_SOL_common_max out_elong_msc_SOL_submax_1 out_elong_msc_SOL_submax_2 out_elong_msc_SOL_max...
                out_elong_leg_trial_max out_elong_leg_ind_max out_elong_leg_common_max out_elong_leg_submax_1 out_elong_leg_submax_2 out_elong_leg_max...
                out_elong_SEE_Fuku_trial_max out_elong_SEE_Fuku_ind_max out_elong_SEE_Fuku_common_max out_elong_SEE_Fuku_submax_1 out_elong_SEE_Fuku_submax_2 out_elong_SEE_Fuku_max...
                out_elong_msc_GM_Fuku_trial_max out_elong_msc_GM_Fuku_ind_max out_elong_msc_GM_Fuku_common_max out_elong_msc_GM_Fuku_submax_1 out_elong_msc_GM_Fuku_submax_2 out_elong_msc_GM_Fuku_max...
                ...
                out_strain_AT_trial_max out_strain_AT_ind_max out_strain_AT_common_max out_strain_AT_submax_1 out_strain_AT_submax_2 out_strain_AT_max...
                out_strain_GMtend_trial_max out_strain_GMtend_ind_max out_strain_GMtend_common_max out_strain_GMtend_submax_1 out_strain_GMtend_submax_2 out_strain_GMtend_max...
                out_strain_GMapo_trial_max out_strain_GMapo_ind_max out_strain_GMapo_common_max out_strain_GMapo_submax_1 out_strain_GMapo_submax_2 out_strain_GMapo_max...
                out_strain_msc_GM_trial_max out_strain_msc_GM_ind_max out_strain_msc_GM_common_max out_strain_msc_GM_submax_1 out_strain_msc_GM_submax_2 out_strain_msc_GM_max...
                out_strain_msc_SOL_trial_max out_strain_msc_SOL_ind_max out_strain_msc_SOL_common_max out_strain_msc_SOL_submax_1 out_strain_msc_SOL_submax_2 out_strain_msc_SOL_max...
                out_strain_leg_trial_max out_strain_leg_ind_max out_strain_leg_common_max out_strain_leg_submax_1 out_strain_leg_submax_2 out_strain_leg_max...
                out_strain_SEE_Fuku_trial_max out_strain_SEE_Fuku_ind_max out_strain_SEE_Fuku_common_max out_strain_SEE_Fuku_submax_1 out_strain_SEE_Fuku_submax_2 out_strain_SEE_Fuku_max...
                out_strain_msc_GM_Fuku_trial_max out_strain_msc_GM_Fuku_ind_max out_strain_msc_GM_Fuku_common_max out_strain_msc_GM_Fuku_submax_1 out_strain_msc_GM_Fuku_submax_2 out_strain_msc_GM_Fuku_max...
                ...                % normalized SEE & muscle
                out_norm_length_leg_trial_max out_norm_length_leg_ind_max out_norm_length_leg_common_max out_norm_length_leg_submax_1 out_norm_length_leg_submax_2 out_norm_length_leg_max...
                out_norm_length_SEE_Fuku_trial_max out_norm_length_SEE_Fuku_ind_max out_norm_length_SEE_Fuku_common_max out_norm_length_SEE_Fuku_submax_1 out_norm_length_SEE_Fuku_submax_2 out_norm_length_SEE_Fuku_max...
                out_norm_length_msc_GM_Fuku_trial_max out_norm_length_msc_GM_Fuku_ind_max out_norm_length_msc_GM_Fuku_common_max  out_norm_length_msc_GM_Fuku_submax_1 out_norm_length_msc_GM_Fuku_submax_2 out_norm_length_msc_GM_Fuku_max...
                ...
                out_norm_elong_leg_trial_max out_norm_elong_leg_ind_max out_norm_elong_leg_common_max out_norm_elong_leg_submax_1 out_norm_elong_leg_submax_2 out_norm_elong_leg_max...
                out_norm_elong_SEE_Fuku_trial_max out_norm_elong_SEE_Fuku_ind_max out_norm_elong_SEE_Fuku_common_max out_norm_elong_SEE_Fuku_submax_1 out_norm_elong_SEE_Fuku_submax_2 out_norm_elong_SEE_Fuku_max...
                out_norm_elong_msc_GM_Fuku_trial_max out_norm_elong_msc_GM_Fuku_ind_max out_norm_elong_msc_GM_Fuku_common_max out_norm_elong_msc_GM_Fuku_submax_1 out_norm_elong_msc_GM_Fuku_submax_2 out_norm_elong_msc_GM_Fuku_max...
                ...
                out_norm_elong_percent_SEE_Fuku_trial_max out_norm_elong_percent_SEE_Fuku_ind_max out_norm_elong_percent_SEE_Fuku_common_max out_norm_elong_percent_SEE_Fuku_submax_1 out_norm_elong_percent_SEE_Fuku_submax_2 out_norm_elong_percent_SEE_Fuku_max...
                out_norm_elong_percent_msc_GM_Fuku_trial_max out_norm_elong_percent_msc_GM_Fuku_ind_max out_norm_elong_percent_msc_GM_Fuku_common_max out_norm_elong_percent_msc_GM_Fuku_submax_1 out_norm_elong_percent_msc_GM_Fuku_submax_2 out_norm_elong_percent_msc_GM_Fuku_max...
                ...                % licht GM and SOL
                out_licht_pennation_GM_trial_max out_licht_pennation_GM_ind_max out_licht_pennation_GM_common_max out_licht_pennation_GM_submax_1 out_licht_pennation_GM_submax_2 out_licht_pennation_GM_max out_licht_pennation_GM_zero...
                out_licht_fas_length_GM_trial_max out_licht_fas_length_GM_ind_max out_licht_fas_length_GM_common_max out_licht_fas_length_GM_submax_1 out_licht_fas_length_GM_submax_2 out_licht_fas_length_GM_max out_licht_fas_length_GM_zero ...
                out_licht_fas_elong_GM_trial_max out_licht_fas_elong_GM_ind_max out_licht_fas_elong_GM_common_max out_licht_fas_elong_GM_submax_1 out_licht_fas_elong_GM_submax_2 out_licht_fas_elong_GM_max out_licht_fas_elong_GM_zero...
                out_licht_fas_strain_GM_trial_max out_licht_fas_strain_GM_ind_max out_licht_fas_strain_GM_common_max out_licht_fas_strain_GM_submax_1 out_licht_fas_strain_GM_submax_2 out_licht_fas_strain_GM_max out_licht_fas_strain_GM_zero...
                ...
                out_licht_pennation_SOL_trial_max out_licht_pennation_SOL_ind_max out_licht_pennation_SOL_common_max out_licht_pennation_SOL_submax_1 out_licht_pennation_SOL_submax_2 out_licht_pennation_SOL_max out_licht_pennation_SOL_zero...
                out_licht_fas_length_SOL_trial_max out_licht_fas_length_SOL_ind_max out_licht_fas_length_SOL_common_max  out_licht_fas_length_SOL_submax_1 out_licht_fas_length_SOL_submax_2 out_licht_fas_length_SOL_max out_licht_fas_length_SOL_zero...
                ...                % normalized GM fascicle length and elong
                out_norm_licht_fas_length_GM_trial_max out_norm_licht_fas_length_GM_ind_max out_norm_licht_fas_length_GM_common_max out_norm_licht_fas_length_GM_submax_1 out_norm_licht_fas_length_GM_submax_2  out_norm_licht_fas_length_GM_max out_norm_licht_fas_length_GM_zero ...
                out_norm_licht_fas_elong_GM_trial_max out_norm_licht_fas_elong_GM_ind_max out_norm_licht_fas_elong_GM_common_max  out_norm_licht_fas_elong_GM_submax_1 out_norm_licht_fas_elong_GM_max out_norm_licht_fas_elong_GM_submax_2 out_norm_licht_fas_elong_GM_zero...
                ...                % EMG
                out_emg_gm_trial_max out_emg_gm_ind_max out_emg_gm_common_max out_emg_gm_submax_1 out_emg_gm_submax_2 out_emg_gm_max out_emg_gm_zero ...
                out_emg_gl_trial_max out_emg_gl_ind_max out_emg_gl_common_max out_emg_gl_submax_1 out_emg_gl_submax_2 out_emg_gl_max out_emg_gl_zero...
                out_emg_sol_trial_max out_emg_sol_ind_max out_emg_sol_common_max out_emg_sol_submax_1 out_emg_sol_submax_2 out_emg_sol_max out_emg_sol_zero...
                out_emg_mean_trial_max out_emg_mean_ind_max out_emg_mean_common_max out_emg_mean_submax_1 out_emg_mean_submax_2 out_emg_mean_max out_emg_mean_zero...
                ];

            clear noraxon_mvc_plantar
            clear GMFAS_* SOL_* GMMTJ_*
            clear data_* MTU_*
            save all_data_passive_inloop
            close all
        end
        
        
        %% LOOP FINISHED /// loop end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        save all_data_passive_endloop
    end
    
    
    %% Truncate angle_vars cells 
        STR_PRE_angle_vars(STR_PRE_count+1:end) = [];
        STR_PRE_prone(STR_PRE_count+1:end,:) = [];
        STR_POST_angle_vars(STR_POST_count+1:end) = [];
        STR_POST_prone(STR_POST_count+1:end,:) = [];
        CON_PRE_angle_vars(CON_PRE_count+1:end) = [];
        CON_PRE_prone(CON_PRE_count+1:end,:) = [];
        CON_POST_angle_vars(CON_POST_count+1:end) = [];
        CON_POST_prone(CON_POST_count+1:end,:) = [];

%TMP
    %if input_resumerun < 2
    
    %% OUTPUT individual trial data XLS 
    
    % write xls
    if ispc
        filename_output = strcat('data_output/all_passive_output_', datestr(now, 'yyyymmdd_HHMM'), '.xlsx');
        
        xlswrite(filename_output, all_passive_output_head, 1, 'A1')
        xlswrite(filename_output, all_passive_output_txt, 1, 'A2')
        xlswrite(filename_output, all_passive_output, 1, 'E2')
    else
        filename_output = strcat('data_output/all_passive_output_', datestr(now, 'yyyymmdd_HHMM'), '.csv');
        csvwrite(filename_output, all_passive_output)
    end
    
    % TMP output fascicle data
%    filename_output = strcat('data_output/all_GMfas_output_', datestr(now, 'yyyymmdd_HHMM'), appendix, '.csv');
%    csvwrite(filename_output, all_GMfas_output)
    
    
    %% OUTPUT individual trial data for Graphpad Prism 
    filename_basis = strcat('data_output/prism_passive/all_passive', '_');   % datestr(now, 'yyyymmdd_HHMM'), '_');

    % GRAPHPAD "BY ROWS" FORMAT - 2 LINES (CON - STR) AND 2 COL (PRE-POST)
    
    % create output file header
    j = 1;
    prism_array_col1 = {'Subj'; all_passive_output_txt{j,4}; all_passive_output_txt{j+2,4} };

    % for each output variable (column):
    loc_relevant_var = 1; % first column of all_passive_data that has data for prism
    no_subjects = size(all_passive_output_txt,1)/4;
    for i = loc_relevant_var:size(all_passive_output,2)
        
        % reset
        prism_array(1:3,1:(2*no_subjects)) = NaN;
        write_col = 1;
        
        % filename with output variable
        filename_variable = strcat(filename_basis, all_passive_output_head{i+4}, '.xlsx');
        
        for j = 1:4:size(all_passive_output_txt,1) % number of lines / trials
            
            % row 1 = subject no
            prism_array(1,write_col) = str2double(all_passive_output_txt{j,1});
            prism_array(1,write_col+no_subjects) = str2double(all_passive_output_txt{j,1});
            
            % rows 2-3 in columns 1 and no_subjects+1 = data
            prism_array(2,write_col) = all_passive_output(j,i); % datamaster 1 = pre con
            prism_array(2,write_col+no_subjects) = all_passive_output(j+1,i); % datamaster 2 = post con
            prism_array(3,write_col) = all_passive_output(j+2,i); % datamaster 3 = pre str
            prism_array(3,write_col+no_subjects) = all_passive_output(j+3,i); % datamaster 4 = post str
            write_col = write_col + 1;
        end
        
        % write first col
        xlswrite(filename_variable, prism_array_col1, 1, 'A1')
        % write data
        xlswrite(filename_variable, prism_array, 1, 'B1')
        
        
        
    end

    
    %% OUTPUT group arrays for CURVE STATS, TO FILE
    
    % variables to export to file (OLD index numbers here):
        %   2 F
        %  18 Torque
        % 36 length msc GM (from Lichtwark/Fukunaga)
        % 28 elong msc GM (from Lichtwark/Fukunaga)
        % 30 strain msc GM (from Lichtwark/Fukunaga)
        % 37 length tend GM (from Lichtwark/Fukunaga)
        % 29 elong tend GM (from Lichtwark/Fukunaga)
        % 31 strain tend GM (from Lichtwark/Fukunaga)
        % 32 length fascicles GM (from Lichtwark)
        % 34 elongation fascicles GM (from Lichtwark)
        % 35 strain fascicles GM (from Lichtwark)
        % 33 pennation GM (from Lichtwark)
        % EMG mean
    out_arrays_input_cols = [col_AV_F col_AV_T col_AV_len_msc_GM_fuku col_AV_elong_msc_GM_fuku col_AV_strain_msc_GM_fuku col_AV_len_SEE_fuku col_AV_elong_SEE_fuku col_AV_strain_SEE_fuku col_AV_len_GMfas_licht col_AV_elong_GMfas_licht col_AV_strain_GMfas_licht col_AV_pennation_GMfas_licht col_AV_EMG_mean];
    out_arrays_input_labels = {'Force' 'Torque' 'GM muscle length' 'GM muscle elong' 'GM muscle strain' 'GM tendon length' 'GM tendon elong' 'GM tendon strain' 'GMfas length' 'GMfas elong' 'GMfas strain' 'GMfas pennation' 'EMG mean 3 msc'};
    
    if CON_PRE_count > 0 && CON_POST_count > 0 && STR_PRE_count > 0 && STR_POST_count > 0
        
        % preallocate
        STR_PRE_angle_vars_STAT{STR_PRE_count} = zeros;
        STR_POST_angle_vars_STAT{STR_POST_count} = zeros;
        CON_PRE_angle_vars_STAT{CON_PRE_count} = zeros;
        CON_POST_angle_vars_STAT{CON_POST_count} = zeros;
        
        % for absolute arrays: select common angle range = the subject/cell item containing the shortest matrix/ROM
        loc_commonROM = min([cellfun('length',STR_PRE_angle_vars) cellfun('length',STR_POST_angle_vars) cellfun('length',CON_PRE_angle_vars) cellfun('length',CON_POST_angle_vars)]); % location of largest common ROM
        
        % resample arrays - ommit columns with NaN only
        
        for i = 1:STR_PRE_count
            locate_NaN_cols = ~all(isnan(STR_PRE_angle_vars{1,i}(1:5,:)),1);
            resample_axis = (0:angle_step_stats_abs:STR_PRE_angle_vars{1,i}(loc_commonROM,1))';
            STR_PRE_angle_vars_STAT{i} = STR_PRE_angle_vars{1,i}(1:numel(resample_axis),:);
            STR_PRE_angle_vars_STAT{i}(:,locate_NaN_cols) = spline(STR_PRE_angle_vars{1,i}(1:loc_commonROM,1)',STR_PRE_angle_vars{1,i}(1:loc_commonROM,locate_NaN_cols).',resample_axis).';
        end
        for i = 1:STR_POST_count
            locate_NaN_cols = ~all(isnan(STR_POST_angle_vars{1,i}(1:5,:)),1);
            resample_axis = (0:angle_step_stats_abs:STR_POST_angle_vars{1,i}(loc_commonROM,1))';
            STR_POST_angle_vars_STAT{i} = STR_POST_angle_vars{1,i}(1:numel(resample_axis),:);
            STR_POST_angle_vars_STAT{i}(:,locate_NaN_cols) = spline(STR_POST_angle_vars{1,i}(1:loc_commonROM,1)',STR_POST_angle_vars{1,i}(1:loc_commonROM,locate_NaN_cols).',resample_axis).';
        end
        for i = 1:CON_PRE_count
            locate_NaN_cols = ~all(isnan(CON_PRE_angle_vars{1,i}(1:5,:)),1);
            resample_axis = (0:angle_step_stats_abs:CON_PRE_angle_vars{1,i}(loc_commonROM,1))';
            CON_PRE_angle_vars_STAT{i} = CON_PRE_angle_vars{1,i}(1:numel(resample_axis),:);
            CON_PRE_angle_vars_STAT{i}(:,locate_NaN_cols) = spline(CON_PRE_angle_vars{1,i}(1:loc_commonROM,1)',CON_PRE_angle_vars{1,i}(1:loc_commonROM,locate_NaN_cols).',resample_axis).';
        end
        for i = 1:CON_POST_count
            locate_NaN_cols = ~all(isnan(CON_POST_angle_vars{1,i}(1:5,:)),1);
            resample_axis = (0:angle_step_stats_abs:CON_POST_angle_vars{1,i}(loc_commonROM,1))';
            CON_POST_angle_vars_STAT{i} = CON_POST_angle_vars{1,i}(1:numel(resample_axis),:);
            CON_POST_angle_vars_STAT{i}(:,locate_NaN_cols) = spline(CON_POST_angle_vars{1,i}(1:loc_commonROM,1)',CON_POST_angle_vars{1,i}(1:loc_commonROM,locate_NaN_cols).',resample_axis).';
        end
        
        % create table headers (subject numbers)
        
        % pre and post values
        out_arrays_headers{1+STR_PRE_count+STR_POST_count+CON_PRE_count+CON_POST_count} = [];
        
        out_arrays_headers{1} = 'Joint_angle';
        for i=1:STR_PRE_count
            out_arrays_headers{i+1} = STR_PRE_ID{i};
        end
        for i=1:CON_PRE_count
            out_arrays_headers{i+1+STR_PRE_count} = CON_PRE_ID{i};
        end
        for i=1:STR_POST_count
            out_arrays_headers{i+1+STR_PRE_count+CON_PRE_count} = STR_POST_ID{i};
        end
        for i=1:CON_POST_count
            out_arrays_headers{i+1+STR_PRE_count+CON_PRE_count+STR_POST_count} = CON_POST_ID{i};
        end
        
        % difference values
        if eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count)
            % output of pre-post DIFFERENCES requires that trials are analysed in systematic order (does not check which subject-numbers to group)
            % - but will not be exported unless equal amount of 4 legs/timepoints - rough coding
            out_arrays_headers_diff{1+STR_PRE_count+CON_PRE_count} = [];
            out_arrays_headers_diff{1} = 'Joint_angle';
            for i=1:STR_PRE_count
                out_arrays_headers_diff{i+1} = strcat('STR_PP_', num2str(STR_PRE_no(i)), '__', num2str(STR_POST_no(i)));
            end
            for i=1:CON_PRE_count
                out_arrays_headers_diff{i+1+STR_PRE_count} = strcat('CON_PP_', num2str(CON_PRE_no(i)), '__',num2str(CON_POST_no(i)));
            end
        end
        
        % preallocate output arrays
        cols_abs = size(STR_PRE_angle_vars_STAT{1},1);
        rows = STR_PRE_count+STR_POST_count+CON_PRE_count+CON_POST_count + 1; % adding 1 for column for joint angles
        rows_diff = STR_PRE_count+CON_PRE_count + 1; % adding 1 for column for joint angles
        out_arrays_abs(cols_abs,rows) = zeros;
        out_arrays_abs_diff(cols_abs,rows_diff) = zeros;

        % organize and output table for each of the selected variables
        for var = 1:length(out_arrays_input_cols)
            % reset output arrays
            out_arrays_abs(cols_abs,rows) = zeros;
            out_arrays_abs_diff(cols_abs,rows_diff) = zeros;
            
            % add as first column, joint angles: abs and normalized angles
            out_arrays_abs(:,1) = STR_PRE_angle_vars_STAT{1}(:,1);
            out_arrays_abs_diff(:,1) = STR_PRE_angle_vars_STAT{1}(:,1);
            
            % add values: pre and post
            
            % add STR PRE first
            for subj = 1:STR_PRE_count
                out_arrays_abs(:,subj+1) = STR_PRE_angle_vars_STAT{subj}(:,out_arrays_input_cols(var));
            end            
            % add STR POST second
            for subj = 1:STR_POST_count
                out_arrays_abs(:,subj+STR_PRE_count+1) = STR_POST_angle_vars_STAT{subj}(:,out_arrays_input_cols(var));
            end            
            % add CON PRE
            for subj = 1:CON_PRE_count
                out_arrays_abs(:,subj+STR_PRE_count+STR_POST_count+1) = CON_PRE_angle_vars_STAT{subj}(:,out_arrays_input_cols(var));
            end            
            % add CON POST
            for subj = 1:CON_POST_count
                out_arrays_abs(:,subj+STR_PRE_count+STR_POST_count+CON_PRE_count+1) = CON_POST_angle_vars_STAT{subj}(:,out_arrays_input_cols(var));
            end
            
            % add values: difference between PRE-POST
            if eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count)
                % add STR first
                for subj = 1:STR_PRE_count
                    out_arrays_abs_diff(:,subj+1) = STR_POST_angle_vars_STAT{subj}(:,out_arrays_input_cols(var)) - STR_PRE_angle_vars_STAT{subj}(:,out_arrays_input_cols(var));
                end
                % add CON second
                for subj = 1:CON_PRE_count
                    out_arrays_abs_diff(:,subj+STR_PRE_count+1) = CON_POST_angle_vars_STAT{subj}(:,out_arrays_input_cols(var)) - CON_PRE_angle_vars_STAT{subj}(:,out_arrays_input_cols(var));
                end
            end
            
            % create tables and save as file
            % pre and post values
            out_arrays_abs_table = array2table(out_arrays_abs,'VariableNames',out_arrays_headers);
            filename_output = strcat('data_output/arrays_passive_acrossangles_', out_arrays_input_labels{var} , '_abs_', datestr(now, 'yyyymmdd_HHMM'));
            writetable(out_arrays_abs_table,filename_output,'Delimiter','\t')
            
            % difference values
            if eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count)
                out_arrays_abs_diff_table = array2table(out_arrays_abs_diff,'VariableNames',out_arrays_headers_diff);
                filename_output = strcat('data_output/arrays_passive_acrossangles_PP_', out_arrays_input_labels{var} , '_abs_', datestr(now, 'yyyymmdd_HHMM'));
                writetable(out_arrays_abs_diff_table,filename_output,'Delimiter','\t')
            end
            clear out_arrays_abs_table out_arrays_abs_diff_table % out_arrays_norm_table out_arrays_norm_diff_table
        end
    end
    
%% TMP  end


    %% GROUP figures - create variables of MEAN + STDAV 
    
    %%%  mean and stdav of each subject's INDIVIDUAL MAX ROM, force, elong, EMG, etc
       
    % prone variables
    STR_PRE_prone_mean = nanmean(STR_PRE_prone);
    STR_POST_prone_mean = nanmean(STR_POST_prone);
    CON_PRE_prone_mean = nanmean(CON_PRE_prone);
    CON_POST_prone_mean = nanmean(CON_POST_prone);

    n_o_array_elements = max( [length(CON_PRE_angle_vars{1,1}(1,:)) length(STR_PRE_angle_vars{1,1}(1,:)) length(CON_POST_angle_vars{1,1}(1,:)) length(STR_POST_angle_vars{1,1}(1,:))] );

    % STR PRE
    if STR_PRE_count > 0
        % preallocate array
        STR_PRE_max(STR_PRE_count,n_o_array_elements) = zeros;
        % collect variables per subject
        for i = 1:STR_PRE_count % per subject
            for j = 1:n_o_array_elements % per element in arrays
                % OLD: Max values
%                STR_PRE_max(i,j) = max(STR_PRE_angle_vars{1,i}(:,j));
                % NEW: end of array values
                STR_PRE_max(i,j) = max(STR_PRE_angle_vars{1,i}(end,j));
            end
        end
        
        % calculate mean and SD of max values across subjects
        
        % misc
        STR_PRE_ROM_mean = nanmean(STR_PRE_max(:,col_AV_angle));
        STR_PRE_ROM_SD = nanstd(STR_PRE_max(:,col_AV_angle));
        STR_PRE_F_mean = nanmean(STR_PRE_max(:,col_AV_F));
        STR_PRE_F_SD = nanstd(STR_PRE_max(:,col_AV_F));
        STR_PRE_EMG_gm_mean = nanmean(STR_PRE_max(:,col_AV_EMG_gm));
        STR_PRE_EMG_gm_SD = nanstd(STR_PRE_max(:,col_AV_EMG_gm));
        STR_PRE_EMG_gl_mean = nanmean(STR_PRE_max(:,col_AV_EMG_gl));
        STR_PRE_EMG_gl_SD = nanstd(STR_PRE_max(:,col_AV_EMG_gl));
        STR_PRE_EMG_sol_mean = nanmean(STR_PRE_max(:,col_AV_EMG_sol));
        STR_PRE_EMG_sol_SD = nanstd(STR_PRE_max(:,col_AV_EMG_sol));
        STR_PRE_EMG_mean_mean = nanmean(STR_PRE_max(:,col_AV_EMG_mean));
        STR_PRE_EMG_mean_SD = nanstd(STR_PRE_max(:,col_AV_EMG_mean));
        
        % elong
        STR_PRE_elong_AT_mean = nanmean(STR_PRE_max(:,col_AV_elong_AT)); % elong AT
        STR_PRE_elong_AT_SD = nanstd(STR_PRE_max(:,col_AV_elong_AT));
        STR_PRE_elong_GMtend_mean = nanmean(STR_PRE_max(:,col_AV_elong_GMtend)); % elong GM tend
        STR_PRE_elong_GMtend_SD = nanstd(STR_PRE_max(:,col_AV_elong_GMtend));
        STR_PRE_elong_MTU_mean = nanmean(STR_PRE_max(:,col_AV_elong_leg)); % elong leg
        STR_PRE_elong_MTU_SD = nanstd(STR_PRE_max(:,col_AV_elong_leg));
        STR_PRE_displ_GMFAS_mean = nanmean(STR_PRE_max(:,col_AV_elong_GMFAS));  % GMFAS displ 
        STR_PRE_displ_GMFAS_SD = nanstd(STR_PRE_max(:,col_AV_elong_GMFAS));
        STR_PRE_elong_GMapo_mean = nanmean(STR_PRE_max(:,col_AV_elong_GMapo)); 
        STR_PRE_elong_GMapo_SD = nanstd(STR_PRE_max(:,col_AV_elong_GMapo));
        STR_PRE_elong_msc_GM_mean = nanmean(STR_PRE_max(:,col_AV_elong_msc_GM)); % GM msc
        STR_PRE_elong_msc_GM_SD = nanstd(STR_PRE_max(:,col_AV_elong_msc_GM));
        STR_PRE_elong_msc_SOL_mean = nanmean(STR_PRE_max(:,col_AV_elong_msc_SOL)); 
        STR_PRE_elong_msc_SOL_SD = nanstd(STR_PRE_max(:,col_AV_elong_msc_SOL));

        % length
        STR_PRE_L_at_SOL_mean = nanmean(STR_PRE_max(:,col_AV_len_AT)); % L AT
        STR_PRE_L_at_SOL_SD = nanstd(STR_PRE_max(:,col_AV_len_AT));
        STR_PRE_L_at_GM_mean = nanmean(STR_PRE_max(:,col_AV_len_GMtend)); % L GM tend
        STR_PRE_L_at_GM_SD = nanstd(STR_PRE_max(:,col_AV_len_GMtend));
        STR_PRE_L_MTU_mean = nanmean(STR_PRE_max(:,col_AV_len_leg)); % L leg
        STR_PRE_L_MTU_SD = nanstd(STR_PRE_max(:,col_AV_len_leg));
        % no L GMFAS
        STR_PRE_L_GMapo_mean = nanmean(STR_PRE_max(:,col_AV_len_GMapo)); 
        STR_PRE_L_GMapo_SD = nanstd(STR_PRE_max(:,col_AV_len_GMapo));
        STR_PRE_L_msc_GM_mean = nanmean(STR_PRE_max(:,col_AV_len_msc_GM)); 
        STR_PRE_L_msc_GM_SD = nanstd(STR_PRE_max(:,col_AV_len_msc_GM));
        STR_PRE_L_msc_SOL_mean = nanmean(STR_PRE_max(:,col_AV_len_msc_SOL)); 
        STR_PRE_L_msc_SOL_SD = nanstd(STR_PRE_max(:,col_AV_len_msc_SOL));

        % strain
        STR_PRE_strain_at_SOL_mean = nanmean(STR_PRE_max(:,col_AV_strain_AT)); % L AT
        STR_PRE_strain_at_SOL_SD = nanstd(STR_PRE_max(:,col_AV_strain_AT));
        STR_PRE_strain_at_GM_mean = nanmean(STR_PRE_max(:,col_AV_strain_GMtend)); % L GM tend
        STR_PRE_strain_at_GM_SD = nanstd(STR_PRE_max(:,col_AV_strain_GMtend));
        STR_PRE_strain_MTU_mean = nanmean(STR_PRE_max(:,col_AV_strain_leg)); % L calf
        STR_PRE_strain_MTU_SD = nanstd(STR_PRE_max(:,col_AV_strain_leg));
        % no strain GMfas
        STR_PRE_strain_GMapo_mean = nanmean(STR_PRE_max(:,col_AV_strain_GMapo)); 
        STR_PRE_strain_GMapo_SD = nanstd(STR_PRE_max(:,col_AV_strain_GMapo));
        STR_PRE_strain_msc_GM_mean = nanmean(STR_PRE_max(:,col_AV_strain_msc_GM)); 
        STR_PRE_strain_msc_GM_SD = nanstd(STR_PRE_max(:,col_AV_strain_msc_GM));
        STR_PRE_strain_msc_SOL_mean = nanmean(STR_PRE_max(:,col_AV_strain_msc_SOL)); 
        STR_PRE_strain_msc_SOL_SD = nanstd(STR_PRE_max(:,col_AV_strain_msc_SOL));

        % elong licht msc/SEE
        STR_PRE_elong_msc_GM_licht_mean = nanmean(STR_PRE_max(:,col_AV_elong_msc_GM_fuku)); 
        STR_PRE_elong_msc_GM_licht_SD = nanstd(STR_PRE_max(:,col_AV_elong_msc_GM_fuku));
        STR_PRE_elong_SEE_licht_mean = nanmean(STR_PRE_max(:,col_AV_elong_SEE_fuku)); 
        STR_PRE_elong_SEE_licht_SD = nanstd(STR_PRE_max(:,col_AV_elong_SEE_fuku));
        % strain licht msc/SEE
        STR_PRE_strain_msc_GM_licht_mean = nanmean(STR_PRE_max(:,col_AV_strain_msc_GM_fuku)); 
        STR_PRE_strain_msc_GM_licht_SD = nanstd(STR_PRE_max(:,col_AV_strain_msc_GM_fuku));
        STR_PRE_strain_SEE_licht_mean = nanmean(STR_PRE_max(:,col_AV_strain_SEE_fuku)); 
        STR_PRE_strain_SEE_licht_SD = nanstd(STR_PRE_max(:,col_AV_strain_SEE_fuku));
        
        %licht fascicle
        STR_PRE_length_GMfas_licht_mean = nanmean(STR_PRE_max(:,col_AV_len_GMfas_licht)); 
        STR_PRE_length_GMfas_licht_SD = nanstd(STR_PRE_max(:,col_AV_len_GMfas_licht)); 
        STR_PRE_pennation_GMfas_licht_mean = nanmean(STR_PRE_max(:,col_AV_pennation_GMfas_licht)); 
        STR_PRE_pennation_GMfas_licht_SD = nanstd(STR_PRE_max(:,col_AV_pennation_GMfas_licht)); 
        STR_PRE_elong_GMfas_licht_mean = nanmean(STR_PRE_max(:,col_AV_elong_GMfas_licht)); 
        STR_PRE_elong_GMfas_licht_SD = nanstd(STR_PRE_max(:,col_AV_elong_GMfas_licht)); 
        STR_PRE_strain_GMfas_licht_mean = nanmean(STR_PRE_max(:,col_AV_strain_GMfas_licht)); 
        STR_PRE_strain_GMfas_licht_SD = nanstd(STR_PRE_max(:,col_AV_strain_GMfas_licht)); 

        % length licht msc/SEE
        STR_PRE_length_msc_GM_licht_mean = nanmean(STR_PRE_max(:,col_AV_len_msc_GM_fuku)); 
        STR_PRE_length_msc_GM_licht_SD = nanstd(STR_PRE_max(:,col_AV_len_msc_GM_fuku));
        STR_PRE_length_SEE_licht_mean = nanmean(STR_PRE_max(:,col_AV_len_SEE_fuku)); 
        STR_PRE_length_SEE_licht_SD = nanstd(STR_PRE_max(:,col_AV_len_SEE_fuku));

        % normalized 10 vars
        STR_PRE_length_norm_MTU_mean = nanmean(STR_PRE_max(:,col_AV_norm_len_leg)); % L leg
        STR_PRE_length_norm_MTU_SD = nanstd(STR_PRE_max(:,col_AV_norm_len_leg));
        STR_PRE_elong_norm_MTU_mean = nanmean(STR_PRE_max(:,col_AV_norm_elong_leg)); % elong leg
        STR_PRE_elong_norm_MTU_SD = nanstd(STR_PRE_max(:,col_AV_norm_elong_leg));

        STR_PRE_length_norm_msc_GM_licht_mean = nanmean(STR_PRE_max(:,col_AV_norm_len_msc_GM_fuku)); 
        STR_PRE_length_norm_msc_GM_licht_SD = nanstd(STR_PRE_max(:,col_AV_norm_len_msc_GM_fuku));
        STR_PRE_elong_norm_msc_GM_licht_mean = nanmean(STR_PRE_max(:,col_AV_norm_elong_msc_GM_fuku)); 
        STR_PRE_elong_norm_msc_GM_licht_SD = nanstd(STR_PRE_max(:,col_AV_norm_elong_msc_GM_fuku));

        STR_PRE_length_norm_SEE_licht_mean = nanmean(STR_PRE_max(:,col_AV_norm_len_SEE_fuku)); 
        STR_PRE_length_norm_SEE_licht_SD = nanstd(STR_PRE_max(:,col_AV_norm_len_SEE_fuku));
        STR_PRE_elong_norm_SEE_licht_mean = nanmean(STR_PRE_max(:,col_AV_norm_elong_SEE_fuku)); 
        STR_PRE_elong_norm_SEE_licht_SD = nanstd(STR_PRE_max(:,col_AV_norm_elong_SEE_fuku));

        STR_PRE_elong_norm_percent_msc_GM_licht_mean = nanmean(STR_PRE_max(:,col_AV_norm_percent_elong_msc_GM_fuku)); 
        STR_PRE_elong_norm_percent_msc_GM_licht_SD = nanstd(STR_PRE_max(:,col_AV_norm_percent_elong_msc_GM_fuku));
        STR_PRE_elong_norm_percent_SEE_licht_mean = nanmean(STR_PRE_max(:,col_AV_norm_percent_elong_SEE_fuku)); 
        STR_PRE_elong_norm_percent_SEE_licht_SD = nanstd(STR_PRE_max(:,col_AV_norm_percent_elong_SEE_fuku));

        STR_PRE_length_norm_GMfas_licht_mean = nanmean(STR_PRE_max(:,col_AV_norm_len_GMfas_licht)); 
        STR_PRE_length_norm_GMfas_licht_SD = nanstd(STR_PRE_max(:,col_AV_norm_len_GMfas_licht)); 
        STR_PRE_elong_norm_GMfas_licht_mean = nanmean(STR_PRE_max(:,col_AV_norm_elong_GMfas_licht)); 
        STR_PRE_elong_norm_GMfas_licht_SD = nanstd(STR_PRE_max(:,col_AV_norm_elong_GMfas_licht)); 
        
        % misc
        STR_PRE_torque_mean = nanmean(STR_PRE_max(:,col_AV_T));
        STR_PRE_torque_SD = nanstd(STR_PRE_max(:,col_AV_T));
        % determine common angle range
        STR_PRE_common_ROM = min(STR_PRE_max(:,col_AV_angle));
    end

    % STR POST
    if STR_POST_count > 0
        % preallocate array
        STR_POST_max(STR_POST_count,n_o_array_elements) = zeros;
 %       STR_POST_max_norm(STR_POST_count,n_o_array_elements) = zeros;
        % collect variables per subject
        for i = 1:STR_POST_count % per subject
            for j = 1:n_o_array_elements % per element in arrays
                % OLD: Max values
%                STR_POST_max(i,j) = max(STR_POST_angle_vars{1,i}(:,j));
%                STR_POST_max_norm(i,j) = max(STR_POST_angle_vars_norm{1,i}(:,j));
                % NEW: end of array values
                STR_POST_max(i,j) = max(STR_POST_angle_vars{1,i}(end,j));
%                STR_POST_max_norm(i,j) = max(STR_POST_angle_vars_norm{1,i}(end,j));
            end
        end
        % calculate mean and SD of max values across subjects
        STR_POST_ROM_mean = nanmean(STR_POST_max(:,col_AV_angle));
        STR_POST_ROM_SD = nanstd(STR_POST_max(:,col_AV_angle));
        STR_POST_F_mean = nanmean(STR_POST_max(:,col_AV_F));
        STR_POST_F_SD = nanstd(STR_POST_max(:,col_AV_F));
        STR_POST_EMG_gm_mean = nanmean(STR_POST_max(:,col_AV_EMG_gm));
        STR_POST_EMG_gm_SD = nanstd(STR_POST_max(:,col_AV_EMG_gm));
        STR_POST_EMG_gl_mean = nanmean(STR_POST_max(:,col_AV_EMG_gl));
        STR_POST_EMG_gl_SD = nanstd(STR_POST_max(:,col_AV_EMG_gl));
        STR_POST_EMG_sol_mean = nanmean(STR_POST_max(:,col_AV_EMG_sol));
        STR_POST_EMG_sol_SD = nanstd(STR_POST_max(:,col_AV_EMG_sol));
        STR_POST_EMG_mean_mean = nanmean(STR_POST_max(:,col_AV_EMG_mean));
        STR_POST_EMG_mean_SD = nanstd(STR_POST_max(:,col_AV_EMG_mean));
        
        STR_POST_elong_AT_mean = nanmean(STR_POST_max(:,col_AV_elong_AT)); % elong AT
        STR_POST_elong_AT_SD = nanstd(STR_POST_max(:,col_AV_elong_AT));
        STR_POST_elong_GMtend_mean = nanmean(STR_POST_max(:,col_AV_elong_GMtend)); % elong GM tend
        STR_POST_elong_GMtend_SD = nanstd(STR_POST_max(:,col_AV_elong_GMtend));
        STR_POST_elong_MTU_mean = nanmean(STR_POST_max(:,col_AV_elong_leg)); % elong leg
        STR_POST_elong_MTU_SD = nanstd(STR_POST_max(:,col_AV_elong_leg));
        STR_POST_displ_GMFAS_mean = nanmean(STR_POST_max(:,col_AV_elong_GMFAS));  % GMFAS displ 
        STR_POST_displ_GMFAS_SD = nanstd(STR_POST_max(:,col_AV_elong_GMFAS));
        STR_POST_elong_GMapo_mean = nanmean(STR_POST_max(:,col_AV_elong_GMapo)); 
        STR_POST_elong_GMapo_SD = nanstd(STR_POST_max(:,col_AV_elong_GMapo));
        STR_POST_elong_msc_GM_mean = nanmean(STR_POST_max(:,col_AV_elong_msc_GM)); % GM msc
        STR_POST_elong_msc_GM_SD = nanstd(STR_POST_max(:,col_AV_elong_msc_GM));
        STR_POST_elong_msc_SOL_mean = nanmean(STR_POST_max(:,col_AV_elong_msc_SOL)); 
        STR_POST_elong_msc_SOL_SD = nanstd(STR_POST_max(:,col_AV_elong_msc_SOL));

        STR_POST_L_at_SOL_mean = nanmean(STR_POST_max(:,col_AV_len_AT)); % L AT
        STR_POST_L_at_SOL_SD = nanstd(STR_POST_max(:,col_AV_len_AT));
        STR_POST_L_at_GM_mean = nanmean(STR_POST_max(:,col_AV_len_GMtend)); % L GM tend
        STR_POST_L_at_GM_SD = nanstd(STR_POST_max(:,col_AV_len_GMtend));
        STR_POST_L_MTU_mean = nanmean(STR_POST_max(:,col_AV_len_leg)); % L leg
        STR_POST_L_MTU_SD = nanstd(STR_POST_max(:,col_AV_len_leg));
        % no L GMFAS
        STR_POST_L_GMapo_mean = nanmean(STR_POST_max(:,col_AV_len_GMapo)); 
        STR_POST_L_GMapo_SD = nanstd(STR_POST_max(:,col_AV_len_GMapo));
        STR_POST_L_msc_GM_mean = nanmean(STR_POST_max(:,col_AV_len_msc_GM)); 
        STR_POST_L_msc_GM_SD = nanstd(STR_POST_max(:,col_AV_len_msc_GM));
        STR_POST_L_msc_SOL_mean = nanmean(STR_POST_max(:,col_AV_len_msc_SOL)); 
        STR_POST_L_msc_SOL_SD = nanstd(STR_POST_max(:,col_AV_len_msc_SOL));

        STR_POST_strain_at_SOL_mean = nanmean(STR_POST_max(:,col_AV_strain_AT)); % L AT
        STR_POST_strain_at_SOL_SD = nanstd(STR_POST_max(:,col_AV_strain_AT));
        STR_POST_strain_at_GM_mean = nanmean(STR_POST_max(:,col_AV_strain_GMtend)); % L GM tend
        STR_POST_strain_at_GM_SD = nanstd(STR_POST_max(:,col_AV_strain_GMtend));
        STR_POST_strain_MTU_mean = nanmean(STR_POST_max(:,col_AV_strain_leg)); % L calf
        STR_POST_strain_MTU_SD = nanstd(STR_POST_max(:,col_AV_strain_leg));
        % no strain GMfas
        STR_POST_strain_GMapo_mean = nanmean(STR_POST_max(:,col_AV_strain_GMapo)); 
        STR_POST_strain_GMapo_SD = nanstd(STR_POST_max(:,col_AV_strain_GMapo));
        STR_POST_strain_msc_GM_mean = nanmean(STR_POST_max(:,col_AV_strain_msc_GM)); 
        STR_POST_strain_msc_GM_SD = nanstd(STR_POST_max(:,col_AV_strain_msc_GM));
        STR_POST_strain_msc_SOL_mean = nanmean(STR_POST_max(:,col_AV_strain_msc_SOL)); 
        STR_POST_strain_msc_SOL_SD = nanstd(STR_POST_max(:,col_AV_strain_msc_SOL));

        STR_POST_elong_msc_GM_licht_mean = nanmean(STR_POST_max(:,col_AV_elong_msc_GM_fuku)); 
        STR_POST_elong_msc_GM_licht_SD = nanstd(STR_POST_max(:,col_AV_elong_msc_GM_fuku));
        STR_POST_elong_tend_GM_licht_mean = nanmean(STR_POST_max(:,col_AV_elong_SEE_fuku)); 
        STR_POST_elong_tend_GM_licht_SD = nanstd(STR_POST_max(:,col_AV_elong_SEE_fuku));
        STR_POST_strain_msc_GM_licht_mean = nanmean(STR_POST_max(:,col_AV_strain_msc_GM_fuku)); 
        STR_POST_strain_msc_GM_licht_SD = nanstd(STR_POST_max(:,col_AV_strain_msc_GM_fuku));
        STR_POST_strain_tend_GM_licht_mean = nanmean(STR_POST_max(:,col_AV_strain_SEE_fuku)); 
        STR_POST_strain_tend_GM_licht_SD = nanstd(STR_POST_max(:,col_AV_strain_SEE_fuku));

        STR_POST_length_GMfas_licht_mean = nanmean(STR_POST_max(:,col_AV_len_GMfas_licht)); 
        STR_POST_length_GMfas_licht_SD = nanstd(STR_POST_max(:,col_AV_len_GMfas_licht)); 
        STR_POST_pennation_GMfas_licht_mean = nanmean(STR_POST_max(:,col_AV_pennation_GMfas_licht)); 
        STR_POST_pennation_GMfas_licht_SD = nanstd(STR_POST_max(:,col_AV_pennation_GMfas_licht)); 
        STR_POST_elong_GMfas_licht_mean = nanmean(STR_POST_max(:,col_AV_elong_GMfas_licht)); 
        STR_POST_elong_GMfas_licht_SD = nanstd(STR_POST_max(:,col_AV_elong_GMfas_licht)); 
        STR_POST_strain_GMfas_licht_mean = nanmean(STR_POST_max(:,col_AV_strain_GMfas_licht)); 
        STR_POST_strain_GMfas_licht_SD = nanstd(STR_POST_max(:,col_AV_strain_GMfas_licht)); 

        STR_POST_length_msc_GM_licht_mean = nanmean(STR_POST_max(:,col_AV_len_msc_GM_fuku)); 
        STR_POST_length_msc_GM_licht_SD = nanstd(STR_POST_max(:,col_AV_len_msc_GM_fuku));
        STR_POST_length_tend_GM_licht_mean = nanmean(STR_POST_max(:,col_AV_len_SEE_fuku)); 
        STR_POST_length_tend_GM_licht_SD = nanstd(STR_POST_max(:,col_AV_len_SEE_fuku));

        STR_POST_length_norm_MTU_mean = nanmean(STR_POST_max(:,col_AV_norm_len_leg)); % L leg
        STR_POST_length_norm_MTU_SD = nanstd(STR_POST_max(:,col_AV_norm_len_leg));
        STR_POST_elong_norm_MTU_mean = nanmean(STR_POST_max(:,col_AV_norm_elong_leg)); % elong leg
        STR_POST_elong_norm_MTU_SD = nanstd(STR_POST_max(:,col_AV_norm_elong_leg));

        STR_POST_length_norm_msc_GM_licht_mean = nanmean(STR_POST_max(:,col_AV_norm_len_msc_GM_fuku)); 
        STR_POST_length_norm_msc_GM_licht_SD = nanstd(STR_POST_max(:,col_AV_norm_len_msc_GM_fuku));
        STR_POST_elong_norm_msc_GM_licht_mean = nanmean(STR_POST_max(:,col_AV_norm_elong_msc_GM_fuku)); 
        STR_POST_elong_norm_msc_GM_licht_SD = nanstd(STR_POST_max(:,col_AV_norm_elong_msc_GM_fuku));

        STR_POST_length_norm_SEE_licht_mean = nanmean(STR_POST_max(:,col_AV_norm_len_SEE_fuku)); 
        STR_POST_length_norm_SEE_licht_SD = nanstd(STR_POST_max(:,col_AV_norm_len_SEE_fuku));
        STR_POST_elong_norm_SEE_licht_mean = nanmean(STR_POST_max(:,col_AV_norm_elong_SEE_fuku)); 
        STR_POST_elong_norm_SEE_licht_SD = nanstd(STR_POST_max(:,col_AV_norm_elong_SEE_fuku));

        STR_POST_elong_norm_percent_msc_GM_licht_mean = nanmean(STR_POST_max(:,col_AV_norm_percent_elong_msc_GM_fuku)); 
        STR_POST_elong_norm_percent_msc_GM_licht_SD = nanstd(STR_POST_max(:,col_AV_norm_percent_elong_msc_GM_fuku));
        STR_POST_elong_norm_percent_SEE_licht_mean = nanmean(STR_POST_max(:,col_AV_norm_percent_elong_SEE_fuku)); 
        STR_POST_elong_norm_percent_SEE_licht_SD = nanstd(STR_POST_max(:,col_AV_norm_percent_elong_SEE_fuku));

        STR_POST_length_norm_GMfas_licht_mean = nanmean(STR_POST_max(:,col_AV_norm_len_GMfas_licht)); 
        STR_POST_length_norm_GMfas_licht_SD = nanstd(STR_POST_max(:,col_AV_norm_len_GMfas_licht)); 
        STR_POST_elong_norm_GMfas_licht_mean = nanmean(STR_POST_max(:,col_AV_norm_elong_GMfas_licht)); 
        STR_POST_elong_norm_GMfas_licht_SD = nanstd(STR_POST_max(:,col_AV_norm_elong_GMfas_licht)); 
        
        STR_POST_torque_mean = nanmean(STR_POST_max(:,col_AV_T));
        STR_POST_torque_SD = nanstd(STR_POST_max(:,col_AV_T));
        % determine common angle range
        STR_POST_common_ROM = min(STR_POST_max(:,col_AV_angle));
    end

    % CON PRE
    if CON_PRE_count > 0
        % preallocate array
        CON_PRE_max(CON_PRE_count,n_o_array_elements) = zeros;
%         CON_PRE_max_norm(CON_PRE_count,n_o_array_elements) = zeros;
        % collect variables per subject
        for i = 1:CON_PRE_count % per subject
            for j = 1:n_o_array_elements % per element in arrays
                % OLD: Max values
%                CON_PRE_max(i,j) = max(CON_PRE_angle_vars{1,i}(:,j));
%                CON_PRE_max_norm(i,j) = max(CON_PRE_angle_vars_norm{1,i}(:,j));
                % NEW: end of array values
                CON_PRE_max(i,j) = max(CON_PRE_angle_vars{1,i}(end,j));
%                CON_PRE_max_norm(i,j) = max(CON_PRE_angle_vars_norm{1,i}(end,j));
            end
        end
        % calculate mean and SD of max values across subjects
        CON_PRE_ROM_mean = nanmean(CON_PRE_max(:,col_AV_angle));
        CON_PRE_ROM_SD = nanstd(CON_PRE_max(:,col_AV_angle));
        CON_PRE_F_mean = nanmean(CON_PRE_max(:,col_AV_F));
        CON_PRE_F_SD = nanstd(CON_PRE_max(:,col_AV_F));
        CON_PRE_EMG_gm_mean = nanmean(CON_PRE_max(:,col_AV_EMG_gm));
        CON_PRE_EMG_gm_SD = nanstd(CON_PRE_max(:,col_AV_EMG_gm));
        CON_PRE_EMG_gl_mean = nanmean(CON_PRE_max(:,col_AV_EMG_gl));
        CON_PRE_EMG_gl_SD = nanstd(CON_PRE_max(:,col_AV_EMG_gl));
        CON_PRE_EMG_sol_mean = nanmean(CON_PRE_max(:,col_AV_EMG_sol));
        CON_PRE_EMG_sol_SD = nanstd(CON_PRE_max(:,col_AV_EMG_sol));
        CON_PRE_EMG_mean_mean = nanmean(CON_PRE_max(:,col_AV_EMG_mean));
        CON_PRE_EMG_mean_SD = nanstd(CON_PRE_max(:,col_AV_EMG_mean));
        
        CON_PRE_elong_AT_mean = nanmean(CON_PRE_max(:,col_AV_elong_AT)); % elong AT
        CON_PRE_elong_AT_SD = nanstd(CON_PRE_max(:,col_AV_elong_AT));
        CON_PRE_elong_GMtend_mean = nanmean(CON_PRE_max(:,col_AV_elong_GMtend)); % elong GM tend
        CON_PRE_elong_GMtend_SD = nanstd(CON_PRE_max(:,col_AV_elong_GMtend));
        CON_PRE_elong_MTU_mean = nanmean(CON_PRE_max(:,col_AV_elong_leg)); % elong leg
        CON_PRE_elong_MTU_SD = nanstd(CON_PRE_max(:,col_AV_elong_leg));
        CON_PRE_displ_GMFAS_mean = nanmean(CON_PRE_max(:,col_AV_elong_GMFAS));  % GMFAS displ 
        CON_PRE_displ_GMFAS_SD = nanstd(CON_PRE_max(:,col_AV_elong_GMFAS));
        CON_PRE_elong_GMapo_mean = nanmean(CON_PRE_max(:,col_AV_elong_GMapo)); 
        CON_PRE_elong_GMapo_SD = nanstd(CON_PRE_max(:,col_AV_elong_GMapo));
        CON_PRE_elong_msc_GM_mean = nanmean(CON_PRE_max(:,col_AV_elong_msc_GM)); % GM msc
        CON_PRE_elong_msc_GM_SD = nanstd(CON_PRE_max(:,col_AV_elong_msc_GM));
        CON_PRE_elong_msc_SOL_mean = nanmean(CON_PRE_max(:,col_AV_elong_msc_SOL)); 
        CON_PRE_elong_msc_SOL_SD = nanstd(CON_PRE_max(:,col_AV_elong_msc_SOL));

        CON_PRE_L_at_SOL_mean = nanmean(CON_PRE_max(:,col_AV_len_AT)); % L AT
        CON_PRE_L_at_SOL_SD = nanstd(CON_PRE_max(:,col_AV_len_AT));
        CON_PRE_L_at_GM_mean = nanmean(CON_PRE_max(:,col_AV_len_GMtend)); % L GM tend
        CON_PRE_L_at_GM_SD = nanstd(CON_PRE_max(:,col_AV_len_GMtend));
        CON_PRE_L_MTU_mean = nanmean(CON_PRE_max(:,col_AV_len_leg)); % L leg
        CON_PRE_L_MTU_SD = nanstd(CON_PRE_max(:,col_AV_len_leg));
        % no L GMFAS
        CON_PRE_L_GMapo_mean = nanmean(CON_PRE_max(:,col_AV_len_GMapo)); 
        CON_PRE_L_GMapo_SD = nanstd(CON_PRE_max(:,col_AV_len_GMapo));
        CON_PRE_L_msc_GM_mean = nanmean(CON_PRE_max(:,col_AV_len_msc_GM)); 
        CON_PRE_L_msc_GM_SD = nanstd(CON_PRE_max(:,col_AV_len_msc_GM));
        CON_PRE_L_msc_SOL_mean = nanmean(CON_PRE_max(:,col_AV_len_msc_SOL)); 
        CON_PRE_L_msc_SOL_SD = nanstd(CON_PRE_max(:,col_AV_len_msc_SOL));

        CON_PRE_strain_at_SOL_mean = nanmean(CON_PRE_max(:,col_AV_strain_AT)); % L AT
        CON_PRE_strain_at_SOL_SD = nanstd(CON_PRE_max(:,col_AV_strain_AT));
        CON_PRE_strain_at_GM_mean = nanmean(CON_PRE_max(:,col_AV_strain_GMtend)); % L GM tend
        CON_PRE_strain_at_GM_SD = nanstd(CON_PRE_max(:,col_AV_strain_GMtend));
        CON_PRE_strain_MTU_mean = nanmean(CON_PRE_max(:,col_AV_strain_leg)); % L calf
        CON_PRE_strain_MTU_SD = nanstd(CON_PRE_max(:,col_AV_strain_leg));
        % no strain GMfas
        CON_PRE_strain_GMapo_mean = nanmean(CON_PRE_max(:,col_AV_strain_GMapo)); 
        CON_PRE_strain_GMapo_SD = nanstd(CON_PRE_max(:,col_AV_strain_GMapo));
        CON_PRE_strain_msc_GM_mean = nanmean(CON_PRE_max(:,col_AV_strain_msc_GM)); 
        CON_PRE_strain_msc_GM_SD = nanstd(CON_PRE_max(:,col_AV_strain_msc_GM));
        CON_PRE_strain_msc_SOL_mean = nanmean(CON_PRE_max(:,col_AV_strain_msc_SOL)); 
        CON_PRE_strain_msc_SOL_SD = nanstd(CON_PRE_max(:,col_AV_strain_msc_SOL));

        CON_PRE_elong_msc_GM_licht_mean = nanmean(CON_PRE_max(:,col_AV_elong_msc_GM_fuku)); 
        CON_PRE_elong_msc_GM_licht_SD = nanstd(CON_PRE_max(:,col_AV_elong_msc_GM_fuku));
        CON_PRE_elong_tend_GM_licht_mean = nanmean(CON_PRE_max(:,col_AV_elong_SEE_fuku)); 
        CON_PRE_elong_tend_GM_licht_SD = nanstd(CON_PRE_max(:,col_AV_elong_SEE_fuku));
        CON_PRE_strain_msc_GM_licht_mean = nanmean(CON_PRE_max(:,col_AV_strain_msc_GM_fuku)); 
        CON_PRE_strain_msc_GM_licht_SD = nanstd(CON_PRE_max(:,col_AV_strain_msc_GM_fuku));
        CON_PRE_strain_tend_GM_licht_mean = nanmean(CON_PRE_max(:,col_AV_strain_SEE_fuku)); 
        CON_PRE_strain_tend_GM_licht_SD = nanstd(CON_PRE_max(:,col_AV_strain_SEE_fuku));

        CON_PRE_length_GMfas_licht_mean = nanmean(CON_PRE_max(:,col_AV_len_GMfas_licht)); 
        CON_PRE_length_GMfas_licht_SD = nanstd(CON_PRE_max(:,col_AV_len_GMfas_licht)); 
        CON_PRE_pennation_GMfas_licht_mean = nanmean(CON_PRE_max(:,col_AV_pennation_GMfas_licht)); 
        CON_PRE_pennation_GMfas_licht_SD = nanstd(CON_PRE_max(:,col_AV_pennation_GMfas_licht)); 
        CON_PRE_elong_GMfas_licht_mean = nanmean(CON_PRE_max(:,col_AV_elong_GMfas_licht)); 
        CON_PRE_elong_GMfas_licht_SD = nanstd(CON_PRE_max(:,col_AV_elong_GMfas_licht)); 
        CON_PRE_strain_GMfas_licht_mean = nanmean(CON_PRE_max(:,col_AV_strain_GMfas_licht)); 
        CON_PRE_strain_GMfas_licht_SD = nanstd(CON_PRE_max(:,col_AV_strain_GMfas_licht)); 

        CON_PRE_length_msc_GM_licht_mean = nanmean(CON_PRE_max(:,col_AV_len_msc_GM_fuku)); 
        CON_PRE_length_msc_GM_licht_SD = nanstd(CON_PRE_max(:,col_AV_len_msc_GM_fuku));
        CON_PRE_length_tend_GM_licht_mean = nanmean(CON_PRE_max(:,col_AV_len_SEE_fuku)); 
        CON_PRE_length_tend_GM_licht_SD = nanstd(CON_PRE_max(:,col_AV_len_SEE_fuku));

        CON_PRE_length_norm_MTU_mean = nanmean(CON_PRE_max(:,col_AV_norm_len_leg)); % L leg
        CON_PRE_length_norm_MTU_SD = nanstd(CON_PRE_max(:,col_AV_norm_len_leg));
        CON_PRE_elong_norm_MTU_mean = nanmean(CON_PRE_max(:,col_AV_norm_elong_leg)); % elong leg
        CON_PRE_elong_norm_MTU_SD = nanstd(CON_PRE_max(:,col_AV_norm_elong_leg));

        CON_PRE_length_norm_msc_GM_licht_mean = nanmean(CON_PRE_max(:,col_AV_norm_len_msc_GM_fuku)); 
        CON_PRE_length_norm_msc_GM_licht_SD = nanstd(CON_PRE_max(:,col_AV_norm_len_msc_GM_fuku));
        CON_PRE_elong_norm_msc_GM_licht_mean = nanmean(CON_PRE_max(:,col_AV_norm_elong_msc_GM_fuku)); 
        CON_PRE_elong_norm_msc_GM_licht_SD = nanstd(CON_PRE_max(:,col_AV_norm_elong_msc_GM_fuku));

        CON_PRE_length_norm_SEE_licht_mean = nanmean(CON_PRE_max(:,col_AV_norm_len_SEE_fuku)); 
        CON_PRE_length_norm_SEE_licht_SD = nanstd(CON_PRE_max(:,col_AV_norm_len_SEE_fuku));
        CON_PRE_elong_norm_SEE_licht_mean = nanmean(CON_PRE_max(:,col_AV_norm_elong_SEE_fuku)); 
        CON_PRE_elong_norm_SEE_licht_SD = nanstd(CON_PRE_max(:,col_AV_norm_elong_SEE_fuku));

        CON_PRE_elong_norm_percent_msc_GM_licht_mean = nanmean(CON_PRE_max(:,col_AV_norm_percent_elong_msc_GM_fuku)); 
        CON_PRE_elong_norm_percent_msc_GM_licht_SD = nanstd(CON_PRE_max(:,col_AV_norm_percent_elong_msc_GM_fuku));
        CON_PRE_elong_norm_percent_SEE_licht_mean = nanmean(CON_PRE_max(:,col_AV_norm_percent_elong_SEE_fuku)); 
        CON_PRE_elong_norm_percent_SEE_licht_SD = nanstd(CON_PRE_max(:,col_AV_norm_percent_elong_SEE_fuku));

        CON_PRE_length_norm_GMfas_licht_mean = nanmean(CON_PRE_max(:,col_AV_norm_len_GMfas_licht)); 
        CON_PRE_length_norm_GMfas_licht_SD = nanstd(CON_PRE_max(:,col_AV_norm_len_GMfas_licht)); 
        CON_PRE_elong_norm_GMfas_licht_mean = nanmean(CON_PRE_max(:,col_AV_norm_elong_GMfas_licht)); 
        CON_PRE_elong_norm_GMfas_licht_SD = nanstd(CON_PRE_max(:,col_AV_norm_elong_GMfas_licht)); 
        
        CON_PRE_torque_mean = nanmean(CON_PRE_max(:,col_AV_T));
        CON_PRE_torque_SD = nanstd(CON_PRE_max(:,col_AV_T));
        % determine common angle range
        CON_PRE_common_ROM = min(CON_PRE_max(:,col_AV_angle));
    end

    % CON POST
    if CON_POST_count > 0
        % preallocate array
        CON_POST_max(CON_POST_count,n_o_array_elements) = zeros;
  %      CON_POST_max_norm(CON_POST_count,n_o_array_elements) = zeros;
        % collect variables per subject
        for i = 1:CON_POST_count % per subject
            for j = 1:n_o_array_elements % per element in arrays
                % OLD: Max values
%                CON_POST_max(i,j) = max(CON_POST_angle_vars{1,i}(:,j));
%                CON_POST_max_norm(i,j) = max(CON_POST_angle_vars_norm{1,i}(:,j));
                % NEW: end of array values
                CON_POST_max(i,j) = max(CON_POST_angle_vars{1,i}(end,j));
%                CON_POST_max_norm(i,j) = max(CON_POST_angle_vars_norm{1,i}(end,j));
            end
        end
        % calculate mean and SD of max values across subjects
        CON_POST_ROM_mean = nanmean(CON_POST_max(:,col_AV_angle));
        CON_POST_ROM_SD = nanstd(CON_POST_max(:,col_AV_angle));
        CON_POST_F_mean = nanmean(CON_POST_max(:,col_AV_F));
        CON_POST_F_SD = nanstd(CON_POST_max(:,col_AV_F));
        CON_POST_EMG_gm_mean = nanmean(CON_POST_max(:,col_AV_EMG_gm));
        CON_POST_EMG_gm_SD = nanstd(CON_POST_max(:,col_AV_EMG_gm));
        CON_POST_EMG_gl_mean = nanmean(CON_POST_max(:,col_AV_EMG_gl));
        CON_POST_EMG_gl_SD = nanstd(CON_POST_max(:,col_AV_EMG_gl));
        CON_POST_EMG_sol_mean = nanmean(CON_POST_max(:,col_AV_EMG_sol));
        CON_POST_EMG_sol_SD = nanstd(CON_POST_max(:,col_AV_EMG_sol));
        CON_POST_EMG_mean_mean = nanmean(CON_POST_max(:,col_AV_EMG_mean));
        CON_POST_EMG_mean_SD = nanstd(CON_POST_max(:,col_AV_EMG_mean));

        CON_POST_elong_AT_mean = nanmean(CON_POST_max(:,col_AV_elong_AT)); % elong AT
        CON_POST_elong_AT_SD = nanstd(CON_POST_max(:,col_AV_elong_AT));
        CON_POST_elong_GMtend_mean = nanmean(CON_POST_max(:,col_AV_elong_GMtend)); % elong GM tend
        CON_POST_elong_GMtend_SD = nanstd(CON_POST_max(:,col_AV_elong_GMtend));
        CON_POST_elong_MTU_mean = nanmean(CON_POST_max(:,col_AV_elong_leg)); % elong leg
        CON_POST_elong_MTU_SD = nanstd(CON_POST_max(:,col_AV_elong_leg));
        CON_POST_displ_GMFAS_mean = nanmean(CON_POST_max(:,col_AV_elong_GMFAS));  % GMFAS displ 
        CON_POST_displ_GMFAS_SD = nanstd(CON_POST_max(:,col_AV_elong_GMFAS));
        CON_POST_elong_GMapo_mean = nanmean(CON_POST_max(:,col_AV_elong_GMapo)); 
        CON_POST_elong_GMapo_SD = nanstd(CON_POST_max(:,col_AV_elong_GMapo));
        CON_POST_elong_msc_GM_mean = nanmean(CON_POST_max(:,col_AV_elong_msc_GM)); % GM msc
        CON_POST_elong_msc_GM_SD = nanstd(CON_POST_max(:,col_AV_elong_msc_GM));
        CON_POST_elong_msc_SOL_mean = nanmean(CON_POST_max(:,col_AV_elong_msc_SOL)); 
        CON_POST_elong_msc_SOL_SD = nanstd(CON_POST_max(:,col_AV_elong_msc_SOL));

        CON_POST_L_at_SOL_mean = nanmean(CON_POST_max(:,col_AV_len_AT)); % L AT
        CON_POST_L_at_SOL_SD = nanstd(CON_POST_max(:,col_AV_len_AT));
        CON_POST_L_at_GM_mean = nanmean(CON_POST_max(:,col_AV_len_GMtend)); % L GM tend
        CON_POST_L_at_GM_SD = nanstd(CON_POST_max(:,col_AV_len_GMtend));
        CON_POST_L_MTU_mean = nanmean(CON_POST_max(:,col_AV_len_leg)); % L leg
        CON_POST_L_MTU_SD = nanstd(CON_POST_max(:,col_AV_len_leg));
        % no L GMFAS
        CON_POST_L_GMapo_mean = nanmean(CON_POST_max(:,col_AV_len_GMapo)); 
        CON_POST_L_GMapo_SD = nanstd(CON_POST_max(:,col_AV_len_GMapo));
        CON_POST_L_msc_GM_mean = nanmean(CON_POST_max(:,col_AV_len_msc_GM)); 
        CON_POST_L_msc_GM_SD = nanstd(CON_POST_max(:,col_AV_len_msc_GM));
        CON_POST_L_msc_SOL_mean = nanmean(CON_POST_max(:,col_AV_len_msc_SOL)); 
        CON_POST_L_msc_SOL_SD = nanstd(CON_POST_max(:,col_AV_len_msc_SOL));

        CON_POST_strain_at_SOL_mean = nanmean(CON_POST_max(:,col_AV_strain_AT)); % L AT
        CON_POST_strain_at_SOL_SD = nanstd(CON_POST_max(:,col_AV_strain_AT));
        CON_POST_strain_at_GM_mean = nanmean(CON_POST_max(:,col_AV_strain_GMtend)); % L GM tend
        CON_POST_strain_at_GM_SD = nanstd(CON_POST_max(:,col_AV_strain_GMtend));
        CON_POST_strain_MTU_mean = nanmean(CON_POST_max(:,col_AV_strain_leg)); % L calf
        CON_POST_strain_MTU_SD = nanstd(CON_POST_max(:,col_AV_strain_leg));
        % no strain GMfas
        CON_POST_strain_GMapo_mean = nanmean(CON_POST_max(:,col_AV_strain_GMapo)); 
        CON_POST_strain_GMapo_SD = nanstd(CON_POST_max(:,col_AV_strain_GMapo));
        CON_POST_strain_msc_GM_mean = nanmean(CON_POST_max(:,col_AV_strain_msc_GM)); 
        CON_POST_strain_msc_GM_SD = nanstd(CON_POST_max(:,col_AV_strain_msc_GM));
        CON_POST_strain_msc_SOL_mean = nanmean(CON_POST_max(:,col_AV_strain_msc_SOL)); 
        CON_POST_strain_msc_SOL_SD = nanstd(CON_POST_max(:,col_AV_strain_msc_SOL));

        CON_POST_elong_msc_GM_licht_mean = nanmean(CON_POST_max(:,col_AV_elong_msc_GM_fuku)); 
        CON_POST_elong_msc_GM_licht_SD = nanstd(CON_POST_max(:,col_AV_elong_msc_GM_fuku));
        CON_POST_elong_tend_GM_licht_mean = nanmean(CON_POST_max(:,col_AV_elong_SEE_fuku)); 
        CON_POST_elong_tend_GM_licht_SD = nanstd(CON_POST_max(:,col_AV_elong_SEE_fuku));
        CON_POST_strain_msc_GM_licht_mean = nanmean(CON_POST_max(:,col_AV_strain_msc_GM_fuku)); 
        CON_POST_strain_msc_GM_licht_SD = nanstd(CON_POST_max(:,col_AV_strain_msc_GM_fuku));
        CON_POST_strain_tend_GM_licht_mean = nanmean(CON_POST_max(:,col_AV_strain_SEE_fuku)); 
        CON_POST_strain_tend_GM_licht_SD = nanstd(CON_POST_max(:,col_AV_strain_SEE_fuku));

        CON_POST_length_GMfas_licht_mean = nanmean(CON_POST_max(:,col_AV_len_GMfas_licht)); 
        CON_POST_length_GMfas_licht_SD = nanstd(CON_POST_max(:,col_AV_len_GMfas_licht)); 
        CON_POST_pennation_GMfas_licht_mean = nanmean(CON_POST_max(:,col_AV_pennation_GMfas_licht)); 
        CON_POST_pennation_GMfas_licht_SD = nanstd(CON_POST_max(:,col_AV_pennation_GMfas_licht)); 
        CON_POST_elong_GMfas_licht_mean = nanmean(CON_POST_max(:,col_AV_elong_GMfas_licht)); 
        CON_POST_elong_GMfas_licht_SD = nanstd(CON_POST_max(:,col_AV_elong_GMfas_licht)); 
        CON_POST_strain_GMfas_licht_mean = nanmean(CON_POST_max(:,col_AV_strain_GMfas_licht)); 
        CON_POST_strain_GMfas_licht_SD = nanstd(CON_POST_max(:,col_AV_strain_GMfas_licht)); 

        CON_POST_length_msc_GM_licht_mean = nanmean(CON_POST_max(:,col_AV_len_msc_GM_fuku)); 
        CON_POST_length_msc_GM_licht_SD = nanstd(CON_POST_max(:,col_AV_len_msc_GM_fuku));
        CON_POST_length_tend_GM_licht_mean = nanmean(CON_POST_max(:,col_AV_len_SEE_fuku)); 
        CON_POST_length_tend_GM_licht_SD = nanstd(CON_POST_max(:,col_AV_len_SEE_fuku));

        CON_POST_length_norm_MTU_mean = nanmean(CON_POST_max(:,col_AV_norm_len_leg)); % L leg
        CON_POST_length_norm_MTU_SD = nanstd(CON_POST_max(:,col_AV_norm_len_leg));
        CON_POST_elong_norm_MTU_mean = nanmean(CON_POST_max(:,col_AV_norm_elong_leg)); % elong leg
        CON_POST_elong_norm_MTU_SD = nanstd(CON_POST_max(:,col_AV_norm_elong_leg));

        CON_POST_length_norm_msc_GM_licht_mean = nanmean(CON_POST_max(:,col_AV_norm_len_msc_GM_fuku)); 
        CON_POST_length_norm_msc_GM_licht_SD = nanstd(CON_POST_max(:,col_AV_norm_len_msc_GM_fuku));
        CON_POST_elong_norm_msc_GM_licht_mean = nanmean(CON_POST_max(:,col_AV_norm_elong_msc_GM_fuku)); 
        CON_POST_elong_norm_msc_GM_licht_SD = nanstd(CON_POST_max(:,col_AV_norm_elong_msc_GM_fuku));

        CON_POST_length_norm_SEE_licht_mean = nanmean(CON_POST_max(:,col_AV_norm_len_SEE_fuku)); 
        CON_POST_length_norm_SEE_licht_SD = nanstd(CON_POST_max(:,col_AV_norm_len_SEE_fuku));
        CON_POST_elong_norm_SEE_licht_mean = nanmean(CON_POST_max(:,col_AV_norm_elong_SEE_fuku)); 
        CON_POST_elong_norm_SEE_licht_SD = nanstd(CON_POST_max(:,col_AV_norm_elong_SEE_fuku));

        CON_POST_elong_norm_percent_msc_GM_licht_mean = nanmean(CON_POST_max(:,col_AV_norm_percent_elong_msc_GM_fuku)); 
        CON_POST_elong_norm_percent_msc_GM_licht_SD = nanstd(CON_POST_max(:,col_AV_norm_percent_elong_msc_GM_fuku));
        CON_POST_elong_norm_percent_SEE_licht_mean = nanmean(CON_POST_max(:,col_AV_norm_percent_elong_SEE_fuku)); 
        CON_POST_elong_norm_percent_SEE_licht_SD = nanstd(CON_POST_max(:,col_AV_norm_percent_elong_SEE_fuku));

        CON_POST_length_norm_GMfas_licht_mean = nanmean(CON_POST_max(:,col_AV_norm_len_GMfas_licht)); 
        CON_POST_length_norm_GMfas_licht_SD = nanstd(CON_POST_max(:,col_AV_norm_len_GMfas_licht)); 
        CON_POST_elong_norm_GMfas_licht_mean = nanmean(CON_POST_max(:,col_AV_norm_elong_GMfas_licht)); 
        CON_POST_elong_norm_GMfas_licht_SD = nanstd(CON_POST_max(:,col_AV_norm_elong_GMfas_licht)); 
        
        CON_POST_torque_mean = nanmean(CON_POST_max(:,col_AV_T));
        CON_POST_torque_SD = nanstd(CON_POST_max(:,col_AV_T));
        % determine common angle range
        CON_POST_common_ROM = min(CON_POST_max(:,col_AV_angle));
    end
    
    
    %% GROUP figures - create AVERAGE ARRAYS 

    if STR_PRE_count > 0
        %%% average ABSOLUTE arrays, up to all subjects' COMMON MAX ROM, for force, elong, EMG

        % preallocate
        len = 10000;
        for i = 1:STR_PRE_count
            if length(STR_PRE_angle_vars{:,i}) < len
                len = length(STR_PRE_angle_vars{:,i});
            end
        end
        STR_PRE_angle_vars_mean_tmp(len,n_o_array_elements,STR_PRE_count) = zeros;

        % STR_PRE_angle_vars has same angles (column 1) for all subjects, so max common ROM will be on the same line for all subjects
        loc_end = find(STR_PRE_angle_vars{1,STR_PRE_count}(:,col_AV_angle) >= (STR_PRE_common_ROM - 0.00001), 1, 'first'); % using last subject (STR_PRE_count) - could use any, angle is the same in all
        for i = 1:STR_PRE_count
            STR_PRE_angle_vars_mean_tmp(:,:,i) = STR_PRE_angle_vars{i}(1:loc_end,:);
        end
        STR_PRE_angle_vars_mean = nanmean(STR_PRE_angle_vars_mean_tmp, 3);
        %STR_PRE_angle_vars_SD = nanstd(STR_PRE_angle_vars_mean_tmp,1,3);

        %%% clean up
        clear STR_PRE_angle_vars_mean_tmp % STR_PRE_angle_vars_norm_mean_tmp
    end


    if STR_POST_count > 0
        %%% average ABSOLUTE arrays, up to all subjects' COMMON MAX ROM, for force, elong, EMG

        % preallocate
        len = 10000;
        for i = 1:STR_POST_count
            if length(STR_POST_angle_vars{:,i}) < len
                len = length(STR_POST_angle_vars{:,i});
            end
        end
        STR_POST_angle_vars_mean_tmp(len,n_o_array_elements,STR_POST_count) = zeros;

        % STR_POST_angle_vars has same angles (column 1) for all subjects, so max common ROM will be on the same line for all subjects
        loc_end = find(STR_POST_angle_vars{1,STR_POST_count}(:,col_AV_angle) >= (STR_POST_common_ROM - 0.00001), 1, 'first'); % using last subject (STR_POST_count) - could use any, angle is the same in all
        for i = 1:STR_POST_count
            STR_POST_angle_vars_mean_tmp(:,:,i) = STR_POST_angle_vars{i}(1:loc_end,:);
        end
        STR_POST_angle_vars_mean = nanmean(STR_POST_angle_vars_mean_tmp, 3);
        %STR_POST_angle_vars_SD = nanstd(STR_POST_angle_vars_mean_tmp,1,3);

        %%% clean up
        clear STR_POST_angle_vars_mean_tmp % STR_POST_angle_vars_norm_mean_tmp
    end


    if CON_PRE_count > 0
        %%% average ABSOLUTE arrays, up to all subjects' COMMON MAX ROM, for force, elong, EMG

        % preallocate
        len = 10000;
        for i = 1:CON_PRE_count
            if length(CON_PRE_angle_vars{:,i}) < len
                len = length(CON_PRE_angle_vars{:,i});
            end
        end
        CON_PRE_angle_vars_mean_tmp(len,n_o_array_elements,CON_PRE_count) = zeros;

        % CON_PRE_angle_vars has same angles (column 1) for all subjects, so max common ROM will be on the same line for all subjects
        loc_end = find(CON_PRE_angle_vars{1,CON_PRE_count}(:,col_AV_angle) >= (CON_PRE_common_ROM - 0.00001), 1, 'first'); % using last subject (CON_PRE_count) - could use any, angle is the same in all
        for i = 1:CON_PRE_count
            CON_PRE_angle_vars_mean_tmp(:,:,i) = CON_PRE_angle_vars{i}(1:loc_end,:);
        end
        CON_PRE_angle_vars_mean = nanmean(CON_PRE_angle_vars_mean_tmp, 3);
        %CON_PRE_angle_vars_SD = nanstd(CON_PRE_angle_vars_mean_tmp,1,3);

        %%% clean up
        clear CON_PRE_angle_vars_mean_tmp % CON_PRE_angle_vars_norm_mean_tmp
    end


    if CON_POST_count > 0
        %%% average ABSOLUTE arrays, up to all subjects' COMMON MAX ROM, for force, elong, EMG

        % preallocate
        len = 10000;
        for i = 1:CON_POST_count
            if length(CON_POST_angle_vars{:,i}) < len
                len = length(CON_POST_angle_vars{:,i});
            end
        end
        CON_POST_angle_vars_mean_tmp(len,n_o_array_elements,CON_POST_count) = zeros;

        % CON_POST_angle_vars has same angles (column 1) for all subjects, so max common ROM will be on the same line for all subjects
        loc_end = find(CON_POST_angle_vars{1,CON_POST_count}(:,col_AV_angle) >= (CON_POST_common_ROM - 0.00001), 1, 'first'); % using last subject (CON_POST_count) - could use any, angle is the same in all
        for i = 1:CON_POST_count
            CON_POST_angle_vars_mean_tmp(:,:,i) = CON_POST_angle_vars{i}(1:loc_end,:);
        end
        CON_POST_angle_vars_mean = nanmean(CON_POST_angle_vars_mean_tmp, 3);
        %CON_POST_angle_vars_SD = nanstd(CON_POST_angle_vars_mean_tmp,1,3);

        %%% clean up
        clear CON_POST_angle_vars_mean_tmp % CON_POST_angle_vars_norm_mean_tmp
    end

    
    %% GROUP figures - CREATE PLOTS /// GROUP PLOTS 
    
        if plot_check
            %% FORCE-angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('force vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_F),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_F),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_F),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_F),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_F_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_F_mean, STR_PRE_F_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_F_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_F_mean, STR_POST_F_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_F_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_F_mean, CON_PRE_F_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_F_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_F_mean, CON_POST_F_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_force)
                xlabel(txt_gonio)
                ylabel('Force (N)')
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('force vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_F),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_F))
                end
                axis(axis_force)
                xlabel(txt_gonio)
                ylabel('Force (N)')
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('force vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_F),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_F))
                end
                axis(axis_force)
                xlabel(txt_gonio)
                ylabel('Force (N)')
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % the following plots require that trials are analysed in systematic order (does not check which subject-numbers to group)
            % - but will not be plot unless equal amount of 4 legs/timepoints - rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND force vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_F),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_F),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_F),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_F),'b','LineStyle','-','LineWidth',1)
                    axis(axis_force)
                    xlabel(txt_gonio)
                    ylabel('Force (N)')
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
           
%             if plot_check
%                 plottitle = horzcat('force vs angle - 4 NORMALIZED');
%                 figure('Name',plottitle)
%                 hold on
%                 plot(STR_PRE_angle_vars_norm_mean(:,col_AV_angle), STR_PRE_angle_vars_norm_mean(:,col_AV_F),'Color',col_lightred,'LineStyle','--','LineWidth',1)
%                 plot(STR_POST_angle_vars_norm_mean(:,col_AV_angle), STR_POST_angle_vars_norm_mean(:,col_AV_F),'r','LineStyle','-','LineWidth',1)
%                 plot(CON_PRE_angle_vars_norm_mean(:,col_AV_angle), CON_PRE_angle_vars_norm_mean(:,col_AV_F),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
%                 plot(CON_POST_angle_vars_norm_mean(:,col_AV_angle), CON_POST_angle_vars_norm_mean(:,col_AV_F),'b','LineStyle','-','LineWidth',1)
%                 
%                 axis(axis_PP)
%                 xlabel('Gonio angle (% of ind max)')
%                 ylabel('Force (% of ind max)')
%                 title(plottitle,'Interpreter', 'none')
%                 legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
%                 print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
%             end

            %% TORQUE-angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('torque vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_T),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_T),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_T),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_T),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_torque_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_torque_mean, STR_PRE_torque_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_torque_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_torque_mean, STR_POST_torque_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_torque_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_torque_mean, CON_PRE_torque_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_torque_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_torque_mean, CON_POST_torque_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_torque)
                xlabel(txt_gonio)
                ylabel('Torque (Nm)')
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('torque vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_T),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_T))
                end
                axis(axis_torque)
                xlabel(txt_gonio)
                ylabel('Torque (Nm)')
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('torque vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_T),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_T))
                end
                axis(axis_torque)
                xlabel(txt_gonio)
                ylabel('Torque (Nm)')
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND torque vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_T),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_T),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_T),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_T),'b','LineStyle','-','LineWidth',1)
                axis(axis_torque)
                xlabel(txt_gonio)
                ylabel('Torque (Nm)')
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% FREE AT: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('free AT length vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_len_AT),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_len_AT),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_len_AT),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_len_AT),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_L_at_SOL_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_L_at_SOL_mean, STR_PRE_L_at_SOL_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_L_at_SOL_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_L_at_SOL_mean, STR_POST_L_at_SOL_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_L_at_SOL_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_L_at_SOL_mean, CON_PRE_L_at_SOL_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_L_at_SOL_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_L_at_SOL_mean, CON_POST_L_at_SOL_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_len_AT)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('free AT length vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_len_AT),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_len_AT))
                end
                axis(axis_len_AT)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('free AT length vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_len_AT),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_len_AT))
                end
                axis(axis_len_AT)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND free AT length vs angle - ',CON_PRE_ID{i}(1:6) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_len_AT),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_len_AT),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_len_AT),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_len_AT),'b','LineStyle','-','LineWidth',1)
                    axis(axis_len_AT)
                    xlabel(txt_gonio)
                    ylabel(txt_length)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% FREE AT: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('free AT elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_elong_AT),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_elong_AT),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_elong_AT),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_elong_AT),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_AT_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_AT_mean, STR_PRE_elong_AT_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_AT_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_AT_mean, STR_POST_elong_AT_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_AT_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_AT_mean, CON_PRE_elong_AT_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_AT_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_AT_mean, CON_POST_elong_AT_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_el_AT)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('free AT elongation vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_elong_AT),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_elong_AT))
                end
                axis(axis_el_AT)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('free AT elongation vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_elong_AT),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_elong_AT))
                end
                axis(axis_el_AT)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND free AT elongation vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_elong_AT),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_elong_AT),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_elong_AT),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_elong_AT),'b','LineStyle','-','LineWidth',1)
                    axis(axis_el_AT)
                    xlabel(txt_gonio)
                    ylabel(txt_elong)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end

            %% FREE AT: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('free AT strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_strain_AT),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_strain_AT),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_strain_AT),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_strain_AT),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_strain_at_SOL_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_strain_at_SOL_mean, STR_PRE_strain_at_SOL_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_strain_at_SOL_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_strain_at_SOL_mean, STR_POST_strain_at_SOL_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_strain_at_SOL_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_strain_at_SOL_mean, CON_PRE_strain_at_SOL_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_strain_at_SOL_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_strain_at_SOL_mean, CON_POST_strain_at_SOL_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_str_AT)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('free AT strain vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_strain_AT),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_strain_AT))
                end
                axis(axis_str_AT)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('free AT strain vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_strain_AT),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_strain_AT))
                end
                axis(axis_str_AT)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND free AT strain vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_strain_AT),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_strain_AT),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_strain_AT),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_strain_AT),'b','LineStyle','-','LineWidth',1)
                axis(axis_str_AT)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            close all
            
            %% GM TENDON: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM tendon (from calc) length vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_len_GMtend),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_len_GMtend),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_len_GMtend),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_len_GMtend),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_L_at_GM_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_L_at_GM_mean, STR_PRE_L_at_GM_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_L_at_GM_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_L_at_GM_mean, STR_POST_L_at_GM_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_L_at_GM_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_L_at_GM_mean, CON_PRE_L_at_GM_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_L_at_GM_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_L_at_GM_mean, CON_POST_L_at_GM_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_len_GMtend)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southeast')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('GM tendon (from calc) length vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_len_GMtend),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_len_GMtend))
                end
                axis(axis_len_GMtend)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('GM tendon (from calc) length vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_len_GMtend),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_len_GMtend))
                end
                axis(axis_len_GMtend)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM tendon (from calc) length vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_len_GMtend),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_len_GMtend),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_len_GMtend),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_len_GMtend),'b','LineStyle','-','LineWidth',1)
                axis(axis_len_GMtend)
                xlabel(txt_gonio)
                ylabel(txt_length)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end

            %% GM TENDON: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM tendon (from calc) elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_elong_GMtend),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_elong_GMtend),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_elong_GMtend),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_elong_GMtend),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_GMtend_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_GMtend_mean, STR_PRE_elong_GMtend_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_GMtend_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_GMtend_mean, STR_POST_elong_GMtend_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_GMtend_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_GMtend_mean, CON_PRE_elong_GMtend_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_GMtend_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_GMtend_mean, CON_POST_elong_GMtend_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_el_GMtend)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('GM tendon (from calc) elongation vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_elong_GMtend),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_elong_GMtend))
                end
                axis(axis_el_GMtend)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('GM tendon (from calc) elongation vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_elong_GMtend),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_elong_GMtend))
                end
                axis(axis_el_GMtend)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM tendon (from calc) elongation vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_elong_GMtend),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_elong_GMtend),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_elong_GMtend),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_elong_GMtend),'b','LineStyle','-','LineWidth',1)
                axis(axis_el_GMtend)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end

            %% GM TENDON: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM tendon (from calc) strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_strain_GMtend),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_strain_GMtend),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_strain_GMtend),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_strain_GMtend),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_strain_at_GM_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_strain_at_GM_mean, STR_PRE_strain_at_GM_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_strain_at_GM_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_strain_at_GM_mean, STR_POST_strain_at_GM_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_strain_at_GM_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_strain_at_GM_mean, CON_PRE_strain_at_GM_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_strain_at_GM_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_strain_at_GM_mean, CON_POST_strain_at_GM_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_str_GMtend)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('GM tendon (from calc) strain vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_strain_GMtend),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_strain_GMtend))
                end
                axis(axis_str_GMtend)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('GM tendon (from calc) strain vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_strain_GMtend),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_strain_GMtend))
                end
                axis(axis_str_GMtend)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM tendon (from calc) strain vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_strain_GMtend),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_strain_GMtend),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_strain_GMtend),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_strain_GMtend),'b','LineStyle','-','LineWidth',1)
                axis(axis_str_GMtend)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% GM APO: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) length vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_len_GMapo),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_len_GMapo),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_len_GMapo),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_len_GMapo),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_L_GMapo_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_L_GMapo_mean, STR_PRE_L_GMapo_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_L_GMapo_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_L_GMapo_mean, STR_POST_L_GMapo_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_L_GMapo_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_L_GMapo_mean, CON_PRE_L_GMapo_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_L_GMapo_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_L_GMapo_mean, CON_POST_L_GMapo_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_len_GMapo)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southeast')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) length vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_len_GMapo),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_len_GMapo))
                end
                axis(axis_len_GMapo)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) length vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_len_GMapo),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_len_GMapo))
                end
                axis(axis_len_GMapo)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM apo (SOL ins-GM ins) length vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_len_GMapo),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_len_GMapo),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_len_GMapo),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_len_GMapo),'b','LineStyle','-','LineWidth',1)
                axis(axis_len_GMapo)
                xlabel(txt_gonio)
                ylabel(txt_length)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% GM APO: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_elong_GMapo),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_elong_GMapo),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_elong_GMapo),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_elong_GMapo),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_GMapo_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_GMapo_mean, STR_PRE_elong_GMapo_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_GMapo_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_GMapo_mean, STR_POST_elong_GMapo_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_GMapo_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_GMapo_mean, CON_PRE_elong_GMapo_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_GMapo_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_GMapo_mean, CON_POST_elong_GMapo_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_el_GMapo)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northeast')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) elongation vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_elong_GMapo),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_elong_GMapo))
                end
                axis(axis_el_GMapo)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) elongation vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_elong_GMapo),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_elong_GMapo))
                end
                axis(axis_el_GMapo)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM apo (SOL ins-GM ins) elongation vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_elong_GMapo),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_elong_GMapo),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_elong_GMapo),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_elong_GMapo),'b','LineStyle','-','LineWidth',1)
                axis(axis_el_GMapo)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% GM APO: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_strain_GMapo),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_strain_GMapo),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_strain_GMapo),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_strain_GMapo),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_strain_GMapo_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_strain_GMapo_mean, STR_PRE_strain_GMapo_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_strain_GMapo_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_strain_GMapo_mean, STR_POST_strain_GMapo_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_strain_GMapo_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_strain_GMapo_mean, CON_PRE_strain_GMapo_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_strain_GMapo_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_strain_GMapo_mean, CON_POST_strain_GMapo_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_str_GMapo)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northeast')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) strain vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_strain_GMapo),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_strain_GMapo))
                end
                axis(axis_str_GMapo)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) strain vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_strain_GMapo),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_strain_GMapo))
                end
                axis(axis_str_GMapo)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM apo (SOL ins-GM ins) strain vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_strain_GMapo),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_strain_GMapo),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_strain_GMapo),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_strain_GMapo),'b','LineStyle','-','LineWidth',1)
                axis(axis_str_GMapo)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            close all
            
            %% SOL MUSCLE portion: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle SOL length vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_len_msc_SOL),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_len_msc_SOL),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_len_msc_SOL),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_len_msc_SOL),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_L_msc_SOL_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_L_msc_SOL_mean, STR_PRE_L_msc_SOL_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_L_msc_SOL_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_L_msc_SOL_mean, STR_POST_L_msc_SOL_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_L_msc_SOL_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_L_msc_SOL_mean, CON_PRE_L_msc_SOL_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_L_msc_SOL_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_L_msc_SOL_mean, CON_POST_L_msc_SOL_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_len_SOL)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southeast')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('muscle SOL length vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_len_msc_SOL),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_len_msc_SOL))
                end
                axis(axis_len_SOL)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('muscle SOL length vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_len_msc_SOL),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_len_msc_SOL))
                end
                axis(axis_len_SOL)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle SOL length vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_len_msc_SOL),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_len_msc_SOL),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_len_msc_SOL),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_len_msc_SOL),'b','LineStyle','-','LineWidth',1)
                axis(axis_len_SOL)
                xlabel(txt_gonio)
                ylabel(txt_length)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Southeast')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% SOL MUSCLE portion: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle SOL elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_elong_msc_SOL),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_elong_msc_SOL),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_elong_msc_SOL),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_elong_msc_SOL),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_msc_SOL_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_msc_SOL_mean, STR_PRE_elong_msc_SOL_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_msc_SOL_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_msc_SOL_mean, STR_POST_elong_msc_SOL_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_msc_SOL_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_msc_SOL_mean, CON_PRE_elong_msc_SOL_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_msc_SOL_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_msc_SOL_mean, CON_POST_elong_msc_SOL_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_el_SOL)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('muscle SOL elongation vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_elong_msc_SOL),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_elong_msc_SOL))
                end
                axis(axis_el_SOL)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('muscle SOL elongation vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_elong_msc_SOL),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_elong_msc_SOL))
                end
                axis(axis_el_SOL)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle SOL elongation vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_elong_msc_SOL),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_elong_msc_SOL),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_elong_msc_SOL),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_elong_msc_SOL),'b','LineStyle','-','LineWidth',1)
                axis(axis_el_SOL)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% SOL MUSCLE portion: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle SOL strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_strain_msc_SOL),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_strain_msc_SOL),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_strain_msc_SOL),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_strain_msc_SOL),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_strain_msc_SOL_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_strain_msc_SOL_mean, STR_PRE_strain_msc_SOL_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_strain_msc_SOL_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_strain_msc_SOL_mean, STR_POST_strain_msc_SOL_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_strain_msc_SOL_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_strain_msc_SOL_mean, CON_PRE_strain_msc_SOL_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_strain_msc_SOL_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_strain_msc_SOL_mean, CON_POST_strain_msc_SOL_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_str_SOL)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('muscle SOL strain vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_strain_msc_SOL),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_strain_msc_SOL))
                end
                axis(axis_str_SOL)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('muscle SOL strain vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_strain_msc_SOL),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_strain_msc_SOL))
                end
                axis(axis_str_SOL)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle SOL strain vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_strain_msc_SOL),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_strain_msc_SOL),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_strain_msc_SOL),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_strain_msc_SOL),'b','LineStyle','-','LineWidth',1)
                axis(axis_str_SOL)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end

            %% GM MUSCLE portion: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_len_msc_GM),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_len_msc_GM),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_len_msc_GM),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_len_msc_GM),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_L_msc_GM_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_L_msc_GM_mean, STR_PRE_L_msc_GM_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_L_msc_GM_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_L_msc_GM_mean, STR_POST_L_msc_GM_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_L_msc_GM_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_L_msc_GM_mean, CON_PRE_L_msc_GM_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_L_msc_GM_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_L_msc_GM_mean, CON_POST_L_msc_GM_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_len_GMmsc)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southeast')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_len_msc_GM),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_len_msc_GM))
                end
                axis(axis_len_GMmsc)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_len_msc_GM),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_len_msc_GM))
                end
                axis(axis_len_GMmsc)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle GM (ins-knee) vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_len_msc_GM),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_len_msc_GM),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_len_msc_GM),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_len_msc_GM),'b','LineStyle','-','LineWidth',1)
                axis(axis_len_GMmsc)
                xlabel(txt_gonio)
                ylabel(txt_length)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% GM MUSCLE portion: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_elong_msc_GM),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_elong_msc_GM),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_elong_msc_GM),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_elong_msc_GM),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_msc_GM_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_msc_GM_mean, STR_PRE_elong_msc_GM_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_msc_GM_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_msc_GM_mean, STR_POST_elong_msc_GM_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_msc_GM_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_msc_GM_mean, CON_PRE_elong_msc_GM_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_msc_GM_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_msc_GM_mean, CON_POST_elong_msc_GM_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_el_GMmsc)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) elongation vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_elong_msc_GM),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_elong_msc_GM))
                end
                axis(axis_el_GMmsc)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) elongation vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_elong_msc_GM),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_elong_msc_GM))
                end
                axis(axis_el_GMmsc)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle GM (ins-knee) elongation vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_elong_msc_GM),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_elong_msc_GM),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_elong_msc_GM),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_elong_msc_GM),'b','LineStyle','-','LineWidth',1)
                axis(axis_el_GMmsc)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% GM MUSCLE portion: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_strain_msc_GM),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_strain_msc_GM),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_strain_msc_GM),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_strain_msc_GM),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_strain_msc_GM_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_strain_msc_GM_mean, STR_PRE_strain_msc_GM_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_strain_msc_GM_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_strain_msc_GM_mean, STR_POST_strain_msc_GM_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_strain_msc_GM_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_strain_msc_GM_mean, CON_PRE_strain_msc_GM_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_strain_msc_GM_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_strain_msc_GM_mean, CON_POST_strain_msc_GM_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_str_GMmsc)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) strain vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_strain_msc_GM),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_strain_msc_GM))
                end
                axis(axis_str_GMmsc)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) strain vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_strain_msc_GM),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_strain_msc_GM))
                end
                axis(axis_str_GMmsc)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle GM (ins-knee) strain vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_strain_msc_GM),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_strain_msc_GM),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_strain_msc_GM),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_strain_msc_GM),'b','LineStyle','-','LineWidth',1)
                axis(axis_str_GMmsc)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
                        
            close all
            
            %% GM MUSCLE Lichtwark/Fukunaga: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle GM (archi) length vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_len_msc_GM_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_len_msc_GM_fuku),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_len_msc_GM_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_len_msc_GM_fuku),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_length_msc_GM_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_length_msc_GM_licht_mean, STR_PRE_length_msc_GM_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_length_msc_GM_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_length_msc_GM_licht_mean, STR_POST_length_msc_GM_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_length_msc_GM_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_length_msc_GM_licht_mean, CON_PRE_length_msc_GM_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_length_msc_GM_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_length_msc_GM_licht_mean, CON_POST_length_msc_GM_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                plot(STR_PRE_prone_mean(col_angle),STR_PRE_prone_mean(col_GMmsc_Fukunaga), 'Color', col_lightblue, 'Marker', 'o', 'MarkerFaceColor', col_lightblue)
                plot(STR_POST_prone_mean(col_angle),STR_POST_prone_mean(col_GMmsc_Fukunaga),'ro')
                plot(CON_PRE_prone_mean(col_angle),CON_PRE_prone_mean(col_GMmsc_Fukunaga), 'Color', col_lightred, 'Marker', 'o', 'MarkerFaceColor', col_lightred)
                plot(CON_POST_prone_mean(col_angle),CON_POST_prone_mean(col_GMmsc_Fukunaga),'bo')
                
                axis(axis_len_GMmsc_arch)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southeast')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('muscle GM (archi) length vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_len_msc_GM_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_len_msc_GM_fuku))
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_PRE_count
                    plot(STR_PRE_prone(i,col_angle),STR_PRE_prone(i,col_GMmsc_Fukunaga),'square')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_prone(i,col_angle),STR_POST_prone(i,col_GMmsc_Fukunaga),'o')
                end

                axis(axis_len_GMmsc_arch)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig(length(fig)/4*3), fig(length(fig)/4*2), fig(length(fig)/4*1)],'PRE', 'POST', 'PRE rest', 'POST rest')

                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('muscle GM (archi) length vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_len_msc_GM_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_len_msc_GM_fuku))
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_PRE_count
                    plot(CON_PRE_prone(i,col_angle),CON_PRE_prone(i,col_GMmsc_Fukunaga),'square')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_prone(i,col_angle),CON_POST_prone(i,col_GMmsc_Fukunaga),'o')
                end
                
                axis(axis_len_GMmsc_arch)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig(length(fig)/4*3), fig(length(fig)/4*2), fig(length(fig)/4*1)],'PRE', 'POST', 'PRE rest', 'POST rest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle GM (archi) length vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_len_msc_GM_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_len_msc_GM_fuku),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_len_msc_GM_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_len_msc_GM_fuku),'b','LineStyle','-','LineWidth',1)
                    axis(axis_len_GMmsc_arch)
                    xlabel(txt_gonio)
                    ylabel(txt_length)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end

            %% GM MUSCLE Lichtwark/Fukunaga: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle GM (archi) elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_elong_msc_GM_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_elong_msc_GM_fuku),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_elong_msc_GM_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_elong_msc_GM_fuku),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_msc_GM_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_msc_GM_licht_mean, STR_PRE_elong_msc_GM_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_msc_GM_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_msc_GM_licht_mean, STR_POST_elong_msc_GM_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_msc_GM_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_msc_GM_licht_mean, CON_PRE_elong_msc_GM_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_msc_GM_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_msc_GM_licht_mean, CON_POST_elong_msc_GM_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_el_GMmsc_arch)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('muscle GM (archi) elongation vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_elong_msc_GM_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_elong_msc_GM_fuku))
                end
                axis(axis_el_GMmsc_arch)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('muscle GM (archi) elongation vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_elong_msc_GM_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_elong_msc_GM_fuku))
                end
                axis(axis_el_GMmsc_arch)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle GM (archi) elongation vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_elong_msc_GM_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_elong_msc_GM_fuku),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_elong_msc_GM_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_elong_msc_GM_fuku),'b','LineStyle','-','LineWidth',1)
                axis(axis_el_GMmsc_arch)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end

            %% GM MUSCLE Lichtwark/Fukunaga: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle GM (archi) strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_strain_msc_GM_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_strain_msc_GM_fuku),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_strain_msc_GM_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_strain_msc_GM_fuku),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_strain_msc_GM_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_strain_msc_GM_licht_mean, STR_PRE_strain_msc_GM_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_strain_msc_GM_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_strain_msc_GM_licht_mean, STR_POST_strain_msc_GM_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_strain_msc_GM_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_strain_msc_GM_licht_mean, CON_PRE_strain_msc_GM_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_strain_msc_GM_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_strain_msc_GM_licht_mean, CON_POST_strain_msc_GM_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_str_GMmsc_arch)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('muscle GM (archi) strain vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_strain_msc_GM_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_strain_msc_GM_fuku))
                end
                axis(axis_str_GMmsc_arch)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('muscle GM (archi) strain vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_strain_msc_GM_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_strain_msc_GM_fuku))
                end
                axis(axis_str_GMmsc_arch)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle GM (archi) strain vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_strain_msc_GM_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_strain_msc_GM_fuku),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_strain_msc_GM_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_strain_msc_GM_fuku),'b','LineStyle','-','LineWidth',1)
                axis(axis_str_GMmsc_arch)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% SEE Lichtwark/Fukunaga: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('SEE (archi) length vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_len_SEE_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_len_SEE_fuku),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_len_SEE_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_len_SEE_fuku),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_length_SEE_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_length_SEE_licht_mean, STR_PRE_length_SEE_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_length_tend_GM_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_length_tend_GM_licht_mean, STR_POST_length_tend_GM_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_length_tend_GM_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_length_tend_GM_licht_mean, CON_PRE_length_tend_GM_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_length_tend_GM_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_length_tend_GM_licht_mean, CON_POST_length_tend_GM_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                plot(STR_PRE_prone_mean(col_angle),STR_PRE_prone_mean(col_SEE_Fukunaga), 'Color', col_lightblue, 'Marker', 'o', 'MarkerFaceColor', col_lightblue)
                plot(STR_POST_prone_mean(col_angle),STR_POST_prone_mean(col_SEE_Fukunaga),'ro')
                plot(CON_PRE_prone_mean(col_angle),CON_PRE_prone_mean(col_SEE_Fukunaga),'Color', col_lightred, 'Marker', 'o', 'MarkerFaceColor', col_lightred)
                plot(CON_POST_prone_mean(col_angle),CON_POST_prone_mean(col_SEE_Fukunaga),'bo')

                axis(axis_len_SEE_arch)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southeast')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('SEE (archi) length vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_len_SEE_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_len_SEE_fuku))
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_PRE_count
                    plot(STR_PRE_prone(i,col_angle),STR_PRE_prone(i,col_SEE_Fukunaga),'square')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_prone(i,col_angle),STR_POST_prone(i,col_SEE_Fukunaga),'o')
                end
                
                axis(axis_len_SEE_arch)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig(length(fig)/4*3), fig(length(fig)/4*2), fig(length(fig)/4*1)],'PRE', 'POST', 'PRE rest', 'POST rest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('SEE (archi) length vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_len_SEE_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_len_SEE_fuku))
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_PRE_count
                    plot(CON_PRE_prone(i,col_angle),CON_PRE_prone(i,col_SEE_Fukunaga),'square')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_prone(i,col_angle),CON_POST_prone(i,col_SEE_Fukunaga),'o')
                end
                
                axis(axis_len_SEE_arch)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig(length(fig)/4*3), fig(length(fig)/4*2), fig(length(fig)/4*1)],'PRE', 'POST', 'PRE rest', 'POST rest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND SEE (archi) length vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_len_SEE_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_len_SEE_fuku),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_len_SEE_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_len_SEE_fuku),'b','LineStyle','-','LineWidth',1)
                    axis(axis_len_SEE_arch)
                    xlabel(txt_gonio)
                    ylabel(txt_elong)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end

            %% SEE Lichtwark/Fukunaga: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('SEE (archi) elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_elong_SEE_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_elong_SEE_fuku),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_elong_SEE_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_elong_SEE_fuku),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_SEE_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_SEE_licht_mean, STR_PRE_elong_SEE_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_tend_GM_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_tend_GM_licht_mean, STR_POST_elong_tend_GM_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_tend_GM_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_tend_GM_licht_mean, CON_PRE_elong_tend_GM_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_tend_GM_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_tend_GM_licht_mean, CON_POST_elong_tend_GM_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_el_SEE_arch)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('SEE (archi) elongation vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_elong_SEE_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_elong_SEE_fuku))
                end
                axis(axis_el_SEE_arch)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('SEE (archi) elongation vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_elong_SEE_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_elong_SEE_fuku))
                end
                axis(axis_el_SEE_arch)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND SEE (archi) elongation vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_elong_SEE_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_elong_SEE_fuku),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_elong_SEE_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_elong_SEE_fuku),'b','LineStyle','-','LineWidth',1)
                axis(axis_el_SEE_arch)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end

            %% SEE Lichtwark/Fukunaga: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('SEE (archi) strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_strain_SEE_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_strain_SEE_fuku),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_strain_SEE_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_strain_SEE_fuku),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_strain_SEE_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_strain_SEE_licht_mean, STR_PRE_strain_SEE_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_strain_tend_GM_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_strain_tend_GM_licht_mean, STR_POST_strain_tend_GM_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_strain_tend_GM_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_strain_tend_GM_licht_mean, CON_PRE_strain_tend_GM_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_strain_tend_GM_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_strain_tend_GM_licht_mean, CON_POST_strain_tend_GM_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_str_SEE_arch)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('SEE (archi) strain vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_strain_SEE_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_strain_SEE_fuku))
                end
                axis(axis_str_SEE_arch)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('SEE (archi) strain vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_strain_SEE_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_strain_SEE_fuku))
                end
                axis(axis_str_SEE_arch)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND SEE (archi) strain vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_strain_SEE_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_strain_SEE_fuku),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_strain_SEE_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_strain_SEE_fuku),'b','LineStyle','-','LineWidth',1)
                axis(axis_str_SEE_arch)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            close all
            
            %% GM fascicle (Lichtwark): Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM fascicle length vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_len_GMfas_licht),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_len_GMfas_licht),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_len_GMfas_licht),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_len_GMfas_licht),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_length_GMfas_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_length_GMfas_licht_mean, STR_PRE_length_GMfas_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_length_GMfas_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_length_GMfas_licht_mean, STR_POST_length_GMfas_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_length_GMfas_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_length_GMfas_licht_mean, CON_PRE_length_GMfas_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_length_GMfas_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_length_GMfas_licht_mean, CON_POST_length_GMfas_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                plot(STR_PRE_prone_mean(col_angle),STR_PRE_prone_mean(col_GMfaslen), 'Color', col_lightblue, 'Marker', 'o', 'MarkerFaceColor', col_lightblue)
                plot(STR_POST_prone_mean(col_angle),STR_POST_prone_mean(col_GMfaslen),'ro')
                plot(CON_PRE_prone_mean(col_angle),CON_PRE_prone_mean(col_GMfaslen), 'Color', col_lightred, 'Marker', 'o', 'MarkerFaceColor', col_lightred)
                plot(CON_POST_prone_mean(col_angle),CON_POST_prone_mean(col_GMfaslen),'bo')

                axis(axis_len_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southeast')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('GM fascicle length vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_len_GMfas_licht),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_len_GMfas_licht))
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_PRE_count
                    plot(STR_PRE_prone(i,col_angle),STR_PRE_prone(i,col_GMfaslen),'square')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_prone(i,col_angle),STR_POST_prone(i,col_GMfaslen),'o')
                end
                
                axis(axis_len_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig(length(fig)/4*3), fig(length(fig)/4*2), fig(length(fig)/4*1)],'PRE', 'POST', 'PRE rest', 'POST rest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('GM fascicle length vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_len_GMfas_licht),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_len_GMfas_licht))
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_PRE_count
                    plot(CON_PRE_prone(i,col_angle),CON_PRE_prone(i,col_GMfaslen),'square')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_prone(i,col_angle),CON_POST_prone(i,col_GMfaslen),'o')
                end
                
                axis(axis_len_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig(length(fig)/4*3), fig(length(fig)/4*2), fig(length(fig)/4*1)],'PRE', 'POST', 'PRE rest', 'POST rest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM fascicle length vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_len_GMfas_licht),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_len_GMfas_licht),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_len_GMfas_licht),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_len_GMfas_licht),'b','LineStyle','-','LineWidth',1)
                axis(axis_len_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_length)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Southeast')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% GM fascicle (Lichtwark): Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM fascicle elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_elong_GMfas_licht),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_elong_GMfas_licht),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_elong_GMfas_licht),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_elong_GMfas_licht),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_GMfas_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_GMfas_licht_mean, STR_PRE_elong_GMfas_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_GMfas_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_GMfas_licht_mean, STR_POST_elong_GMfas_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_GMfas_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_GMfas_licht_mean, CON_PRE_elong_GMfas_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_GMfas_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_GMfas_licht_mean, CON_POST_elong_GMfas_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_el_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('GM fascicle elongation vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_elong_GMfas_licht),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_elong_GMfas_licht))
                end
                axis(axis_el_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('GM fascicle elongation vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_elong_GMfas_licht),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_elong_GMfas_licht))
                end
                axis(axis_el_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM fascicle elongation vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_elong_GMfas_licht),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_elong_GMfas_licht),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_elong_GMfas_licht),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_elong_GMfas_licht),'b','LineStyle','-','LineWidth',1)
                axis(axis_el_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Southeast')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% GM fascicle (Lichtwark): Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM fascicle strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_strain_GMfas_licht),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_strain_GMfas_licht),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_strain_GMfas_licht),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_strain_GMfas_licht),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_strain_GMfas_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_strain_GMfas_licht_mean, STR_PRE_strain_GMfas_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_strain_GMfas_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_strain_GMfas_licht_mean, STR_POST_strain_GMfas_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_strain_GMfas_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_strain_GMfas_licht_mean, CON_PRE_strain_GMfas_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_strain_GMfas_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_strain_GMfas_licht_mean, CON_POST_strain_GMfas_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_str_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('GM fascicle strain vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_strain_GMfas_licht),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_strain_GMfas_licht))
                end
                axis(axis_str_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('GM fascicle strain vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_strain_GMfas_licht),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_strain_GMfas_licht))
                end
                axis(axis_str_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM fascicle strain vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_strain_GMfas_licht),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_strain_GMfas_licht),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_strain_GMfas_licht),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_strain_GMfas_licht),'b','LineStyle','-','LineWidth',1)
                axis(axis_str_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Southeast')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% GM fascicle (Lichtwark): Pennation angle vs ankle angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM pennation angle vs ankle angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_pennation_GMfas_licht),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_pennation_GMfas_licht),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_pennation_GMfas_licht),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_pennation_GMfas_licht),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_pennation_GMfas_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_pennation_GMfas_licht_mean, STR_PRE_pennation_GMfas_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_pennation_GMfas_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_pennation_GMfas_licht_mean, STR_POST_pennation_GMfas_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_pennation_GMfas_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_pennation_GMfas_licht_mean, CON_PRE_pennation_GMfas_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_pennation_GMfas_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_pennation_GMfas_licht_mean, CON_POST_pennation_GMfas_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                plot(STR_PRE_prone_mean(col_angle),STR_PRE_prone_mean(col_penn_ang), 'Color', col_lightblue, 'Marker', 'o', 'MarkerFaceColor', col_lightblue)
                plot(STR_POST_prone_mean(col_angle),STR_POST_prone_mean(col_penn_ang),'ro')
                plot(CON_PRE_prone_mean(col_angle),CON_PRE_prone_mean(col_penn_ang), 'Color', col_lightred, 'Marker', 'o', 'MarkerFaceColor', col_lightred)
                plot(CON_POST_prone_mean(col_angle),CON_POST_prone_mean(col_penn_ang),'bo')
                
                axis(axis_penn_GMFAS)
                xlabel(txt_gonio)
                ylabel('Pennation angle (°)')
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('GM pennation angle vs ankle angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_pennation_GMfas_licht),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_pennation_GMfas_licht))
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_PRE_count
                    plot(STR_PRE_prone(i,col_angle),STR_PRE_prone(i,col_penn_ang),'square')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_prone(i,col_angle),STR_POST_prone(i,col_penn_ang),'o')
                end
                
                axis(axis_penn_GMFAS)
                xlabel(txt_gonio)
                ylabel('Pennation angle (°)')
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig(length(fig)/4*3), fig(length(fig)/4*2), fig(length(fig)/4*1)],'PRE', 'POST', 'PRE rest', 'POST rest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('GM pennation angle vs ankle angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_pennation_GMfas_licht),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_pennation_GMfas_licht))
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_PRE_count
                    plot(CON_PRE_prone(i,col_angle),CON_PRE_prone(i,col_penn_ang),'square')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_prone(i,col_angle),CON_POST_prone(i,col_penn_ang),'o')
                end
                
                axis(axis_penn_GMFAS)
                xlabel(txt_gonio)
                ylabel('Pennation angle (°)')
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig(length(fig)/4*3), fig(length(fig)/4*2), fig(length(fig)/4*1)],'PRE', 'POST', 'PRE rest', 'POST rest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM pennation angle vs ankle angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_pennation_GMfas_licht),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_pennation_GMfas_licht),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_pennation_GMfas_licht),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_pennation_GMfas_licht),'b','LineStyle','-','LineWidth',1)
                axis(axis_penn_GMFAS)
                xlabel(txt_gonio)
                ylabel('Pennation angle (°)')
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
                    
            close all
            
            %% Full MTU: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('MTU length vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_len_leg),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_len_leg),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_len_leg),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_len_leg),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_L_MTU_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_L_MTU_mean, STR_PRE_L_MTU_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_L_MTU_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_L_MTU_mean, STR_POST_L_MTU_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_L_MTU_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_L_MTU_mean, CON_PRE_L_MTU_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_L_MTU_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_L_MTU_mean, CON_POST_L_MTU_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                plot(STR_PRE_prone_mean(col_angle),STR_PRE_prone_mean(col_leg), 'Color', col_lightblue, 'Marker', 'o', 'MarkerFaceColor', col_lightblue)
                plot(STR_POST_prone_mean(col_angle),STR_POST_prone_mean(col_leg),'ro')
                plot(CON_PRE_prone_mean(col_angle),CON_PRE_prone_mean(col_leg), 'Color', col_lightred, 'Marker', 'o', 'MarkerFaceColor', col_lightred)
                plot(CON_POST_prone_mean(col_angle),CON_POST_prone_mean(col_leg),'bo')

                axis(axis_len_MTU)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southeast')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('MTU length vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_len_leg),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_len_leg))
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_PRE_count
                    plot(STR_PRE_prone(i,col_angle),STR_PRE_prone(i,col_leg),'square')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_prone(i,col_angle),STR_POST_prone(i,col_leg),'o')
                end
                
                axis(axis_len_MTU)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig(length(fig)/4*3), fig(length(fig)/4*2), fig(length(fig)/4*1)],'PRE', 'POST', 'PRE rest', 'POST rest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('MTU length vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_len_leg),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_len_leg))
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_PRE_count
                    plot(CON_PRE_prone(i,col_angle),CON_PRE_prone(i,col_leg),'square')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_prone(i,col_angle),CON_POST_prone(i,col_leg),'o')
                end
                
                axis(axis_len_MTU)
                xlabel(txt_gonio)
                ylabel(txt_length)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig(length(fig)/4*3), fig(length(fig)/4*2), fig(length(fig)/4*1)],'PRE', 'POST', 'PRE rest', 'POST rest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND MTU length vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_len_leg),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_len_leg),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_len_leg),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_len_leg),'b','LineStyle','-','LineWidth',1)
                axis(axis_len_MTU)
                xlabel(txt_gonio)
                ylabel(txt_length)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Southeast')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% Full MTU: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('MTU elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_elong_leg),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_elong_leg),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_elong_leg),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_elong_leg),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_MTU_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_MTU_mean, STR_PRE_elong_MTU_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_MTU_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_MTU_mean, STR_POST_elong_MTU_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_MTU_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_MTU_mean, CON_PRE_elong_MTU_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_MTU_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_MTU_mean, CON_POST_elong_MTU_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_el_MTU)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('MTU elongation vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_elong_leg),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_elong_leg))
                end
                axis(axis_el_MTU)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('MTU elongation vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_elong_leg),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_elong_leg))
                end
                axis(axis_el_MTU)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND MTU elongation vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_elong_leg),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_elong_leg),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_elong_leg),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_elong_leg),'b','LineStyle','-','LineWidth',1)
                axis(axis_el_MTU)
                xlabel(txt_gonio)
                ylabel(txt_elong)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% Full MTU: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % MTU strain plots simply show the Grieve calculation factor
            % removed from output through toggle "plot_old"
            if plot_check && plot_old
                plottitle = horzcat('MTU strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_strain_leg),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_strain_leg),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_strain_leg),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_strain_leg),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_strain_MTU_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_strain_MTU_mean, STR_PRE_strain_MTU_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_strain_MTU_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_strain_MTU_mean, STR_POST_strain_MTU_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_strain_MTU_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_strain_MTU_mean, CON_PRE_strain_MTU_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_strain_MTU_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_strain_MTU_mean, CON_POST_strain_MTU_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_str_MTU)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check && plot_old
                plottitle = horzcat('MTU strain vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_strain_leg),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_strain_leg))
                end
                axis(axis_str_MTU)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check && plot_old
                plottitle = horzcat('MTU strain vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_strain_leg),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_strain_leg))
                end
                axis(axis_str_MTU)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND MTU strain vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_strain_leg),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_strain_leg),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_strain_leg),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_strain_leg),'b','LineStyle','-','LineWidth',1)
                axis(axis_str_MTU)
                xlabel(txt_gonio)
                ylabel(txt_strain)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% GMFAS scans: Displacement vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('displacement GMFAS vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_elong_GMFAS),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_elong_GMFAS),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_elong_GMFAS),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_elong_GMFAS),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_displ_GMFAS_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_displ_GMFAS_mean, STR_PRE_displ_GMFAS_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_displ_GMFAS_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_displ_GMFAS_mean, STR_POST_displ_GMFAS_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_displ_GMFAS_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_displ_GMFAS_mean, CON_PRE_displ_GMFAS_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_displ_GMFAS_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_displ_GMFAS_mean, CON_POST_displ_GMFAS_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_displ_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_displ)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('displacement GMFAS vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_elong_GMFAS),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_elong_GMFAS))
                end
                axis(axis_displ_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_displ)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('displacement GMFAS vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_elong_GMFAS),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_elong_GMFAS))
                end
                axis(axis_displ_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_displ)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual && plot_old
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND displacement GMFAS vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_elong_GMFAS),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_elong_GMFAS),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_elong_GMFAS),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_elong_GMFAS),'b','LineStyle','-','LineWidth',1)
                axis(axis_displ_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_displ)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            
            
            close all

            
            %% EMG vs angle GM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('EMG gas med vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_EMG_gm),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_EMG_gm),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_EMG_gm),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_EMG_gm),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_EMG_gm_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_EMG_gm_mean, STR_PRE_EMG_gm_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_EMG_gm_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_EMG_gm_mean, STR_POST_EMG_gm_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_EMG_gm_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_EMG_gm_mean, CON_PRE_EMG_gm_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_EMG_gm_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_EMG_gm_mean, CON_POST_EMG_gm_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_EMG_narrow)
                xlabel(txt_gonio)
                ylabel(txt_emg)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
                
            end
            
            if plot_check
                plottitle = horzcat('EMG gas med vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_EMG_gm),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_EMG_gm))
                end
                axis(axis_EMG_wide)
                xlabel(txt_gonio)
                ylabel(txt_emg)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('EMG gas med vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_EMG_gm),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_EMG_gm))
                end
                axis(axis_EMG_wide)
                xlabel(txt_gonio)
                ylabel(txt_emg)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND EMG gas med vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_EMG_gm),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_EMG_gm),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_EMG_gm),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_EMG_gm),'b','LineStyle','-','LineWidth',1)
                axis(axis_EMG_wide)
                xlabel(txt_gonio)
                ylabel(txt_emg)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% EMG vs angle GL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('EMG gas lat vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_EMG_gl),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_EMG_gl),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_EMG_gl),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_EMG_gl),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_EMG_gl_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_EMG_gl_mean, STR_PRE_EMG_gl_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_EMG_gl_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_EMG_gl_mean, STR_POST_EMG_gl_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_EMG_gl_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_EMG_gl_mean, CON_PRE_EMG_gl_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_EMG_gl_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_EMG_gl_mean, CON_POST_EMG_gl_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_EMG_narrow)
                xlabel(txt_gonio)
                ylabel(txt_emg)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('EMG gas lat vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_EMG_gl),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_EMG_gl))
                end
                axis(axis_EMG_narrow)
                xlabel(txt_gonio)
                ylabel(txt_emg)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('EMG gas lat vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_EMG_gl),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_EMG_gl))
                end
                axis(axis_EMG_narrow)
                xlabel(txt_gonio)
                ylabel(txt_emg)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND EMG gas lat vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_EMG_gl),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_EMG_gl),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_EMG_gl),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_EMG_gl),'b','LineStyle','-','LineWidth',1)
                    axis(axis_EMG_narrow)
                    xlabel(txt_gonio)
                    ylabel(txt_emg)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% EMG vs angle SOL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('EMG soleus vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_EMG_sol),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_EMG_sol),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_EMG_sol),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_EMG_sol),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_EMG_sol_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_EMG_sol_mean, STR_PRE_EMG_sol_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_EMG_sol_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_EMG_sol_mean, STR_POST_EMG_sol_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_EMG_sol_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_EMG_sol_mean, CON_PRE_EMG_sol_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_EMG_sol_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_EMG_sol_mean, CON_POST_EMG_sol_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_EMG_narrow)
                xlabel(txt_gonio)
                ylabel(txt_emg)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('EMG soleus vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_EMG_sol),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_EMG_sol))
                end
                axis(axis_EMG_wide)
                xlabel(txt_gonio)
                ylabel(txt_emg)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('EMG soleus vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_EMG_sol),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_EMG_sol))
                end
                axis(axis_EMG_wide)
                xlabel(txt_gonio)
                ylabel(txt_emg)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND EMG soleus vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_EMG_sol),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_EMG_sol),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_EMG_sol),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_EMG_sol),'b','LineStyle','-','LineWidth',1)
                axis(axis_EMG_wide)
                xlabel(txt_gonio)
                ylabel(txt_emg)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            
            
            %% EMG vs angle, average 3 muscles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('EMG 3msc vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_EMG_mean),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_EMG_mean),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_EMG_mean),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_EMG_mean),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_EMG_mean_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_EMG_mean_mean, STR_PRE_EMG_mean_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_EMG_mean_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_EMG_mean_mean, STR_POST_EMG_mean_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_EMG_mean_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_EMG_mean_mean, CON_PRE_EMG_mean_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_EMG_mean_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_EMG_mean_mean, CON_POST_EMG_mean_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_EMG_narrow)
                xlabel(txt_gonio)
                ylabel(txt_emg)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('EMG 3msc vs angle - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_EMG_mean),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_EMG_mean))
                end
                axis(axis_EMG_narrow)
                xlabel(txt_gonio)
                ylabel(txt_emg)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('EMG 3msc vs angle - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_EMG_mean),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_EMG_mean))
                end
                axis(axis_EMG_narrow)
                xlabel(txt_gonio)
                ylabel(txt_emg)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND EMG 3msc vs angle - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_EMG_mean),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_EMG_mean),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_EMG_mean),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_EMG_mean),'b','LineStyle','-','LineWidth',1)
                axis(axis_EMG_narrow)
                xlabel(txt_gonio)
                ylabel(txt_emg)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            
            close all
            
                        
            %% NORMALIZED GM MUSCLE Lichtwark/Fukunaga: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle GM (archi) length vs angle NORM - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_norm_len_msc_GM_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_norm_len_msc_GM_fuku),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_norm_len_msc_GM_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_norm_len_msc_GM_fuku),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_length_norm_msc_GM_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_length_norm_msc_GM_licht_mean, STR_PRE_length_norm_msc_GM_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_length_norm_msc_GM_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_length_norm_msc_GM_licht_mean, STR_POST_length_norm_msc_GM_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_length_norm_msc_GM_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_length_norm_msc_GM_licht_mean, CON_PRE_length_norm_msc_GM_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_length_norm_msc_GM_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_length_norm_msc_GM_licht_mean, CON_POST_length_norm_msc_GM_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')

                %axis(axis_len_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_length_norm)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southeast')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('muscle GM (archi) length vs angle NORM - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_norm_len_msc_GM_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_norm_len_msc_GM_fuku))
                end
                
                %axis(axis_len_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_length_norm)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('muscle GM (archi) length vs angle NORM - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_norm_len_msc_GM_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_norm_len_msc_GM_fuku))
                end
                
                %axis(axis_len_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_length_norm)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle GM (archi) length vs angle NORM - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_norm_len_msc_GM_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_norm_len_msc_GM_fuku),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_norm_len_msc_GM_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_norm_len_msc_GM_fuku),'b','LineStyle','-','LineWidth',1)
                    %axis(axis_len_GMFAS)
                    xlabel(txt_gonio)
                    ylabel(txt_length_norm)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
                        
            %% NORMALIZED GM MUSCLE Lichtwark/Fukunaga: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle GM (archi) elongation vs angle NORM - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_norm_elong_msc_GM_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_norm_elong_msc_GM_fuku),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_norm_elong_msc_GM_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_norm_elong_msc_GM_fuku),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_norm_msc_GM_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_norm_msc_GM_licht_mean, STR_PRE_elong_norm_msc_GM_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_norm_msc_GM_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_norm_msc_GM_licht_mean, STR_POST_elong_norm_msc_GM_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_norm_msc_GM_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_norm_msc_GM_licht_mean, CON_PRE_elong_norm_msc_GM_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_norm_msc_GM_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_norm_msc_GM_licht_mean, CON_POST_elong_norm_msc_GM_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                %axis(axis_el_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('muscle GM (archi) elongation vs angle NORM - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_norm_elong_msc_GM_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_norm_elong_msc_GM_fuku))
                end
                %axis(axis_el_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('muscle GM (archi) elongation vs angle NORM - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_norm_elong_msc_GM_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_norm_elong_msc_GM_fuku))
                end
                %axis(axis_el_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle GM (archi) elongation vs angle NORM - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_norm_elong_msc_GM_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_norm_elong_msc_GM_fuku),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_norm_elong_msc_GM_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_norm_elong_msc_GM_fuku),'b','LineStyle','-','LineWidth',1)
                    %axis(axis_el_GMFAS)
                    xlabel(txt_gonio)
                    ylabel(txt_elong_norm)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
                                    
            %% NORMALIZED SEE Lichtwark/Fukunaga: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('SEE (archi) length vs angle NORM - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_norm_len_SEE_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_norm_len_SEE_fuku),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_norm_len_SEE_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_norm_len_SEE_fuku),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_length_norm_SEE_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_length_norm_SEE_licht_mean, STR_PRE_length_norm_SEE_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_length_norm_SEE_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_length_norm_SEE_licht_mean, STR_POST_length_norm_SEE_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_length_norm_SEE_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_length_norm_SEE_licht_mean, CON_PRE_length_norm_SEE_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_length_norm_SEE_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_length_norm_SEE_licht_mean, CON_POST_length_norm_SEE_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')

                %axis(axis_len_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_length_norm)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southeast')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('SEE (archi) length vs angle NORM - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_norm_len_SEE_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_norm_len_SEE_fuku))
                end
                
                %axis(axis_len_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_length_norm)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('SEE (archi) length vs angle NORM - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_norm_len_SEE_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_norm_len_SEE_fuku))
                end
                
                %axis(axis_len_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_length_norm)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND SEE (archi) length vs angle NORM - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_norm_len_SEE_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_norm_len_SEE_fuku),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_norm_len_SEE_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_norm_len_SEE_fuku),'b','LineStyle','-','LineWidth',1)
                    %axis(axis_len_GMFAS)
                    xlabel(txt_gonio)
                    ylabel(txt_length_norm)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% NORMALIZED SEE Lichtwark/Fukunaga: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('SEE (archi) elongation vs angle NORM - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_norm_elong_SEE_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_norm_elong_SEE_fuku),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_norm_elong_SEE_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_norm_elong_SEE_fuku),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_norm_SEE_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_norm_SEE_licht_mean, STR_PRE_elong_norm_SEE_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_norm_SEE_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_norm_SEE_licht_mean, STR_POST_elong_norm_SEE_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_norm_SEE_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_norm_SEE_licht_mean, CON_PRE_elong_norm_SEE_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_norm_SEE_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_norm_SEE_licht_mean, CON_POST_elong_norm_SEE_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                %axis(axis_el_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('SEE (archi) elongation vs angle NORM - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_norm_elong_SEE_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_norm_elong_SEE_fuku))
                end
                %axis(axis_el_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('SEE (archi) elongation vs angle NORM - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_norm_elong_SEE_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_norm_elong_SEE_fuku))
                end
                %axis(axis_el_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND SEE (archi) elongation vs angle NORM - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_norm_elong_SEE_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_norm_elong_SEE_fuku),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_norm_elong_SEE_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_norm_elong_SEE_fuku),'b','LineStyle','-','LineWidth',1)
                    %axis(axis_el_GMFAS)
                    xlabel(txt_gonio)
                    ylabel(txt_elong_norm)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            close all
            
            
            %% NORMALIZED TO ELONG, GM MUSCLE Lichtwark/Fukunaga: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle GM (archi) elongation vs angle NORM% - 1 % contrib');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_norm_percent_elong_msc_GM_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_norm_percent_elong_msc_GM_fuku),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_norm_percent_elong_msc_GM_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_norm_percent_elong_msc_GM_fuku),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_norm_percent_msc_GM_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_norm_percent_msc_GM_licht_mean, STR_PRE_elong_norm_percent_msc_GM_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_norm_percent_msc_GM_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_norm_percent_msc_GM_licht_mean, STR_POST_elong_norm_percent_msc_GM_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_norm_percent_msc_GM_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_norm_percent_msc_GM_licht_mean, CON_PRE_elong_norm_percent_msc_GM_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_norm_percent_msc_GM_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_norm_percent_msc_GM_licht_mean, CON_POST_elong_norm_percent_msc_GM_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_el_percent)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm_perc)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('muscle GM (archi) elongation vs angle NORM% - 2 STRETCHERS ind PRE+POST % contrib');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_norm_percent_elong_msc_GM_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_norm_percent_elong_msc_GM_fuku))
                end
                axis(axis_el_percent)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm_perc)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('muscle GM (archi) elongation vs angle NORM% - 3 CONTROLS ind PRE+POST % contrib');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_norm_percent_elong_msc_GM_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_norm_percent_elong_msc_GM_fuku))
                end
                axis(axis_el_percent)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm_perc)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle GM (archi) elongation vs angle NORM% - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST % contrib');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_norm_percent_elong_msc_GM_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_norm_percent_elong_msc_GM_fuku),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_norm_percent_elong_msc_GM_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_norm_percent_elong_msc_GM_fuku),'b','LineStyle','-','LineWidth',1)
                    axis(axis_el_percent)
                    xlabel(txt_gonio)
                    ylabel(txt_elong_norm_perc)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% NORMALIZED TO ELONG, SEE Lichtwark/Fukunaga: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('SEE (archi) elongation vs angle NORM% - 1 % contrib');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_norm_percent_elong_SEE_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_norm_percent_elong_SEE_fuku),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_norm_percent_elong_SEE_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_norm_percent_elong_SEE_fuku),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_norm_percent_SEE_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_norm_percent_SEE_licht_mean, STR_PRE_elong_norm_percent_SEE_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_norm_percent_SEE_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_norm_percent_SEE_licht_mean, STR_POST_elong_norm_percent_SEE_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_norm_percent_SEE_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_norm_percent_SEE_licht_mean, CON_PRE_elong_norm_percent_SEE_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_norm_percent_SEE_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_norm_percent_SEE_licht_mean, CON_POST_elong_norm_percent_SEE_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                axis(axis_el_percent)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm_perc)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('SEE (archi) elongation vs angle NORM% - 2 STRETCHERS ind PRE+POST % contrib');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_norm_percent_elong_SEE_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_norm_percent_elong_SEE_fuku))
                end
                axis(axis_el_percent)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm_perc)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('SEE (archi) elongation vs angle NORM% - 3 CONTROLS ind PRE+POST % contrib');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_norm_percent_elong_SEE_fuku),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_norm_percent_elong_SEE_fuku))
                end
                axis(axis_el_percent)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm_perc)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND SEE (archi) elongation vs angle NORM% - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST % contrib');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_norm_percent_elong_SEE_fuku),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_norm_percent_elong_SEE_fuku),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_norm_percent_elong_SEE_fuku),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_norm_percent_elong_SEE_fuku),'b','LineStyle','-','LineWidth',1)
                    axis(axis_el_percent)
                    xlabel(txt_gonio)
                    ylabel(txt_elong_norm_perc)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            %% NORMALIZED GM fascicle (Lichtwark): Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM fascicle length vs angle NORM - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_norm_len_GMfas_licht),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_norm_len_GMfas_licht),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_norm_len_GMfas_licht),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_norm_len_GMfas_licht),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_length_norm_GMfas_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_length_norm_GMfas_licht_mean, STR_PRE_length_norm_GMfas_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_length_norm_GMfas_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_length_norm_GMfas_licht_mean, STR_POST_length_norm_GMfas_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_length_norm_GMfas_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_length_norm_GMfas_licht_mean, CON_PRE_length_norm_GMfas_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_length_norm_GMfas_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_length_norm_GMfas_licht_mean, CON_POST_length_norm_GMfas_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                %axis(axis_len_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_length_norm)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southeast')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('GM fascicle length vs angle NORM - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_norm_len_GMfas_licht),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_norm_len_GMfas_licht))
                end
                
                %axis(axis_len_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_length_norm)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('GM fascicle length vs angle NORM - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_norm_len_GMfas_licht),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_norm_len_GMfas_licht))
                end
                
                %axis(axis_len_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_length_norm)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM fascicle length vs angle NORM - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_norm_len_GMfas_licht),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_norm_len_GMfas_licht),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_norm_len_GMfas_licht),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_norm_len_GMfas_licht),'b','LineStyle','-','LineWidth',1)
                    %axis(axis_len_GMFAS)
                    xlabel(txt_gonio)
                    ylabel(txt_length_norm)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Southeast')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
                        
            %% NORMALIZED GM fascicle (Lichtwark): Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM fascicle elongation vs angle NORM - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_norm_elong_GMfas_licht),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_norm_elong_GMfas_licht),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_norm_elong_GMfas_licht),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_norm_elong_GMfas_licht),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_norm_GMfas_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_norm_GMfas_licht_mean, STR_PRE_elong_norm_GMfas_licht_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_norm_GMfas_licht_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_norm_GMfas_licht_mean, STR_POST_elong_norm_GMfas_licht_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_norm_GMfas_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_norm_GMfas_licht_mean, CON_PRE_elong_norm_GMfas_licht_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_norm_GMfas_licht_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_norm_GMfas_licht_mean, CON_POST_elong_norm_GMfas_licht_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                %axis(axis_el_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('GM fascicle elongation vs angle NORM - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_norm_elong_GMfas_licht),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_norm_elong_GMfas_licht))
                end
                %axis(axis_el_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('GM fascicle elongation vs angle NORM - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_norm_elong_GMfas_licht),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_norm_elong_GMfas_licht))
                end
                %axis(axis_el_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM fascicle elongation vs angle NORM - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_norm_elong_GMfas_licht),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_norm_elong_GMfas_licht),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_norm_elong_GMfas_licht),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_norm_elong_GMfas_licht),'b','LineStyle','-','LineWidth',1)
                    %axis(axis_el_GMFAS)
                    xlabel(txt_gonio)
                    ylabel(txt_elong_norm)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Southeast')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
                                               
            %% NORMALIZED leg/MTU length: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('MTU length vs angle NORM - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_norm_len_leg),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_norm_len_leg),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_norm_len_leg),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_norm_len_leg),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_length_norm_MTU_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_length_norm_MTU_mean, STR_PRE_length_norm_MTU_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_length_norm_MTU_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_length_norm_MTU_mean, STR_POST_length_norm_MTU_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_length_norm_MTU_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_length_norm_MTU_mean, CON_PRE_length_norm_MTU_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_length_norm_MTU_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_length_norm_MTU_mean, CON_POST_length_norm_MTU_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')

                %axis(axis_len_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_length_norm)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southeast')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('MTU length vs angle NORM - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_norm_len_leg),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_norm_len_leg))
                end
                
                %axis(axis_len_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_length_norm)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('MTU length vs angle NORM - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_norm_len_leg),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_norm_len_leg))
                end
                
                %axis(axis_len_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_length_norm)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND MTU length vs angle NORM - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_norm_len_leg),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_norm_len_leg),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_norm_len_leg),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_norm_len_leg),'b','LineStyle','-','LineWidth',1)
                    %axis(axis_len_GMFAS)
                    xlabel(txt_gonio)
                    ylabel(txt_length_norm)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
                        
            %% NORMALIZED leg/MTU length: Elongation vs angle NORM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('MTU elongation vs angle NORM - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,col_AV_angle), STR_PRE_angle_vars_mean(:,col_AV_norm_elong_leg),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,col_AV_angle), STR_POST_angle_vars_mean(:,col_AV_norm_elong_leg),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,col_AV_angle), CON_PRE_angle_vars_mean(:,col_AV_norm_elong_leg),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,col_AV_angle), CON_POST_angle_vars_mean(:,col_AV_norm_elong_leg),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_norm_MTU_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_norm_MTU_mean, STR_PRE_elong_norm_MTU_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_norm_MTU_mean, STR_POST_ROM_SD, 'r.')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_norm_MTU_mean, STR_POST_elong_norm_MTU_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_norm_MTU_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_norm_MTU_mean, CON_PRE_elong_norm_MTU_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_norm_MTU_mean, CON_POST_ROM_SD, 'b.')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_norm_MTU_mean, CON_POST_elong_norm_MTU_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
                
                %axis(axis_el_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm)
                title(plottitle,'Interpreter', 'none')
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            if plot_check
                plottitle = horzcat('MTU elongation vs angle NORM - 2 STRETCHERS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle),STR_PRE_angle_vars{1,i}(:,col_AV_norm_elong_leg),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle),STR_POST_angle_vars{1,i}(:,col_AV_norm_elong_leg))
                end
                %axis(axis_el_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            if plot_check
                plottitle = horzcat('MTU elongation vs angle NORM - 3 CONTROLS ind PRE+POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle),CON_PRE_angle_vars{1,i}(:,col_AV_norm_elong_leg),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle),CON_POST_angle_vars{1,i}(:,col_AV_norm_elong_leg))
                end
                %axis(axis_el_GMFAS)
                xlabel(txt_gonio)
                ylabel(txt_elong_norm)
                title(plottitle,'Interpreter', 'none')
                fig=get(gca,'Children');
                legend([fig(end), fig((length(fig)/2))], 'PRE', 'POST')
                print(horzcat('data_plots/GRP_INT ',plottitle),'-dpng')
            end
            
            % rough coding
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND MTU elongation vs angle NORM - ', CON_PRE_ID{i}(1:6) ,' 2legs PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,col_AV_angle), STR_PRE_angle_vars{1,i}(:,col_AV_norm_elong_leg),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,col_AV_angle), STR_POST_angle_vars{1,i}(:,col_AV_norm_elong_leg),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,col_AV_angle), CON_PRE_angle_vars{1,i}(:,col_AV_norm_elong_leg),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,col_AV_angle), CON_POST_angle_vars{1,i}(:,col_AV_norm_elong_leg),'b','LineStyle','-','LineWidth',1)
                    %axis(axis_el_GMFAS)
                    xlabel(txt_gonio)
                    ylabel(txt_elong_norm)
                    title(plottitle,'Interpreter', 'none')
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    print(horzcat('data_plots/',plottitle),'-dpng')
                end
            end
            
            close all
            
        end
        
end