%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main file for analysis of passive dorsiflexion with US
% Marie Moltubakk 4.2.2015
% 
% Note 1:
% The scripts assume a certain structure for input files from Tracker and
% Noraxon. I.e. number of and order of EMG and other channels. If these
% scripts are to be used for other projects, some of the code must be
% modified. These lines are marked % PROJECTSPECIFIC
% 
% 16.09.15: adapted to read datamaster including "Lichtwark" file names
% 03.08.16: INPUT ARGUMENT: pass 0/1/2 for amount of plots (see below)
% 13.03.17: adapted for pre-post data comparison
%
% input argument 1 = project selection (1 = BD, 2 = intervent)
% input argument 2 = plot selection (0 = none, 1 = group plots, 2 = ind plots)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






function [] = passiveUS(input_project, input_plot)
    close all
    
    
    
    %% PLOTS - determine which plots to display
    global plot_achilles plot_norm plot_emg plot_check plot_us subject_id plot_licht plot_individual plot_conversion

    if input_plot >= 1 
        plot_check = 1; % LEVEL 1: group plots
    else
        plot_check = 0;
    end
    if input_plot >= 2
        plot_individual = 1; % LEVEL 1B: individual plots
    else
        plot_individual = 0;
    end
    
    toggle_normalization = 1; % 0 = GM muscle/tendon/fascicle length/elong in absolute values, 1 = % of leg length

    % LEVEL 2: main checkpoint plots
    plot_norm = 0;
    
    % LEVEL 3:
    plot_conversion = 0; % turn on/off plots for data conversion Norm
    plot_us = 0; % tracked feature vs external marker 
    plot_emg = 0; % RMS 3 EMG channels per trial
    plot_achilles = 0; % turn on/off Achilles machine plots
    plot_licht = 0; % plot averaging of trials from Lichtwark US fascicle tracking
    %% PLOTS - determine which plots to display



    %% Set constants and globals % PROJECTSPECIFIC

    % sampling frequencies
    global us_zerodispframes noraxonfreq freq_default
    us_zerodispframes = 1; % No of US frames to average as zero displacement
    noraxonfreq = 1500; % sampling frequency of noraxon data
    freq_default = 100; % output frequency for noraxon data without US video (MVC, etc)
    angle_step = 0.05; % reshaped, averaged data extracted every x degrees

    % variables for NORM conversion factors calculated from actual data
    global convert_norm_angle_a convert_norm_angle_b convert_norm_torque_a convert_norm_torque_b convert_norm_velocity_a convert_norm_velocity_b convert_norm_direction_b
    % default conversion factors from Norm manuals and Achilles calibration
    global convert_achilles norm_volt_per_degree norm_volt_per_velocity norm_volt_per_nm_a norm_volt_per_nm_b
    convert_achilles = -81.9909;  % conversion factor ACHILLES torque, V -> Nm
    norm_volt_per_degree = (2048*((1*16)+0)/1024)*(10000/32768); %   (Gain * (Sampled Position + Offset) / 1024) * (10v/32768) = Volt-value, Sampled Position is in units of 1/16 degree.
    norm_volt_per_velocity = (1024*((1*16)+0)/1024)*(10000/32768); %   (Gain * (Sampled Velocity + Offset) / 1024) * (10v/32768) = Volt-value. Sampled Velocity is in units of 1/16 degree.
    % norm_volt_per_nm = (1024*((1*(1.355818*32768/500))+0)/1024)*(10000/32768); % 27.116 µV/Nm   Sampled Torque is in units of Foot-Pounds * 32768 / 500
    norm_volt_per_nm_a = 0.07490071; % from file "M M M conversion volt-torque 2014-DES NEW DATA.xlsx"
    norm_volt_per_nm_b = 0.69283161;
    
    % cutoff frequencies for filters
    global emg_bandpass emg_rms_ms mvc_window_ms 
    emg_bandpass = [10/(noraxonfreq/2) 500/(noraxonfreq/2)]; % cutoff frequencies for EMG butterworth bandpass filter
    emg_rms_ms = 500; % milliseconds RMS window, must be divisible by 4 - ref Basmajian 1985 = 50 ms, ref Aagaard paper Passive tensile stress = 200 ms
    mvc_window_ms = 500; % milliseconds window for determining MVC torque and EMG
    global angle_cutoff velocity_cutoff torque_cutoff_bandstop torque_cutoff_active angle_cutoff_active velocity_cutoff_active
    angle_cutoff = 10/(noraxonfreq/2); % 9/(noraxonfreq/2); % cutoff freq, Norm angle PASSIVE - ref Winter 1990 = 15 hz. Kongsgaard = 8hz
    angle_cutoff_active = 20/(noraxonfreq/2); %  20/(noraxonfreq/2);
    velocity_cutoff = 20/(noraxonfreq/2); %  12/(noraxonfreq/2);
    velocity_cutoff_active = 20/(noraxonfreq/2); %  12/(noraxonfreq/2);
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
    %% Set constants and globals 
    
           
    
    %% set AXES etc for plots
    col_lightblue = [0.6 0.8 1];
    col_lightred = [1 0.6 0.8];

    if input_project == 1 
        %% dancer study
        axis_ind_elong = [-1 40 -6 20];
        axis_ind_strain = [-1 40 -6 20];

        axis_EMG = [-1 40 0 30];
        
        axis_faslen = [-1 40 40 85];
        axis_pennation = [-1 40 8 18];
        axis_faslen_elong = [-1 40 -inf inf];
        axis_faslen_str = [-1 40 -inf inf];
        
        axis_len_MTU = [-1 35 -0 550];
        axis_el_MTU = [-1 35 -1 35];
        %axis_str_MTU = [-1 35 0 8];
        %axis_str_MTUP = [-5 100 0 8];

        axis_displ_GMFAS = [-1 35 -2.5 2.5];
        axis_displ_SOL = [-1 35 -3 11];

        axis_el_SOL = [-1 35 -1 10];
        axis_str_SOL = [-1 35 -0.5 5];
        axis_str_SOLP = [-5 100 -0.5 5];
        axis_len_SOL = [-1 35 0 350];

        axis_el_GM = [-1 35 -1 22];
        axis_str_GM = [-1 35 -1 10];
        axis_str_GMP = [-5 100 -1 10];
        axis_len_GM = [-1 35 0 350];

        axis_el_GMapo = [-1 35 -7 3];
        axis_str_GMapo = [-1 35 -5 2];
        axis_str_GMapoP = [-5 100 -5 2];
        axis_len_GMapo = [-1 35 0 180];

        axis_el_GMtend = [-1 35 -2 20];
        axis_str_GMtend = [-1 35 -1 10];
        axis_str_GMtendP = [-5 100 -1 10];
        axis_len_GMtend = [-1 35 0 300];

        axis_el_AT = [-1 35 -1 25];
        axis_str_AT = [-1 35 -1 30];
        axis_str_ATP = [-5 100 -1 30];
        axis_len_AT = [-1 35 0 140];

        axis_force = [-1 35 0 1500];
        axis_torque = [-1 35 0 80];
        axis_PP = [-5 100 0 105];
    else
        %% intervention study

        axis_ind_elong = [-1 37 -6 30];
        axis_ind_strain = [-1 37 -6 30];

        axis_EMG = [-1 37 0 30];

        axis_faslen = [-1 37 -inf inf];
        axis_pennation = [-1 37 -inf inf];
        
        axis_len_MTU = [-1 37 -0 550];
        axis_el_MTU = [-1 37 -1 35];
        %axis_str_MTU = [-1 37 0 8];
        %axis_str_MTUP = [-5 100 0 8];

        axis_displ_GMFAS = [-1 37 -2 4.5];
        axis_displ_SOL = [-1 37 -3 11];

        axis_el_SOL = [-1 37 -1 11];
        axis_str_SOL = [-1 37 -1 6];
        axis_str_SOLP = [-5 100 -1 4];
        axis_len_SOL = [-1 37 0 350];

        axis_el_GM = [-1 37 -1 15];
        axis_str_GM = [-1 37 -0.5 6];
        axis_str_GMP = [-5 100 -0.5 6];
        axis_len_GM = [-1 37 0 350];

        axis_el_GMapo = [-1 37 -5 2];
        axis_str_GMapo = [-1 37 -5 2];
        axis_str_GMapoP = [-5 100 -4 2];
        axis_len_GMapo = [-1 37 0 150];

        axis_el_GMtend = [-1 37 -1 28];
        axis_str_GMtend = [-1 37 -1 14];
        axis_str_GMtendP = [-5 100 -1 14];
        axis_len_GMtend = [-1 37 0 300];

        axis_el_AT = [-1 37 -1 32];
        axis_str_AT = [-1 37 -1 45];
        axis_str_ATP = [-5 100 -1 70];
        axis_len_AT = [-1 37 0 140];

        axis_force = [-1 37 0 1900];
        axis_torque = [-1 37 0 100];
        axis_PP = [-5 100 0 105];
    end
    %% set AXES for plots
    
    
    %% Read max angles and forces from "create_angles_passive.m"
    %%% To extract max angles and forces per trial/subject/common
    %%% Produces arrays with angles and forces, to be retrieved later

%    global ang_subjectno 
%    global ang_pre_r ang_pre_l ang_post_r ang_post_l ang_ind_max ang_common_max % ang = norm ankle angle, probably not to be used
    global input_gon_pre_r input_gon_pre_l input_gon_post_r input_gon_post_l input_gon_ind_max input_gon_common_max % gon = goniometer ankle angle
    global input_for_pre_r input_for_pre_l input_for_post_r input_for_post_l input_for_ind_max input_for_common_max input_for_ind_rmax input_for_ind_lmax
    dm_filename = 'angles_output.tsv';
    if exist(dm_filename, 'file') == 2
        read_angles_passive(dm_filename);
    else
        cprintf('*red', 'ERROR: Subject angle data file does not exist. Generate by running "create_angles_passive".\n')
        return
    end
    dm_filename = 'forces_output.tsv';
    if exist(dm_filename, 'file') == 2
        read_forces_passive(dm_filename);
    else
        cprintf('*red', 'ERROR: Subject force data file does not exist. Generate by running "create_angles_passive".\n')
        return
    end
    %% Read ... from "create_angles_passive.m"



    %% Read datamaster file, to connect corresponding data files
    %%% Produces arrays with file names and variables per trial, to be retrieved later

    global dm_subjectno dm_timepoint dm_side dm_trial 
    global dm_ROM_sol1_NX dm_ROM_sol1_US dm_ROM_sol1_US_frame dm_ROM_sol2_NX dm_ROM_sol2_US dm_ROM_sol2_US_frame
    global dm_ROM_gmmtj1_NX dm_ROM_gmmtj1_US dm_ROM_gmmtj1_US_frame dm_ROM_gmmtj2_NX dm_ROM_gmmtj2_US dm_ROM_gmmtj2_US_frame
    global dm_ROM_gmfas1_NX dm_ROM_gmfas1_US dm_ROM_gmfas1_US_frame dm_ROM_gmfas2_NX dm_ROM_gmfas2_US dm_ROM_gmfas2_US_frame  dm_ROM_gmfas1_licht dm_ROM_gmfas2_licht
    global dm_MVC_PF dm_MVC_DF %dm_CPM_calc_NX dm_CPM_calc_US dm_CPM_calc_US_frame dm_CPM_sol_NX dm_CPM_sol_US dm_CPM_sol_US_frame 
    global dm_leg_length dm_at_SOL_length dm_at_GM_length
    global at_momentarm
    global filepath
    dm_filename = 'data/datamaster_passive.tsv';
    dm_columns = 35; % number of data columns entered per subject % PROJECTSPECIFIC
    linestotal = read_datamaster_passive(dm_filename,dm_columns);
    %% Read datamaster file, to connect corresponding data files
        
    
    
    %% preallocate output arrays
    % common arrays for all subjects:
        all_passive_output_head = {'Subject', 'Time', 'Side', 'Trial', ...
            'ROM trial (°)', 'ROM subject (L/R/PRE/POST)', 'ROM common (all subjects)', '10°', '67% ROM', ...
            'force @ trial ROM (N)', 'force @ subject max ROM', 'force @ subject max force', 'force @ common max ROM', 'force @ 0 deg', 'force @ 10°', 'force @ 67% ROM', ...
            'torque @ trial ROM (N)', 'torque @ subject max ROM', 'torque @ subject max torque', 'torque @ common max ROM', 'torque @ 0 deg', 'torque @ 10°', 'torque @ 67% ROM', ...
            'angle @ trial max force (°)', 'angle @ subject max force', 'angle @ common max force', 'angle @ ind R max force', 'angle @ ind L max force', ...            
            'Stiffness @ trial ROM (Nm/°)', 'Stiffness @ subject max ROM', 'Stiffness @ common max ROM', 'Stiffness @ 15 deg', 'Stiffness @ 10°', 'Stiffness @ 67% ROM',...
            'Stiffness index (Nm/°^2)', ...
            'displ SOL @ trial ROM (mm)', 'displ SOL @ subject max ROM', 'displ SOL @ common max ROM', 'displ SOL @ 10°', 'displ SOL @ 67% ROM', ...
            'displ GMMTJ @ trial ROM', 'displ GMMTJ @ subject max ROM', 'displ GMMTJ @ common max ROM', 'displ GMMTJ @ 10°', 'displ GMMTJ @ 67% ROM', ...
            'displ GMfas @ trial ROM', 'displ GMfas @ subject max ROM', 'displ GMfas @ common max ROM', 'displ GMFAS @ 10°', 'displ GMFAS @ 67% ROM', ...
            'length AT trial max', 'length GMtend trial max', 'length leg trial max', 'length GMapo trial max', 'length msc GM trial max', 'length msc SOL trial max', ...
            'length AT ind max', 'length GMtend ind max', 'length leg ind max', 'length GMapo ind max', 'length msc GM ind max', 'length msc SOL ind max', ...
            'length AT common max', 'length GMtend common max', 'length leg common max', 'length GMapo common max', 'length msc GM common max', 'length msc SOL common max', ...
            'length AT submax 1', 'length GMtend submax 1', 'length leg submax 1', 'length GMapo submax 1', 'length msc GM submax 1', 'length msc SOL submax 1', ...
            'length AT submax 2', 'length GMtend submax 2', 'length leg submax 2', 'length GMapo submax 2', 'length msc GM submax 2', 'length msc SOL submax 2', ...
            'elong AT trial max', 'elong GMtend trial max', 'elong leg trial max', 'elong GMapo trial max', 'elong msc GM trial max', 'elong msc SOL trial max',...
            'elong AT ind max', 'elong GMtend ind max', 'elong leg ind max', 'elong GMapo ind max', 'elong msc GM ind max', 'elong msc SOL ind max', ...
            'elong AT common max', 'elong GMtend common max', 'elong leg common max', 'elong GMapo common max', 'elong msc GM common max', 'elong msc SOL common max', ...
            'elong AT submax 1', 'elong GMtend submax 1', 'elong leg submax 1', 'elong GMapo submax 1', 'elong msc GM submax 1', 'elong msc SOL submax 1', ...
            'elong AT submax 2', 'elong GMtend submax 2', 'elong leg submax 2', 'elong GMapo submax 2', 'elong msc GM submax 2', 'elong msc SOL submax 2', ...
            'strain AT trial max', 'strain GMtend trial max', 'strain leg trial max', 'strain GMapo trial max', 'strain msc GM trial max', 'strain msc SOL trial max', ...
            'strain AT ind max', 'strain GMtend ind max', 'strain leg ind max', 'strain GMapo ind max', 'strain msc GM ind max', 'strain msc SOL ind max', ...
            'strain AT common max', 'strain GMtend common max', 'strain leg common max', 'strain GMapo common max', 'strain msc GM common max', 'strain msc SOL common max', ...
            'strain AT submax 1', 'strain GMtend submax 1', 'strain leg submax 1', 'strain GMapo submax 1', 'strain msc GM submax 1', 'strain msc SOL submax 1', ...
            'strain AT submax 2', 'strain GMtend submax 2', 'strain leg submax 2', 'strain GMapo submax 2', 'strain msc GM submax 2', 'strain msc SOL submax 2', ...
            'length GMtend Fuku trial max', 'length GMtend Fuku ind max', 'length GMtend Fuku common max', 'length GMtend Fuku submax 1', 'length GMtend Fuku submax 2' ...
            'elong GMtend Fuku trial max', 'elong GMtend Fuku ind max', 'elong GMtend Fuku common max', 'elong GMtend Fuku submax 1', 'elong GMtend Fuku submax 2' ...
            'strain GMtend Fuku trial max', 'strain GMtend Fuku ind max', 'strain GMtend Fuku common max', 'strain GMtend Fuku submax 1', 'strain GMtend Fuku submax 2' ...
            'length msc GM Fuku trial max', 'length msc GM Fuku ind max', 'length msc GM Fuku common max', 'length msc GM Fuku submax 1', 'length msc GM Fuku submax 2' ...
            'elong msc GM Fuku trial max', 'elong msc GM Fuku ind max', 'elong msc GM Fuku common max', 'elong msc GM Fuku submax 1', 'elong msc GM Fuku submax 2' ...
            'strain msc GM Fuku trial max', 'strain msc GM Fuku ind max', 'strain msc GM Fuku common max', 'strain msc GM Fuku submax 1', 'strain msc GM Fuku submax 2' ...
            'faslen GM @ trial ROM (mm)', 'faslen GM @ subject max ROM', 'faslen GM @ common max ROM', 'faslen GM @ 0 deg', 'faslen GM @ 10°', 'faslen GM @ 67% ROM', ...
            'pennation GM @ trial ROM (°)', 'pennation GM @ subject max ROM', 'pennation GM @ common max ROM', 'pennation GM @ 0 deg', 'pennation GM @ 10°', 'pennation GM @ 67% ROM', ...
            'faslen SOL @ trial ROM', 'faslen SOL @ subject max ROM', 'faslen SOL @ common max ROM', 'faslen SOL @ 0 deg', 'faslen SOL @ 10°', 'faslen SOL @ 67% ROM', ...
            'pennation SOL @ trial ROM', 'pennation SOL @ subject max ROM', 'pennation SOL @ common max ROM', 'pennation SOL @ 0 deg', 'pennation SOL @ 10°', 'pennation SOL @ 67% ROM', ...
            'GM elong contribution@ trial ROM (%)', 'GM elong contribution@ subject max ROM', 'GM elong contribution@ common max ROM', 'GM elong contribution@ 10°', 'GM elong contribution@ 67% ROM', ...
            'EMG GM @ trial ROM (%)', 'EMG GM @ subject max ROM', 'EMG GM @ common max ROM', 'EMG GM @ 10°', 'EMG GM @ 67% ROM', ...
            'EMG GL @ trial ROM', 'EMG GL @ subject max ROM', 'EMG GL @ common max ROM', 'EMG GL @ 10°', 'EMG GL @ 67% ROM', ...
            'EMG SOL @ trial ROM', 'EMG SOL @ subject max ROM', 'EMG SOL @ common max ROM', 'EMG SOL @ 10°', 'EMG SOL @ 67% ROM', ...
            }; % PROJECTSPECIFIC
    all_passive_output = zeros(ceil(linestotal),length(all_passive_output_head)-4); 
    all_passive_output_txt = cell(ceil(linestotal),4);

    if input_project == 1 % BD study
        BD_count = 0; % # of ballet dancer subjects
        CON_count = 0; % # of controls = intervention study subjects

        BD_angle_vars{ceil(linestotal)} = zeros;
        BD_angle_vars_mean{ceil(linestotal)} = zeros;
        BD_angle_vars_norm{ceil(linestotal)} = zeros;
        BD_angle_vars_norm_indlength{ceil(linestotal)} = zeros;
        BD_angle_vars_norm_mean{ceil(linestotal)} = zeros;

        CON_angle_vars{ceil(linestotal)} = zeros;
        CON_angle_vars_mean{ceil(linestotal)} = zeros;
        CON_angle_vars_norm{ceil(linestotal)} = zeros;
        CON_angle_vars_norm_indlength{ceil(linestotal)} = zeros;
        CON_angle_vars_norm_mean{ceil(linestotal)} = zeros;
    else % intervention study
        CON_PRE_count = 0;
        CON_POST_count = 0;
        STR_PRE_count = 0;
        STR_POST_count = 0;
        CON_PRE_subject_ID(ceil(linestotal/4)) = zeros; % holds subject numbers for intervention study - rough coding

        CON_PRE_angle_vars{ceil(linestotal)} = zeros;
        CON_PRE_angle_vars_mean{ceil(linestotal)} = zeros;
        CON_PRE_angle_vars_norm{ceil(linestotal)} = zeros;
        CON_PRE_angle_vars_norm_indlength{ceil(linestotal)} = zeros;
        CON_PRE_angle_vars_norm_mean{ceil(linestotal)} = zeros;

        CON_POST_angle_vars{ceil(linestotal)} = zeros;
        CON_POST_angle_vars_mean{ceil(linestotal)} = zeros;
        CON_POST_angle_vars_norm{ceil(linestotal)} = zeros;
        CON_POST_angle_vars_norm_indlength{ceil(linestotal)} = zeros;
        CON_POST_angle_vars_norm_mean{ceil(linestotal)} = zeros;

        STR_PRE_angle_vars{ceil(linestotal)} = zeros;
        STR_PRE_angle_vars_mean{ceil(linestotal)} = zeros;
        STR_PRE_angle_vars_norm{ceil(linestotal)} = zeros;
        STR_PRE_angle_vars_norm_indlength{ceil(linestotal)} = zeros;
        STR_PRE_angle_vars_norm_mean{ceil(linestotal)} = zeros;

        STR_POST_angle_vars{ceil(linestotal)} = zeros;
        STR_POST_angle_vars_mean{ceil(linestotal)} = zeros;
        STR_POST_angle_vars_norm{ceil(linestotal)} = zeros;
        STR_POST_angle_vars_norm_indlength{ceil(linestotal)} = zeros;
        STR_POST_angle_vars_norm_mean{ceil(linestotal)} = zeros;
    end
    % preallocate output arrays
    
    
    
    
    
    
    
    
    
    
    for line = 1:linestotal %% LOOP through all lines in datamaster file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear emg_all;
        emg_all{3,6} = zeros; % 6 EMG channels for 1 subject, to be reused across subjects
        
        
        
        %% subject/trial identifier
        trial_subjectno = str2double(dm_subjectno{line});
        trial_timepoint = strcmp(dm_timepoint{line},'POST'); % 0 = PRE, 1 = POST
        trial_leg = strcmp(dm_trial{line},'STR'); % 0 = CON, 1 = STR
        
        if input_project == 1 % BD study
            if trial_subjectno > 100
                filepath = 'data\BD\';
                subject_id = horzcat('Dancer ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
                BD_count = BD_count + 1;
            else
                filepath = 'data\';
                subject_id = horzcat('Control ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
                CON_count = CON_count + 1;
            end
        elseif input_project == 2 % intervention
            filepath = 'data\';
            subject_id = horzcat('Intervent ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
            if trial_timepoint == 0 && trial_leg == 0 % PRE, CON
                CON_PRE_count = CON_PRE_count + 1;
                CON_PRE_subject_ID(CON_PRE_count) = trial_subjectno;
            elseif trial_timepoint == 0 && trial_leg == 1 % PRE, STR
                STR_PRE_count = STR_PRE_count + 1;
            elseif trial_timepoint == 1 && trial_leg == 0 % POST, CON
                CON_POST_count = CON_POST_count + 1;
            elseif trial_timepoint == 1 && trial_leg == 1 % POST, STR
                STR_POST_count = STR_POST_count + 1;
            end
        end
        cprintf('*black', horzcat('----------------', subject_id, '------------------\n'))
        %%
        
        
        
        %% calculate EMG data 

        % prepare column placement
        if strcmpi(dm_side{line},'R') == 1
            column_tibant = column_r_tibant;
            column_gm = column_r_gm;
            column_gl = column_r_gl;
            column_sol = column_r_sol;
        else % left
            column_tibant = column_l_tibant;
            column_gm = column_l_gm;
            column_gl = column_l_gl;
            column_sol = column_l_sol;
        end

        % Read noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample - DORSIFLEXION
        % Produce a new noraxon data array
        noraxon_mvc_dorsi = read_noraxon_stiffness(strcat(filepath, dm_MVC_DF{line}), freq_default, dm_side{line}, 'MVC dorsi');

        % Calculate co-activation constants
        % Read complete, prepared noraxon array + number of frames to average (freq * time)
        % Produce max torque, max EMG constants
        [~,EMG_max_TA] = calculate_EMG_max(noraxon_mvc_dorsi, freq_default*(mvc_window_ms/1000), column_tibant, 1); % 1 = invert torque for dorsiflexion

        % Read noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample - PLANTAR FLEXION
        % Produce a new noraxon data array
        noraxon_mvc_plantar = read_noraxon_stiffness(strcat(filepath, dm_MVC_PF{line}), freq_default, dm_side{line}, 'MVC plantar');

        % Calculate co-activation constants
        % Read complete, prepared noraxon array + number of frames to average (freq * time)
        % Produce max torque, max EMG constants
        [~,EMG_max_gm] = calculate_EMG_max(noraxon_mvc_plantar, freq_default*(mvc_window_ms/1000), column_gm, 0);
        [~,EMG_max_gl] = calculate_EMG_max(noraxon_mvc_plantar, freq_default*(mvc_window_ms/1000), column_gl, 0);
        [~,EMG_max_sol] = calculate_EMG_max(noraxon_mvc_plantar, freq_default*(mvc_window_ms/1000), column_sol, 0);
        %%


        
        %% calculate NORM conversion factors
        % retrieve conversion constants for Norm data
        [convert_norm_angle_a, convert_norm_angle_b, convert_norm_torque_a, convert_norm_torque_b, convert_norm_velocity_a, convert_norm_velocity_b, convert_norm_direction_b] = calculate_norm_constants_complete(dm_subjectno{line}, dm_side{line}, dm_timepoint{line}, line, 'passive');
        %%



        %% calculate ACHILLES TENDON MOMENT ARM
        at_momentarm = calculate_momentarm(0, 0, dm_leg_length{line});
        %%
        
        

        %% calculations for 2x SOL trials
        
        % extract force, gonio, angle, displacement for EACH TRIAL
        % NB: extract_force_displ_singletrial_passive_EMG is where torque is converted to force
        if(strcmpi(dm_ROM_sol1_NX{line}, 'null'))
            % allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case in final_stiffness
            SOL_force_1 = zeros;
            SOL_gonio_1 = zeros;
            SOL_angle_1 = zeros;
            SOL_displacement_1 = zeros;
        else 
            [SOL_force_1, SOL_gonio_1, SOL_angle_1, SOL_displacement_1, SOL_emg_gm_1, SOL_emg_gl_1, SOL_emg_sol_1, SOL_time_1] = extract_force_displ_singletrial_passive_EMG(dm_ROM_sol1_NX{line}, dm_ROM_sol1_US{line}, dm_ROM_sol1_US_frame{line}, EMG_max_TA, EMG_max_gm, EMG_max_gl, EMG_max_sol, dm_leg_length{line}, dm_side{line}, line, 'SOL1');
            emg_all{1,1} = [SOL_gonio_1 SOL_emg_gm_1];
            emg_all{2,1} = [SOL_gonio_1 SOL_emg_gl_1];
            emg_all{3,1} = [SOL_gonio_1 SOL_emg_sol_1];
        end
        if(strcmpi(dm_ROM_sol2_NX{line}, 'null'))
            SOL_force_2 = zeros;
            SOL_gonio_2 = zeros;
            SOL_angle_2 = zeros;
            SOL_displacement_2 = zeros;
        else 
            [SOL_force_2, SOL_gonio_2, SOL_angle_2, SOL_displacement_2, SOL_emg_gm_2, SOL_emg_gl_2, SOL_emg_sol_2, SOL_time_2] = extract_force_displ_singletrial_passive_EMG(dm_ROM_sol2_NX{line}, dm_ROM_sol2_US{line}, dm_ROM_sol2_US_frame{line}, EMG_max_TA, EMG_max_gm, EMG_max_gl, EMG_max_sol, dm_leg_length{line}, dm_side{line}, line, 'SOL2');    
            emg_all{1,2} = [SOL_gonio_2 SOL_emg_gm_2];
            emg_all{2,2} = [SOL_gonio_2 SOL_emg_gl_2];
            emg_all{3,2} = [SOL_gonio_2 SOL_emg_sol_2];
        end

        % average two passive trials
        if(strcmpi(dm_ROM_sol1_NX{line}, 'null')) % trial 1 not existing
            data_SOL = average_passive_trials_EMG(SOL_force_2, SOL_gonio_2, SOL_angle_2, SOL_displacement_2, SOL_emg_gm_2, SOL_emg_gl_2, SOL_emg_sol_2, SOL_time_2);
        elseif(strcmpi(dm_ROM_sol2_NX{line}, 'null'))
            data_SOL = average_passive_trials_EMG(SOL_force_1, SOL_gonio_1, SOL_angle_1, SOL_displacement_1, SOL_emg_gm_1, SOL_emg_gl_1, SOL_emg_sol_1, SOL_time_1);
        else % if 2 trials exist (none of variables are 'null')
            data_SOL = average_passive_trials_EMG(SOL_force_1, SOL_gonio_1, SOL_angle_1, SOL_displacement_1, SOL_emg_gm_1, SOL_emg_gl_1, SOL_emg_sol_1, SOL_time_1, SOL_force_2, SOL_gonio_2, SOL_angle_2, SOL_displacement_2, SOL_emg_gm_2, SOL_emg_gl_2, SOL_emg_sol_2, SOL_time_2);
        end

        % plot summary of 2 trials: Displacement-angle
        if plot_check && plot_individual
            plottitle = horzcat('IND displacement SOL vs angle, ', subject_id);
            figure('Name',plottitle);
            plot(SOL_gonio_1,SOL_displacement_1,'LineWidth',2)
            hold on
            plot(SOL_gonio_2,SOL_displacement_2,'LineWidth',2)
            axis(axis_displ_SOL)
            ylabel('Displacement (mm)')
            xlabel('Gonio angle (°)')
            title(plottitle)
            legend('Trial 1','Trial 2','Location','Southeast')
            saveas(gcf, horzcat('data_plots/', plottitle,'.jpg'))
        end
        %%

        %% calculations for 2x GM MTJ trials

        % extract force, gonio, angle, displacement for EACH TRIAL
        if(strcmpi(dm_ROM_gmmtj1_NX{line}, 'null'))
            % allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case in final_stiffness
            GMMTJ_force_1 = zeros;
            GMMTJ_angle_1 = zeros;
            GMMTJ_gonio_1 = zeros;
            GMMTJ_displacement_1 = zeros;
        else 
            [GMMTJ_force_1, GMMTJ_gonio_1, GMMTJ_angle_1, GMMTJ_displacement_1, GMMTJ_emg_gm_1, GMMTJ_emg_gl_1, GMMTJ_emg_sol_1, GMMTJ_time_1] = extract_force_displ_singletrial_passive_EMG(dm_ROM_gmmtj1_NX{line}, dm_ROM_gmmtj1_US{line}, dm_ROM_gmmtj1_US_frame{line}, EMG_max_TA, EMG_max_gm, EMG_max_gl, EMG_max_sol, dm_leg_length{line}, dm_side{line}, line, 'GMMTJ1');
            emg_all{1,3} = [GMMTJ_gonio_1 GMMTJ_emg_gm_1];
            emg_all{2,3} = [GMMTJ_gonio_1 GMMTJ_emg_gl_1];
            emg_all{3,3} = [GMMTJ_gonio_1 GMMTJ_emg_sol_1];
        end
        if(strcmpi(dm_ROM_gmmtj2_NX{line}, 'null'))
            GMMTJ_force_2 = zeros;
            GMMTJ_angle_2 = zeros;
            GMMTJ_gonio_2 = zeros;
            GMMTJ_displacement_2 = zeros;
        else 
            [GMMTJ_force_2, GMMTJ_gonio_2, GMMTJ_angle_2, GMMTJ_displacement_2, GMMTJ_emg_gm_2, GMMTJ_emg_gl_2, GMMTJ_emg_sol_2, GMMTJ_time_2] = extract_force_displ_singletrial_passive_EMG(dm_ROM_gmmtj2_NX{line}, dm_ROM_gmmtj2_US{line}, dm_ROM_gmmtj2_US_frame{line}, EMG_max_TA, EMG_max_gm, EMG_max_gl, EMG_max_sol, dm_leg_length{line}, dm_side{line}, line, 'GMMTJ2');    
            emg_all{1,4} = [GMMTJ_gonio_2 GMMTJ_emg_gm_2];
            emg_all{2,4} = [GMMTJ_gonio_2 GMMTJ_emg_gl_2];
            emg_all{3,4} = [GMMTJ_gonio_2 GMMTJ_emg_sol_2];
        end

        % average two passive trials
        if(strcmpi(dm_ROM_gmmtj1_NX{line}, 'null'))
            data_GMMTJ = average_passive_trials_EMG(GMMTJ_force_2, GMMTJ_gonio_2, GMMTJ_angle_2, GMMTJ_displacement_2, GMMTJ_emg_gm_2, GMMTJ_emg_gl_2, GMMTJ_emg_sol_2, GMMTJ_time_2);
        elseif(strcmpi(dm_ROM_gmmtj2_NX{line}, 'null'))
            data_GMMTJ = average_passive_trials_EMG(GMMTJ_force_1, GMMTJ_gonio_1, GMMTJ_angle_1, GMMTJ_displacement_1, GMMTJ_emg_gm_1, GMMTJ_emg_gl_1, GMMTJ_emg_sol_1, GMMTJ_time_1);
        else % if 2 trials exist (none of variables are 'null')
            data_GMMTJ = average_passive_trials_EMG(GMMTJ_force_1, GMMTJ_gonio_1, GMMTJ_angle_1, GMMTJ_displacement_1, GMMTJ_emg_gm_1, GMMTJ_emg_gl_1, GMMTJ_emg_sol_1, GMMTJ_time_1, GMMTJ_force_2, GMMTJ_gonio_2, GMMTJ_angle_2, GMMTJ_displacement_2, GMMTJ_emg_gm_2, GMMTJ_emg_gl_2, GMMTJ_emg_sol_2, GMMTJ_time_2);
        end
    
        % plot summary of 2 trials: Displacement-angle
        if plot_check && plot_individual
            plottitle = horzcat('IND displacement GMMTJ vs angle, ', subject_id);
            figure('Name',plottitle);
            plot(GMMTJ_gonio_1,GMMTJ_displacement_1,'LineWidth',2)
            hold on
            plot(GMMTJ_gonio_2,GMMTJ_displacement_2,'LineWidth',2)
            axis(axis_el_GM)
            ylabel('Displacement (mm)')
            xlabel('Gonio angle (°)')
            title(plottitle)
            legend('Trial 1','Trial 2','Location','Southeast')
            saveas(gcf, horzcat('data_plots/', plottitle,'.jpg'))
        end
        %%
        
        %% calculations for 2x GM fascicle trials - ORIGINAL calculations, as for GMMTJ and SOL

        % extract force, gonio, angle, displacement for EACH TRIAL
        if(strcmpi(dm_ROM_gmfas1_NX{line}, 'null'))
            % allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case in final_stiffness
            GMFAS_force_1 = zeros;
            GMFAS_angle_1 = zeros;
            GMFAS_gonio_1 = zeros;
            GMFAS_displacement_1 = zeros;
        else 
            [GMFAS_force_1, GMFAS_gonio_1, GMFAS_angle_1, GMFAS_displacement_1, GMFAS_emg_gm_1, GMFAS_emg_gl_1, GMFAS_emg_sol_1, GMFAS_time_1] = extract_force_displ_singletrial_passive_EMG(dm_ROM_gmfas1_NX{line}, dm_ROM_gmfas1_US{line}, dm_ROM_gmfas1_US_frame{line}, EMG_max_TA, EMG_max_gm, EMG_max_gl, EMG_max_sol, dm_leg_length{line}, dm_side{line}, line, 'GMFAS1');
            emg_all{1,5} = [GMFAS_gonio_1 GMFAS_emg_gm_1];
            emg_all{2,5} = [GMFAS_gonio_1 GMFAS_emg_gl_1];
            emg_all{3,5} = [GMFAS_gonio_1 GMFAS_emg_sol_1];
        end
        if(strcmpi(dm_ROM_gmfas2_NX{line}, 'null'))
            GMFAS_force_2 = zeros;
            GMFAS_angle_2 = zeros;
            GMFAS_gonio_2 = zeros;
            GMFAS_displacement_2 = zeros;
        else 
            [GMFAS_force_2, GMFAS_gonio_2, GMFAS_angle_2, GMFAS_displacement_2, GMFAS_emg_gm_2, GMFAS_emg_gl_2, GMFAS_emg_sol_2, GMFAS_time_2] = extract_force_displ_singletrial_passive_EMG(dm_ROM_gmfas2_NX{line}, dm_ROM_gmfas2_US{line}, dm_ROM_gmfas2_US_frame{line}, EMG_max_TA, EMG_max_gm, EMG_max_gl, EMG_max_sol, dm_leg_length{line}, dm_side{line}, line, 'GMFAS2');    
            emg_all{1,6} = [GMFAS_gonio_2 GMFAS_emg_gm_2];
            emg_all{2,6} = [GMFAS_gonio_2 GMFAS_emg_gl_2];
            emg_all{3,6} = [GMFAS_gonio_2 GMFAS_emg_sol_2];
        end

        % average two passive trials
        if(strcmpi(dm_ROM_gmfas1_NX{line}, 'null'))
            data_GMFAS = average_passive_trials_EMG(GMFAS_force_2, GMFAS_gonio_2, GMFAS_angle_2, GMFAS_displacement_2, GMFAS_emg_gm_2, GMFAS_emg_gl_2, GMFAS_emg_sol_2, GMFAS_time_2);
        elseif(strcmpi(dm_ROM_gmfas2_NX{line}, 'null'))
            data_GMFAS = average_passive_trials_EMG(GMFAS_force_1, GMFAS_gonio_1, GMFAS_angle_1, GMFAS_displacement_1, GMFAS_emg_gm_1, GMFAS_emg_gl_1, GMFAS_emg_sol_1, GMFAS_time_1);
        else % if 2 trials exist (none of variables are 'null')
            data_GMFAS = average_passive_trials_EMG(GMFAS_force_1, GMFAS_gonio_1, GMFAS_angle_1, GMFAS_displacement_1, GMFAS_emg_gm_1, GMFAS_emg_gl_1, GMFAS_emg_sol_1, GMFAS_time_1, GMFAS_force_2, GMFAS_gonio_2, GMFAS_angle_2, GMFAS_displacement_2, GMFAS_emg_gm_2, GMFAS_emg_gl_2, GMFAS_emg_sol_2, GMFAS_time_2);
        end

        % plot summary of 2 trials: Displacement-angle
        if plot_check && plot_individual
            plottitle = horzcat('IND displacement GMFAS vs angle, ', subject_id);
            figure('Name',plottitle);
            plot(GMFAS_gonio_1,GMFAS_displacement_1,'LineWidth',2)
            hold on
            plot(GMFAS_gonio_2,GMFAS_displacement_2,'LineWidth',2)
            axis(axis_displ_GMFAS)
            ylabel('Displacement (mm)')
            xlabel('Gonio angle (°)')
            title(plottitle)
            legend('Trial 1','Trial 2','Location','Southeast')
            saveas(gcf, horzcat('data_plots/', plottitle,'.jpg'))
        end
        %%
                
        %% calculations for 2x GM fascicle trials - LICHTWARK analysis
        
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
        %   averaged fasicle length
        %   averaged pennation angle
        % OR containing zeros (if nonexistent)
        if(strcmpi(dm_ROM_gmfas1_licht{line}, 'null')==0 && strcmpi(dm_ROM_gmfas2_licht{line}, 'null')==0)
            
            % if 2 GM trials exist (none of variables are 'null')

            % GM Lichtwark: perform averaging of trial 1 and trial 2
            data_GMFAS_licht_GM = average_passive_trials_licht(GMFAS_gonio_1, GMFAS_angle_1, GMFAS_licht_GM_faslen_1, GMFAS_licht_GM_pennation_1, GMFAS_time_1, GMFAS_gonio_2, GMFAS_angle_2, GMFAS_licht_GM_faslen_2, GMFAS_licht_GM_pennation_2, GMFAS_time_2, 'GM fascicles');

            % SOL Lichtwark: check for existence of SOL data
            if GMFAS_licht_SOL_1_exists && GMFAS_licht_SOL_2_exists
                % average two trials:
                data_GMFAS_licht_SOL = average_passive_trials_licht(GMFAS_gonio_1, GMFAS_angle_1, GMFAS_licht_SOL_faslen_1, GMFAS_licht_SOL_pennation_1, GMFAS_time_1, GMFAS_gonio_2, GMFAS_angle_2, GMFAS_licht_SOL_faslen_2, GMFAS_licht_SOL_pennation_2, GMFAS_time_2, 'SOL fascicles');
            elseif GMFAS_licht_SOL_1_exists
                data_GMFAS_licht_SOL = average_passive_trials_licht(GMFAS_gonio_1, GMFAS_angle_1, GMFAS_licht_SOL_faslen_1, GMFAS_licht_SOL_pennation_1, GMFAS_time_1, 'SOL fascicles');
            elseif GMFAS_licht_SOL_2_exists
                data_GMFAS_licht_SOL = average_passive_trials_licht(GMFAS_gonio_2, GMFAS_angle_2, GMFAS_licht_SOL_faslen_2, GMFAS_licht_SOL_pennation_2, GMFAS_time_2, 'SOL fascicles');
            else % no SOL exists
                data_GMFAS_licht_SOL = zeros(1,3);
            end
            
        elseif(strcmpi(dm_ROM_gmfas1_licht{line}, 'null')==0) % only trial 1 exists
            % keep GM trial 1
            data_GMFAS_licht_GM = average_passive_trials_licht(GMFAS_gonio_1, GMFAS_angle_1, GMFAS_licht_GM_faslen_1, GMFAS_licht_GM_pennation_1, GMFAS_time_1, 'GM fascicles');
            % keep eventual SOL trial 1
            if GMFAS_licht_SOL_1_exists
                data_GMFAS_licht_SOL = average_passive_trials_licht(GMFAS_gonio_1, GMFAS_angle_1, GMFAS_licht_SOL_faslen_1, GMFAS_licht_SOL_pennation_1, GMFAS_time_1, 'SOL fascicles');
            else % no SOL exists
                data_GMFAS_licht_SOL = zeros(1,3);
            end    
            
        elseif(strcmpi(dm_ROM_gmfas2_licht{line}, 'null')==0) % only trial 2 exists
            % keep GM trial 2
            data_GMFAS_licht_GM = average_passive_trials_licht(GMFAS_gonio_2, GMFAS_angle_2, GMFAS_licht_GM_faslen_2, GMFAS_licht_GM_pennation_2, GMFAS_time_2, 'GM fascicles');
            % keep eventual SOL trial 2
            if GMFAS_licht_SOL_2_exists
                data_GMFAS_licht_SOL = average_passive_trials_licht(GMFAS_gonio_2, GMFAS_angle_2, GMFAS_licht_SOL_faslen_2, GMFAS_licht_SOL_pennation_2, GMFAS_time_2, 'SOL fascicles');
            else % no SOL exists
                data_GMFAS_licht_SOL = zeros(1,3);
            end    
            
        else % no trials exist
            data_GMFAS_licht_GM = zeros(1,3);
            data_GMFAS_licht_SOL = zeros(1,3);
        end
        %%
        
        
        
        %% average SOL + GMMTJ + GMFAS trials for force, gonio, angle, EMG
        
        % average 3 force arrays into one
        data_force_gonio = average_passive_forces_EMG(data_SOL, data_GMMTJ, data_GMFAS);
        
        % plot all 6 trials separately: Force-angle
        if plot_check && plot_individual
            plottitle = horzcat('IND force vs angle, ', subject_id);
            figure('Name',plottitle);
            plot(SOL_gonio_1,SOL_force_1,'LineWidth',2, 'Color',[1 0 0])
            hold on
            plot(SOL_gonio_2,SOL_force_2,'LineWidth',2, 'Color',[1 0.6 0])
            plot(GMMTJ_gonio_1,GMMTJ_force_1,'LineWidth',2, 'Color',[1 1 0])
            plot(GMMTJ_gonio_2,GMMTJ_force_2,'LineWidth',2, 'Color',[0 1 0])
            plot(GMFAS_gonio_1,GMFAS_force_1,'LineWidth',2, 'Color',[0 0 1])
            plot(GMFAS_gonio_2,GMFAS_force_2,'LineWidth',2, 'Color',[1 0 1])
            axis(axis_force)
            ylabel('Force (N)')
            xlabel('Gonio angle (°)')
            title(plottitle)
            legend('SOL1','SOL2','GMMTJ1','GMMTJ2','GMFAS1','GMFAS2','Location','Southeast')
            saveas(gcf, horzcat('data_plots/', plottitle,'.jpg'))
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
            plot(data_force_gonio(:,2),data_force_gonio(:,3),'y','LineWidth',2);
            
            for i = 1:6 % 6 trials GL
                if ~isempty(emg_all{1,i})
                    plot(emg_all{2,i}(:,1),emg_all{2,i}(:,2),'m');
                end
            end
            plot(data_force_gonio(:,2),data_force_gonio(:,4),'m','LineWidth',2);
            
            for i = 1:6 % 6 trials SOL
                if ~isempty(emg_all{1,i})
                    plot(emg_all{3,i}(:,1),emg_all{3,i}(:,2),'c');
                end
            end
            plot(data_force_gonio(:,2),data_force_gonio(:,5),'c','LineWidth',2);
            
            axis(axis_EMG)
            ylabel('Muscle activation (% of MVC)')
            xlabel('Gonio angle (°)')
            title(plottitle)
            legend('Gastr.med.','Gastr.lat.','Soleus','Location','Northwest')
            saveas(gcf, horzcat('data_plots/', plottitle,'.jpg'))
        end
        %%
        
        
        
        %% check conformation of GONIOMETER to norm angle
        if plot_check && plot_individual
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
                plot(SOL_angle_2, SOL_angle_2, 'LineWidth',2, 'Color',[0.3 0.3 0.3], 'LineStyle',':') % conformation line
            else
                plot(SOL_angle_1, SOL_angle_1, 'LineWidth',2, 'Color',[0.3 0.3 0.3], 'LineStyle',':') % conformation line
            end
            plot([0 0],[0 0], 'Marker','o', 'MarkerSize',10, 'MarkerFaceColor',[0.3 0.3 0.3])% zero point
            xlabel('Norm angle (°)')
            ylabel('Gonio angle (°)')
            title(plottitle)
            legend('SOL1','SOL2','GMMTJ1','GMMTJ2','GMFAS1','GMFAS2','Norm/Norm','Location','Southeast')
            saveas(gcf, horzcat('data_plots/', plottitle,'.jpg'))
        end
        %%
        
        
        
        %% GM MUSCLE LENGTH CHANGE from penn.ang. & fas.len. (Lichtwark/Fukunaga)
        % Fukunaga 2001: In vivo behavor of human ...

        %   data_GMFAS_licht_GM
        % containing:
        %   1 averaged angle (currently calculated from gonio)
        %   2 averaged fasicle length - in mm
        %   3 averaged pennation angle - in radians?
        % OR containing 3x zeros (if nonexistent)
        
        if length(data_GMFAS_licht_GM) > 3
            % fascicle longitudinal component = faslen * cos penn_ang, at each joint angle
            horiz_l = data_GMFAS_licht_GM(:,2) .* cosd(data_GMFAS_licht_GM(:,3));

            % elongation: longitudinal component at angle zero = zero elongation
            loc_frame = find(data_GMFAS_licht_GM(:,1)>=0,1,'first');
            data_GMFAS_elong_Lichtwark = [data_GMFAS_licht_GM(:,1) horiz_l - horiz_l(loc_frame)];

       %     x axis - elong and pennation, y axis - fas len
            if plot_check && plot_individual
                plottitle = horzcat('IND Lichtwark elong vs angle, ', subject_id);
                figure('Name',plottitle);
                hold on
                yyaxis left
                plot(data_GMFAS_elong_Lichtwark(:,1), data_GMFAS_elong_Lichtwark(:,2),'LineWidth',2)
                plot(data_GMFAS_licht_GM(:,1), data_GMFAS_licht_GM(:,3),':b','LineWidth',1) % pennation
                ylabel('Elongation (mm) / Pennation angle (rad?)')
                yyaxis right
                plot(data_GMFAS_licht_GM(:,1), data_GMFAS_licht_GM(:,2),'--r','LineWidth',1) % faslen
                ylabel('Fascicle length (mm)')
                xlabel('Gonio angle (°)')
                title(plottitle)
                legend('Elongation', 'Pennation angle','Fascicle length','Location','Southeast')
                saveas(gcf, horzcat('data_plots/', plottitle,'.jpg'))
            end
        else % GM licht data don't exist
            data_GMFAS_elong_Lichtwark = 0;
        end
        %% 

        
         
        %% calculate and plot arrays of LENGTH, ELONGATION, STRAIN

        % MTU_length_array =                            MTU_elong_array =      MTU_strain_array =
        %1       angle
        %2       calc to SOL ins = free tendon                 
        %3       calc go GM ins = GM tendon
        %4       calc to knee = entire calf
        %5       DISPLACEMENT from GMFAS tracking       GM muscle elongation from faslen + penn.ang. (Lichtwark/Fukunaga)
        %6       SOL to GM = GM apo
        %7       GM to knee = GM msc
        %8       SOL length (H&H) minus free AT = SOL msc
        %9       GM tendon Fukunaga
        %10      GM muscle Fukunaga
        
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
        
        % calculate MTU lengths
        [MTU_length_array, MTU_elong_array, MTU_strain_array] = calculate_mtu_length(data_SOL(:,2:3), data_GMMTJ(:,2:3), data_GMFAS(:,2:3), data_GMFAS_elong_Lichtwark, dm_at_SOL_length{line}, dm_at_GM_length{line}, dm_leg_length{line}, out_ROM_trial_max);
        
        if plot_check && plot_individual
            % raw lengths (mm)
            plottitle = horzcat('IND MTU length vs angle, ', subject_id);
            figure('Name',plottitle);
            hold on
            plot(MTU_length_array(:,1),MTU_length_array(:,4),'LineWidth',1.5,'Color','blue') % AT + GM apo + GM msc = full GM MTU / calf length
            plot(MTU_length_array(:,1),MTU_length_array(:,8)+MTU_length_array(:,2),'LineWidth',1.5,'Color','black') % AT + SOL msc = SOL MTU
            plot(MTU_length_array(:,1),MTU_length_array(:,3),'LineWidth',1.5,'Color','red') % AT + GM apo = GM tendon (linear)
            plot(MTU_length_array(:,1),MTU_length_array(:,9),'LineWidth',1.5,'LineStyle','--','Color','red') % GM tendon from anthropometry Lichtwark/Fukunaga
            plot(MTU_length_array(:,1),MTU_length_array(:,2),'LineWidth',1.5,'Color','green') % AT
            % plot vertical lines to show color coding of subsequent figs
            plot([0,0], [MTU_length_array(1,4), MTU_length_array(1,3)],'Color','cyan')
            plot([1,1], [MTU_length_array(1,8)+MTU_length_array(1,2), MTU_length_array(1,3)],'Color','black')
            plot([0,0], [MTU_length_array(1,3), MTU_length_array(1,2)],'Color',[1 0.75 0])
            plot([0,0], [MTU_length_array(1,2), 0],'Color','green')
            axis(axis_len_MTU)
            ylabel('Length (mm)')
            xlabel('Gonio angle (°)')
            title(plottitle)
            legend('Full GM MTU', 'SOL MTU', 'GM tendon (linear)', 'GM tendon (anthro)', 'AT free', 'GM msc.', 'SOL msc.', 'GM apo.', 'Location','Southeast');
            saveas(gcf, horzcat('data_plots/', plottitle,'.jpg'))
            
            % raw elongation (mm)
            plottitle = horzcat('IND MTU elongation vs angle, ', subject_id);
            figure('Name',plottitle);
            hold on
            plot(MTU_elong_array(:,1),MTU_elong_array(:,4),'LineWidth',1.5,'Color','blue') % MTU = AT + GM apo + msc
            plot(MTU_elong_array(:,1),MTU_elong_array(:,2),'LineWidth',1.5,'Color','green') % AT
            plot(MTU_elong_array(:,1),MTU_elong_array(:,3),'LineWidth',1.5,'Color','red') % GM tendon = AT + GM apo
            plot(MTU_elong_array(:,1),MTU_elong_array(:,9),'LineWidth',1.5,'Color','red','LineStyle','--') % GM tendon from anthropometry Lichtwark/Fukunaga
            plot(MTU_elong_array(:,1),MTU_elong_array(:,7),'LineWidth',1.5,'Color','cyan') % GM msc
            plot(MTU_elong_array(:,1),MTU_elong_array(:,10),'LineWidth',1.5,'Color','cyan','LineStyle','--') % GM msc from anthropometry Lichtwark/Fukunaga
            plot(MTU_elong_array(:,1),MTU_elong_array(:,8),'LineWidth',1.5,'Color','black') % SOL msc
            plot(MTU_elong_array(:,1),MTU_elong_array(:,5),'LineWidth',1.5,'Color','yellow') % GM FAS displacement
            plot(MTU_elong_array(:,1),MTU_elong_array(:,6),'LineWidth',1.5,'Color',[1 0.75 0]) % GM apo
            axis(axis_ind_elong)
            ylabel('Elongation (mm)')
            xlabel('Gonio angle (°)')
            title(plottitle)
            legend('GM MTU', 'AT free', 'GM tendon (linear)', 'GM tendon (anthro)', 'GM msc. (linear)', 'GM msc. (anthro)', 'SOL msc.', 'GM fasc.DISPL.', 'GM apo.', 'Location','Northwest');
            saveas(gcf, horzcat('data_plots/', plottitle,'.jpg'))
            xlim([-0.5 inf]) 
            ylim([-inf inf]) 
            saveas(gcf, horzcat('data_plots/', plottitle,' ZOOM.jpg'))
            
            % strain (percent of initial length)
            plottitle = horzcat('IND MTU strain vs angle, ', subject_id);
            figure('Name',plottitle);
            hold on
            plot(MTU_strain_array(:,1),MTU_strain_array(:,4),'LineWidth',1.5,'Color','blue') % MTU = AT + GM apo + msc
            plot(MTU_strain_array(:,1),MTU_strain_array(:,2),'LineWidth',1.5,'Color','green') % AT
            plot(MTU_strain_array(:,1),MTU_strain_array(:,3),'LineWidth',1.5,'Color','red') % GM tendon = AT + GM apo
            plot(MTU_strain_array(:,1),MTU_strain_array(:,9),'LineWidth',1.5,'Color','red','LineStyle','--') % GM tendon from anthropometry Lichtwark/Fukunaga
            plot(MTU_strain_array(:,1),MTU_strain_array(:,7),'LineWidth',1.5,'Color','cyan') % GM msc
            plot(MTU_strain_array(:,1),MTU_strain_array(:,10),'LineWidth',1.5,'Color','cyan','LineStyle','--') % GM msc from anthropometry Lichtwark/Fukunaga
            plot(MTU_strain_array(:,1),MTU_strain_array(:,8),'LineWidth',1.5,'Color','black') % SOL msc
            plot(MTU_strain_array(:,1),MTU_strain_array(:,6),'LineWidth',1.5,'Color',[1 0.75 0]) % GM apo
            if max(MTU_strain_array(:,2)) > axis_ind_strain(4)
                text(axis_ind_strain(2)-10,axis_ind_strain(4)-1, horzcat('Max AT str = ', num2str(round(max(MTU_strain_array(:,2)),1))),'Color','green') % TEXT: max AT strain
            end
            axis(axis_ind_strain)
            ylabel('Strain (% of initial length)')
            xlabel('Gonio angle (°)')
            title(plottitle)
            legend('GM MTU', 'AT free', 'GM tendon (linear)', 'GM tendon (anthro)', 'GM msc. (linear)', 'GM msc. (anthro)', 'SOL msc.', 'GM apo.', 'Location','Northwest');
            saveas(gcf, horzcat('data_plots/', plottitle,'.jpg'))
            xlim([-0.5 inf]) 
            ylim([-inf inf]) 
            saveas(gcf, horzcat('data_plots/', plottitle,' ZOOM.jpg'))
        end
        %%
        
        
        
        %% extract FORCE, ANGLE, LENGTH/ELONGATION/STRAIN, EMG at various joint angles
        %   force, angle, EMG from all 3 scan locations / 6 trials, averaged
        %   displacement of SOL, GMMTJ, GMFAS from 2 trials per scan location
        %   length, elongation, strain of AT, GM apo, msc etc - from displacement data
        
        %   joint angles = 
        %       out_ROM_trial_max = trial max (different PRE, POST, L, R)
        %       out_ROM_ind_max = subject ind max (lowest PRE, POST, L, R)
        %       out_ROM_common_max = common max (lowest of all subjects)
        %       out_ROM_submax_1 OLD = additional predetermined angle: 1/3 of trial max ROM
        %       out_ROM_submax_2 = additional predetermined angle: 2/3 of trial max ROM
        %       out_ROM_submax_1 = 10 degrees for those subjects which have this ROM available 

        
        
        % column choices force, torque, displacement
        col_force = 1;      % from data_force_gonio / data_SOL etc
        col_angle = 2;
        col_displ = 3;
        
        % print error if max angles do not exist --> create_angles_passive needs to be run
        if str2double(input_gon_ind_max{trial_subjectno}) == 100
            cprintf('*red', 'ERROR: Max ROM values are not calculated for subject. Run create_angles_passive first.\n')
        end

        % find relevant goniometer angles 
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
        out_ROM_ind_max = str2double(input_gon_ind_max{trial_subjectno}) - 0.00000001;
        out_ROM_common_max = str2double(input_gon_common_max{trial_subjectno}) - 0.00000001;
       % out_ROM_submax_1 = out_ROM_trial_max * 1/3; %VAR --- old, replaced by 10 degrees (reused submax_1 for simplicity)
        out_ROM_submax_2 =  out_ROM_trial_max * 2/3; %VAR
        out_ROM_submax_1 = 10; %VAR
        
        % forces (using data_force_gonio = averaged data from 3 scan locations / 6 trials)
        loc_frame = find(data_force_gonio(:,col_angle)>=out_ROM_trial_max,1,'first'); % when averaging angles across trials, intervals of 0.05 degrees are used. This means that any required angle will exist in all data series
        out_F_trial_max_ROM = data_force_gonio(loc_frame,col_force); % force at highest angle in array
        out_F_trial_max_F = max(data_force_gonio(:,col_force)); % highest force in array
        loc_frame = find(data_force_gonio(:,col_angle)>=out_ROM_ind_max,1,'first'); 
        out_F_ind_max = data_force_gonio(loc_frame,col_force);
        loc_frame = find(data_force_gonio(:,col_angle)>=out_ROM_common_max,1,'first'); 
        out_F_common_max = data_force_gonio(loc_frame,col_force);
        loc_frame = find(data_force_gonio(:,col_angle)>=0,1,'first'); % zero angle
        out_F_zero = data_force_gonio(loc_frame,col_force);
        loc_frame = find(data_force_gonio(:,col_angle)>=out_ROM_submax_1,1,'first'); 
        if isempty(loc_frame) == 0
            out_F_submax_1 = data_force_gonio(loc_frame,col_force);
            out_T_submax_1 = out_F_submax_1 *at_momentarm;
        else
            out_F_submax_1 = 100;
            out_T_submax_1 = 100;
        end
        loc_frame = find(data_force_gonio(:,col_angle)>=out_ROM_submax_2,1,'first'); 
        out_F_submax_2 = data_force_gonio(loc_frame,col_force);
        
        % torques
        out_T_trial_max_ROM = out_F_trial_max_ROM *at_momentarm;
        out_T_trial_max_F = out_F_trial_max_F *at_momentarm;
        out_T_ind_max = out_F_ind_max *at_momentarm;
        out_T_common_max = out_F_common_max *at_momentarm;
        out_T_zero = out_F_zero *at_momentarm;
        % submex 2 - with force above
        out_T_submax_2 = out_F_submax_2 *at_momentarm;

        % displacements SOL
        loc_frame = find(data_SOL(:,col_angle)>=out_ROM_trial_max,1,'first'); 
        out_displ_SOL_trial_max = data_SOL(loc_frame,col_displ);
        loc_frame = find(data_SOL(:,col_angle)>=out_ROM_ind_max,1,'first'); 
        out_displ_SOL_ind_max = data_SOL(loc_frame,col_displ); 
        loc_frame = find(data_SOL(:,col_angle)>=out_ROM_common_max,1,'first'); 
        out_displ_SOL_common_max = data_SOL(loc_frame,col_displ); 
        loc_frame = find(data_SOL(:,col_angle)>=out_ROM_submax_1,1,'first'); 
        if isempty(loc_frame) == 0
            out_displ_SOL_submax_1 = data_SOL(loc_frame,col_displ); 
        else
            out_displ_SOL_submax_1 = 100;
        end
        loc_frame = find(data_SOL(:,col_angle)>=out_ROM_submax_2,1,'first'); 
        out_displ_SOL_submax_2 = data_SOL(loc_frame,col_displ); 

        % displacements GMMTJ
        loc_frame = find(data_GMMTJ(:,col_angle)>=out_ROM_trial_max,1,'first'); 
        out_displ_GMMTJ_trial_max = data_GMMTJ(loc_frame,col_displ); 
        loc_frame = find(data_GMMTJ(:,col_angle)>=out_ROM_ind_max,1,'first'); 
        out_displ_GMMTJ_ind_max = data_GMMTJ(loc_frame,col_displ); 
        loc_frame = find(data_GMMTJ(:,col_angle)>=out_ROM_common_max,1,'first'); 
        out_displ_GMMTJ_common_max = data_GMMTJ(loc_frame,col_displ); 
        loc_frame = find(data_GMMTJ(:,col_angle)>=out_ROM_submax_1,1,'first'); 
        if isempty(loc_frame) == 0
            out_displ_GMMTJ_submax_1 = data_GMMTJ(loc_frame,col_displ); 
        else
            out_displ_GMMTJ_submax_1 = 100;
        end
        loc_frame = find(data_GMMTJ(:,col_angle)>=out_ROM_submax_2,1,'first'); 
        out_displ_GMMTJ_submax_2 = data_GMMTJ(loc_frame,col_displ); 

        % displacements GMFAS
        loc_frame = find(data_GMFAS(:,col_angle)>=out_ROM_trial_max,1,'first'); 
        out_displ_GMFAS_trial_max = data_GMFAS(loc_frame,col_displ); 
        loc_frame = find(data_GMFAS(:,col_angle)>=out_ROM_ind_max,1,'first'); 
        out_displ_GMFAS_ind_max = data_GMFAS(loc_frame,col_displ); 
        loc_frame = find(data_GMFAS(:,col_angle)>=out_ROM_common_max,1,'first'); 
        out_displ_GMFAS_common_max = data_GMFAS(loc_frame,col_displ); 
        loc_frame = find(data_GMFAS(:,col_angle)>=out_ROM_submax_1,1,'first'); 
        if isempty(loc_frame) == 0
            out_displ_GMFAS_submax_1 = data_GMFAS(loc_frame,col_displ); 
        else
            out_displ_GMFAS_submax_1 = 100;
        end
        loc_frame = find(data_GMFAS(:,col_angle)>=out_ROM_submax_2,1,'first'); 
        out_displ_GMFAS_submax_2 = data_GMFAS(loc_frame,col_displ); 

        % column choices ELONGATION / LENGTH / STRAIN
        col_angle_elong = 1;
        col_AT = 2;
        col_GM_tend = 3;
        col_leg = 4;
        %col_GMFAS = 5; % in elongation: GMFAS displacement
        col_GM_apo = 6;
        col_msc_GM = 7;
        col_msc_SOL = 8;
        col_GM_tend_Fukunaga = 9;
        col_GM_msc_Fukunaga = 10;
        
        loc_frame = find(MTU_elong_array(:,col_angle_elong)>=out_ROM_trial_max,1,'first'); 
        out_elong_AT_trial_max = MTU_elong_array(loc_frame,col_AT);
        out_elong_GMtend_trial_max = MTU_elong_array(loc_frame,col_GM_tend);
        out_elong_leg_trial_max = MTU_elong_array(loc_frame,col_leg);
        out_elong_GMapo_trial_max = MTU_elong_array(loc_frame,col_GM_apo);
        out_elong_msc_GM_trial_max = MTU_elong_array(loc_frame,col_msc_GM);
        out_elong_msc_SOL_trial_max = MTU_elong_array(loc_frame,col_msc_SOL);
        out_strain_AT_trial_max = MTU_strain_array(loc_frame,col_AT);
        out_strain_GMtend_trial_max = MTU_strain_array(loc_frame,col_GM_tend);
        out_strain_leg_trial_max = MTU_strain_array(loc_frame,col_leg);
        out_strain_GMapo_trial_max = MTU_strain_array(loc_frame,col_GM_apo);
        out_strain_msc_GM_trial_max = MTU_strain_array(loc_frame,col_msc_GM);
        out_strain_msc_SOL_trial_max = MTU_strain_array(loc_frame,col_msc_SOL);
        out_length_AT_trial_max = MTU_length_array(loc_frame,col_AT);
        out_length_GMtend_trial_max = MTU_length_array(loc_frame,col_GM_tend);
        out_length_leg_trial_max = MTU_length_array(loc_frame,col_leg);
        out_length_GMapo_trial_max = MTU_length_array(loc_frame,col_GM_apo);
        out_length_msc_GM_trial_max = MTU_length_array(loc_frame,col_msc_GM);
        out_length_msc_SOL_trial_max = MTU_length_array(loc_frame,col_msc_SOL);
        out_contrib_GM_trial_max = 100 * MTU_elong_array(loc_frame,col_GM_msc_Fukunaga) / MTU_elong_array(loc_frame,col_leg);
        out_elong_GMtend_Fuku_trial_max = MTU_elong_array(loc_frame,col_GM_tend_Fukunaga);
        out_elong_msc_GM_Fuku_trial_max = MTU_elong_array(loc_frame,col_GM_msc_Fukunaga);
        out_strain_GMtend_Fuku_trial_max = MTU_strain_array(loc_frame,col_GM_tend_Fukunaga);
        out_strain_msc_GM_Fuku_trial_max = MTU_strain_array(loc_frame,col_GM_msc_Fukunaga);
        out_length_GMtend_Fuku_trial_max = MTU_length_array(loc_frame,col_GM_tend_Fukunaga);
        out_length_msc_GM_Fuku_trial_max = MTU_length_array(loc_frame,col_GM_msc_Fukunaga);
        
        loc_frame = find(MTU_elong_array(:,col_angle_elong)>=out_ROM_ind_max,1,'first'); 
        out_elong_AT_ind_max = MTU_elong_array(loc_frame,col_AT); 
        out_elong_GMtend_ind_max = MTU_elong_array(loc_frame,col_GM_tend);
        out_elong_leg_ind_max = MTU_elong_array(loc_frame,col_leg);
        out_elong_GMapo_ind_max = MTU_elong_array(loc_frame,col_GM_apo); 
        out_elong_msc_GM_ind_max = MTU_elong_array(loc_frame,col_msc_GM); 
        out_elong_msc_SOL_ind_max = MTU_elong_array(loc_frame,col_msc_SOL); 
        out_strain_AT_ind_max = MTU_strain_array(loc_frame,col_AT); 
        out_strain_GMtend_ind_max = MTU_strain_array(loc_frame,col_GM_tend);
        out_strain_leg_ind_max = MTU_strain_array(loc_frame,col_leg);
        out_strain_GMapo_ind_max = MTU_strain_array(loc_frame,col_GM_apo); 
        out_strain_msc_GM_ind_max = MTU_strain_array(loc_frame,col_msc_GM); 
        out_strain_msc_SOL_ind_max = MTU_strain_array(loc_frame,col_msc_SOL); 
        out_length_AT_ind_max = MTU_length_array(loc_frame,col_AT); 
        out_length_GMtend_ind_max = MTU_length_array(loc_frame,col_GM_tend);
        out_length_leg_ind_max = MTU_length_array(loc_frame,col_leg);
        out_length_GMapo_ind_max = MTU_length_array(loc_frame,col_GM_apo); 
        out_length_msc_GM_ind_max = MTU_length_array(loc_frame,col_msc_GM); 
        out_length_msc_SOL_ind_max = MTU_length_array(loc_frame,col_msc_SOL); 
        out_contrib_GM_ind_max = 100 * MTU_elong_array(loc_frame,col_GM_msc_Fukunaga) / MTU_elong_array(loc_frame,col_leg);
        out_elong_GMtend_Fuku_ind_max = MTU_elong_array(loc_frame,col_GM_tend_Fukunaga);
        out_elong_msc_GM_Fuku_ind_max = MTU_elong_array(loc_frame,col_GM_msc_Fukunaga);
        out_strain_GMtend_Fuku_ind_max = MTU_strain_array(loc_frame,col_GM_tend_Fukunaga);
        out_strain_msc_GM_Fuku_ind_max = MTU_strain_array(loc_frame,col_GM_msc_Fukunaga);
        out_length_GMtend_Fuku_ind_max = MTU_length_array(loc_frame,col_GM_tend_Fukunaga);
        out_length_msc_GM_Fuku_ind_max = MTU_length_array(loc_frame,col_GM_msc_Fukunaga);
        
        loc_frame = find(MTU_elong_array(:,col_angle_elong)>=out_ROM_common_max,1,'first'); 
        out_elong_AT_common_max = MTU_elong_array(loc_frame,col_AT); 
        out_elong_GMtend_common_max = MTU_elong_array(loc_frame,col_GM_tend);
        out_elong_leg_common_max = MTU_elong_array(loc_frame,col_leg);
        out_elong_GMapo_common_max = MTU_elong_array(loc_frame,col_GM_apo); 
        out_elong_msc_GM_common_max = MTU_elong_array(loc_frame,col_msc_GM); 
        out_elong_msc_SOL_common_max = MTU_elong_array(loc_frame,col_msc_SOL); 
        out_strain_AT_common_max = MTU_strain_array(loc_frame,col_AT); 
        out_strain_GMtend_common_max = MTU_strain_array(loc_frame,col_GM_tend);
        out_strain_leg_common_max = MTU_strain_array(loc_frame,col_leg);
        out_strain_GMapo_common_max = MTU_strain_array(loc_frame,col_GM_apo); 
        out_strain_msc_GM_common_max = MTU_strain_array(loc_frame,col_msc_GM); 
        out_strain_msc_SOL_common_max = MTU_strain_array(loc_frame,col_msc_SOL); 
        out_length_AT_common_max = MTU_length_array(loc_frame,col_AT); 
        out_length_GMtend_common_max = MTU_length_array(loc_frame,col_GM_tend);
        out_length_leg_common_max = MTU_length_array(loc_frame,col_leg);
        out_length_GMapo_common_max = MTU_length_array(loc_frame,col_GM_apo); 
        out_length_msc_GM_common_max = MTU_length_array(loc_frame,col_msc_GM); 
        out_length_msc_SOL_common_max = MTU_length_array(loc_frame,col_msc_SOL); 
        out_contrib_GM_common_max = 100 * MTU_elong_array(loc_frame,col_GM_msc_Fukunaga) / MTU_elong_array(loc_frame,col_leg);
        out_elong_GMtend_Fuku_common_max = MTU_elong_array(loc_frame,col_GM_tend_Fukunaga);
        out_elong_msc_GM_Fuku_common_max = MTU_elong_array(loc_frame,col_GM_msc_Fukunaga);
        out_strain_GMtend_Fuku_common_max = MTU_strain_array(loc_frame,col_GM_tend_Fukunaga);
        out_strain_msc_GM_Fuku_common_max = MTU_strain_array(loc_frame,col_GM_msc_Fukunaga);
        out_length_GMtend_Fuku_common_max = MTU_length_array(loc_frame,col_GM_tend_Fukunaga);
        out_length_msc_GM_Fuku_common_max = MTU_length_array(loc_frame,col_GM_msc_Fukunaga);
        
        loc_frame = find(MTU_elong_array(:,col_angle_elong)>=out_ROM_submax_1,1,'first'); 
        if isempty(loc_frame)
            out_elong_AT_submax_1 = 100;
            out_elong_GMtend_submax_1 = 100;
            out_elong_leg_submax_1 = 100;
            out_elong_GMapo_submax_1 = 100;
            out_elong_msc_GM_submax_1 = 100;
            out_elong_msc_SOL_submax_1 = 100;
            out_strain_AT_submax_1 = 100;
            out_strain_GMtend_submax_1 = 100;
            out_strain_leg_submax_1 = 100;
            out_strain_GMapo_submax_1 = 100;
            out_strain_msc_GM_submax_1 = 100;
            out_strain_msc_SOL_submax_1 = 100;
            out_length_AT_submax_1 = 100;
            out_length_GMtend_submax_1 = 100;
            out_length_leg_submax_1 = 100;
            out_length_GMapo_submax_1 = 100;
            out_length_msc_GM_submax_1 = 100;
            out_length_msc_SOL_submax_1 = 100;
            out_contrib_GM_submax_1 = 100;
            out_elong_GMtend_Fuku_submax_1 = 100;
            out_elong_msc_GM_Fuku_submax_1 = 100;
            out_strain_GMtend_Fuku_submax_1 = 100;
            out_strain_msc_GM_Fuku_submax_1 = 100;
            out_length_GMtend_Fuku_submax_1 = 100;
            out_length_msc_GM_Fuku_submax_1 = 100;
        else
            out_elong_AT_submax_1 = MTU_elong_array(loc_frame,col_AT); 
            out_elong_GMtend_submax_1 = MTU_elong_array(loc_frame,col_GM_tend);
            out_elong_leg_submax_1 = MTU_elong_array(loc_frame,col_leg);
            out_elong_GMapo_submax_1 = MTU_elong_array(loc_frame,col_GM_apo); 
            out_elong_msc_GM_submax_1 = MTU_elong_array(loc_frame,col_msc_GM); 
            out_elong_msc_SOL_submax_1 = MTU_elong_array(loc_frame,col_msc_SOL); 
            out_strain_AT_submax_1 = MTU_strain_array(loc_frame,col_AT); 
            out_strain_GMtend_submax_1 = MTU_strain_array(loc_frame,col_GM_tend);
            out_strain_leg_submax_1 = MTU_strain_array(loc_frame,col_leg);
            out_strain_GMapo_submax_1 = MTU_strain_array(loc_frame,col_GM_apo); 
            out_strain_msc_GM_submax_1 = MTU_strain_array(loc_frame,col_msc_GM); 
            out_strain_msc_SOL_submax_1 = MTU_strain_array(loc_frame,col_msc_SOL); 
            out_length_AT_submax_1 = MTU_length_array(loc_frame,col_AT); 
            out_length_GMtend_submax_1 = MTU_length_array(loc_frame,col_GM_tend);
            out_length_leg_submax_1 = MTU_length_array(loc_frame,col_leg);
            out_length_GMapo_submax_1 = MTU_length_array(loc_frame,col_GM_apo); 
            out_length_msc_GM_submax_1 = MTU_length_array(loc_frame,col_msc_GM); 
            out_length_msc_SOL_submax_1 = MTU_length_array(loc_frame,col_msc_SOL); 
            out_contrib_GM_submax_1 = 100 * MTU_elong_array(loc_frame,col_GM_msc_Fukunaga) / MTU_elong_array(loc_frame,col_leg);
            out_elong_GMtend_Fuku_submax_1 = MTU_elong_array(loc_frame,col_GM_tend_Fukunaga);
            out_elong_msc_GM_Fuku_submax_1 = MTU_elong_array(loc_frame,col_GM_msc_Fukunaga);
            out_strain_GMtend_Fuku_submax_1 = MTU_strain_array(loc_frame,col_GM_tend_Fukunaga);
            out_strain_msc_GM_Fuku_submax_1 = MTU_strain_array(loc_frame,col_GM_msc_Fukunaga);
            out_length_GMtend_Fuku_submax_1 = MTU_length_array(loc_frame,col_GM_tend_Fukunaga);
            out_length_msc_GM_Fuku_submax_1 = MTU_length_array(loc_frame,col_GM_msc_Fukunaga);
        end
        
        loc_frame = find(MTU_elong_array(:,col_angle_elong)>=out_ROM_submax_2,1,'first');
        out_elong_AT_submax_2 = MTU_elong_array(loc_frame,col_AT); 
        out_elong_GMtend_submax_2 = MTU_elong_array(loc_frame,col_GM_tend);
        out_elong_leg_submax_2 = MTU_elong_array(loc_frame,col_leg);
        out_elong_GMapo_submax_2 = MTU_elong_array(loc_frame,col_GM_apo); 
        out_elong_msc_GM_submax_2 = MTU_elong_array(loc_frame,col_msc_GM); 
        out_elong_msc_SOL_submax_2 = MTU_elong_array(loc_frame,col_msc_SOL); 
        out_strain_AT_submax_2 = MTU_strain_array(loc_frame,col_AT); 
        out_strain_GMtend_submax_2 = MTU_strain_array(loc_frame,col_GM_tend);
        out_strain_leg_submax_2 = MTU_strain_array(loc_frame,col_leg);
        out_strain_GMapo_submax_2 = MTU_strain_array(loc_frame,col_GM_apo); 
        out_strain_msc_GM_submax_2 = MTU_strain_array(loc_frame,col_msc_GM); 
        out_strain_msc_SOL_submax_2 = MTU_strain_array(loc_frame,col_msc_SOL); 
        out_length_AT_submax_2 = MTU_length_array(loc_frame,col_AT); 
        out_length_GMtend_submax_2 = MTU_length_array(loc_frame,col_GM_tend);
        out_length_leg_submax_2 = MTU_length_array(loc_frame,col_leg);
        out_length_GMapo_submax_2 = MTU_length_array(loc_frame,col_GM_apo); 
        out_length_msc_GM_submax_2 = MTU_length_array(loc_frame,col_msc_GM); 
        out_length_msc_SOL_submax_2 = MTU_length_array(loc_frame,col_msc_SOL); 
        out_contrib_GM_submax_2 = 100 * MTU_elong_array(loc_frame,col_GM_msc_Fukunaga) / MTU_elong_array(loc_frame,col_leg);
        out_elong_GMtend_Fuku_submax_2 = MTU_elong_array(loc_frame,col_GM_tend_Fukunaga);
        out_elong_msc_GM_Fuku_submax_2 = MTU_elong_array(loc_frame,col_GM_msc_Fukunaga);
        out_strain_GMtend_Fuku_submax_2 = MTU_strain_array(loc_frame,col_GM_tend_Fukunaga);
        out_strain_msc_GM_Fuku_submax_2 = MTU_strain_array(loc_frame,col_GM_msc_Fukunaga);
        out_length_GMtend_Fuku_submax_2 = MTU_length_array(loc_frame,col_GM_tend_Fukunaga);
        out_length_msc_GM_Fuku_submax_2 = MTU_length_array(loc_frame,col_GM_msc_Fukunaga);
        
        
        
        % EMG (from all 3 scan locations / 6 trials, averaged)
        
        emg_step = 9; %VAR - number of EMG values BEFORE relevant angle, to include in average. 9 values = a span of 0.5 degrees.
        % cannot average values AROUND relevant angle, because for a few (least flexible) trials, we want the data AT the last available angle.

        % identify locations of the relevant goniometer angles, in the array with all 2+2+2 trials averaged:
        loc_frame_trial_max = find(data_force_gonio(:,col_angle)>=out_ROM_trial_max,1,'first'); % when averaging angles across trials, intervals of 0.05 degrees are used. This means that any required angle will exist in all data series
        loc_frame_ind_max = find(data_force_gonio(:,col_angle)>=out_ROM_ind_max,1,'first'); 
        loc_frame_common_max = find(data_force_gonio(:,col_angle)>=out_ROM_common_max,1,'first'); 
        loc_frame_submax_1 = find(data_force_gonio(:,col_angle)>=out_ROM_submax_1,1,'first'); 
        loc_frame_submax_2 = find(data_force_gonio(:,col_angle)>=out_ROM_submax_2,1,'first'); 
        loc_frame_zero = find(data_force_gonio(:,col_angle)>=0,1,'first'); 

        % break if frame not found
        if isempty(loc_frame_ind_max)
            % cprintf('error', horzcat('ERROR: Computed trial max ROM (', num2str(out_ROM_trial_max), ') does not exist in the current data series (max = ', num2str(max(data_force_gonio(:,col_angle))), '). Check "create_angles_passive" vs current data.\n' ));
            error(horzcat('ERROR: Computed trial max ROM (', num2str(out_ROM_trial_max), ') does not exist in the current data series (max = ', num2str(max(data_force_gonio(:,col_angle))), '). Check "create_angles_passive" vs current data.\n' ))
        end
        
        % EMG gm
        out_emg_gm_trial_max = mean(data_force_gonio(loc_frame_trial_max-emg_step:loc_frame_trial_max,3)); % EMG gm = column 3
        out_emg_gm_ind_max = mean(data_force_gonio(loc_frame_ind_max-emg_step:loc_frame_ind_max,3));
        out_emg_gm_common_max = mean(data_force_gonio(loc_frame_common_max-emg_step:loc_frame_common_max,3));
        out_emg_gm_submax_2 = mean(data_force_gonio(loc_frame_submax_2-emg_step:loc_frame_submax_2,3));

        % EMG gl
        out_emg_gl_trial_max = mean(data_force_gonio(loc_frame_trial_max-emg_step:loc_frame_trial_max,4)); % EMG gl = column 4
        out_emg_gl_ind_max = mean(data_force_gonio(loc_frame_ind_max-emg_step:loc_frame_ind_max,4)); 
        out_emg_gl_common_max = mean(data_force_gonio(loc_frame_common_max-emg_step:loc_frame_common_max,4));
        out_emg_gl_submax_2 = mean(data_force_gonio(loc_frame_submax_2-emg_step:loc_frame_submax_2,4));

        % EMG sol
        out_emg_sol_trial_max = mean(data_force_gonio(loc_frame_trial_max-emg_step:loc_frame_trial_max,5)); % EMG sol = column 5
        out_emg_sol_ind_max = mean(data_force_gonio(loc_frame_ind_max-emg_step:loc_frame_ind_max,5));
        out_emg_sol_common_max = mean(data_force_gonio(loc_frame_common_max-emg_step:loc_frame_common_max,5));
        out_emg_sol_submax_2 = mean(data_force_gonio(loc_frame_submax_2-emg_step:loc_frame_submax_2,5));

        if isempty(loc_frame_submax_1)
            out_emg_gm_submax_1 = 100;
            out_emg_gl_submax_1 = 100;
            out_emg_sol_submax_1 = 100;
        else
            out_emg_gm_submax_1 = mean(data_force_gonio(loc_frame_submax_1-emg_step:loc_frame_submax_1,3));
            out_emg_gl_submax_1 = mean(data_force_gonio(loc_frame_submax_1-emg_step:loc_frame_submax_1,4));
            out_emg_sol_submax_1 = mean(data_force_gonio(loc_frame_submax_1-emg_step:loc_frame_submax_1,5));
        end
        %%
        
        
        
        %% extract FASCICLE LENGTH, PENNATION ANGLE, GM ELONG at various joint angles
        %     data_GMFAS_licht_GM
        %     data_GMFAS_licht_SOL
        % containing:
        %   averaged angle (currently calculated from gonio)
        %   averaged fasicle length
        %   averaged pennation angle
        % OR containing zeros (if nonexistent)
        
        %    data_GM_elong_Lichtwark
        % containing:
        %  angle
        %  elongation (Lichtwark/Fukunaga)
        
        % set columns:
        col_lichtangle = 1;
        col_lichtfas = 2;   % from data_GMFAS_licht
        col_lichtpenn = 3;
        % col_lichtelong = 2; % from data_GM_elong_Lichtwark
        
        % preallocate / set zero values:
        
        % --- GM pennation and fascicle length
        out_licht_faslen_GM_trial_max = 100;
        out_licht_pennation_GM_trial_max = 100;
        out_licht_faslen_GM_common_max = 100;
        out_licht_pennation_GM_common_max = 100;
        out_licht_faslen_GM_ind_max = 100;
        out_licht_pennation_GM_ind_max = 100;
        out_licht_faslen_GM_zero = 100;
        out_licht_pennation_GM_zero = 100;
        out_licht_faslen_GM_submax_1 = 100;
        out_licht_pennation_GM_submax_1 = 100;
        out_licht_faslen_GM_submax_2 = 100;
        out_licht_pennation_GM_submax_2 = 100;
        % --- GM muscle elongation from Lichtwark/Fukunaga - below vars are
        %       not used from this part of analysis - instead extracted from MTU
        %       length arrays
%         out_licht_elong_GM_trial_max = 100;
%         out_licht_elong_GM_common_max = 100;
%         out_licht_elong_GM_ind_max = 100;
%         out_licht_elong_GM_zero = 100;
%         out_licht_elong_GM_submax_1 = 100;
%         out_licht_elong_GM_submax_2 = 100;
        % --- SOL pennation and fascicle length
        out_licht_faslen_SOL_trial_max = 100;
        out_licht_pennation_SOL_trial_max = 100;
        out_licht_faslen_SOL_common_max = 100;
        out_licht_pennation_SOL_common_max = 100;
        out_licht_faslen_SOL_ind_max = 100;
        out_licht_pennation_SOL_ind_max = 100;
        out_licht_faslen_SOL_zero = 100;
        out_licht_pennation_SOL_zero = 100;
        out_licht_faslen_SOL_submax_1 = 100;
        out_licht_pennation_SOL_submax_1 = 100;
        out_licht_faslen_SOL_submax_2 = 100;
        out_licht_pennation_SOL_submax_2 = 100;
            
        if data_GMFAS_licht_GM == 0
            % no licht data existing
        else
            % at trial max angle:
            loc_frame = find(data_GMFAS(:,col_angle)>=out_ROM_trial_max,1,'first'); 
            out_licht_faslen_GM_trial_max = data_GMFAS_licht_GM(loc_frame,col_lichtfas);
            out_licht_pennation_GM_trial_max = data_GMFAS_licht_GM(loc_frame,col_lichtpenn);
            %out_licht_elong_GM_trial_max = data_GM_elong_Lichtwark(loc_frame,col_lichtelong);
            if (length((data_GMFAS_licht_SOL)) == 3) == 0
                % SOL exists
                out_licht_faslen_SOL_trial_max = data_GMFAS_licht_SOL(loc_frame,col_lichtfas); 
                out_licht_pennation_SOL_trial_max = data_GMFAS_licht_SOL(loc_frame,col_lichtpenn); 
            end

            % at individual max angle (across sides/timepoints):
            loc_frame = find(data_GMFAS(:,col_angle)>=out_ROM_ind_max,1,'first'); 
            out_licht_faslen_GM_ind_max = data_GMFAS_licht_GM(loc_frame,col_lichtfas); 
            out_licht_pennation_GM_ind_max = data_GMFAS_licht_GM(loc_frame,col_lichtpenn); 
            %out_licht_elong_GM_ind_max = data_GM_elong_Lichtwark(loc_frame,col_lichtelong);
            if (length((data_GMFAS_licht_SOL)) == 3) == 0
                % SOL exists
                out_licht_faslen_SOL_ind_max = data_GMFAS_licht_SOL(loc_frame,col_lichtfas); 
                out_licht_pennation_SOL_ind_max = data_GMFAS_licht_SOL(loc_frame,col_lichtpenn); 
            end

            % at subject common max angle:
            loc_frame = find(data_GMFAS(:,col_angle)>=out_ROM_common_max,1,'first'); 
            out_licht_faslen_GM_common_max = data_GMFAS_licht_GM(loc_frame,col_lichtfas); 
            out_licht_pennation_GM_common_max = data_GMFAS_licht_GM(loc_frame,col_lichtpenn); 
            %out_licht_elong_GM_common_max = data_GM_elong_Lichtwark(loc_frame,col_lichtelong);
            if (length((data_GMFAS_licht_SOL)) == 3) == 0
                % SOL exists
                out_licht_faslen_SOL_common_max = data_GMFAS_licht_SOL(loc_frame,col_lichtfas); 
                out_licht_pennation_SOL_common_max = data_GMFAS_licht_SOL(loc_frame,col_lichtpenn); 
            end

            % at zero angle:
            loc_frame = find(data_GMFAS(:,2)>=0,1,'first'); 
            out_licht_faslen_GM_zero = data_GMFAS_licht_GM(loc_frame,col_lichtfas); 
            out_licht_pennation_GM_zero = data_GMFAS_licht_GM(loc_frame,col_lichtpenn); 
            %out_licht_elong_GM_zero = data_GM_elong_Lichtwark(loc_frame,col_lichtelong);
            if (length((data_GMFAS_licht_SOL)) == 3) == 0
                % SOL exists
                out_licht_faslen_SOL_zero = data_GMFAS_licht_SOL(loc_frame,col_lichtfas); 
                out_licht_pennation_SOL_zero = data_GMFAS_licht_SOL(loc_frame,col_lichtpenn); 
            end

            % at submax_1 angle:
            loc_frame = find(data_GMFAS(:,col_angle)>=out_ROM_submax_1,1,'first'); 
            % does angle exist in data?
            if isempty(loc_frame)
                %
            else
                % GM data exist:
                out_licht_faslen_GM_submax_1 = data_GMFAS_licht_GM(loc_frame,col_lichtfas); 
                out_licht_pennation_GM_submax_1 = data_GMFAS_licht_GM(loc_frame,col_lichtpenn); 
                %out_licht_elong_GM_submax_1 = data_GM_elong_Lichtwark(loc_frame,col_lichtelong);
                if (length((data_GMFAS_licht_SOL)) == 3) == 0
                    % SOL exists
                    out_licht_faslen_SOL_submax_1 = data_GMFAS_licht_SOL(loc_frame,col_lichtfas); 
                    out_licht_pennation_SOL_submax_1 = data_GMFAS_licht_SOL(loc_frame,col_lichtpenn); 
                end
            end
            
            % at submax_2 angle:
            loc_frame = find(data_GMFAS(:,col_angle)>=out_ROM_submax_2,1,'first'); 
            out_licht_faslen_GM_submax_2 = data_GMFAS_licht_GM(loc_frame,col_lichtfas); 
            out_licht_pennation_GM_submax_2 = data_GMFAS_licht_GM(loc_frame,col_lichtpenn); 
            %out_licht_elong_GM_submax_2 = data_GM_elong_Lichtwark(loc_frame,col_lichtelong);
            if (length((data_GMFAS_licht_SOL)) == 3) == 0
                % SOL exists
                out_licht_faslen_SOL_submax_2 = data_GMFAS_licht_SOL(loc_frame,col_lichtfas); 
                out_licht_pennation_SOL_submax_2 = data_GMFAS_licht_SOL(loc_frame,col_lichtpenn); 
            end
        end
        %% 
        
        
        


        %% PASSIVE STIFFNESS and STIFFNESS INDEX (Nordez 2006)
        
        % gonio angle = data_force_gonio(:,col_angle)
        % force = data_force_gonio(:,col_force)
        % multiplying force with at_momentarm to convert to torque
        
        
        
        %%% PASSIVE STIFFNESS:
        
        % passive stiffness = delta torque / delta angle, at various angles
        % fit 4th order polynomial to averaged torque-angle curve, using data from zero angle to subject's individual max ROM (across both legs)
        fit_ind_max = polyfit(data_force_gonio(loc_frame_zero:loc_frame_ind_max,col_angle), at_momentarm*data_force_gonio(loc_frame_zero:loc_frame_ind_max,col_force), 4);

        % extract passive stiffness (derivate of 4th order poly) at:
        out_pstiff_trial_max = (4 * fit_ind_max(1) * out_ROM_trial_max^3) + (3 * fit_ind_max(2) * out_ROM_trial_max^2) + (2 * fit_ind_max(3) * out_ROM_trial_max) + fit_ind_max(4);
        out_pstiff_ind_max = (4 * fit_ind_max(1) * out_ROM_ind_max^3) + (3 * fit_ind_max(2) * out_ROM_ind_max^2) + (2 * fit_ind_max(3) * out_ROM_ind_max) + fit_ind_max(4);
        out_pstiff_common_max = (4 * fit_ind_max(1) * out_ROM_common_max^3) + (3 * fit_ind_max(2) * out_ROM_common_max^2) + (2 * fit_ind_max(3) * out_ROM_common_max) + fit_ind_max(4);
        if out_ROM_trial_max > out_ROM_submax_1 
            out_pstiff_submax_1 = (4 * fit_ind_max(1) * out_ROM_submax_1^3) + (3 * fit_ind_max(2) * out_ROM_submax_1^2) + (2 * fit_ind_max(3) * out_ROM_submax_1) + fit_ind_max(4);
        else
            out_pstiff_submax_1 = 100;
        end
        out_pstiff_submax_2 = (4 * fit_ind_max(1) * out_ROM_submax_2^3) + (3 * fit_ind_max(2) * out_ROM_submax_2^2) + (2 * fit_ind_max(3) * out_ROM_submax_2) + fit_ind_max(4);
        out_pstiff_angle = 15; %VAR
        if out_ROM_trial_max > out_pstiff_angle
            out_pstiff_15 = (4 * fit_ind_max(1) * out_pstiff_angle^3) + (3 * fit_ind_max(2) * out_pstiff_angle^2) + (2 * fit_ind_max(3) * out_pstiff_angle) + fit_ind_max(4);
        else
            out_pstiff_15 = 100;
        end
        
%         % plotting various methods of curve fit equations:
%         fit_ind_max2 = polyfit(data_force_gonio(loc_frame_zero:loc_frame_ind_max,col_angle), at_momentarm*data_force_gonio(loc_frame_zero:loc_frame_ind_max,col_force), 3);
%         fit_ind_max3 = polyfit(data_force_gonio(loc_frame_zero:loc_frame_ind_max,col_angle), at_momentarm*data_force_gonio(loc_frame_zero:loc_frame_ind_max,col_force), 2);
%                 plottitle = horzcat('IND torque-angle fit for stiffness, ', subject_id);
%                 figure('Name',plottitle)
%                 hold on
%                 plot(data_force_gonio(1:loc_frame_ind_max,col_angle), at_momentarm*data_force_gonio(1:loc_frame_ind_max,col_force))
%                 plot(0:0.2:max([20 max(data_force_gonio(loc_frame_zero:loc_frame_ind_max,col_angle))]), polyval(fit_ind_max,0:0.2:max([20 max(data_force_gonio(loc_frame_zero:loc_frame_ind_max,col_angle))])),':')
%                 plot(0:0.2:max([20 max(data_force_gonio(loc_frame_zero:loc_frame_ind_max,col_angle))]), polyval(fit_ind_max2,0:0.2:max([20 max(data_force_gonio(loc_frame_zero:loc_frame_ind_max,col_angle))])),'--')
%                 plot(0:0.2:max([20 max(data_force_gonio(loc_frame_zero:loc_frame_ind_max,col_angle))]), polyval(fit_ind_max3,0:0.2:max([20 max(data_force_gonio(loc_frame_zero:loc_frame_ind_max,col_angle))])),'-.')
%                 axis([0 max([20 max(data_force_gonio(loc_frame_zero:loc_frame_ind_max,col_angle))]) 0 1.5*max(at_momentarm*data_force_gonio(loc_frame_zero:loc_frame_ind_max,col_force))])
%                 xlabel('Gonio angle (°)')
%                 ylabel('Force (N)')
%                 title(plottitle)
%                 legend('Raw data', 'Fit 4th order','Fit 3rd order','Fit 2nd order','Location','Northwest')
%                 saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))

                
        
        %%% STIFFNESS INDEX:
        
        % fit 2nd order polynomial to averaged torque-angle curve, using data from zero angle to ind max ROM:
        fit_ind_max = polyfit(data_force_gonio(loc_frame_zero:loc_frame_ind_max,col_angle), data_force_gonio(loc_frame_zero:loc_frame_ind_max,col_force), 2);

        % extract stiffness index as 2 * a
        out_pstiff_index = 2 * fit_ind_max(1);
        
        
        
        %%% AT STIFFNESS:     % MMM LATER
        
        % stiffness = delta force / delta length for AT ? check literature
        %%
        
        
        
        %% extract ANGLES at specific FORCE levels
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
            else % L
                loc_F_trial_max = str2double(input_for_pre_l{trial_subjectno}) - 0.0001;
            end
        else % POST
            if strcmpi(dm_side{line}, 'R') == 1
                loc_F_trial_max = str2double(input_for_post_r{trial_subjectno}) - 0.0001;
            else % L
                loc_F_trial_max = str2double(input_for_post_l{trial_subjectno}) - 0.0001;
            end
        end
        loc_F_ind_max = str2double(input_for_ind_max{trial_subjectno}) - 0.0001;
        loc_F_common_max = str2double(input_for_common_max{trial_subjectno}) - 0.0001;
        loc_F_ind_rmax = str2double(input_for_ind_rmax{trial_subjectno}) - 0.0001;
        loc_F_ind_lmax = str2double(input_for_ind_lmax{trial_subjectno}) - 0.0001;

        % find goniometer angles (from all 3 scan locations / 6 trials, averaged)
        loc_frame = find(data_force_gonio(:,col_force)>=loc_F_trial_max,1,'first'); % when averaging angles across trials, intervals of 0.05 degrees are used. This means that any required angle will exist in all data series
        if isempty(loc_frame) % catch error of FORCE not found:
            loc_frame = find(data_force_gonio(:,col_force)>=max(data_force_gonio(:,col_force)),1,'first'); 
            cprintf('*red', 'ERROR: Recorded ind max FORCE not found in data. Run create_angles_passive?\n')
        end
        out_angle_trial_max = data_force_gonio(loc_frame,col_angle); 
        loc_frame = find(data_force_gonio(:,col_force)>=loc_F_ind_max,1,'first'); 
        if isempty(loc_frame) % catch error of FORCE not found:
            loc_frame = find(data_force_gonio(:,col_force)>=max(data_force_gonio(:,col_force)),1,'first'); 
            cprintf('*red', 'ERROR: Recorded ind max FORCE not found in data. Run create_angles_passive?\n')
        end
        out_angle_ind_max = data_force_gonio(loc_frame,col_angle);
        loc_frame = find(data_force_gonio(:,col_force)>=loc_F_common_max,1,'first');
        out_angle_common_max = data_force_gonio(loc_frame,col_angle);
        
        % for current trial (L or R), find the angle at the max force for
        % the other leg, if the other leg has a lower max force as lowest
        % PRE/POST (not really sure what these data should be used for...)
        if str2double(input_for_ind_rmax(trial_subjectno)) > 9000
            % data do not exist for the right side (array preloaded with "empty" values of 10000)
            out_angle_ind_rmax = 100;
        else
            loc_frame = find(data_force_gonio(:,col_force)>=loc_F_ind_rmax,1,'first'); 
            if isempty(loc_frame) % force level does not exist
                out_angle_ind_rmax = 100;
            else
                out_angle_ind_rmax = data_force_gonio(loc_frame,col_angle);
            end
        end
        
        if str2double(input_for_ind_lmax(trial_subjectno)) > 9000
            % data do not exist for the left side (array preloaded with "empty" values of 10000)
            out_angle_ind_lmax = 100;
        else
            loc_frame = find(data_force_gonio(:,col_force)>=loc_F_ind_lmax,1,'first'); 
            if isempty(loc_frame) % force level does not exist
                out_angle_ind_lmax = 100;
            else
                out_angle_ind_lmax = data_force_gonio(loc_frame,col_angle);
            end
        end
        %%
        
               
        
        %% prepare arrays for group plots 

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
        %   --- OR containing 3x zeros (if nonexistent)

        
        % extract angle range common to all trials for current subject
        angle_start = 0 - 0.00000001; % tweak for computer storage of numbers, where 8.3000 is stored as 8.2999999999999999 %VAR
        angle_stop = out_ROM_trial_max;

        % identify locations of start/stop angles in above mentioned arrays
        loc_angle_start = find(data_force_gonio(:,col_angle)>=angle_start,1,'first');
        loc_angle_stop = find(data_force_gonio(:,col_angle)>=angle_stop,1,'first');
        loc_angle_licht_start = find(data_GMFAS_licht_GM(:,col_lichtangle)>=angle_start,1,'first');
        loc_angle_licht_stop = find(data_GMFAS_licht_GM(:,col_lichtangle)>=angle_stop,1,'first');
        
        
        
        %%% contents of below angle_vars arrays:
                %   1 angle
                %   2 F 
                %   3 EMG_gm 
                %   4 EMG_gl 
                %   5 EMG_sol 

                %   6 elong AT
                %   7 elong GMtend
                %   8 elong MTU/calf        
                %   9 displ*** GMFAS     
                %  10 elong GMapo       
                %  11 elong msc GM         

                %  12 L_AT              
                %  13 L_GMtend          
                %  14 L_MTU/calf
                %  15 --- (GMFAS)
                %  16 L_GMapo           
                %  17 L_msc_GM             

                %  18 Torque            

                %  19 elong msc SOL
                %  20 L msc SOL

                % 21 strain AT
                % 22 strain GMtend
                % 23 strain MTJ
                % 24 ------ (GMFAS)
                % 25 strain GMapo
                % 26 strain msc GM
                
                % 27 strain msc SOL
        
                % 28 elong msc GM (from Lichtwark/Fukunaga)
                % 29 elong tend GM (from Lichtwark/Fukunaga)
                % 30 strain msc GM (from Lichtwark/Fukunaga)
                % 31 strain tend GM (from Lichtwark/Fukunaga)
                
                % 32 length fascicles GM (from Lichtwark)
                % 33 pennation GM (from Lichtwark)
                % 34 elongation fascicles GM (from Lichtwark)
                % 35 strain fascicles GM (from Lichtwark)
                
        if input_project == 1 % BD study
         %% BD study
            if trial_subjectno > 100 % BD subject
                %% BD subjects
                % all data in ONE cell, up to each subject's max angle, RAW data:
                BD_angle_vars{BD_count} = [ ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_angle) ...  1
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_force) ...  2
                    data_force_gonio(loc_angle_start:loc_angle_stop,3) ...          3
                    data_force_gonio(loc_angle_start:loc_angle_stop,4)...           4
                    data_force_gonio(loc_angle_start:loc_angle_stop,5) ...          5
                    MTU_elong_array(:,2) ...                                        6
                    MTU_elong_array(:,3) ...                                        7
                    MTU_elong_array(:,4) ...                                        8
                    MTU_elong_array(:,5) ...                                        9     GMFAS displacement
                    MTU_elong_array(:,6) ...                                        10
                    MTU_elong_array(:,7) ...                                        11
                    MTU_length_array(:,2) ...                                       12
                    MTU_length_array(:,3) ...                                       13
                    MTU_length_array(:,4) ...                                       14
                    MTU_length_array(:,5) ...                                       15     empty - GMFAS
                    MTU_length_array(:,6) ...                                       16
                    MTU_length_array(:,7) ...                                       17
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_force)*at_momentarm ... % 18 Torque
                    MTU_elong_array(:,8) ...                                        19
                    MTU_length_array(:,8) ...                                       20
                    MTU_strain_array(:,2) ...                                       21
                    MTU_strain_array(:,3) ...                                       22
                    MTU_strain_array(:,4) ...                                       23
                    MTU_strain_array(:,5) ...                                       24     empty - GMFAS
                    MTU_strain_array(:,6) ...                                       25
                    MTU_strain_array(:,7) ...                                       26
                    MTU_strain_array(:,8) ...                                       27
                    MTU_elong_array(:,9) ...                                        28
                    MTU_elong_array(:,10) ...                                       29
                    MTU_strain_array(:,9) ...                                       30
                    MTU_strain_array(:,10) ...                                      31
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas) ...   32
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtpenn) ...  33
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas)-data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas) ...   34
                    (data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas)-data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas))/data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas)*100 ...   35
                    ];

                % all data in ONE cell, NORMALIZED data:
                %    - EMG only normalized to %MVC, not to max EMG @ trial
                %    - elongations (6-11) not normalized to anything
                %    NB: in #12-17, normalizing to initial length in array, instead of using constants: str2double(dm_at_SOL_length{line}), dm_at_GM_length, (cm to mm) 10*dm_leg_length. Tested and values are the same except rounding at 2nd or 3rd decimal.
                % 
                % 2017-03-23: After strain was added to BD_angle_vars,
                % additional array with normalized data could be used only for
                % force, torque etc - strain could be taken from the end 
                % of non-normalized array ... but the script stays 
                % unchanged for simplicity

                
                BD_angle_vars_norm_indlength{BD_count} = [ ...
                    BD_angle_vars{1,BD_count}(:,1)*100/out_ROM_trial_max ...                            1 angle - normalized to trial max ROM 
                    BD_angle_vars{1,BD_count}(:,2)*100/max(BD_angle_vars{1,BD_count}(:,2)) ...          2 force - to maximal force in trial
                    BD_angle_vars{1,BD_count}(:,3) ...                                                  3 EMG - NOT normalized
                    BD_angle_vars{1,BD_count}(:,4) ...                                                  4 EMG - NOT normalized
                    BD_angle_vars{1,BD_count}(:,5) ...                                                  5 EMG - NOT normalized
                    BD_angle_vars{1,BD_count}(:,6) ...                                                  6 displ - NOT normalized
                    BD_angle_vars{1,BD_count}(:,7) ...                                                  7 displ - NOT normalized
                    BD_angle_vars{1,BD_count}(:,8) ...                                                  8 displ - NOT normalized
                    BD_angle_vars{1,BD_count}(:,9) ...                                                  9 displ - NOT normalized
                    BD_angle_vars{1,BD_count}(:,10) ...                                                 10 displ - NOT normalized
                    BD_angle_vars{1,BD_count}(:,11) ...                                                 11 displ - NOT normalized
                    (BD_angle_vars{1,BD_count}(:,12)-BD_angle_vars{1,BD_count}(1,12)) *100/BD_angle_vars{1,BD_count}(1,12) ...   12 length - to initial length of free AT
                    (BD_angle_vars{1,BD_count}(:,13)-BD_angle_vars{1,BD_count}(1,13)) *100/BD_angle_vars{1,BD_count}(1,13) ...    13 length - to initial length of GM tend
                    (BD_angle_vars{1,BD_count}(:,14)-BD_angle_vars{1,BD_count}(1,14)) *100/BD_angle_vars{1,BD_count}(1,14) ... 14 leg length - to initial leg length
                    BD_angle_vars{1,BD_count}(:,15) ...                                                 15 GMfas - no length exists
                    (BD_angle_vars{1,BD_count}(:,16)-BD_angle_vars{1,BD_count}(1,16)) *100/BD_angle_vars{1,BD_count}(1,16) ... 16 GM apo length - normalized to initial length of apo
                    (BD_angle_vars{1,BD_count}(:,17)-BD_angle_vars{1,BD_count}(1,17)) *100/BD_angle_vars{1,BD_count}(1,17) ... 17 GM msc length - normalized to initial msc length
                    BD_angle_vars{1,BD_count}(:,18)*100/max(BD_angle_vars{1,BD_count}(:,18)) ...        18 torque - to max torque in trial
                    BD_angle_vars{1,BD_count}(:,19) ...                                                 19 displ - NOT normalized
                    (BD_angle_vars{1,BD_count}(:,20)-BD_angle_vars{1,BD_count}(1,20)) *100/BD_angle_vars{1,BD_count}(1,20) ... 20 SOL msc length - normalized to initial msc length
                    BD_angle_vars{1,BD_count}(:,21) ...                                                 21 strain - NOT normalized
                    BD_angle_vars{1,BD_count}(:,22) ...                                                 22 strain - NOT normalized
                    BD_angle_vars{1,BD_count}(:,23) ...                                                 23 strain - NOT normalized
                    BD_angle_vars{1,BD_count}(:,24) ...                                                 24 strain - NOT normalized
                    BD_angle_vars{1,BD_count}(:,25) ...                                                 25 strain - NOT normalized
                    BD_angle_vars{1,BD_count}(:,26) ...                                                 26 strain - NOT normalized
                    BD_angle_vars{1,BD_count}(:,27) ...                                                 27 strain - NOT normalized
                    BD_angle_vars{1,BD_count}(:,28) ...                                                 28 NOT normalized
                    BD_angle_vars{1,BD_count}(:,29) ...                                                 29 NOT normalized
                    BD_angle_vars{1,BD_count}(:,30) ...                                                 30 NOT normalized
                    BD_angle_vars{1,BD_count}(:,31) ...                                                 31 NOT normalized
                    BD_angle_vars{1,BD_count}(:,32) ...                                                 32 GM fas len (Licht) - norm to initial leg length --- TODO for all?
                    BD_angle_vars{1,BD_count}(:,33) ...                                                 33 NOT normalized
                    BD_angle_vars{1,BD_count}(:,34) ...                                                 34 NOT normalized
                    BD_angle_vars{1,BD_count}(:,35) ...                                                 35 NOT normalized
                    ];

                % reshape
                BD_angle_vars_norm{BD_count} = [ (0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,2), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,3), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,4), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,5), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,6), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,7), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,8), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,9), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,10), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,11), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,12), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,13), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,14), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,15), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,16), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,17), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,18), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,19), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,20), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,21), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,22), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,23), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,24), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,25), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,26), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,27), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,28), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,29), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,30), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,31), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,32), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,33), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,34), 0:angle_step:100)', ...
                    spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,35), 0:angle_step:100)', ...
                    ];
                %% END BD subjects
            else
                %% CON subjects

                % all data in ONE cell, common angles, RAW data:
                CON_angle_vars{CON_count} = [ ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_angle) ...  1
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_force) ...  2
                    data_force_gonio(loc_angle_start:loc_angle_stop,3) ...          3
                    data_force_gonio(loc_angle_start:loc_angle_stop,4)...           4
                    data_force_gonio(loc_angle_start:loc_angle_stop,5) ...          5
                    MTU_elong_array(:,2) ...                                        6
                    MTU_elong_array(:,3) ...                                        7
                    MTU_elong_array(:,4) ...                                        8
                    MTU_elong_array(:,5) ...                                        9
                    MTU_elong_array(:,6) ...                                        10
                    MTU_elong_array(:,7) ...                                        11
                    MTU_length_array(:,2) ...                                       12
                    MTU_length_array(:,3) ...                                       13
                    MTU_length_array(:,4) ...                                       14
                    MTU_length_array(:,5) ...                                       15
                    MTU_length_array(:,6) ...                                       16
                    MTU_length_array(:,7) ...                                       17
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_force)*at_momentarm ... % 18
                    MTU_elong_array(:,8) ...                                        19
                    MTU_length_array(:,8) ...                                       20
                    MTU_strain_array(:,2) ...                                       21
                    MTU_strain_array(:,3) ...                                       22
                    MTU_strain_array(:,4) ...                                       23
                    MTU_strain_array(:,5) ...                                       24
                    MTU_strain_array(:,6) ...                                       25
                    MTU_strain_array(:,7) ...                                       26
                    MTU_strain_array(:,8) ...                                       27
                    MTU_elong_array(:,9) ...                                        28
                    MTU_elong_array(:,10) ...                                       29
                    MTU_strain_array(:,9) ...                                       30
                    MTU_strain_array(:,10) ...                                      31
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas) ...   32
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtpenn) ...  33
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas)-data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas) ...   34
                    (data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas)-data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas))/data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas)*100 ...   35
                    ];

                % all data in ONE cell, NORMALIZED data:
                %    - EMG only normalized to %MVC, not to max EMG @ trial
                %    - elongations (6-11) not normalized to anything
                %    NB: in #12-17, normalizing to initial length in array, instead of using constants: str2double(dm_at_SOL_length{line}), dm_at_GM_length, (cm to mm) 10*dm_leg_length. Tested and values are the same except rounding at 2nd or 3rd decimal.
                CON_angle_vars_norm_indlength{CON_count} = [ ...
                    CON_angle_vars{1,CON_count}(:,1)*100/out_ROM_trial_max ...                            1 angle - normalized to trial max ROM
                    CON_angle_vars{1,CON_count}(:,2)*100/max(CON_angle_vars{1,CON_count}(:,2)) ...          2 force - to maximal force in trial
                    CON_angle_vars{1,CON_count}(:,3) ...                                                  3 EMG - NOT normalized
                    CON_angle_vars{1,CON_count}(:,4) ...                                                  4 EMG - NOT normalized
                    CON_angle_vars{1,CON_count}(:,5) ...                                                  5 EMG - NOT normalized
                    CON_angle_vars{1,CON_count}(:,6) ...                                                  6 displ - NOT normalized
                    CON_angle_vars{1,CON_count}(:,7) ...                                                  7 displ - NOT normalized
                    CON_angle_vars{1,CON_count}(:,8) ...                                                  8 displ - NOT normalized
                    CON_angle_vars{1,CON_count}(:,9) ...                                                  9 displ - NOT normalized
                    CON_angle_vars{1,CON_count}(:,10) ...                                                 10 displ - NOT normalized
                    CON_angle_vars{1,CON_count}(:,11) ...                                                 11 displ - NOT normalized
                    (CON_angle_vars{1,CON_count}(:,12)-CON_angle_vars{1,CON_count}(1,12)) *100/CON_angle_vars{1,CON_count}(1,12) ...   12 length - to initial length of free AT
                    (CON_angle_vars{1,CON_count}(:,13)-CON_angle_vars{1,CON_count}(1,13)) *100/CON_angle_vars{1,CON_count}(1,13) ...    13 length - to initial length of GM tend
                    (CON_angle_vars{1,CON_count}(:,14)-CON_angle_vars{1,CON_count}(1,14)) *100/CON_angle_vars{1,CON_count}(1,14) ... 14 leg length - to initial leg length
                    CON_angle_vars{1,CON_count}(:,15) ...                                                 15 GMfas - no length exists
                    (CON_angle_vars{1,CON_count}(:,16)-CON_angle_vars{1,CON_count}(1,16)) *100/CON_angle_vars{1,CON_count}(1,16) ... 16 GM apo length - normalized to initial length of apo
                    (CON_angle_vars{1,CON_count}(:,17)-CON_angle_vars{1,CON_count}(1,17)) *100/CON_angle_vars{1,CON_count}(1,17) ... 17 GM msc length - normalized to initial msc length
                    CON_angle_vars{1,CON_count}(:,18)*100/max(CON_angle_vars{1,CON_count}(:,18)) ...       18 torque - to max torque in trial
                    CON_angle_vars{1,CON_count}(:,19) ...                                                 19 displ - NOT normalized
                    (CON_angle_vars{1,CON_count}(:,20)-CON_angle_vars{1,CON_count}(1,20)) *100/CON_angle_vars{1,CON_count}(1,20) ... 20 SOL msc length - normalized to initial msc length
                    CON_angle_vars{1,CON_count}(:,21) ...                                                 21 strain - NOT normalized
                    CON_angle_vars{1,CON_count}(:,22) ...                                                 22 strain - NOT normalized
                    CON_angle_vars{1,CON_count}(:,23) ...                                                 23 strain - NOT normalized
                    CON_angle_vars{1,CON_count}(:,24) ...                                                 24 strain - NOT normalized
                    CON_angle_vars{1,CON_count}(:,25) ...                                                 25 strain - NOT normalized
                    CON_angle_vars{1,CON_count}(:,26) ...                                                 26 strain - NOT normalized
                    CON_angle_vars{1,CON_count}(:,27) ...                                                 27 strain - NOT normalized
                    CON_angle_vars{1,CON_count}(:,28) ...                                                 28 NOT normalized
                    CON_angle_vars{1,CON_count}(:,29) ...                                                 29 NOT normalized
                    CON_angle_vars{1,CON_count}(:,30) ...                                                 30 NOT normalized
                    CON_angle_vars{1,CON_count}(:,31) ...                                                 31 NOT normalized
                    CON_angle_vars{1,CON_count}(:,32) ...                                                 32 NOT normalized
                    CON_angle_vars{1,CON_count}(:,33) ...                                                 33 NOT normalized
                    CON_angle_vars{1,CON_count}(:,34) ...                                                 34 NOT normalized
                    CON_angle_vars{1,CON_count}(:,35) ...                                                 35 NOT normalized
                    ];

                % reshape
                CON_angle_vars_norm{CON_count} = [ (0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,2), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,3), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,4), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,5), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,6), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,7), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,8), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,9), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,10), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,11), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,12), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,13), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,14), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,15), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,16), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,17), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,18), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,19), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,20), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,21), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,22), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,23), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,24), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,25), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,26), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,27), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,28), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,29), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,30), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,31), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,32), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,33), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,34), 0:angle_step:100)', ...
                    spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,35), 0:angle_step:100)', ...
                    ];
                %% END CON subjects
            end
         %%
        else % project == 2 == intervention
         %% intervention study
            if trial_timepoint == 0 && trial_leg == 0 % PRE, CON
                %% PRE CON
                % all data in ONE cell, up to each subject's max angle, RAW data:
                CON_PRE_angle_vars{CON_PRE_count} = [ ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_angle) ...  1
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_force) ...  2
                    data_force_gonio(loc_angle_start:loc_angle_stop,3) ...          3
                    data_force_gonio(loc_angle_start:loc_angle_stop,4)...           4
                    data_force_gonio(loc_angle_start:loc_angle_stop,5) ...          5
                    MTU_elong_array(:,2) ...                                        6
                    MTU_elong_array(:,3) ...                                        7
                    MTU_elong_array(:,4) ...                                        8
                    MTU_elong_array(:,5) ...                                        9
                    MTU_elong_array(:,6) ...                                        10
                    MTU_elong_array(:,7) ...                                        11
                    MTU_length_array(:,2) ...                                       12
                    MTU_length_array(:,3) ...                                       13
                    MTU_length_array(:,4) ...                                       14
                    MTU_length_array(:,5) ...                                       15
                    MTU_length_array(:,6) ...                                       16
                    MTU_length_array(:,7) ...                                       17
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_force)*at_momentarm ... % 18
                    MTU_elong_array(:,8) ...                                        19
                    MTU_length_array(:,8) ...                                       20
                    MTU_strain_array(:,2) ...                                       21
                    MTU_strain_array(:,3) ...                                       22
                    MTU_strain_array(:,4) ...                                       23
                    MTU_strain_array(:,5) ...                                       24
                    MTU_strain_array(:,6) ...                                       25
                    MTU_strain_array(:,7) ...                                       26
                    MTU_strain_array(:,8) ...                                       27
                    MTU_elong_array(:,9) ...                                        28
                    MTU_elong_array(:,10) ...                                       29
                    MTU_strain_array(:,9) ...                                       30
                    MTU_strain_array(:,10) ...                                      31
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas) ...   32
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtpenn) ...  33
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas)-data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas) ...   34
                    (data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas)-data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas))/data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas)*100 ...   35
                    ];

                % all data in ONE cell, NORMALIZED data:
                %    - EMG only normalized to %MVC, not to max EMG @ trial
                %    - elongations (6-11) not normalized to anything
                %    NB: in #12-17, normalizing to initial length in array, instead of using constants: str2double(dm_at_SOL_length{line}), dm_at_GM_length, (cm to mm) 10*dm_leg_length. Tested and values are the same except rounding at 2nd or 3rd decimal.
                CON_PRE_angle_vars_norm_indlength{CON_PRE_count} = [ ...
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,1)*100/out_ROM_trial_max ...                            1 angle - normalized to trial max ROM 
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,2)*100/max(CON_PRE_angle_vars{1,CON_PRE_count}(:,2)) ...          2 force - to maximal force in trial
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,3) ...                                                  3 EMG - NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,4) ...                                                  4 EMG - NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,5) ...                                                  5 EMG - NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,6) ...                                                  6 displ - NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,7) ...                                                  7 displ - NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,8) ...                                                  8 displ - NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,9) ...                                                  9 displ - NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,10) ...                                                 10 displ - NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,11) ...                                                 11 displ - NOT normalized
                    (CON_PRE_angle_vars{1,CON_PRE_count}(:,12)-CON_PRE_angle_vars{1,CON_PRE_count}(1,12)) *100/CON_PRE_angle_vars{1,CON_PRE_count}(1,12) ...   12 length - to initial length of free AT
                    (CON_PRE_angle_vars{1,CON_PRE_count}(:,13)-CON_PRE_angle_vars{1,CON_PRE_count}(1,13)) *100/CON_PRE_angle_vars{1,CON_PRE_count}(1,13) ...    13 length - to initial length of GM tend
                    (CON_PRE_angle_vars{1,CON_PRE_count}(:,14)-CON_PRE_angle_vars{1,CON_PRE_count}(1,14)) *100/CON_PRE_angle_vars{1,CON_PRE_count}(1,14) ... 14 leg length - to initial leg length
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,15) ...                                                 15 GMfas - no length exists
                    (CON_PRE_angle_vars{1,CON_PRE_count}(:,16)-CON_PRE_angle_vars{1,CON_PRE_count}(1,16)) *100/CON_PRE_angle_vars{1,CON_PRE_count}(1,16) ... 16 GM apo length - normalized to initial length of apo
                    (CON_PRE_angle_vars{1,CON_PRE_count}(:,17)-CON_PRE_angle_vars{1,CON_PRE_count}(1,17)) *100/CON_PRE_angle_vars{1,CON_PRE_count}(1,17) ... 17 GM msc length - normalized to initial msc length
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,18)*100/max(CON_PRE_angle_vars{1,CON_PRE_count}(:,18)) ...        18 torque - to max torque in trial
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,19) ...                                                 19 displ - NOT normalized
                    (CON_PRE_angle_vars{1,CON_PRE_count}(:,20)-CON_PRE_angle_vars{1,CON_PRE_count}(1,20)) *100/CON_PRE_angle_vars{1,CON_PRE_count}(1,20) ... 20 SOL msc length - normalized to initial msc length
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,21) ...                                                 21 strain - NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,22) ...                                                 22 strain - NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,23) ...                                                 23 strain - NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,24) ...                                                 24 strain - NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,25) ...                                                 25 strain - NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,26) ...                                                 26 strain - NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,27) ...                                                 27 strain - NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,28) ...                                                 28 NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,29) ...                                                 29 NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,30) ...                                                 30 NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,31) ...                                                 31 NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,32) ...                                                 32 NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,33) ...                                                 33 NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,34) ...                                                 34 NOT normalized
                    CON_PRE_angle_vars{1,CON_PRE_count}(:,35) ...                                                 35 NOT normalized
                    ];

                % reshape
                CON_PRE_angle_vars_norm{CON_PRE_count} = [ (0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,2), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,3), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,4), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,5), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,6), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,7), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,8), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,9), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,10), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,11), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,12), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,13), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,14), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,15), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,16), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,17), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,18), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,19), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,20), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,21), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,22), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,23), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,24), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,25), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,26), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,27), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,28), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,29), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,30), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,31), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,32), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,33), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,34), 0:angle_step:100)', ...
                    spline(CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,1), CON_PRE_angle_vars_norm_indlength{1,CON_PRE_count}(:,35), 0:angle_step:100)', ...
                    ];
                %% END PRE CON
            elseif trial_timepoint == 0 && trial_leg == 1 % PRE, STR
                %% PRE STR
                % all data in ONE cell, common angles, RAW data:
                STR_PRE_angle_vars{STR_PRE_count} = [ ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_angle) ...  1
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_force) ...  2
                    data_force_gonio(loc_angle_start:loc_angle_stop,3) ...          3
                    data_force_gonio(loc_angle_start:loc_angle_stop,4)...           4
                    data_force_gonio(loc_angle_start:loc_angle_stop,5) ...          5
                    MTU_elong_array(:,2) ...                                        6
                    MTU_elong_array(:,3) ...                                        7
                    MTU_elong_array(:,4) ...                                        8
                    MTU_elong_array(:,5) ...                                        9
                    MTU_elong_array(:,6) ...                                        10
                    MTU_elong_array(:,7) ...                                        11
                    MTU_length_array(:,2) ...                                       12
                    MTU_length_array(:,3) ...                                       13
                    MTU_length_array(:,4) ...                                       14
                    MTU_length_array(:,5) ...                                       15
                    MTU_length_array(:,6) ...                                       16
                    MTU_length_array(:,7) ...                                       17
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_force)*at_momentarm ... % 18
                    MTU_elong_array(:,8) ...                                        19
                    MTU_length_array(:,8) ...                                       20
                    MTU_strain_array(:,2) ...                                       21
                    MTU_strain_array(:,3) ...                                       22
                    MTU_strain_array(:,4) ...                                       23
                    MTU_strain_array(:,5) ...                                       24
                    MTU_strain_array(:,6) ...                                       25
                    MTU_strain_array(:,7) ...                                       26
                    MTU_strain_array(:,8) ...                                       27
                    MTU_elong_array(:,9) ...                                        28
                    MTU_elong_array(:,10) ...                                       29
                    MTU_strain_array(:,9) ...                                       30
                    MTU_strain_array(:,10) ...                                      31
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas) ...   32
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtpenn) ...  33
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas)-data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas) ...   34
                    (data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas)-data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas))/data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas)*100 ...   35
                    ];

                % all data in ONE cell, NORMALIZED data:
                %    - EMG only normalized to %MVC, not to max EMG @ trial
                %    - elongations (6-11) not normalized to anything
                %    NB: in #12-17, normalizing to initial length in array, instead of using constants: str2double(dm_at_SOL_length{line}), dm_at_GM_length, (cm to mm) 10*dm_leg_length. Tested and values are the same except rounding at 2nd or 3rd decimal.
                STR_PRE_angle_vars_norm_indlength{STR_PRE_count} = [ ...
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,1)*100/out_ROM_trial_max ...                            1 angle - normalized to trial max ROM
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,2)*100/max(STR_PRE_angle_vars{1,STR_PRE_count}(:,2)) ...          2 force - to maximal force in trial
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,3) ...                                                  3 EMG - NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,4) ...                                                  4 EMG - NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,5) ...                                                  5 EMG - NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,6) ...                                                  6 displ - NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,7) ...                                                  7 displ - NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,8) ...                                                  8 displ - NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,9) ...                                                  9 displ - NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,10) ...                                                 10 displ - NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,11) ...                                                 11 displ - NOT normalized
                    (STR_PRE_angle_vars{1,STR_PRE_count}(:,12)-STR_PRE_angle_vars{1,STR_PRE_count}(1,12)) *100/STR_PRE_angle_vars{1,STR_PRE_count}(1,12) ...   12 length - to initial length of free AT
                    (STR_PRE_angle_vars{1,STR_PRE_count}(:,13)-STR_PRE_angle_vars{1,STR_PRE_count}(1,13)) *100/STR_PRE_angle_vars{1,STR_PRE_count}(1,13) ...    13 length - to initial length of GM tend
                    (STR_PRE_angle_vars{1,STR_PRE_count}(:,14)-STR_PRE_angle_vars{1,STR_PRE_count}(1,14)) *100/STR_PRE_angle_vars{1,STR_PRE_count}(1,14) ... 14 leg length - to initial leg length
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,15) ...                                                 15 GMfas - no length exists
                    (STR_PRE_angle_vars{1,STR_PRE_count}(:,16)-STR_PRE_angle_vars{1,STR_PRE_count}(1,16)) *100/STR_PRE_angle_vars{1,STR_PRE_count}(1,16) ... 16 GM apo length - normalized to initial length of apo
                    (STR_PRE_angle_vars{1,STR_PRE_count}(:,17)-STR_PRE_angle_vars{1,STR_PRE_count}(1,17)) *100/STR_PRE_angle_vars{1,STR_PRE_count}(1,17) ... 17 GM msc length - normalized to initial msc length
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,18)*100/max(STR_PRE_angle_vars{1,STR_PRE_count}(:,18)) ...       18 torque - to max torque in trial
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,19) ...                                                 19 displ - NOT normalized
                    (STR_PRE_angle_vars{1,STR_PRE_count}(:,20)-STR_PRE_angle_vars{1,STR_PRE_count}(1,20)) *100/STR_PRE_angle_vars{1,STR_PRE_count}(1,20) ... 20 SOL msc length - normalized to initial msc length
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,21) ...                                                 21 strain - NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,22) ...                                                 22 strain - NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,23) ...                                                 23 strain - NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,24) ...                                                 24 strain - NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,25) ...                                                 25 strain - NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,26) ...                                                 26 strain - NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,27) ...                                                 27 strain - NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,28) ...                                                 28 NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,29) ...                                                 29 NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,30) ...                                                 30 NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,31) ...                                                 31 NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,32) ...                                                 32 NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,33) ...                                                 33 NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,34) ...                                                 34 NOT normalized
                    STR_PRE_angle_vars{1,STR_PRE_count}(:,35) ...                                                 35 NOT normalized
                    ];

                % reshape
                STR_PRE_angle_vars_norm{STR_PRE_count} = [ (0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,2), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,3), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,4), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,5), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,6), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,7), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,8), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,9), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,10), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,11), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,12), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,13), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,14), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,15), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,16), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,17), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,18), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,19), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,20), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,21), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,22), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,23), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,24), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,25), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,26), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,27), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,28), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,29), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,30), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,31), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,32), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,33), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,34), 0:angle_step:100)', ...
                    spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1), STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,35), 0:angle_step:100)', ...
                    ];
                %% END PRE STR
            elseif trial_timepoint == 1 && trial_leg == 0 % POST, CON
                %% POST CON
                % all data in ONE cell, up to each subject's max angle, RAW data:
                CON_POST_angle_vars{CON_POST_count} = [ ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_angle) ...  1
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_force) ...  2
                    data_force_gonio(loc_angle_start:loc_angle_stop,3) ...          3
                    data_force_gonio(loc_angle_start:loc_angle_stop,4)...           4
                    data_force_gonio(loc_angle_start:loc_angle_stop,5) ...          5
                    MTU_elong_array(:,2) ...                                        6
                    MTU_elong_array(:,3) ...                                        7
                    MTU_elong_array(:,4) ...                                        8
                    MTU_elong_array(:,5) ...                                        9
                    MTU_elong_array(:,6) ...                                        10
                    MTU_elong_array(:,7) ...                                        11
                    MTU_length_array(:,2) ...                                       12
                    MTU_length_array(:,3) ...                                       13
                    MTU_length_array(:,4) ...                                       14
                    MTU_length_array(:,5) ...                                       15
                    MTU_length_array(:,6) ...                                       16
                    MTU_length_array(:,7) ...                                       17
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_force)*at_momentarm ... % 18
                    MTU_elong_array(:,8) ...                                        19
                    MTU_length_array(:,8) ...                                       20
                    MTU_strain_array(:,2) ...                                       21
                    MTU_strain_array(:,3) ...                                       22
                    MTU_strain_array(:,4) ...                                       23
                    MTU_strain_array(:,5) ...                                       24
                    MTU_strain_array(:,6) ...                                       25
                    MTU_strain_array(:,7) ...                                       26
                    MTU_strain_array(:,8) ...                                       27
                    MTU_elong_array(:,9) ...                                        28
                    MTU_elong_array(:,10) ...                                       29
                    MTU_strain_array(:,9) ...                                       30
                    MTU_strain_array(:,10) ...                                      31
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas) ...   32
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtpenn) ...  33
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas)-data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas) ...   34
                    (data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas)-data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas))/data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas)*100 ...   35
                    ];

                % all data in ONE cell, NORMALIZED data:
                %    - EMG only normalized to %MVC, not to max EMG @ trial
                %    - elongations (6-11) not normalized to anything
                %    NB: in #12-17, normalizing to initial length in array, instead of using constants: str2double(dm_at_SOL_length{line}), dm_at_GM_length, (cm to mm) 10*dm_leg_length. Tested and values are the same except rounding at 2nd or 3rd decimal.
                CON_POST_angle_vars_norm_indlength{CON_POST_count} = [ ...
                    CON_POST_angle_vars{1,CON_POST_count}(:,1)*100/out_ROM_trial_max ...                            1 angle - normalized to trial max ROM 
                    CON_POST_angle_vars{1,CON_POST_count}(:,2)*100/max(CON_POST_angle_vars{1,CON_POST_count}(:,2)) ...          2 force - to maximal force in trial
                    CON_POST_angle_vars{1,CON_POST_count}(:,3) ...                                                  3 EMG - NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,4) ...                                                  4 EMG - NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,5) ...                                                  5 EMG - NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,6) ...                                                  6 displ - NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,7) ...                                                  7 displ - NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,8) ...                                                  8 displ - NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,9) ...                                                  9 displ - NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,10) ...                                                 10 displ - NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,11) ...                                                 11 displ - NOT normalized
                    (CON_POST_angle_vars{1,CON_POST_count}(:,12)-CON_POST_angle_vars{1,CON_POST_count}(1,12)) *100/CON_POST_angle_vars{1,CON_POST_count}(1,12) ...   12 length - to initial length of free AT
                    (CON_POST_angle_vars{1,CON_POST_count}(:,13)-CON_POST_angle_vars{1,CON_POST_count}(1,13)) *100/CON_POST_angle_vars{1,CON_POST_count}(1,13) ...    13 length - to initial length of GM tend
                    (CON_POST_angle_vars{1,CON_POST_count}(:,14)-CON_POST_angle_vars{1,CON_POST_count}(1,14)) *100/CON_POST_angle_vars{1,CON_POST_count}(1,14) ... 14 leg length - to initial leg length
                    CON_POST_angle_vars{1,CON_POST_count}(:,15) ...                                                 15 GMfas - no length exists
                    (CON_POST_angle_vars{1,CON_POST_count}(:,16)-CON_POST_angle_vars{1,CON_POST_count}(1,16)) *100/CON_POST_angle_vars{1,CON_POST_count}(1,16) ... 16 GM apo length - normalized to initial length of apo
                    (CON_POST_angle_vars{1,CON_POST_count}(:,17)-CON_POST_angle_vars{1,CON_POST_count}(1,17)) *100/CON_POST_angle_vars{1,CON_POST_count}(1,17) ... 17 GM msc length - normalized to initial msc length
                    CON_POST_angle_vars{1,CON_POST_count}(:,18)*100/max(CON_POST_angle_vars{1,CON_POST_count}(:,18)) ...        18 torque - to max torque in trial
                    CON_POST_angle_vars{1,CON_POST_count}(:,19) ...                                                 19 displ - NOT normalized
                    (CON_POST_angle_vars{1,CON_POST_count}(:,20)-CON_POST_angle_vars{1,CON_POST_count}(1,20)) *100/CON_POST_angle_vars{1,CON_POST_count}(1,20) ... 20 SOL msc length - normalized to initial msc length
                    CON_POST_angle_vars{1,CON_POST_count}(:,21) ...                                                 21 strain - NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,22) ...                                                 22 strain - NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,23) ...                                                 23 strain - NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,24) ...                                                 24 strain - NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,25) ...                                                 25 strain - NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,26) ...                                                 26 strain - NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,27) ...                                                 27 strain - NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,28) ...                                                 28 NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,29) ...                                                 29 NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,30) ...                                                 30 NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,31) ...                                                 31 NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,32) ...                                                 32 NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,33) ...                                                 33 NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,34) ...                                                 34 NOT normalized
                    CON_POST_angle_vars{1,CON_POST_count}(:,35) ...                                                 35 NOT normalized
                    ];

                % reshape
                CON_POST_angle_vars_norm{CON_POST_count} = [ (0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,2), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,3), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,4), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,5), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,6), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,7), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,8), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,9), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,10), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,11), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,12), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,13), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,14), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,15), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,16), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,17), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,18), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,19), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,20), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,21), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,22), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,23), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,24), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,25), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,26), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,27), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,28), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,29), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,30), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,31), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,32), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,33), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,34), 0:angle_step:100)', ...
                    spline(CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,1), CON_POST_angle_vars_norm_indlength{1,CON_POST_count}(:,35), 0:angle_step:100)', ...
                    ];
                %% END POST CON
            elseif trial_timepoint == 1 && trial_leg == 1 % POST, STR
                %% POST STR
                % all data in ONE cell, up to each subject's max angle, RAW data:
                STR_POST_angle_vars{STR_POST_count} = [ ...
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_angle) ...  1
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_force) ...  2
                    data_force_gonio(loc_angle_start:loc_angle_stop,3) ...          3
                    data_force_gonio(loc_angle_start:loc_angle_stop,4)...           4
                    data_force_gonio(loc_angle_start:loc_angle_stop,5) ...          5
                    MTU_elong_array(:,2) ...                                        6
                    MTU_elong_array(:,3) ...                                        7
                    MTU_elong_array(:,4) ...                                        8
                    MTU_elong_array(:,5) ...                                        9
                    MTU_elong_array(:,6) ...                                        10
                    MTU_elong_array(:,7) ...                                        11
                    MTU_length_array(:,2) ...                                       12
                    MTU_length_array(:,3) ...                                       13
                    MTU_length_array(:,4) ...                                       14
                    MTU_length_array(:,5) ...                                       15
                    MTU_length_array(:,6) ...                                       16
                    MTU_length_array(:,7) ...                                       17
                    data_force_gonio(loc_angle_start:loc_angle_stop,col_force)*at_momentarm ... % 18
                    MTU_elong_array(:,8) ...                                        19
                    MTU_length_array(:,8) ...                                       20
                    MTU_strain_array(:,2) ...                                       21
                    MTU_strain_array(:,3) ...                                       22
                    MTU_strain_array(:,4) ...                                       23
                    MTU_strain_array(:,5) ...                                       24
                    MTU_strain_array(:,6) ...                                       25
                    MTU_strain_array(:,7) ...                                       26
                    MTU_strain_array(:,8) ...                                       27
                    MTU_elong_array(:,9) ...                                        28
                    MTU_elong_array(:,10) ...                                       29
                    MTU_strain_array(:,9) ...                                       30
                    MTU_strain_array(:,10) ...                                      31
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas) ...   32
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtpenn) ...  33
                    data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas)-data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas) ...   34
                    (data_GMFAS_licht_GM(loc_angle_licht_start:loc_angle_licht_stop,col_lichtfas)-data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas))/data_GMFAS_licht_GM(loc_angle_licht_start,col_lichtfas)*100 ...   35
                    ];

                % all data in ONE cell, NORMALIZED data:
                %    - EMG only normalized to %MVC, not to max EMG @ trial
                %    - elongations (6-11) not normalized to anything
                %    NB: in #12-17, normalizing to initial length in array, instead of using constants: str2double(dm_at_SOL_length{line}), dm_at_GM_length, (cm to mm) 10*dm_leg_length. Tested and values are the same except rounding at 2nd or 3rd decimal.
                STR_POST_angle_vars_norm_indlength{STR_POST_count} = [ ...
                    STR_POST_angle_vars{1,STR_POST_count}(:,1)*100/out_ROM_trial_max ...                            1 angle - normalized to trial max ROM 
                    STR_POST_angle_vars{1,STR_POST_count}(:,2)*100/max(STR_POST_angle_vars{1,STR_POST_count}(:,2)) ...          2 force - to maximal force in trial
                    STR_POST_angle_vars{1,STR_POST_count}(:,3) ...                                                  3 EMG - NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,4) ...                                                  4 EMG - NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,5) ...                                                  5 EMG - NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,6) ...                                                  6 displ - NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,7) ...                                                  7 displ - NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,8) ...                                                  8 displ - NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,9) ...                                                  9 displ - NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,10) ...                                                 10 displ - NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,11) ...                                                 11 displ - NOT normalized
                    (STR_POST_angle_vars{1,STR_POST_count}(:,12)-STR_POST_angle_vars{1,STR_POST_count}(1,12)) *100/STR_POST_angle_vars{1,STR_POST_count}(1,12) ...   12 length - to initial length of free AT
                    (STR_POST_angle_vars{1,STR_POST_count}(:,13)-STR_POST_angle_vars{1,STR_POST_count}(1,13)) *100/STR_POST_angle_vars{1,STR_POST_count}(1,13) ...    13 length - to initial length of GM tend
                    (STR_POST_angle_vars{1,STR_POST_count}(:,14)-STR_POST_angle_vars{1,STR_POST_count}(1,14)) *100/STR_POST_angle_vars{1,STR_POST_count}(1,14) ... 14 leg length - to initial leg length
                    STR_POST_angle_vars{1,STR_POST_count}(:,15) ...                                                 15 GMfas - no length exists
                    (STR_POST_angle_vars{1,STR_POST_count}(:,16)-STR_POST_angle_vars{1,STR_POST_count}(1,16)) *100/STR_POST_angle_vars{1,STR_POST_count}(1,16) ... 16 GM apo length - normalized to initial length of apo
                    (STR_POST_angle_vars{1,STR_POST_count}(:,17)-STR_POST_angle_vars{1,STR_POST_count}(1,17)) *100/STR_POST_angle_vars{1,STR_POST_count}(1,17) ... 17 GM msc length - normalized to initial msc length
                    STR_POST_angle_vars{1,STR_POST_count}(:,18)*100/max(STR_POST_angle_vars{1,STR_POST_count}(:,18)) ...        18 torque - to max torque in trial
                    STR_POST_angle_vars{1,STR_POST_count}(:,19) ...                                                 19 displ - NOT normalized
                    (STR_POST_angle_vars{1,STR_POST_count}(:,20)-STR_POST_angle_vars{1,STR_POST_count}(1,20)) *100/STR_POST_angle_vars{1,STR_POST_count}(1,20) ... 20 SOL msc length - normalized to initial msc length
                    STR_POST_angle_vars{1,STR_POST_count}(:,21) ...                                                 21 strain - NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,22) ...                                                 22 strain - NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,23) ...                                                 23 strain - NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,24) ...                                                 24 strain - NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,25) ...                                                 25 strain - NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,26) ...                                                 26 strain - NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,27) ...                                                 27 strain - NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,28) ...                                                 28 NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,29) ...                                                 29 NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,30) ...                                                 30 NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,31) ...                                                 31 NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,32) ...                                                 32 NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,33) ...                                                 33 NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,34) ...                                                 34 NOT normalized
                    STR_POST_angle_vars{1,STR_POST_count}(:,35) ...                                                 35 NOT normalized
                    ];

                % reshape
                STR_POST_angle_vars_norm{STR_POST_count} = [ (0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,2), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,3), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,4), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,5), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,6), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,7), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,8), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,9), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,10), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,11), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,12), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,13), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,14), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,15), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,16), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,17), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,18), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,19), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,20), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,21), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,22), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,23), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,24), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,25), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,26), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,27), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,28), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,29), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,30), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,31), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,32), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,33), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,34), 0:angle_step:100)', ...
                    spline(STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,1), STR_POST_angle_vars_norm_indlength{1,STR_POST_count}(:,35), 0:angle_step:100)', ...
                    ];
                %% END POST STR
            end
         %% 
        end
        %%
        
        
        
        %% prepare arrays for individual trial data to file

        % txt trial ID
        all_passive_output_txt(line,:) = [dm_subjectno(line) dm_timepoint(line) dm_side(line) dm_trial(line)];
        
        % add data to a common array for all subjects    
        all_passive_output(line,:) = [out_ROM_trial_max out_ROM_ind_max out_ROM_common_max out_ROM_submax_1 out_ROM_submax_2...
        	out_F_trial_max_ROM out_F_trial_max_F out_F_ind_max out_F_common_max out_F_zero out_F_submax_1 out_F_submax_2...
        	out_T_trial_max_ROM out_T_trial_max_F out_T_ind_max out_T_common_max out_T_zero out_T_submax_1 out_T_submax_2...
        	out_angle_trial_max out_angle_ind_max out_angle_common_max out_angle_ind_rmax out_angle_ind_lmax...
        	out_pstiff_trial_max out_pstiff_ind_max out_pstiff_common_max out_pstiff_15 out_pstiff_submax_1 out_pstiff_submax_2 out_pstiff_index...
        	out_displ_SOL_trial_max out_displ_SOL_ind_max out_displ_SOL_common_max out_displ_SOL_submax_1 out_displ_SOL_submax_2...
            out_displ_GMMTJ_trial_max out_displ_GMMTJ_ind_max out_displ_GMMTJ_common_max out_displ_GMMTJ_submax_1 out_displ_GMMTJ_submax_2...
            out_displ_GMFAS_trial_max out_displ_GMFAS_ind_max out_displ_GMFAS_common_max out_displ_GMFAS_submax_1 out_displ_GMFAS_submax_2...
        	out_length_AT_trial_max out_length_GMtend_trial_max out_length_leg_trial_max out_length_GMapo_trial_max out_length_msc_GM_trial_max out_length_msc_SOL_trial_max...
        	out_length_AT_ind_max out_length_GMtend_ind_max out_length_leg_ind_max out_length_GMapo_ind_max out_length_msc_GM_ind_max out_length_msc_SOL_ind_max...
        	out_length_AT_common_max out_length_GMtend_common_max out_length_leg_common_max out_length_GMapo_common_max out_length_msc_GM_common_max out_length_msc_SOL_common_max...
        	out_length_AT_submax_1 out_length_GMtend_submax_1 out_length_leg_submax_1 out_length_GMapo_submax_1 out_length_msc_GM_submax_1 out_length_msc_SOL_submax_1...
        	out_length_AT_submax_2 out_length_GMtend_submax_2 out_length_leg_submax_2 out_length_GMapo_submax_2 out_length_msc_GM_submax_2 out_length_msc_SOL_submax_2...
        	out_elong_AT_trial_max out_elong_GMtend_trial_max out_elong_leg_trial_max out_elong_GMapo_trial_max out_elong_msc_GM_trial_max out_elong_msc_SOL_trial_max...
        	out_elong_AT_ind_max out_elong_GMtend_ind_max out_elong_leg_ind_max out_elong_GMapo_ind_max out_elong_msc_GM_ind_max out_elong_msc_SOL_ind_max...
        	out_elong_AT_common_max out_elong_GMtend_common_max out_elong_leg_common_max out_elong_GMapo_common_max out_elong_msc_GM_common_max out_elong_msc_SOL_common_max...
        	out_elong_AT_submax_1 out_elong_GMtend_submax_1 out_elong_leg_submax_1 out_elong_GMapo_submax_1 out_elong_msc_GM_submax_1 out_elong_msc_SOL_submax_1...
        	out_elong_AT_submax_2 out_elong_GMtend_submax_2 out_elong_leg_submax_2 out_elong_GMapo_submax_2 out_elong_msc_GM_submax_2 out_elong_msc_SOL_submax_2...
        	out_strain_AT_trial_max out_strain_GMtend_trial_max out_strain_leg_trial_max out_strain_GMapo_trial_max out_strain_msc_GM_trial_max out_strain_msc_SOL_trial_max...
        	out_strain_AT_ind_max out_strain_GMtend_ind_max out_strain_leg_ind_max out_strain_GMapo_ind_max out_strain_msc_GM_ind_max out_strain_msc_SOL_ind_max...
        	out_strain_AT_common_max out_strain_GMtend_common_max out_strain_leg_common_max out_strain_GMapo_common_max out_strain_msc_GM_common_max out_strain_msc_SOL_common_max...
        	out_strain_AT_submax_1 out_strain_GMtend_submax_1 out_strain_leg_submax_1 out_strain_GMapo_submax_1 out_strain_msc_GM_submax_1 out_strain_msc_SOL_submax_1...
        	out_strain_AT_submax_2 out_strain_GMtend_submax_2 out_strain_leg_submax_2 out_strain_GMapo_submax_2 out_strain_msc_GM_submax_2 out_strain_msc_SOL_submax_2...
            out_length_GMtend_Fuku_trial_max out_length_GMtend_Fuku_ind_max out_length_GMtend_Fuku_common_max out_length_GMtend_Fuku_submax_1 out_length_GMtend_Fuku_submax_2 ...
            out_elong_GMtend_Fuku_trial_max out_elong_GMtend_Fuku_ind_max out_elong_GMtend_Fuku_common_max out_elong_GMtend_Fuku_submax_1 out_elong_GMtend_Fuku_submax_2 ...
            out_strain_GMtend_Fuku_trial_max out_strain_GMtend_Fuku_ind_max out_strain_GMtend_Fuku_common_max out_strain_GMtend_Fuku_submax_1 out_strain_GMtend_Fuku_submax_2 ...
            out_length_msc_GM_Fuku_trial_max out_length_msc_GM_Fuku_ind_max out_length_msc_GM_Fuku_common_max out_length_msc_GM_Fuku_submax_1 out_length_msc_GM_Fuku_submax_2 ...
            out_elong_msc_GM_Fuku_trial_max out_elong_msc_GM_Fuku_ind_max out_elong_msc_GM_Fuku_common_max out_elong_msc_GM_Fuku_submax_1 out_elong_msc_GM_Fuku_submax_2 ...
            out_strain_msc_GM_Fuku_trial_max out_strain_msc_GM_Fuku_ind_max out_strain_msc_GM_Fuku_common_max out_strain_msc_GM_Fuku_submax_1 out_strain_msc_GM_Fuku_submax_2 ...
        	out_licht_faslen_GM_trial_max out_licht_faslen_GM_ind_max out_licht_faslen_GM_common_max out_licht_faslen_GM_zero out_licht_faslen_GM_submax_1 out_licht_faslen_GM_submax_2...
            out_licht_pennation_GM_trial_max out_licht_pennation_GM_ind_max out_licht_pennation_GM_common_max out_licht_pennation_GM_zero out_licht_pennation_GM_submax_1 out_licht_pennation_GM_submax_2...
        	out_licht_faslen_SOL_trial_max out_licht_faslen_SOL_ind_max out_licht_faslen_SOL_common_max out_licht_faslen_SOL_zero out_licht_faslen_SOL_submax_1 out_licht_faslen_SOL_submax_2...
            out_licht_pennation_SOL_trial_max out_licht_pennation_SOL_ind_max out_licht_pennation_SOL_common_max out_licht_pennation_SOL_zero out_licht_pennation_SOL_submax_1 out_licht_pennation_SOL_submax_2...
            out_contrib_GM_trial_max out_contrib_GM_ind_max out_contrib_GM_common_max out_contrib_GM_submax_1 out_contrib_GM_submax_2...
        	out_emg_gm_trial_max out_emg_gm_ind_max out_emg_gm_common_max out_emg_gm_submax_1 out_emg_gm_submax_2...
            out_emg_gl_trial_max out_emg_gl_ind_max out_emg_gl_common_max out_emg_gl_submax_1 out_emg_gl_submax_2...
            out_emg_sol_trial_max out_emg_sol_ind_max out_emg_sol_common_max out_emg_sol_submax_1 out_emg_sol_submax_2];
        %%

        
        
    end
    %% LOOP FINISHED  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    
    
    %% OUTPUT individual trial data TO FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write xls
    if ispc
        filename_output = strcat('data_output/all_passive_output_', datestr(now, 'yyyy-mm-dd HH-MM'));
        
        xlswrite(filename_output, all_passive_output_head, 1, 'A1')
        xlswrite(filename_output, all_passive_output_txt, 1, 'A2')
        xlswrite(filename_output, all_passive_output, 1, 'E2')
    else
        filename_output = strcat('data_output/all_passive_output_', datestr(now, 'yyyy-mm-dd HH-MM'), '.csv');
        csvwrite(filename_output, all_passive_output)
    end
    %% OUTPUT individual trial data TO FILE FINISHED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    

    %% GROUP CALCULATIONS - MEAN + STDAV FOR PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%  mean and stdav of each subject's INDIVIDUAL MAX ROM, force, elong, EMG, etc
    if input_project == 1
        %% BD study
        n_o_array_elements = length(CON_angle_vars{1,1}(1,:));
        if BD_count > 0
            % preallocate array
            BD_max(BD_count,n_o_array_elements) = zeros;
     %       BD_max_norm(BD_count,n_o_array_elements) = zeros;
            % collect variables per subject
            for i = 1:BD_count % per subject
                for j = 1:n_o_array_elements % per element in arrays
                     % OLD: max value
    %                BD_max(i,j) = max(BD_angle_vars{1,i}(:,j));
    %                BD_max_norm(i,j) = max(BD_angle_vars_norm{1,i}(:,j));
                     % NEW: END of array value
                    BD_max(i,j) = max(BD_angle_vars{1,i}(end,j));
    %                BD_max_norm(i,j) = max(BD_angle_vars_norm{1,i}(end,j));
                end
            end
            % calculate mean and SD of max values across subjects
            BD_ROM_mean = mean(BD_max(:,1));
            BD_ROM_SD = std(BD_max(:,1));
            BD_F_mean = mean(BD_max(:,2));
            BD_F_SD = std(BD_max(:,2));
            BD_EMG_gm_mean = mean(BD_max(:,3));
            BD_EMG_gm_SD = std(BD_max(:,3));
            BD_EMG_gl_mean = mean(BD_max(:,4));
            BD_EMG_gl_SD = std(BD_max(:,4));
            BD_EMG_sol_mean = mean(BD_max(:,5));
            BD_EMG_sol_SD = std(BD_max(:,5));

            BD_elong_AT_mean = mean(BD_max(:,6));
            BD_elong_AT_SD = std(BD_max(:,6));
            BD_elong_GMtend_mean = mean(BD_max(:,7));
            BD_elong_GMtend_SD = std(BD_max(:,7));
            BD_elong_MTU_mean = mean(BD_max(:,8)); % elong calf
            BD_elong_MTU_SD = std(BD_max(:,8));
            BD_displ_GMFAS_mean = mean(BD_max(:,9));  % GMFAS displ
            BD_displ_GMFAS_SD = std(BD_max(:,9));
            BD_elong_GMapo_mean = mean(BD_max(:,10)); 
            BD_elong_GMapo_SD = std(BD_max(:,10));
            BD_elong_msc_GM_mean = mean(BD_max(:,11)); % GM msc
            BD_elong_msc_GM_SD = std(BD_max(:,11));
            BD_elong_msc_SOL_mean = mean(BD_max(:,19)); 
            BD_elong_msc_SOL_SD = std(BD_max(:,19));

            BD_L_at_SOL_mean = mean(BD_max(:,12)); % L AT
            BD_L_at_SOL_SD = std(BD_max(:,12));
            BD_L_at_GM_mean = mean(BD_max(:,13)); % L GM tend
            BD_L_at_GM_SD = std(BD_max(:,13));
            BD_L_MTU_mean = mean(BD_max(:,14)); % L calf
            BD_L_MTU_SD = std(BD_max(:,14));
            % 15 - no Length GMfas
            BD_L_GMapo_mean = mean(BD_max(:,16)); 
            BD_L_GMapo_SD = std(BD_max(:,16));
            BD_L_msc_GM_mean = mean(BD_max(:,17)); 
            BD_L_msc_GM_SD = std(BD_max(:,17));
            BD_L_msc_SOL_mean = mean(BD_max(:,20)); 
            BD_L_msc_SOL_SD = std(BD_max(:,20));

            BD_strain_at_SOL_mean = mean(BD_max(:,21)); % L AT
            BD_strain_at_SOL_SD = std(BD_max(:,21));
            BD_strain_at_GM_mean = mean(BD_max(:,22)); % L GM tend
            BD_strain_at_GM_SD = std(BD_max(:,22));
            %BD_strain_MTU_mean = mean(BD_max(:,23)); % L calf
            %BD_strain_MTU_SD = std(BD_max(:,23));
            % 24 - no strain GMfas
            BD_strain_GMapo_mean = mean(BD_max(:,25)); 
            BD_strain_GMapo_SD = std(BD_max(:,25));
            BD_strain_msc_GM_mean = mean(BD_max(:,26)); 
            BD_strain_msc_GM_SD = std(BD_max(:,26));
            BD_strain_msc_SOL_mean = mean(BD_max(:,27)); 
            BD_strain_msc_SOL_SD = std(BD_max(:,27));
            
            BD_elong_msc_GM_licht_mean = mean(BD_max(:,28)); 
            BD_elong_msc_GM_licht_SD = std(BD_max(:,28));
            BD_elong_tend_GM_licht_mean = mean(BD_max(:,29)); 
            BD_elong_tend_GM_licht_SD = std(BD_max(:,29));
            BD_strain_msc_GM_licht_mean = mean(BD_max(:,30)); 
            BD_strain_msc_GM_licht_SD = std(BD_max(:,30));
            BD_strain_tend_GM_licht_mean = mean(BD_max(:,31)); 
            BD_strain_tend_GM_licht_SD = std(BD_max(:,31));
            
            BD_GM_length_faslen_licht_mean = mean(BD_max(:,32)); 
            BD_GM_length_faslen_licht_SD = std(BD_max(:,32)); 
            BD_GM_pennation_licht_mean = mean(BD_max(:,33)); 
            BD_GM_pennation_licht_SD = std(BD_max(:,33)); 
            BD_GM_elong_faslen_licht_mean = mean(BD_max(:,34)); 
            BD_GM_elong_faslen_licht_SD = std(BD_max(:,34)); 
            BD_GM_strain_faslen_licht_mean = mean(BD_max(:,35)); 
            BD_GM_strain_faslen_licht_SD = std(BD_max(:,35)); 
            
            BD_torque_mean = mean(BD_max(:,18));
            BD_torque_SD = std(BD_max(:,18));
            % determine common angle range
            BD_common_ROM = min(BD_max(:,1));
        end

        if CON_count > 0
            % preallocate array
            CON_max(CON_count,n_o_array_elements) = zeros;
    %        CON_max_norm(CON_count,n_o_array_elements) = zeros;
            % collect variables per subject
            for i = 1:CON_count % per subject
                for j = 1:n_o_array_elements % per element in arrays
                    % OLD: Max values
    %                CON_max(i,j) = max(CON_angle_vars{1,i}(:,j));
    %                CON_max_norm(i,j) = max(CON_angle_vars_norm{1,i}(:,j));
                    % NEW: end of array values
                    CON_max(i,j) = max(CON_angle_vars{1,i}(end,j));
    %                CON_max_norm(i,j) = max(CON_angle_vars_norm{1,i}(end,j));
                end
            end
            % calculate mean and SD of max values across subjects
            CON_ROM_mean = mean(CON_max(:,1));
            CON_ROM_SD = std(CON_max(:,1));
            CON_F_mean = mean(CON_max(:,2));
            CON_F_SD = std(CON_max(:,2));
            CON_EMG_gm_mean = mean(CON_max(:,3));
            CON_EMG_gm_SD = std(CON_max(:,3));
            CON_EMG_gl_mean = mean(CON_max(:,4));
            CON_EMG_gl_SD = std(CON_max(:,4));
            CON_EMG_sol_mean = mean(CON_max(:,5));
            CON_EMG_sol_SD = std(CON_max(:,5));

            CON_elong_AT_mean = mean(CON_max(:,6)); % elong AT
            CON_elong_AT_SD = std(CON_max(:,6));
            CON_elong_GMtend_mean = mean(CON_max(:,7)); % elong GM tend
            CON_elong_GMtend_SD = std(CON_max(:,7));
            CON_elong_MTU_mean = mean(CON_max(:,8)); % elong leg
            CON_elong_MTU_SD = std(CON_max(:,8));
            CON_displ_GMFAS_mean = mean(CON_max(:,9));  % GMFAS displ 
            CON_displ_GMFAS_SD = std(CON_max(:,9));
            CON_elong_GMapo_mean = mean(CON_max(:,10)); 
            CON_elong_GMapo_SD = std(CON_max(:,10));
            CON_elong_msc_GM_mean = mean(CON_max(:,11)); % GM msc
            CON_elong_msc_GM_SD = std(CON_max(:,11));
            CON_elong_msc_SOL_mean = mean(CON_max(:,19)); 
            CON_elong_msc_SOL_SD = std(CON_max(:,19));

            CON_L_at_SOL_mean = mean(CON_max(:,12)); % L AT
            CON_L_at_SOL_SD = std(CON_max(:,12));
            CON_L_at_GM_mean = mean(CON_max(:,13)); % L GM tend
            CON_L_at_GM_SD = std(CON_max(:,13));
            CON_L_MTU_mean = mean(CON_max(:,14)); % L leg
            CON_L_MTU_SD = std(CON_max(:,14));
            % no L GMFAS
            CON_L_GMapo_mean = mean(CON_max(:,16)); 
            CON_L_GMapo_SD = std(CON_max(:,16));
            CON_L_msc_GM_mean = mean(CON_max(:,17)); 
            CON_L_msc_GM_SD = std(CON_max(:,17));
            CON_L_msc_SOL_mean = mean(CON_max(:,20)); 
            CON_L_msc_SOL_SD = std(CON_max(:,20));

            CON_strain_at_SOL_mean = mean(CON_max(:,21)); % L AT
            CON_strain_at_SOL_SD = std(CON_max(:,21));
            CON_strain_at_GM_mean = mean(CON_max(:,22)); % L GM tend
            CON_strain_at_GM_SD = std(CON_max(:,22));
            %CON_strain_MTU_mean = mean(CON_max(:,23)); % L calf
            %CON_strain_MTU_SD = std(CON_max(:,23));
            % no strain GMfas
            CON_strain_GMapo_mean = mean(CON_max(:,25)); 
            CON_strain_GMapo_SD = std(CON_max(:,25));
            CON_strain_msc_GM_mean = mean(CON_max(:,26)); 
            CON_strain_msc_GM_SD = std(CON_max(:,26));
            CON_strain_msc_SOL_mean = mean(CON_max(:,27)); 
            CON_strain_msc_SOL_SD = std(CON_max(:,27));

            CON_elong_msc_GM_licht_mean = mean(CON_max(:,28)); 
            CON_elong_msc_GM_licht_SD = std(CON_max(:,28));
            CON_elong_tend_GM_licht_mean = mean(CON_max(:,29)); 
            CON_elong_tend_GM_licht_SD = std(CON_max(:,29));
            CON_strain_msc_GM_licht_mean = mean(CON_max(:,30)); 
            CON_strain_msc_GM_licht_SD = std(CON_max(:,30));
            CON_strain_tend_GM_licht_mean = mean(CON_max(:,31)); 
            CON_strain_tend_GM_licht_SD = std(CON_max(:,31));
            
            CON_GM_length_faslen_licht_mean = mean(CON_max(:,32)); 
            CON_GM_length_faslen_licht_SD = std(CON_max(:,32)); 
            CON_GM_pennation_licht_mean = mean(CON_max(:,33)); 
            CON_GM_pennation_licht_SD = std(CON_max(:,33)); 
            CON_GM_elong_faslen_licht_mean = mean(CON_max(:,34)); 
            CON_GM_elong_faslen_licht_SD = std(CON_max(:,34)); 
            CON_GM_strain_faslen_licht_mean = mean(CON_max(:,35)); 
            CON_GM_strain_faslen_licht_SD = std(CON_max(:,35)); 

            CON_torque_mean = mean(CON_max(:,18));
            CON_torque_SD = std(CON_max(:,18));
            % determine common angle range
            CON_common_ROM = min(CON_max(:,1));
        end
        %%
            
    elseif input_project == 2
        %% intervention study
        n_o_array_elements = max( [length(CON_PRE_angle_vars{1,1}(1,:)) length(STR_PRE_angle_vars{1,1}(1,:)) length(CON_POST_angle_vars{1,1}(1,:)) length(STR_POST_angle_vars{1,1}(1,:))] );
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
            CON_PRE_ROM_mean = mean(CON_PRE_max(:,1));
            CON_PRE_ROM_SD = std(CON_PRE_max(:,1));
            CON_PRE_F_mean = mean(CON_PRE_max(:,2));
            CON_PRE_F_SD = std(CON_PRE_max(:,2));
            CON_PRE_EMG_gm_mean = mean(CON_PRE_max(:,3));
            CON_PRE_EMG_gm_SD = std(CON_PRE_max(:,3));
            CON_PRE_EMG_gl_mean = mean(CON_PRE_max(:,4));
            CON_PRE_EMG_gl_SD = std(CON_PRE_max(:,4));
            CON_PRE_EMG_sol_mean = mean(CON_PRE_max(:,5));
            CON_PRE_EMG_sol_SD = std(CON_PRE_max(:,5));

            CON_PRE_elong_AT_mean = mean(CON_PRE_max(:,6)); % elong AT
            CON_PRE_elong_AT_SD = std(CON_PRE_max(:,6));
            CON_PRE_elong_GMtend_mean = mean(CON_PRE_max(:,7)); % elong GM tend
            CON_PRE_elong_GMtend_SD = std(CON_PRE_max(:,7));
            CON_PRE_elong_MTU_mean = mean(CON_PRE_max(:,8)); % elong leg
            CON_PRE_elong_MTU_SD = std(CON_PRE_max(:,8));
            CON_PRE_displ_GMFAS_mean = mean(CON_PRE_max(:,9));  % GMFAS displ 
            CON_PRE_displ_GMFAS_SD = std(CON_PRE_max(:,9));
            CON_PRE_elong_GMapo_mean = mean(CON_PRE_max(:,10)); 
            CON_PRE_elong_GMapo_SD = std(CON_PRE_max(:,10));
            CON_PRE_elong_msc_GM_mean = mean(CON_PRE_max(:,11)); % GM msc
            CON_PRE_elong_msc_GM_SD = std(CON_PRE_max(:,11));
            CON_PRE_elong_msc_SOL_mean = mean(CON_PRE_max(:,19)); 
            CON_PRE_elong_msc_SOL_SD = std(CON_PRE_max(:,19));

            CON_PRE_L_at_SOL_mean = mean(CON_PRE_max(:,12)); % L AT
            CON_PRE_L_at_SOL_SD = std(CON_PRE_max(:,12));
            CON_PRE_L_at_GM_mean = mean(CON_PRE_max(:,13)); % L GM tend
            CON_PRE_L_at_GM_SD = std(CON_PRE_max(:,13));
            CON_PRE_L_MTU_mean = mean(CON_PRE_max(:,14)); % L leg
            CON_PRE_L_MTU_SD = std(CON_PRE_max(:,14));
            % no L GMFAS
            CON_PRE_L_GMapo_mean = mean(CON_PRE_max(:,16)); 
            CON_PRE_L_GMapo_SD = std(CON_PRE_max(:,16));
            CON_PRE_L_msc_GM_mean = mean(CON_PRE_max(:,17)); 
            CON_PRE_L_msc_GM_SD = std(CON_PRE_max(:,17));
            CON_PRE_L_msc_SOL_mean = mean(CON_PRE_max(:,20)); 
            CON_PRE_L_msc_SOL_SD = std(CON_PRE_max(:,20));

            CON_PRE_strain_at_SOL_mean = mean(CON_PRE_max(:,21)); % L AT
            CON_PRE_strain_at_SOL_SD = std(CON_PRE_max(:,21));
            CON_PRE_strain_at_GM_mean = mean(CON_PRE_max(:,22)); % L GM tend
            CON_PRE_strain_at_GM_SD = std(CON_PRE_max(:,22));
            %CON_PRE_strain_MTU_mean = mean(CON_PRE_max(:,23)); % L calf
            %CON_PRE_strain_MTU_SD = std(CON_PRE_max(:,23));
            % no strain GMfas
            CON_PRE_strain_GMapo_mean = mean(CON_PRE_max(:,25)); 
            CON_PRE_strain_GMapo_SD = std(CON_PRE_max(:,25));
            CON_PRE_strain_msc_GM_mean = mean(CON_PRE_max(:,26)); 
            CON_PRE_strain_msc_GM_SD = std(CON_PRE_max(:,26));
            CON_PRE_strain_msc_SOL_mean = mean(CON_PRE_max(:,27)); 
            CON_PRE_strain_msc_SOL_SD = std(CON_PRE_max(:,27));
            
            CON_PRE_elong_msc_GM_licht_mean = mean(CON_PRE_max(:,28)); 
            CON_PRE_elong_msc_GM_licht_SD = std(CON_PRE_max(:,28));
            CON_PRE_elong_tend_GM_licht_mean = mean(CON_PRE_max(:,29)); 
            CON_PRE_elong_tend_GM_licht_SD = std(CON_PRE_max(:,29));
            CON_PRE_strain_msc_GM_licht_mean = mean(CON_PRE_max(:,30)); 
            CON_PRE_strain_msc_GM_licht_SD = std(CON_PRE_max(:,30));
            CON_PRE_strain_tend_GM_licht_mean = mean(CON_PRE_max(:,31)); 
            CON_PRE_strain_tend_GM_licht_SD = std(CON_PRE_max(:,31));

            CON_PRE_GM_length_faslen_licht_mean = mean(CON_PRE_max(:,32)); 
            CON_PRE_GM_length_faslen_licht_SD = std(CON_PRE_max(:,32)); 
            CON_PRE_GM_pennation_licht_mean = mean(CON_PRE_max(:,33)); 
            CON_PRE_GM_pennation_licht_SD = std(CON_PRE_max(:,33)); 
            CON_PRE_GM_elong_faslen_licht_mean = mean(CON_PRE_max(:,34)); 
            CON_PRE_GM_elong_faslen_licht_SD = std(CON_PRE_max(:,34)); 
            CON_PRE_GM_strain_faslen_licht_mean = mean(CON_PRE_max(:,35)); 
            CON_PRE_GM_strain_faslen_licht_SD = std(CON_PRE_max(:,35)); 
            
            CON_PRE_torque_mean = mean(CON_PRE_max(:,18));
            CON_PRE_torque_SD = std(CON_PRE_max(:,18));
            % determine common angle range
            CON_PRE_common_ROM = min(CON_PRE_max(:,1));
        end
        % STR PRE
        if STR_PRE_count > 0
            % preallocate array
            STR_PRE_max(STR_PRE_count,n_o_array_elements) = zeros;
    %        STR_PRE_max_norm(STR_PRE_count,n_o_array_elements) = zeros;
            % collect variables per subject
            for i = 1:STR_PRE_count % per subject
                for j = 1:n_o_array_elements % per element in arrays
                    % OLD: Max values
    %                STR_PRE_max(i,j) = max(STR_PRE_angle_vars{1,i}(:,j));
    %                STR_PRE_max_norm(i,j) = max(STR_PRE_angle_vars_norm{1,i}(:,j));
                    % NEW: end of array values
                    STR_PRE_max(i,j) = max(STR_PRE_angle_vars{1,i}(end,j));
     %               STR_PRE_max_norm(i,j) = max(STR_PRE_angle_vars_norm{1,i}(end,j));
                end
            end
            % calculate mean and SD of max values across subjects
            STR_PRE_ROM_mean = mean(STR_PRE_max(:,1));
            STR_PRE_ROM_SD = std(STR_PRE_max(:,1));
            STR_PRE_F_mean = mean(STR_PRE_max(:,2));
            STR_PRE_F_SD = std(STR_PRE_max(:,2));
            STR_PRE_EMG_gm_mean = mean(STR_PRE_max(:,3));
            STR_PRE_EMG_gm_SD = std(STR_PRE_max(:,3));
            STR_PRE_EMG_gl_mean = mean(STR_PRE_max(:,4));
            STR_PRE_EMG_gl_SD = std(STR_PRE_max(:,4));
            STR_PRE_EMG_sol_mean = mean(STR_PRE_max(:,5));
            STR_PRE_EMG_sol_SD = std(STR_PRE_max(:,5));

            STR_PRE_elong_AT_mean = mean(STR_PRE_max(:,6)); % elong AT
            STR_PRE_elong_AT_SD = std(STR_PRE_max(:,6));
            STR_PRE_elong_GMtend_mean = mean(STR_PRE_max(:,7)); % elong GM tend
            STR_PRE_elong_GMtend_SD = std(STR_PRE_max(:,7));
            STR_PRE_elong_MTU_mean = mean(STR_PRE_max(:,8)); % elong leg
            STR_PRE_elong_MTU_SD = std(STR_PRE_max(:,8));
            STR_PRE_displ_GMFAS_mean = mean(STR_PRE_max(:,9));  % GMFAS displ 
            STR_PRE_displ_GMFAS_SD = std(STR_PRE_max(:,9));
            STR_PRE_elong_GMapo_mean = mean(STR_PRE_max(:,10)); 
            STR_PRE_elong_GMapo_SD = std(STR_PRE_max(:,10));
            STR_PRE_elong_msc_GM_mean = mean(STR_PRE_max(:,11)); % GM msc
            STR_PRE_elong_msc_GM_SD = std(STR_PRE_max(:,11));
            STR_PRE_elong_msc_SOL_mean = mean(STR_PRE_max(:,19)); 
            STR_PRE_elong_msc_SOL_SD = std(STR_PRE_max(:,19));

            STR_PRE_L_at_SOL_mean = mean(STR_PRE_max(:,12)); % L AT
            STR_PRE_L_at_SOL_SD = std(STR_PRE_max(:,12));
            STR_PRE_L_at_GM_mean = mean(STR_PRE_max(:,13)); % L GM tend
            STR_PRE_L_at_GM_SD = std(STR_PRE_max(:,13));
            STR_PRE_L_MTU_mean = mean(STR_PRE_max(:,14)); % L leg
            STR_PRE_L_MTU_SD = std(STR_PRE_max(:,14));
            % no L GMFAS
            STR_PRE_L_GMapo_mean = mean(STR_PRE_max(:,16)); 
            STR_PRE_L_GMapo_SD = std(STR_PRE_max(:,16));
            STR_PRE_L_msc_GM_mean = mean(STR_PRE_max(:,17)); 
            STR_PRE_L_msc_GM_SD = std(STR_PRE_max(:,17));
            STR_PRE_L_msc_SOL_mean = mean(STR_PRE_max(:,20)); 
            STR_PRE_L_msc_SOL_SD = std(STR_PRE_max(:,20));

            STR_PRE_strain_at_SOL_mean = mean(STR_PRE_max(:,21)); % L AT
            STR_PRE_strain_at_SOL_SD = std(STR_PRE_max(:,21));
            STR_PRE_strain_at_GM_mean = mean(STR_PRE_max(:,22)); % L GM tend
            STR_PRE_strain_at_GM_SD = std(STR_PRE_max(:,22));
            %STR_PRE_strain_MTU_mean = mean(STR_PRE_max(:,23)); % L calf
            %STR_PRE_strain_MTU_SD = std(STR_PRE_max(:,23));
            % no strain GMfas
            STR_PRE_strain_GMapo_mean = mean(STR_PRE_max(:,25)); 
            STR_PRE_strain_GMapo_SD = std(STR_PRE_max(:,25));
            STR_PRE_strain_msc_GM_mean = mean(STR_PRE_max(:,26)); 
            STR_PRE_strain_msc_GM_SD = std(STR_PRE_max(:,26));
            STR_PRE_strain_msc_SOL_mean = mean(STR_PRE_max(:,27)); 
            STR_PRE_strain_msc_SOL_SD = std(STR_PRE_max(:,27));

            STR_PRE_elong_msc_GM_licht_mean = mean(STR_PRE_max(:,28)); 
            STR_PRE_elong_msc_GM_licht_SD = std(STR_PRE_max(:,28));
            STR_PRE_elong_tend_GM_licht_mean = mean(STR_PRE_max(:,29)); 
            STR_PRE_elong_tend_GM_licht_SD = std(STR_PRE_max(:,29));
            STR_PRE_strain_msc_GM_licht_mean = mean(STR_PRE_max(:,30)); 
            STR_PRE_strain_msc_GM_licht_SD = std(STR_PRE_max(:,30));
            STR_PRE_strain_tend_GM_licht_mean = mean(STR_PRE_max(:,31)); 
            STR_PRE_strain_tend_GM_licht_SD = std(STR_PRE_max(:,31));
            
            STR_PRE_GM_length_faslen_licht_mean = mean(STR_PRE_max(:,32)); 
            STR_PRE_GM_length_faslen_licht_SD = std(STR_PRE_max(:,32)); 
            STR_PRE_GM_pennation_licht_mean = mean(STR_PRE_max(:,33)); 
            STR_PRE_GM_pennation_licht_SD = std(STR_PRE_max(:,33)); 
            STR_PRE_GM_elong_faslen_licht_mean = mean(STR_PRE_max(:,34)); 
            STR_PRE_GM_elong_faslen_licht_SD = std(STR_PRE_max(:,34)); 
            STR_PRE_GM_strain_faslen_licht_mean = mean(STR_PRE_max(:,35)); 
            STR_PRE_GM_strain_faslen_licht_SD = std(STR_PRE_max(:,35)); 

            STR_PRE_torque_mean = mean(STR_PRE_max(:,18));
            STR_PRE_torque_SD = std(STR_PRE_max(:,18));
            % determine common angle range
            STR_PRE_common_ROM = min(STR_PRE_max(:,1));
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
            CON_POST_ROM_mean = mean(CON_POST_max(:,1));
            CON_POST_ROM_SD = std(CON_POST_max(:,1));
            CON_POST_F_mean = mean(CON_POST_max(:,2));
            CON_POST_F_SD = std(CON_POST_max(:,2));
            CON_POST_EMG_gm_mean = mean(CON_POST_max(:,3));
            CON_POST_EMG_gm_SD = std(CON_POST_max(:,3));
            CON_POST_EMG_gl_mean = mean(CON_POST_max(:,4));
            CON_POST_EMG_gl_SD = std(CON_POST_max(:,4));
            CON_POST_EMG_sol_mean = mean(CON_POST_max(:,5));
            CON_POST_EMG_sol_SD = std(CON_POST_max(:,5));

            CON_POST_elong_AT_mean = mean(CON_POST_max(:,6)); % elong AT
            CON_POST_elong_AT_SD = std(CON_POST_max(:,6));
            CON_POST_elong_GMtend_mean = mean(CON_POST_max(:,7)); % elong GM tend
            CON_POST_elong_GMtend_SD = std(CON_POST_max(:,7));
            CON_POST_elong_MTU_mean = mean(CON_POST_max(:,8)); % elong leg
            CON_POST_elong_MTU_SD = std(CON_POST_max(:,8));
            CON_POST_displ_GMFAS_mean = mean(CON_POST_max(:,9));  % GMFAS displ 
            CON_POST_displ_GMFAS_SD = std(CON_POST_max(:,9));
            CON_POST_elong_GMapo_mean = mean(CON_POST_max(:,10)); 
            CON_POST_elong_GMapo_SD = std(CON_POST_max(:,10));
            CON_POST_elong_msc_GM_mean = mean(CON_POST_max(:,11)); % GM msc
            CON_POST_elong_msc_GM_SD = std(CON_POST_max(:,11));
            CON_POST_elong_msc_SOL_mean = mean(CON_POST_max(:,19)); 
            CON_POST_elong_msc_SOL_SD = std(CON_POST_max(:,19));

            CON_POST_L_at_SOL_mean = mean(CON_POST_max(:,12)); % L AT
            CON_POST_L_at_SOL_SD = std(CON_POST_max(:,12));
            CON_POST_L_at_GM_mean = mean(CON_POST_max(:,13)); % L GM tend
            CON_POST_L_at_GM_SD = std(CON_POST_max(:,13));
            CON_POST_L_MTU_mean = mean(CON_POST_max(:,14)); % L leg
            CON_POST_L_MTU_SD = std(CON_POST_max(:,14));
            % no L GMFAS
            CON_POST_L_GMapo_mean = mean(CON_POST_max(:,16)); 
            CON_POST_L_GMapo_SD = std(CON_POST_max(:,16));
            CON_POST_L_msc_GM_mean = mean(CON_POST_max(:,17)); 
            CON_POST_L_msc_GM_SD = std(CON_POST_max(:,17));
            CON_POST_L_msc_SOL_mean = mean(CON_POST_max(:,20)); 
            CON_POST_L_msc_SOL_SD = std(CON_POST_max(:,20));

            CON_POST_strain_at_SOL_mean = mean(CON_POST_max(:,21)); % L AT
            CON_POST_strain_at_SOL_SD = std(CON_POST_max(:,21));
            CON_POST_strain_at_GM_mean = mean(CON_POST_max(:,22)); % L GM tend
            CON_POST_strain_at_GM_SD = std(CON_POST_max(:,22));
            %CON_POST_strain_MTU_mean = mean(CON_POST_max(:,23)); % L calf
            %CON_POST_strain_MTU_SD = std(CON_POST_max(:,23));
            % no strain GMfas
            CON_POST_strain_GMapo_mean = mean(CON_POST_max(:,25)); 
            CON_POST_strain_GMapo_SD = std(CON_POST_max(:,25));
            CON_POST_strain_msc_GM_mean = mean(CON_POST_max(:,26)); 
            CON_POST_strain_msc_GM_SD = std(CON_POST_max(:,26));
            CON_POST_strain_msc_SOL_mean = mean(CON_POST_max(:,27)); 
            CON_POST_strain_msc_SOL_SD = std(CON_POST_max(:,27));

            CON_POST_elong_msc_GM_licht_mean = mean(CON_POST_max(:,28)); 
            CON_POST_elong_msc_GM_licht_SD = std(CON_POST_max(:,28));
            CON_POST_elong_tend_GM_licht_mean = mean(CON_POST_max(:,29)); 
            CON_POST_elong_tend_GM_licht_SD = std(CON_POST_max(:,29));
            CON_POST_strain_msc_GM_licht_mean = mean(CON_POST_max(:,30)); 
            CON_POST_strain_msc_GM_licht_SD = std(CON_POST_max(:,30));
            CON_POST_strain_tend_GM_licht_mean = mean(CON_POST_max(:,31)); 
            CON_POST_strain_tend_GM_licht_SD = std(CON_POST_max(:,31));
            
            CON_POST_GM_length_faslen_licht_mean = mean(CON_POST_max(:,32)); 
            CON_POST_GM_length_faslen_licht_SD = std(CON_POST_max(:,32)); 
            CON_POST_GM_pennation_licht_mean = mean(CON_POST_max(:,33)); 
            CON_POST_GM_pennation_licht_SD = std(CON_POST_max(:,33)); 
            CON_POST_GM_elong_faslen_licht_mean = mean(CON_POST_max(:,34)); 
            CON_POST_GM_elong_faslen_licht_SD = std(CON_POST_max(:,34)); 
            CON_POST_GM_strain_faslen_licht_mean = mean(CON_POST_max(:,35)); 
            CON_POST_GM_strain_faslen_licht_SD = std(CON_POST_max(:,35)); 

            CON_POST_torque_mean = mean(CON_POST_max(:,18));
            CON_POST_torque_SD = std(CON_POST_max(:,18));
            % determine common angle range
            CON_POST_common_ROM = min(CON_POST_max(:,1));
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
            STR_POST_ROM_mean = mean(STR_POST_max(:,1));
            STR_POST_ROM_SD = std(STR_POST_max(:,1));
            STR_POST_F_mean = mean(STR_POST_max(:,2));
            STR_POST_F_SD = std(STR_POST_max(:,2));
            STR_POST_EMG_gm_mean = mean(STR_POST_max(:,3));
            STR_POST_EMG_gm_SD = std(STR_POST_max(:,3));
            STR_POST_EMG_gl_mean = mean(STR_POST_max(:,4));
            STR_POST_EMG_gl_SD = std(STR_POST_max(:,4));
            STR_POST_EMG_sol_mean = mean(STR_POST_max(:,5));
            STR_POST_EMG_sol_SD = std(STR_POST_max(:,5));

            STR_POST_elong_AT_mean = mean(STR_POST_max(:,6)); % elong AT
            STR_POST_elong_AT_SD = std(STR_POST_max(:,6));
            STR_POST_elong_GMtend_mean = mean(STR_POST_max(:,7)); % elong GM tend
            STR_POST_elong_GMtend_SD = std(STR_POST_max(:,7));
            STR_POST_elong_MTU_mean = mean(STR_POST_max(:,8)); % elong leg
            STR_POST_elong_MTU_SD = std(STR_POST_max(:,8));
            STR_POST_displ_GMFAS_mean = mean(STR_POST_max(:,9));  % GMFAS displ 
            STR_POST_displ_GMFAS_SD = std(STR_POST_max(:,9));
            STR_POST_elong_GMapo_mean = mean(STR_POST_max(:,10)); 
            STR_POST_elong_GMapo_SD = std(STR_POST_max(:,10));
            STR_POST_elong_msc_GM_mean = mean(STR_POST_max(:,11)); % GM msc
            STR_POST_elong_msc_GM_SD = std(STR_POST_max(:,11));
            STR_POST_elong_msc_SOL_mean = mean(STR_POST_max(:,19)); 
            STR_POST_elong_msc_SOL_SD = std(STR_POST_max(:,19));

            STR_POST_L_at_SOL_mean = mean(STR_POST_max(:,12)); % L AT
            STR_POST_L_at_SOL_SD = std(STR_POST_max(:,12));
            STR_POST_L_at_GM_mean = mean(STR_POST_max(:,13)); % L GM tend
            STR_POST_L_at_GM_SD = std(STR_POST_max(:,13));
            STR_POST_L_MTU_mean = mean(STR_POST_max(:,14)); % L leg
            STR_POST_L_MTU_SD = std(STR_POST_max(:,14));
            % no L GMFAS
            STR_POST_L_GMapo_mean = mean(STR_POST_max(:,16)); 
            STR_POST_L_GMapo_SD = std(STR_POST_max(:,16));
            STR_POST_L_msc_GM_mean = mean(STR_POST_max(:,17)); 
            STR_POST_L_msc_GM_SD = std(STR_POST_max(:,17));
            STR_POST_L_msc_SOL_mean = mean(STR_POST_max(:,20)); 
            STR_POST_L_msc_SOL_SD = std(STR_POST_max(:,20));

            STR_POST_strain_at_SOL_mean = mean(STR_POST_max(:,21)); % L AT
            STR_POST_strain_at_SOL_SD = std(STR_POST_max(:,21));
            STR_POST_strain_at_GM_mean = mean(STR_POST_max(:,22)); % L GM tend
            STR_POST_strain_at_GM_SD = std(STR_POST_max(:,22));
            %STR_POST_strain_MTU_mean = mean(STR_POST_max(:,23)); % L calf
            %STR_POST_strain_MTU_SD = std(STR_POST_max(:,23));
            % no strain GMfas
            STR_POST_strain_GMapo_mean = mean(STR_POST_max(:,25)); 
            STR_POST_strain_GMapo_SD = std(STR_POST_max(:,25));
            STR_POST_strain_msc_GM_mean = mean(STR_POST_max(:,26)); 
            STR_POST_strain_msc_GM_SD = std(STR_POST_max(:,26));
            STR_POST_strain_msc_SOL_mean = mean(STR_POST_max(:,27)); 
            STR_POST_strain_msc_SOL_SD = std(STR_POST_max(:,27));

            STR_POST_elong_msc_GM_licht_mean = mean(STR_POST_max(:,28)); 
            STR_POST_elong_msc_GM_licht_SD = std(STR_POST_max(:,28));
            STR_POST_elong_tend_GM_licht_mean = mean(STR_POST_max(:,29)); 
            STR_POST_elong_tend_GM_licht_SD = std(STR_POST_max(:,29));
            STR_POST_strain_msc_GM_licht_mean = mean(STR_POST_max(:,30)); 
            STR_POST_strain_msc_GM_licht_SD = std(STR_POST_max(:,30));
            STR_POST_strain_tend_GM_licht_mean = mean(STR_POST_max(:,31)); 
            STR_POST_strain_tend_GM_licht_SD = std(STR_POST_max(:,31));

            STR_POST_GM_length_faslen_licht_mean = mean(STR_POST_max(:,32)); 
            STR_POST_GM_length_faslen_licht_SD = std(STR_POST_max(:,32)); 
            STR_POST_GM_pennation_licht_mean = mean(STR_POST_max(:,33)); 
            STR_POST_GM_pennation_licht_SD = std(STR_POST_max(:,33)); 
            STR_POST_GM_elong_faslen_licht_mean = mean(STR_POST_max(:,34)); 
            STR_POST_GM_elong_faslen_licht_SD = std(STR_POST_max(:,34)); 
            STR_POST_GM_strain_faslen_licht_mean = mean(STR_POST_max(:,35)); 
            STR_POST_GM_strain_faslen_licht_SD = std(STR_POST_max(:,35)); 

            STR_POST_torque_mean = mean(STR_POST_max(:,18));
            STR_POST_torque_SD = std(STR_POST_max(:,18));
            % determine common angle range
            STR_POST_common_ROM = min(STR_POST_max(:,1));
        end
        %%    
    end
    %% GROUP CALCULATIONS - NUMBERS FINISHED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    
    %% GROUP CALCULATIONS - % of leg lengths for BD study %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if input_project == 1 && toggle_normalization == 1

        %%% in angle_vars arrays, replace absolute values with
        %%% normalization to initial leg length:
        
        % --- angle_vars arrays:
        % 14 L_MTU/calf
        
        % 13 L_GMtend
        % 29 elong tend GM (from Lichtwark/Fukunaga)
        % 31 strain tend GM (from Lichtwark/Fukunaga)
        
        % 17 L_msc_GM
        % 28 elong msc GM (from Lichtwark/Fukunaga)
        % 30 strain msc GM (from Lichtwark/Fukunaga)
        
        % 32 LENGTH faslen GM (from Lichtwark)
        for i = 1:length(BD_count)
            % GM tendon length, elong
            BD_angle_vars{i}(:,13) = BD_angle_vars{i}(:,13)/BD_angle_vars{i}(1,14);
            BD_angle_vars{i}(:,29) = BD_angle_vars{i}(:,29)/BD_angle_vars{i}(1,14);
            % GM muscle length, elong
            BD_angle_vars{i}(:,17) = BD_angle_vars{i}(:,17)/BD_angle_vars{i}(1,14);
            BD_angle_vars{i}(:,28) = BD_angle_vars{i}(:,28)/BD_angle_vars{i}(1,14);
            % GM fascicle length, elong
            BD_angle_vars{i}(:,32) = BD_angle_vars{i}(:,32)/BD_angle_vars{i}(1,14);
            BD_angle_vars{i}(:,34) = BD_angle_vars{i}(:,34)/BD_angle_vars{i}(1,14);
        end
        for i = 1:length(CON_count)
            % GM tendon length, elong
            CON_angle_vars{i}(:,13) = CON_angle_vars{i}(:,13)/CON_angle_vars{i}(1,14);
            CON_angle_vars{i}(:,29) = CON_angle_vars{i}(:,29)/CON_angle_vars{i}(1,14);
            % GM muscle length, elong
            CON_angle_vars{i}(:,17) = CON_angle_vars{i}(:,17)/CON_angle_vars{i}(1,14);
            CON_angle_vars{i}(:,28) = CON_angle_vars{i}(:,28)/CON_angle_vars{i}(1,14);
            % GM fascicle length, elong
            CON_angle_vars{i}(:,32) = CON_angle_vars{i}(:,32)/CON_angle_vars{i}(1,14);
            CON_angle_vars{i}(:,34) = CON_angle_vars{i}(:,34)/CON_angle_vars{i}(1,14);
        end

        %%% adjust relevant axes: % MMM TODO values
        axis_faslen = [-1 40 -inf inf];
        axis_faslen_elong = [-1 40 -inf inf];
        axis_faslen_str = [-1 40 -inf inf];
        
        axis_el_GM = [-1 40 -inf inf];
        axis_str_GM = [-1 40 -inf inf];
        axis_len_GM = [-1 40 -inf inf];

        axis_el_GMtend = [-1 40 -inf inf];
        axis_str_GMtend = [-1 40 -inf inf];
        axis_len_GMtend = [-1 40 -inf inf];

        
    % MMM GOON TODO - what about mean values?

    % old:
%         % preallocate
%         BD_norm_len_GM{BD_count} = zeros;
%         CON_norm_len_GM{CON_count} = zeros;
%         
%         % GM lengths, normalized to leg length:
%         %   1 = angle
%         %   2 = tendon     - LENGTHS
%         %   3 = muscle
%         %   4 = fascicle
%         %   5 = tendon     - ELONGATIONS
%         %   6 = muscle
%         %   7 = fascicle
%         for i = 1:length(BD_count)
%             BD_norm_len_GM{i}(:,1) = BD_angle_vars{i}(:,1);
%             BD_norm_len_GM{i}(:,2) = (BD_angle_vars{i}(:,29)+BD_angle_vars{i}(1,13))/BD_angle_vars{i}(1,14); % elong->len by initial tendon length. Divide by initial leg length.
%             BD_norm_len_GM{i}(:,3) = (BD_angle_vars{i}(:,28)+BD_angle_vars{i}(1,17))/BD_angle_vars{i}(1,14); % elong->len by initial muscle length. divide by initial leg length
%             BD_norm_len_GM{i}(:,4) = BD_angle_vars{i}(:,32)/BD_angle_vars{i}(1,14); % divide by initial leg length
%             BD_norm_len_GM{i}(:,5) = BD_angle_vars{i}(:,29)/BD_angle_vars{i}(1,14); % Divide by initial leg length
%             BD_norm_len_GM{i}(:,6) = BD_angle_vars{i}(:,28)/BD_angle_vars{i}(1,14); % divide by initial leg length
%             BD_norm_len_GM{i}(:,7) = (BD_angle_vars{i}(:,32)-BD_angle_vars{i}(1,32))/BD_angle_vars{i}(1,14); % len->elong by initial faslen. divide by initial leg length
%         end
%         for i = 1:length(CON_count)
%             CON_norm_len_GM{i}(:,1) = CON_angle_vars{i}(:,1);
%             CON_norm_len_GM{i}(:,2) = (CON_angle_vars{i}(:,29)+CON_angle_vars{i}(1,13))/CON_angle_vars{i}(1,14); % elong->len by initial tendon length. Divide by initial leg length.
%             CON_norm_len_GM{i}(:,3) = (CON_angle_vars{i}(:,28)+CON_angle_vars{i}(1,17))/CON_angle_vars{i}(1,14); % elong->len by initial muscle length. divide by initial leg length
%             CON_norm_len_GM{i}(:,4) = CON_angle_vars{i}(:,32)/CON_angle_vars{i}(1,14); % divide by initial leg length
%             CON_norm_len_GM{i}(:,5) = CON_angle_vars{i}(:,29)/CON_angle_vars{i}(1,14); % Divide by initial leg length
%             CON_norm_len_GM{i}(:,6) = CON_angle_vars{i}(:,28)/CON_angle_vars{i}(1,14); % divide by initial leg length
%             CON_norm_len_GM{i}(:,7) = (CON_angle_vars{i}(:,32)-CON_angle_vars{i}(1,32))/CON_angle_vars{i}(1,14); % len->elong by initial faslen. divide by initial leg length
%         end
        
    end
    
    
    
    
    
    
    %% GROUP CALCULATIONS - MEAN ARRAYS FOR PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if input_project == 1
        %% BD study
        %%% average ABSOLUTE arrays, up to all subjects' COMMON MAX ROM, for force, elong, EMG
        % BD_angle_vars_mean, CON_angle_vars_mean
        if BD_count > 0
            % preallocate
            len = 10000;
            for i = 1:BD_count
                if length(BD_angle_vars{:,i}) < len
                    len = length(BD_angle_vars{:,i});
                end
            end
            BD_angle_vars_mean_tmp(len,n_o_array_elements,BD_count) = zeros;

            % BD_angle_vars has same angles (column 1) for all subjects, so max common ROM will be on the same line for all subjects
            loc_end = find(BD_angle_vars{1,BD_count}(:,1) >= (BD_common_ROM - 0.00001), 1, 'first'); % using last subject (BD_count) - could use any, angle is the same in all
            for i = 1:BD_count
                BD_angle_vars_mean_tmp(:,:,i) = BD_angle_vars{i}(1:loc_end,:);
            end
            BD_angle_vars_mean = nanmean(BD_angle_vars_mean_tmp, 3);
            BD_angle_vars_SD = nanstd(BD_angle_vars_mean_tmp,1,3);
        end
        if CON_count > 0
            % preallocate
            len = 10000;
            for i = 1:CON_count
                if length(CON_angle_vars{:,i}) < len
                    len = length(CON_angle_vars{:,i});
                end
            end
            CON_angle_vars_mean_tmp(len,n_o_array_elements,CON_count) = zeros;

            % BD_angle_vars has same angles (column 1) for all subjects, so max common ROM will be on the same line for all subjects
            loc_end = find(CON_angle_vars{1,CON_count}(:,1) >= (CON_common_ROM - 0.00001), 1, 'first'); % using last subject (CON_count) - could use any, angle is the same in all
            for i = 1:CON_count
                CON_angle_vars_mean_tmp(:,:,i) = CON_angle_vars{i}(1:loc_end,:);
            end
            CON_angle_vars_mean = nanmean(CON_angle_vars_mean_tmp, 3);
            CON_angle_vars_SD = nanstd(CON_angle_vars_mean_tmp,1,3);
        end



        %%% average NORMALIZED arrays, up to all subjects' INDIVIDUAL ROM, for force, elong, EMG
        % BD_angle_vars_norm_mean CON_angle_vars_norm_mean
        if BD_count > 0
            % preallocate
            BD_angle_vars_norm_mean_tmp(length(BD_angle_vars_norm{:,1}),n_o_array_elements,BD_count) = zeros;

            % BD_angle_vars_norm has same angles (column 1) for all subjects, from 0 to 100% for all subjects
            for i = 1:BD_count
                BD_angle_vars_norm_mean_tmp(:,:,i) = BD_angle_vars_norm{i}(:,:);
            end
            BD_angle_vars_norm_mean = nanmean(BD_angle_vars_norm_mean_tmp, 3);
        end
        if CON_count > 0
            % preallocate
            CON_angle_vars_norm_mean_tmp(length(CON_angle_vars_norm{:,1}),n_o_array_elements,CON_count) = zeros;

            for i = 1:CON_count
                CON_angle_vars_norm_mean_tmp(:,:,i) = CON_angle_vars_norm{i}(:,:);
            end
            CON_angle_vars_norm_mean = nanmean(CON_angle_vars_norm_mean_tmp, 3);
        end

        clear BD_angle_vars_mean_tmp BD_angle_vars_norm_mean_tmp CON_angle_vars_mean_tmp CON_angle_vars_norm_mean_tmp
        %%
    elseif input_project == 2
        %% intervention
        
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
            loc_end = find(CON_PRE_angle_vars{1,CON_PRE_count}(:,1) >= (CON_PRE_common_ROM - 0.00001), 1, 'first'); % using last subject (CON_PRE_count) - could use any, angle is the same in all
            for i = 1:CON_PRE_count
                CON_PRE_angle_vars_mean_tmp(:,:,i) = CON_PRE_angle_vars{i}(1:loc_end,:);
            end
            CON_PRE_angle_vars_mean = nanmean(CON_PRE_angle_vars_mean_tmp, 3);
            %CON_PRE_angle_vars_SD = nanstd(CON_PRE_angle_vars_mean_tmp,1,3);
            
            %%% average NORMALIZED arrays, up to all subjects' INDIVIDUAL ROM, for force, elong, EMG
            
            % preallocate
            CON_PRE_angle_vars_norm_mean_tmp(length(CON_PRE_angle_vars_norm{:,1}),n_o_array_elements,CON_PRE_count) = zeros;

            % CON_PRE_angle_vars_norm has same angles (column 1) for all subjects, from 0 to 100% for all subjects
            for i = 1:CON_PRE_count
                CON_PRE_angle_vars_norm_mean_tmp(:,:,i) = CON_PRE_angle_vars_norm{i}(:,:);
            end
            CON_PRE_angle_vars_norm_mean = nanmean(CON_PRE_angle_vars_norm_mean_tmp, 3);
            
            %%% clean up
            clear CON_PRE_angle_vars_mean_tmp CON_PRE_angle_vars_norm_mean_tmp
        end

        
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
            loc_end = find(STR_PRE_angle_vars{1,STR_PRE_count}(:,1) >= (STR_PRE_common_ROM - 0.00001), 1, 'first'); % using last subject (STR_PRE_count) - could use any, angle is the same in all
            for i = 1:STR_PRE_count
                STR_PRE_angle_vars_mean_tmp(:,:,i) = STR_PRE_angle_vars{i}(1:loc_end,:);
            end
            STR_PRE_angle_vars_mean = nanmean(STR_PRE_angle_vars_mean_tmp, 3);
            %STR_PRE_angle_vars_SD = nanstd(STR_PRE_angle_vars_mean_tmp,1,3);

            %%% average NORMALIZED arrays, up to all subjects' INDIVIDUAL ROM, for force, elong, EMG

            % preallocate
            STR_PRE_angle_vars_norm_mean_tmp(length(STR_PRE_angle_vars_norm{:,1}),n_o_array_elements,STR_PRE_count) = zeros;

            % STR_PRE_angle_vars_norm has same angles (column 1) for all subjects, from 0 to 100% for all subjects
            for i = 1:STR_PRE_count
                STR_PRE_angle_vars_norm_mean_tmp(:,:,i) = STR_PRE_angle_vars_norm{i}(:,:);
            end
            STR_PRE_angle_vars_norm_mean = nanmean(STR_PRE_angle_vars_norm_mean_tmp, 3);

            %%% clean up

            clear STR_PRE_angle_vars_mean_tmp STR_PRE_angle_vars_norm_mean_tmp
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
            loc_end = find(CON_POST_angle_vars{1,CON_POST_count}(:,1) >= (CON_POST_common_ROM - 0.00001), 1, 'first'); % using last subject (CON_POST_count) - could use any, angle is the same in all
            for i = 1:CON_POST_count
                CON_POST_angle_vars_mean_tmp(:,:,i) = CON_POST_angle_vars{i}(1:loc_end,:);
            end
            CON_POST_angle_vars_mean = nanmean(CON_POST_angle_vars_mean_tmp, 3);
            %CON_POST_angle_vars_SD = nanstd(CON_POST_angle_vars_mean_tmp,1,3);

            %%% average NORMALIZED arrays, up to all subjects' INDIVIDUAL ROM, for force, elong, EMG

            % preallocate
            CON_POST_angle_vars_norm_mean_tmp(length(CON_POST_angle_vars_norm{:,1}),n_o_array_elements,CON_POST_count) = zeros;

            % CON_POST_angle_vars_norm has same angles (column 1) for all subjects, from 0 to 100% for all subjects
            for i = 1:CON_POST_count
                CON_POST_angle_vars_norm_mean_tmp(:,:,i) = CON_POST_angle_vars_norm{i}(:,:);
            end
            CON_POST_angle_vars_norm_mean = nanmean(CON_POST_angle_vars_norm_mean_tmp, 3);

            %%% clean up

            clear CON_POST_angle_vars_mean_tmp CON_POST_angle_vars_norm_mean_tmp
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
            loc_end = find(STR_POST_angle_vars{1,STR_POST_count}(:,1) >= (STR_POST_common_ROM - 0.00001), 1, 'first'); % using last subject (STR_POST_count) - could use any, angle is the same in all
            for i = 1:STR_POST_count
                STR_POST_angle_vars_mean_tmp(:,:,i) = STR_POST_angle_vars{i}(1:loc_end,:);
            end
            STR_POST_angle_vars_mean = nanmean(STR_POST_angle_vars_mean_tmp, 3);
            %STR_POST_angle_vars_SD = nanstd(STR_POST_angle_vars_mean_tmp,1,3);

            %%% average NORMALIZED arrays, up to all subjects' INDIVIDUAL ROM, for force, elong, EMG

            % preallocate
            STR_POST_angle_vars_norm_mean_tmp(length(STR_POST_angle_vars_norm{:,1}),n_o_array_elements,STR_POST_count) = zeros;

            % STR_POST_angle_vars_norm has same angles (column 1) for all subjects, from 0 to 100% for all subjects
            for i = 1:STR_POST_count
                STR_POST_angle_vars_norm_mean_tmp(:,:,i) = STR_POST_angle_vars_norm{i}(:,:);
            end
            STR_POST_angle_vars_norm_mean = nanmean(STR_POST_angle_vars_norm_mean_tmp, 3);

            %%% clean up

            clear STR_POST_angle_vars_mean_tmp STR_POST_angle_vars_norm_mean_tmp
        end
        %%
    end
    %% GROUP CALCULATIONS - ARRAYS FINISHED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    
    
    %% GROUP CALCULATIONS - % of leg lengths for BD study %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %% PLOT GROUP FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if input_project == 1
        %% BD study
            %% FORCE-angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('force vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,2),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,2),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_F_mean, BD_F_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_F_mean, CON_F_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_F_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_F_mean, CON_ROM_SD, '*b')
                axis(axis_force)
                xlabel('Gonio angle (°)')
                ylabel('Force (N)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('force vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,2))
                end
                axis(axis_force)
                xlabel('Gonio angle (°)')
                ylabel('Force (N)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('force vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,2))
                end
                axis(axis_force)
                xlabel('Gonio angle (°)')
                ylabel('Force (N)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('force vs angle - 4 NORMALIZED');
                figure('Name',plottitle)
                plot(BD_angle_vars_norm_mean(:,1), BD_angle_vars_norm_mean(:,2),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_norm_mean(:,1), CON_angle_vars_norm_mean(:,2),'b','LineWidth',2)
                axis(axis_PP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Force (% of ind max)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('force vs angle - 5 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars_norm{1,i}(:,1),BD_angle_vars_norm{1,i}(:,2))
                end
                axis(axis_PP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Force (% of ind max)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('force vs angle - 6 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars_norm{1,i}(:,1),CON_angle_vars_norm{1,i}(:,2))
                end
                axis(axis_PP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Force (% of ind max)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            
            %% TORQUE-angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('torque vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,18),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,18),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_torque_mean, BD_torque_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_torque_mean, CON_torque_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_torque_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_torque_mean, CON_ROM_SD, '*b')
                axis(axis_torque)
                xlabel('Gonio angle (°)')
                ylabel('Torque (Nm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('torque vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,18))
                end
                axis(axis_torque)
                xlabel('Gonio angle (°)')
                ylabel('Torque (Nm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('torque vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,18))
                end
                axis(axis_torque)
                xlabel('Gonio angle (°)')
                ylabel('Torque (Nm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('torque vs angle - 4 NORMALIZED');
                figure('Name',plottitle)
                plot(BD_angle_vars_norm_mean(:,1), BD_angle_vars_norm_mean(:,18),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_norm_mean(:,1), CON_angle_vars_norm_mean(:,18),'b','LineWidth',2)
                axis(axis_PP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Torque (% of ind max)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('torque vs angle - 5 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars_norm{1,i}(:,1),BD_angle_vars_norm{1,i}(:,18))
                end
                axis(axis_PP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Torque (% of ind max)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('torque vs angle - 6 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars_norm{1,i}(:,1),CON_angle_vars_norm{1,i}(:,18))
                end
                axis(axis_PP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Torque (% of ind max)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            
            %% FREE AT: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('free AT length vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,12),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,12),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_L_at_SOL_mean, BD_L_at_SOL_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_L_at_SOL_mean, CON_L_at_SOL_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_L_at_SOL_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_L_at_SOL_mean, CON_ROM_SD, '*b')
                axis(axis_len_AT)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Southeast')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('free AT length vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,12))
                end
                axis(axis_len_AT)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('free AT length vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,12))
                end
                axis(axis_len_AT)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('free AT strain vs angle - 4 NORMALIZED');
                figure('Name',plottitle)
                plot(BD_angle_vars_norm_mean(:,1), BD_angle_vars_norm_mean(:,12),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_norm_mean(:,1), CON_angle_vars_norm_mean(:,12),'b','LineWidth',2)
                axis(axis_str_ATP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('free AT strain vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,21),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,21),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_strain_at_SOL_mean, BD_strain_at_SOL_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_strain_at_SOL_mean, CON_strain_at_SOL_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_strain_at_SOL_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_strain_at_SOL_mean, CON_ROM_SD, '*b')
                axis(axis_str_AT)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('free AT strain vs angle - 5 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars_norm{1,i}(:,1),BD_angle_vars_norm{1,i}(:,12))
                end
                axis(axis_str_ATP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('free AT strain vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars_norm_indlength{1,i}(:,12))
                end
                axis(axis_str_AT)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('free AT strain vs angle - 6 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars_norm{1,i}(:,1),CON_angle_vars_norm{1,i}(:,12))
                end
                axis(axis_str_ATP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('free AT strain vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars_norm_indlength{1,i}(:,12))
                end
                axis(axis_str_AT)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            
            %% FREE AT: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('free AT elongation vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,6),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,6),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_elong_AT_mean, BD_elong_AT_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_elong_AT_mean, CON_elong_AT_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_elong_AT_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_elong_AT_mean, CON_ROM_SD, '*b')
                axis(axis_el_AT)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('free AT elongation vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,6))
                end
                axis(axis_el_AT)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('free AT elongation vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,6))
                end
                axis(axis_el_AT)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            
            %% GM TENDON: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('GM tendon (from calc) length vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,13),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,13),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_L_at_GM_mean, BD_L_at_GM_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_L_at_GM_mean, CON_L_at_GM_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_L_at_GM_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_L_at_GM_mean, CON_ROM_SD, '*b')
                axis(axis_len_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Southeast')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('GM tendon (from calc) length vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,13))
                end
                axis(axis_len_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('GM tendon (from calc) length vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,13))
                end
                axis(axis_len_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('GM tendon (from calc) strain vs angle - 4 NORMALIZED');
                figure('Name',plottitle)
                plot(BD_angle_vars_norm_mean(:,1), BD_angle_vars_norm_mean(:,13),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_norm_mean(:,1), CON_angle_vars_norm_mean(:,13),'b','LineWidth',2)
                axis(axis_str_GMtendP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('GM tendon (from calc) strain vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,22),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,22),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_strain_at_GM_mean, BD_strain_at_GM_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_strain_at_GM_mean, CON_strain_at_GM_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_strain_at_GM_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_strain_at_GM_mean, CON_ROM_SD, '*b')
                axis(axis_str_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('GM tendon (from calc) strain vs angle - 5 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars_norm{1,i}(:,1),BD_angle_vars_norm{1,i}(:,13))
                end
                axis(axis_str_GMtendP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('GM tendon (from calc) strain vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars_norm_indlength{1,i}(:,13))
                end
                axis(axis_str_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('GM tendon (from calc) strain vs angle - 6 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars_norm{1,i}(:,1),CON_angle_vars_norm{1,i}(:,13))
                end
                axis(axis_str_GMtendP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('GM tendon (from calc) strain vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars_norm_indlength{1,i}(:,13))
                end
                axis(axis_str_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
                       
            %% GM TENDON: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('GM tendon (from calc) elongation vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,7),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,7),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_elong_GMtend_mean, BD_elong_GMtend_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_elong_GMtend_mean, CON_elong_GMtend_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_elong_GMtend_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_elong_GMtend_mean, CON_ROM_SD, '*b')
                axis(axis_el_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('GM tendon (from calc) elongation vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,7))
                end
                axis(axis_el_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('GM tendon (from calc) elongation vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,7))
                end
                axis(axis_el_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            
            %% GM APO: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) length vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,16),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,16),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_L_GMapo_mean, BD_L_GMapo_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_L_GMapo_mean, CON_L_GMapo_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_L_GMapo_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_L_GMapo_mean, CON_ROM_SD, '*b')
                axis(axis_len_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Southeast')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) length vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,16))
                end
                axis(axis_len_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) length vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,16))
                end
                axis(axis_len_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) strain vs angle - 4 NORMALIZED');
                figure('Name',plottitle)
                plot(BD_angle_vars_norm_mean(:,1), BD_angle_vars_norm_mean(:,16),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_norm_mean(:,1), CON_angle_vars_norm_mean(:,16),'b','LineWidth',2)
                axis(axis_str_GMapoP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) strain vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,25),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,25),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_strain_GMapo_mean, BD_strain_GMapo_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_strain_GMapo_mean, CON_strain_GMapo_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_strain_GMapo_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_strain_GMapo_mean, CON_ROM_SD, '*b')
                axis(axis_str_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            if BD_count > 1 && plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) strain vs angle - 5 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars_norm{1,i}(:,1),BD_angle_vars_norm{1,i}(:,16))
                end
                axis(axis_str_GMapoP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) strain vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars_norm_indlength{1,i}(:,16))
                end
                axis(axis_str_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) strain vs angle - 6 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars_norm{1,i}(:,1),CON_angle_vars_norm{1,i}(:,16))
                end
                axis(axis_str_GMapoP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) strain vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars_norm_indlength{1,i}(:,16))
                end
                axis(axis_str_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            
            %% GM APO: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) elongation vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,10),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,10),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_elong_GMapo_mean, BD_elong_GMapo_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_elong_GMapo_mean, CON_elong_GMapo_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_elong_GMapo_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_elong_GMapo_mean, CON_ROM_SD, '*b')
                axis(axis_el_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) elongation vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,10))
                end
                axis(axis_el_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) elongation vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,10))
                end
                axis(axis_el_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
                        
            %% GM MUSCLE portion: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('muscle GM (ins-knee) length vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,17),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,17),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_L_msc_GM_mean, BD_L_msc_GM_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_L_msc_GM_mean, CON_L_msc_GM_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_L_msc_GM_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_L_msc_GM_mean, CON_ROM_SD, '*b')
                axis(axis_len_GM)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Southeast')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('muscle GM (ins-knee) length vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,17))
                end
                axis(axis_len_GM)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('muscle GM (ins-knee) length vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,17))
                end
                axis(axis_len_GM)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('muscle GM (ins-knee) strain vs angle - 4 NORMALIZED');
                figure('Name',plottitle)
                plot(BD_angle_vars_norm_mean(:,1), BD_angle_vars_norm_mean(:,17),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_norm_mean(:,1), CON_angle_vars_norm_mean(:,17),'b','LineWidth',2)
                axis(axis_str_GMP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('muscle GM (ins-knee) strain vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,26),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,26),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_strain_msc_GM_mean, BD_strain_msc_GM_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_strain_msc_GM_mean, CON_strain_msc_GM_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_strain_msc_GM_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_strain_msc_GM_mean, CON_ROM_SD, '*b')
                axis(axis_str_GM)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            if BD_count > 1 && plot_check
                plottitle = horzcat('muscle GM (ins-knee) strain vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars_norm{1,i}(:,1),BD_angle_vars_norm{1,i}(:,17))
                end
                axis(axis_str_GMP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('muscle GM (ins-knee) strain vs angle - 5 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars_norm_indlength{1,i}(:,17))
                end
                axis(axis_str_GM)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('muscle GM (ins-knee) strain vs angle - 6 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars_norm{1,i}(:,1),CON_angle_vars_norm{1,i}(:,17))
                end
                axis(axis_str_GMP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('muscle GM (ins-knee) strain vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars_norm_indlength{1,i}(:,17))
                end
                axis(axis_str_GM)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            
            %% GM MUSCLE portion: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('muscle GM (ins-knee) elongation vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,11),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,11),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_elong_msc_GM_mean, BD_elong_msc_GM_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_elong_msc_GM_mean, CON_elong_msc_GM_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_elong_msc_GM_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_elong_msc_GM_mean, CON_ROM_SD, '*b')
                axis(axis_el_GM)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('muscle GM (ins-knee) elongation vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,11))
                end
                axis(axis_el_GM)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('muscle GM (ins-knee) elongation vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,11))
                end
                axis(axis_el_GM)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            

            %% GM MUSCLE Lichtwark/Fukunaga: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('muscle GM (architecture) elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,28),'r','LineStyle','-','LineWidth',2)
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,28),'b','LineStyle','-','LineWidth',2)
                herrorbar(BD_ROM_mean, BD_elong_msc_GM_licht_mean, BD_ROM_SD, '*r')
                errorbar(BD_ROM_mean, BD_elong_msc_GM_licht_mean, BD_elong_msc_GM_licht_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                herrorbar(CON_ROM_mean, CON_elong_msc_GM_licht_mean, CON_ROM_SD, '*b')
                errorbar(CON_ROM_mean, CON_elong_msc_GM_licht_mean, CON_elong_msc_GM_licht_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                axis(axis_el_GM)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('muscle GM (architecture) elongation vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,28))
                end
                axis(axis_el_GM)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('muscle GM (architecture) elongation vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,28))
                end
                axis(axis_el_GM)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            
            %% GM TENDON Lichtwark/Fukunaga: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('GM tendon (architecture) elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,29),'r','LineStyle','-','LineWidth',2)
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,29),'b','LineStyle','-','LineWidth',2)
                herrorbar(BD_ROM_mean, BD_elong_tend_GM_licht_mean, BD_ROM_SD, '*r')
                errorbar(BD_ROM_mean, BD_elong_tend_GM_licht_mean, BD_elong_tend_GM_licht_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                herrorbar(CON_ROM_mean, CON_elong_tend_GM_licht_mean, CON_ROM_SD, '*b')
                errorbar(CON_ROM_mean, CON_elong_tend_GM_licht_mean, CON_elong_tend_GM_licht_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                axis(axis_el_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('GM tendon (architecture) elongation vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,29))
                end
                axis(axis_el_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('GM tendon (architecture) elongation vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,29))
                end
                axis(axis_el_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end

            
            %% GM MUSCLE Lichtwark/Fukunaga: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('muscle GM (architecture) strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,30),'r','LineStyle','-','LineWidth',2)
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,30),'b','LineStyle','-','LineWidth',2)
                herrorbar(BD_ROM_mean, BD_strain_msc_GM_licht_mean, BD_ROM_SD, '*r')
                errorbar(BD_ROM_mean, BD_strain_msc_GM_licht_mean, BD_strain_msc_GM_licht_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                herrorbar(CON_ROM_mean, CON_strain_msc_GM_licht_mean, CON_ROM_SD, '*b')
                errorbar(CON_ROM_mean, CON_strain_msc_GM_licht_mean, CON_strain_msc_GM_licht_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                axis(axis_str_GM)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('muscle GM (architecture) strain vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,30))
                end
                axis(axis_str_GM)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('muscle GM (architecture) strain vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,30))
                end
                axis(axis_str_GM)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            
            %% GM TENDON Lichtwark/Fukunaga: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('GM tendon (architecture) strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,31),'r','LineStyle','-','LineWidth',2)
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,31),'b','LineStyle','-','LineWidth',2)
                herrorbar(BD_ROM_mean, BD_strain_tend_GM_licht_mean, BD_ROM_SD, '*r')
                errorbar(BD_ROM_mean, BD_strain_tend_GM_licht_mean, BD_strain_tend_GM_licht_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                herrorbar(CON_ROM_mean, CON_strain_tend_GM_licht_mean, CON_ROM_SD, '*b')
                errorbar(CON_ROM_mean, CON_strain_tend_GM_licht_mean, CON_strain_tend_GM_licht_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                axis(axis_str_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('GM tendon (architecture) strain vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,31))
                end
                axis(axis_str_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('GM tendon (architecture) strain vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,31))
                end
                axis(axis_str_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end

            
            %% GM fascicle length (Lichtwark): Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('GM fascicle length vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,32),'r','LineStyle','-','LineWidth',2)
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,32),'b','LineStyle','-','LineWidth',2)
                herrorbar(BD_ROM_mean, BD_GM_length_faslen_licht_mean, BD_ROM_SD, '*r')
                errorbar(BD_ROM_mean, BD_GM_length_faslen_licht_mean, BD_GM_length_faslen_licht_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                herrorbar(CON_ROM_mean, CON_GM_length_faslen_licht_mean, CON_ROM_SD, '*b')
                errorbar(CON_ROM_mean, CON_GM_length_faslen_licht_mean, CON_GM_length_faslen_licht_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                axis(axis_faslen)
                xlabel('Gonio angle (°)')
                ylabel('Fascicle length (mm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('GM fascicle length vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,32))
                end
                axis(axis_faslen)
                xlabel('Gonio angle (°)')
                ylabel('Fascicle length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('GM fascicle length vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,32))
                end
                axis(axis_faslen)
                xlabel('Gonio angle (°)')
                ylabel('Fascicle length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end

            
            %% GM fascicle length (Lichtwark): Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('GM fascicle elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,34),'r','LineStyle','-','LineWidth',2)
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,34),'b','LineStyle','-','LineWidth',2)
                herrorbar(BD_ROM_mean, BD_GM_elong_faslen_licht_mean, BD_ROM_SD, '*r')
                errorbar(BD_ROM_mean, BD_GM_elong_faslen_licht_mean, BD_GM_elong_faslen_licht_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                herrorbar(CON_ROM_mean, CON_GM_elong_faslen_licht_mean, CON_ROM_SD, '*b')
                errorbar(CON_ROM_mean, CON_GM_elong_faslen_licht_mean, CON_GM_elong_faslen_licht_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                axis(axis_faslen_elong)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('GM fascicle elongation vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,34))
                end
                axis(axis_faslen_elong)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('GM fascicle elongation vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,34))
                end
                axis(axis_faslen_elong)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end

            
            %% GM fascicle length (Lichtwark): Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('GM fascicle strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,35),'r','LineStyle','-','LineWidth',2)
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,35),'b','LineStyle','-','LineWidth',2)
                herrorbar(BD_ROM_mean, BD_GM_strain_faslen_licht_mean, BD_ROM_SD, '*r')
                errorbar(BD_ROM_mean, BD_GM_strain_faslen_licht_mean, BD_GM_strain_faslen_licht_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                herrorbar(CON_ROM_mean, CON_GM_strain_faslen_licht_mean, CON_ROM_SD, '*b')
                errorbar(CON_ROM_mean, CON_GM_strain_faslen_licht_mean, CON_GM_strain_faslen_licht_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                axis(axis_faslen_str)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('GM fascicle strain vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,35))
                end
                axis(axis_faslen_str)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('GM fascicle strain vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,35))
                end
                axis(axis_faslen_str)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end


            %% GM pennation angle (Lichtwark) vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('GM pennation angle vs ankle angle - 1');
                figure('Name',plottitle)
                hold on
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,33),'r','LineStyle','-','LineWidth',2)
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,33),'b','LineStyle','-','LineWidth',2)
                herrorbar(BD_ROM_mean, BD_GM_pennation_licht_mean, BD_ROM_SD, '*r')
                errorbar(BD_ROM_mean, BD_GM_pennation_licht_mean, BD_GM_pennation_licht_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                herrorbar(CON_ROM_mean, CON_GM_pennation_licht_mean, CON_ROM_SD, '*b')
                errorbar(CON_ROM_mean, CON_GM_pennation_licht_mean, CON_GM_pennation_licht_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                axis(axis_pennation)
                xlabel('Gonio angle (°)')
                ylabel('Pennation angle (°)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('GM pennation angle vs ankle angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,33))
                end
                axis(axis_pennation)
                xlabel('Gonio angle (°)')
                ylabel('Pennation angle (°)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('GM pennation angle vs ankle angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,33))
                end
                axis(axis_pennation)
                xlabel('Gonio angle (°)')
                ylabel('Pennation angle (°)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            
            %% SOL MUSCLE portion: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('muscle SOL length vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,20),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,20),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_L_msc_SOL_mean, BD_L_msc_SOL_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_L_msc_SOL_mean, CON_L_msc_SOL_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_L_msc_SOL_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_L_msc_SOL_mean, CON_ROM_SD, '*b')
                axis(axis_len_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Southeast')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('muscle SOL length vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,20))
                end
                axis(axis_len_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('muscle SOL length vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,20))
                end
                axis(axis_len_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('muscle SOL strain vs angle - 4 NORMALIZED');
                figure('Name',plottitle)
                plot(BD_angle_vars_norm_mean(:,1), BD_angle_vars_norm_mean(:,20),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_norm_mean(:,1), CON_angle_vars_norm_mean(:,20),'b','LineWidth',2)
                axis(axis_str_SOLP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('muscle SOL strain vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,27),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,27),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_strain_msc_SOL_mean, BD_strain_msc_SOL_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_strain_msc_SOL_mean, CON_strain_msc_SOL_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_strain_msc_SOL_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_strain_msc_SOL_mean, CON_ROM_SD, '*b')
                axis(axis_str_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('muscle SOL strain vs angle - 5 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars_norm{1,i}(:,1),BD_angle_vars_norm{1,i}(:,20))
                end
                axis(axis_str_SOLP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('muscle SOL strain vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars_norm_indlength{1,i}(:,20))
                end
                axis(axis_str_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('muscle SOL strain vs angle - 6 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars_norm{1,i}(:,1),CON_angle_vars_norm{1,i}(:,20))
                end
                axis(axis_str_SOLP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('muscle SOL strain vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars_norm_indlength{1,i}(:,20))
                end
                axis(axis_str_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            
            %% SOL MUSCLE portion: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('muscle SOL elongation vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,19),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,19),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_elong_msc_SOL_mean, BD_elong_msc_SOL_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_elong_msc_SOL_mean, CON_elong_msc_SOL_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_elong_msc_SOL_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_elong_msc_SOL_mean, CON_ROM_SD, '*b')
                axis(axis_el_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('muscle SOL elongation vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,19))
                end
                axis(axis_el_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('muscle SOL elongation vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,19))
                end
                axis(axis_el_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            
            %% Full MTU: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('MTU length vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,14),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,14),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_L_MTU_mean, BD_L_MTU_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_L_MTU_mean, CON_L_MTU_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_L_MTU_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_L_MTU_mean, CON_ROM_SD, '*b')
                axis(axis_len_MTU)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Southeast')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('MTU length vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,14))
                end
                axis(axis_len_MTU)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('MTU length vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,14))
                end
                axis(axis_len_MTU)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            % commenting out MTU strain plots since they simply show the Grieve calculation factor
%             if BD_count > 1 && CON_count > 1 && plot_check
%                 plottitle = horzcat('MTU strain vs angle - 4 NORMALIZED');
%                 figure('Name',plottitle)
%                 plot(BD_angle_vars_norm_mean(:,1), BD_angle_vars_norm_mean(:,14),'r','LineWidth',2)
%                 hold on
%                 plot(CON_angle_vars_norm_mean(:,1), CON_angle_vars_norm_mean(:,14),'b','LineWidth',2)
%                 axis(axis_str_MTUP)
%                 xlabel('Gonio angle (% of ind max)')
%                 ylabel('Strain (% of initial length)')
%                 title(plottitle)
%                 legend('Dancer avg', 'Control avg','Location','Northwest')
%                 saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
%             end
%             if BD_count > 1 && CON_count > 1 && plot_check
%                 plottitle = horzcat('MTU strain vs angle - 1');
%                 figure('Name',plottitle)
%                 plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,23),'r','LineWidth',2)
%                 hold on
%                 plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,23),'b','LineWidth',2)
%                 errorbar(BD_ROM_mean, BD_strain_MTU_mean, BD_strain_MTU_SD, '*r', 'MarkerFaceColor', 'r')
%                 errorbar(CON_ROM_mean, CON_strain_MTU_mean, CON_strain_MTU_SD, '*b', 'MarkerFaceColor', 'b')
%                 herrorbar(BD_ROM_mean, BD_strain_MTU_mean, BD_ROM_SD, '*r')
%                 herrorbar(CON_ROM_mean, CON_strain_MTU_mean, CON_ROM_SD, '*b')
%                 axis(axis_str_MTU)
%                 xlabel('Gonio angle (°)')
%                 ylabel('Strain (% of initial length)')
%                 title(plottitle)
%                 legend('Dancer avg', 'Control avg','Location','Northwest')
%                 saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
%             end
%             if BD_count > 1 && plot_check
%                 plottitle = horzcat('MTU strain vs angle - 5 ind curves dancers');
%                 figure('Name',plottitle)
%                 hold on
%                 for i = 1:BD_count
%                     plot(BD_angle_vars_norm{1,i}(:,1),BD_angle_vars_norm{1,i}(:,14))
%                 end
%                 axis(axis_str_MTUP)
%                 xlabel('Gonio angle (% of ind max)')
%                 ylabel('Strain (% of initial length)')
%                 title(plottitle)
%                 %legend
%                 saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
%             end
%             if BD_count > 1 && plot_check
%                 plottitle = horzcat('MTU strain vs angle - 2 ind curves dancers');
%                 figure('Name',plottitle)
%                 hold on
%                 for i = 1:BD_count
%                     plot(BD_angle_vars{1,i}(:,1),BD_angle_vars_norm_indlength{1,i}(:,14))
%                 end
%                 axis(axis_str_MTU)
%                 xlabel('Gonio angle (°)')
%                 ylabel('Strain (% of initial length)')
%                 title(plottitle)
%                 %legend
%                 saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
%             end
%             if CON_count > 1 && plot_check
%                 plottitle = horzcat('MTU strain vs angle - 6 ind curves controls');
%                 figure('Name',plottitle)
%                 hold on
%                 for i = 1:CON_count
%                     plot(CON_angle_vars_norm{1,i}(:,1),CON_angle_vars_norm{1,i}(:,14))
%                 end
%                 axis(axis_str_MTUP)
%                 xlabel('Gonio angle (% of ind max)')
%                 ylabel('Strain (% of initial length)')
%                 title(plottitle)
%                 %legend
%                 saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
%             end
%             if CON_count > 1 && plot_check
%                 plottitle = horzcat('MTU strain vs angle - 3 ind curves controls');
%                 figure('Name',plottitle)
%                 hold on
%                 for i = 1:CON_count
%                     plot(CON_angle_vars{1,i}(:,1),CON_angle_vars_norm_indlength{1,i}(:,14))
%                 end
%                 axis(axis_str_MTU)
%                 xlabel('Gonio angle (°)')
%                 ylabel('Strain (% of initial length)')
%                 title(plottitle)
%                 %legend
%                 saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
%             end
            
            
            %% Full MTU: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('MTU elongation vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,8),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,8),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_elong_MTU_mean, BD_elong_MTU_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_elong_MTU_mean, CON_elong_MTU_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_elong_MTU_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_elong_MTU_mean, CON_ROM_SD, '*b')
                axis(axis_el_MTU)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('MTU elongation vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,8))
                end
                axis(axis_el_MTU)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('MTU elongation vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,8))
                end
                axis(axis_el_MTU)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            
            %% GMFAS scans: Displacement vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('displacement GMFAS vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,9),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,9),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_displ_GMFAS_mean, BD_displ_GMFAS_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_displ_GMFAS_mean, CON_displ_GMFAS_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_displ_GMFAS_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_displ_GMFAS_mean, CON_ROM_SD, '*b')
                axis(axis_displ_GMFAS)
                xlabel('Gonio angle (°)')
                ylabel('Displacement (mm)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('displacement GMFAS vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,9))
                end
                axis(axis_displ_GMFAS)
                xlabel('Gonio angle (°)')
                ylabel('Displacement (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('displacement GMFAS vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,9))
                end
                axis(axis_displ_GMFAS)
                xlabel('Gonio angle (°)')
                ylabel('Displacement (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end

            
            %% Length +- SD, 3 MTU components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('length ALL MTU components vs angle');
                figure('Name',plottitle)
                hold on
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,14),'r','LineWidth',2) % full MTU - calc-knee
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,13),'r','LineWidth',2) % GM tend  - calc-GM ins
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,12),'r','LineWidth',2) % free AT  - calc-SOL ins
                
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,14),'b','LineWidth',2)
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,13),'b','LineWidth',2)
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,12),'b','LineWidth',2)
                
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,14)+BD_angle_vars_SD(:,14),'r','LineWidth',0.25)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,14)-BD_angle_vars_SD(:,14),'r','LineWidth',0.25)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,13)+BD_angle_vars_SD(:,13),'r','LineWidth',0.25)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,13)-BD_angle_vars_SD(:,13),'r','LineWidth',0.25)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,12)+BD_angle_vars_SD(:,12),'r','LineWidth',0.25)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,12)-BD_angle_vars_SD(:,12),'r','LineWidth',0.25)
                
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,14)+CON_angle_vars_SD(:,14),'b','LineWidth',0.25)
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,14)-CON_angle_vars_SD(:,14),'b','LineWidth',0.25)
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,13)+CON_angle_vars_SD(:,13),'b','LineWidth',0.25)
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,13)-CON_angle_vars_SD(:,13),'b','LineWidth',0.25)
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,12)+CON_angle_vars_SD(:,12),'b','LineWidth',0.25)
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,12)-CON_angle_vars_SD(:,12),'b','LineWidth',0.25)
                
                axis(axis_len_MTU)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                legend('Dancer MTU', 'Dancer GM tend','Dancer free AT','Control MTU','Control GM tend','Control free AT','Location','Southeast')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            
            %% EMG vs angle GM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('EMG gas.med. vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,3),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,3),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_EMG_gm_mean, BD_EMG_gm_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_EMG_gm_mean, CON_EMG_gm_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_EMG_gm_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_EMG_gm_mean, CON_ROM_SD, '*b')
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('EMG gas.med. vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,3))
                end
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('EMG gas.med. vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,3))
                end
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            
            %% EMG vs angle GL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('EMG gas.lat. vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,4),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,4),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_EMG_gl_mean, BD_EMG_gl_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_EMG_gl_mean, CON_EMG_gl_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_EMG_gl_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_EMG_gl_mean, CON_ROM_SD, '*b')
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('EMG gas.lat. vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,4))
                end
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('EMG gas.lat. vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,4))
                end
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
            
            %% EMG vs angle SOL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if BD_count > 1 && CON_count > 1 && plot_check
                plottitle = horzcat('EMG soleus vs angle - 1');
                figure('Name',plottitle)
                plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,5),'r','LineWidth',2)
                hold on
                plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,5),'b','LineWidth',2)
                errorbar(BD_ROM_mean, BD_EMG_sol_mean, BD_EMG_sol_SD, '*r', 'MarkerFaceColor', 'r')
                errorbar(CON_ROM_mean, CON_EMG_sol_mean, CON_EMG_sol_SD, '*b', 'MarkerFaceColor', 'b')
                herrorbar(BD_ROM_mean, BD_EMG_sol_mean, BD_ROM_SD, '*r')
                herrorbar(CON_ROM_mean, CON_EMG_sol_mean, CON_ROM_SD, '*b')
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if BD_count > 1 && plot_check
                plottitle = horzcat('EMG soleus vs angle - 2 ind curves dancers');
                figure('Name',plottitle)
                hold on
                for i = 1:BD_count
                    plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,5))
                end
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            if CON_count > 1 && plot_check
                plottitle = horzcat('EMG soleus vs angle - 3 ind curves controls');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_count
                    plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,5))
                end
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_BD ',plottitle,'.jpg'))
            end
            
    elseif input_project == 2
        %% intervention study
        % MMM LATER? add figures showing US displacements
            %% FORCE-angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('force vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,2),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,2),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,2),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,2),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_F_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_F_mean, STR_PRE_F_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_F_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_F_mean, STR_POST_F_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_F_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_F_mean, CON_PRE_F_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_F_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_F_mean, CON_POST_F_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_force)
                xlabel('Gonio angle (°)')
                ylabel('Force (N)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('force vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,2),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,2))
                end
                axis(axis_force)
                xlabel('Gonio angle (°)')
                ylabel('Force (N)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('force vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,2),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,2))
                end
                axis(axis_force)
                xlabel('Gonio angle (°)')
                ylabel('Force (N)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('force vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,2),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,2),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,2),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,2),'b','LineStyle','-','LineWidth',1)
                    axis(axis_force)
                    xlabel('Gonio angle (°)')
                    ylabel('Force (N)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end

            
            
            % MMM LATER? other type of normalization - normalize to what? add plots with individual data normalised (same as #2 and #3)?
            if plot_check
                plottitle = horzcat('force vs angle - 4 NORMALIZED');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_norm_mean(:,1), STR_PRE_angle_vars_norm_mean(:,2),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_norm_mean(:,1), STR_POST_angle_vars_norm_mean(:,2),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_norm_mean(:,1), CON_PRE_angle_vars_norm_mean(:,2),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_norm_mean(:,1), CON_POST_angle_vars_norm_mean(:,2),'b','LineStyle','-','LineWidth',1)
                
                axis(axis_PP)
                xlabel('Gonio angle (% of ind max)')
                ylabel('Force (% of ind max)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            %% TORQUE-angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('torque vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,18),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,18),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,18),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,18),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_torque_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_torque_mean, STR_PRE_torque_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_torque_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_torque_mean, STR_POST_torque_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_torque_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_torque_mean, CON_PRE_torque_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_torque_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_torque_mean, CON_POST_torque_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_torque)
                xlabel('Gonio angle (°)')
                ylabel('Torque (Nm)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('torque vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,18),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,18))
                end
                axis(axis_torque)
                xlabel('Gonio angle (°)')
                ylabel('Torque (Nm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('torque vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,18),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,18))
                end
                axis(axis_torque)
                xlabel('Gonio angle (°)')
                ylabel('Torque (Nm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND torque vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,18),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,18),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,18),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,18),'b','LineStyle','-','LineWidth',1)
                axis(axis_torque)
                xlabel('Gonio angle (°)')
                ylabel('Torque (Nm)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% FREE AT: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('free AT length vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,12),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,12),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,12),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,12),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_L_at_SOL_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_L_at_SOL_mean, STR_PRE_L_at_SOL_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_L_at_SOL_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_L_at_SOL_mean, STR_POST_L_at_SOL_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_L_at_SOL_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_L_at_SOL_mean, CON_PRE_L_at_SOL_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_L_at_SOL_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_L_at_SOL_mean, CON_POST_L_at_SOL_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_len_AT)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('free AT length vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,12),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,12))
                end
                axis(axis_len_AT)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('free AT length vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,12),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,12))
                end
                axis(axis_len_AT)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND free AT length vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,12),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,12),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,12),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,12),'b','LineStyle','-','LineWidth',1)
                    axis(axis_len_AT)
                    xlabel('Gonio angle (°)')
                    ylabel('Length (mm)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% FREE AT: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('free AT elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,6),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,6),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,6),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,6),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_AT_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_AT_mean, STR_PRE_elong_AT_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_AT_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_AT_mean, STR_POST_elong_AT_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_AT_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_AT_mean, CON_PRE_elong_AT_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_AT_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_AT_mean, CON_POST_elong_AT_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_el_AT)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('free AT elongation vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,6),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,6))
                end
                axis(axis_el_AT)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('free AT elongation vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,6),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,6))
                end
                axis(axis_el_AT)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND free AT elongation vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,6),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,6),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,6),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,6),'b','LineStyle','-','LineWidth',1)
                    axis(axis_el_AT)
                    xlabel('Gonio angle (°)')
                    ylabel('Elongation (mm)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end

            %% FREE AT: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('free AT strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,21),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,21),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,21),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,21),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_strain_at_SOL_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_strain_at_SOL_mean, STR_PRE_strain_at_SOL_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_strain_at_SOL_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_strain_at_SOL_mean, STR_POST_strain_at_SOL_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_strain_at_SOL_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_strain_at_SOL_mean, CON_PRE_strain_at_SOL_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_strain_at_SOL_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_strain_at_SOL_mean, CON_POST_strain_at_SOL_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_str_AT)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('free AT strain vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,21),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,21))
                end
                axis(axis_str_AT)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('free AT strain vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,21),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,21))
                end
                axis(axis_str_AT)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND free AT strain vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,21),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,21),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,21),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,21),'b','LineStyle','-','LineWidth',1)
                axis(axis_str_AT)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% GM TENDON: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM tendon (from calc) length vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,13),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,13),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,13),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,13),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_L_at_GM_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_L_at_GM_mean, STR_PRE_L_at_GM_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_L_at_GM_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_L_at_GM_mean, STR_POST_L_at_GM_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_L_at_GM_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_L_at_GM_mean, CON_PRE_L_at_GM_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_L_at_GM_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_L_at_GM_mean, CON_POST_L_at_GM_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_len_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southeast')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('GM tendon (from calc) length vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,13),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,13))
                end
                axis(axis_len_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('GM tendon (from calc) length vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,13),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,13))
                end
                axis(axis_len_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM tendon (from calc) length vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,13),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,13),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,13),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,13),'b','LineStyle','-','LineWidth',1)
                axis(axis_len_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end

            %% GM TENDON: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM tendon (from calc) elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,7),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,7),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,7),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,7),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_GMtend_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_GMtend_mean, STR_PRE_elong_GMtend_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_GMtend_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_GMtend_mean, STR_POST_elong_GMtend_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_GMtend_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_GMtend_mean, CON_PRE_elong_GMtend_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_GMtend_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_GMtend_mean, CON_POST_elong_GMtend_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_el_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('GM tendon (from calc) elongation vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,7),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,7))
                end
                axis(axis_el_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('GM tendon (from calc) elongation vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,7),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,7))
                end
                axis(axis_el_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM tendon (from calc) elongation vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,7),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,7),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,7),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,7),'b','LineStyle','-','LineWidth',1)
                axis(axis_el_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end

            %% GM TENDON: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM tendon (from calc) strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,22),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,22),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,22),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,22),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_strain_at_GM_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_strain_at_GM_mean, STR_PRE_strain_at_GM_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_strain_at_GM_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_strain_at_GM_mean, STR_POST_strain_at_GM_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_strain_at_GM_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_strain_at_GM_mean, CON_PRE_strain_at_GM_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_strain_at_GM_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_strain_at_GM_mean, CON_POST_strain_at_GM_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_str_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('GM tendon (from calc) strain vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,22),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,22))
                end
                axis(axis_str_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('GM tendon (from calc) strain vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,22),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,22))
                end
                axis(axis_str_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM tendon (from calc) strain vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,22),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,22),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,22),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,22),'b','LineStyle','-','LineWidth',1)
                axis(axis_str_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% GM APO: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) length vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,16),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,16),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,16),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,16),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_L_GMapo_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_L_GMapo_mean, STR_PRE_L_GMapo_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_L_GMapo_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_L_GMapo_mean, STR_POST_L_GMapo_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_L_GMapo_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_L_GMapo_mean, CON_PRE_L_GMapo_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_L_GMapo_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_L_GMapo_mean, CON_POST_L_GMapo_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_len_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southeast')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) length vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,16),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,16))
                end
                axis(axis_len_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) length vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,16),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,16))
                end
                axis(axis_len_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM apo (SOL ins-GM ins) length vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,16),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,16),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,16),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,16),'b','LineStyle','-','LineWidth',1)
                axis(axis_len_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% GM APO: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,10),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,10),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,10),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,10),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_GMapo_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_GMapo_mean, STR_PRE_elong_GMapo_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_GMapo_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_GMapo_mean, STR_POST_elong_GMapo_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_GMapo_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_GMapo_mean, CON_PRE_elong_GMapo_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_GMapo_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_GMapo_mean, CON_POST_elong_GMapo_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_el_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northeast')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) elongation vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,10),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,10))
                end
                axis(axis_el_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) elongation vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,10),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,10))
                end
                axis(axis_el_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM apo (SOL ins-GM ins) elongation vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,10),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,10),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,10),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,10),'b','LineStyle','-','LineWidth',1)
                axis(axis_el_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% GM APO: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,25),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,25),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,25),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,25),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_strain_GMapo_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_strain_GMapo_mean, STR_PRE_strain_GMapo_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_strain_GMapo_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_strain_GMapo_mean, STR_POST_strain_GMapo_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_strain_GMapo_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_strain_GMapo_mean, CON_PRE_strain_GMapo_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_strain_GMapo_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_strain_GMapo_mean, CON_POST_strain_GMapo_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_str_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northeast')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) strain vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,25),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,25))
                end
                axis(axis_str_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('GM apo (SOL ins-GM ins) strain vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,25),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,25))
                end
                axis(axis_str_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM apo (SOL ins-GM ins) strain vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,25),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,25),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,25),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,25),'b','LineStyle','-','LineWidth',1)
                axis(axis_str_GMapo)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
                        
            %% GM MUSCLE portion: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,17),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,17),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,17),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,17),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_L_msc_GM_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_L_msc_GM_mean, STR_PRE_L_msc_GM_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_L_msc_GM_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_L_msc_GM_mean, STR_POST_L_msc_GM_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_L_msc_GM_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_L_msc_GM_mean, CON_PRE_L_msc_GM_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_L_msc_GM_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_L_msc_GM_mean, CON_POST_L_msc_GM_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_len_GM)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southeast')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,17),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,17))
                end
                axis(axis_len_GM)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,17),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,17))
                end
                axis(axis_len_GM)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle GM (ins-knee) vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,17),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,17),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,17),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,17),'b','LineStyle','-','LineWidth',1)
                axis(axis_len_GM)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% GM MUSCLE portion: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,11),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,11),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,11),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,11),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_msc_GM_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_msc_GM_mean, STR_PRE_elong_msc_GM_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_msc_GM_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_msc_GM_mean, STR_POST_elong_msc_GM_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_msc_GM_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_msc_GM_mean, CON_PRE_elong_msc_GM_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_msc_GM_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_msc_GM_mean, CON_POST_elong_msc_GM_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_el_GM)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) elongation vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,11),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,11))
                end
                axis(axis_el_GM)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) elongation vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,11),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,11))
                end
                axis(axis_el_GM)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle GM (ins-knee) elongation vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,11),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,11),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,11),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,11),'b','LineStyle','-','LineWidth',1)
                axis(axis_el_GM)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% GM MUSCLE portion: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,26),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,26),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,26),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,26),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_strain_msc_GM_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_strain_msc_GM_mean, STR_PRE_strain_msc_GM_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_strain_msc_GM_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_strain_msc_GM_mean, STR_POST_strain_msc_GM_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_strain_msc_GM_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_strain_msc_GM_mean, CON_PRE_strain_msc_GM_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_strain_msc_GM_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_strain_msc_GM_mean, CON_POST_strain_msc_GM_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_str_GM)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) strain vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,26),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,26))
                end
                axis(axis_str_GM)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('muscle GM (ins-knee) strain vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,26),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,26))
                end
                axis(axis_str_GM)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle GM (ins-knee) strain vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,26),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,26),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,26),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,26),'b','LineStyle','-','LineWidth',1)
                axis(axis_str_GM)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% GM MUSCLE Lichtwark/Fukunaga: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle GM (architecture) elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,28),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,28),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,28),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,28),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_msc_GM_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_msc_GM_licht_mean, STR_PRE_elong_msc_GM_licht_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_msc_GM_licht_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_msc_GM_licht_mean, STR_POST_elong_msc_GM_licht_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_msc_GM_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_msc_GM_licht_mean, CON_PRE_elong_msc_GM_licht_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_msc_GM_licht_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_msc_GM_licht_mean, CON_POST_elong_msc_GM_licht_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_el_GM)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('muscle GM (architecture) elongation vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,28),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,28))
                end
                axis(axis_el_GM)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('muscle GM (architecture) elongation vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,28),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,28))
                end
                axis(axis_el_GM)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle GM (architecture) elongation vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,28),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,28),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,28),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,28),'b','LineStyle','-','LineWidth',1)
                axis(axis_el_GM)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end

            %% GM TENDON Lichtwark/Fukunaga: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM tendon (architecture) elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,29),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,29),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,29),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,29),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_tend_GM_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_tend_GM_licht_mean, STR_PRE_elong_tend_GM_licht_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_tend_GM_licht_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_tend_GM_licht_mean, STR_POST_elong_tend_GM_licht_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_tend_GM_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_tend_GM_licht_mean, CON_PRE_elong_tend_GM_licht_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_tend_GM_licht_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_tend_GM_licht_mean, CON_POST_elong_tend_GM_licht_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_el_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('GM tendon (architecture) elongation vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,29),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,29))
                end
                axis(axis_el_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('GM tendon (architecture) elongation vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,29),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,29))
                end
                axis(axis_el_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM tendon (architecture) elongation vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,29),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,29),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,29),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,29),'b','LineStyle','-','LineWidth',1)
                axis(axis_el_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end

            %% GM MUSCLE Lichtwark/Fukunaga: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle GM (architecture) strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,30),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,30),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,30),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,30),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_strain_msc_GM_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_strain_msc_GM_licht_mean, STR_PRE_strain_msc_GM_licht_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_strain_msc_GM_licht_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_strain_msc_GM_licht_mean, STR_POST_strain_msc_GM_licht_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_strain_msc_GM_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_strain_msc_GM_licht_mean, CON_PRE_strain_msc_GM_licht_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_strain_msc_GM_licht_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_strain_msc_GM_licht_mean, CON_POST_strain_msc_GM_licht_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_str_GM)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('muscle GM (architecture) strain vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,30),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,30))
                end
                axis(axis_str_GM)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('muscle GM (architecture) strain vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,30),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,30))
                end
                axis(axis_str_GM)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle GM (architecture) strain vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,30),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,30),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,30),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,30),'b','LineStyle','-','LineWidth',1)
                axis(axis_str_GM)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end

            %% GM TENDON Lichtwark/Fukunaga: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM tendon (architecture) strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,31),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,31),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,31),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,31),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_strain_tend_GM_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_strain_tend_GM_licht_mean, STR_PRE_strain_tend_GM_licht_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_strain_tend_GM_licht_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_strain_tend_GM_licht_mean, STR_POST_strain_tend_GM_licht_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_strain_tend_GM_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_strain_tend_GM_licht_mean, CON_PRE_strain_tend_GM_licht_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_strain_tend_GM_licht_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_strain_tend_GM_licht_mean, CON_POST_strain_tend_GM_licht_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_str_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('GM tendon (architecture) strain vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,31),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,31))
                end
                axis(axis_str_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('GM tendon (architecture) strain vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,31),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,31))
                end
                axis(axis_str_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM tendon (architecture) strain vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,31),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,31),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,31),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,31),'b','LineStyle','-','LineWidth',1)
                axis(axis_str_GMtend)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end

            %% GM fascicle length (Lichtwark): length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM fascicle length vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,32),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,32),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,32),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,32),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_GM_length_faslen_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_GM_length_faslen_licht_mean, STR_PRE_GM_length_faslen_licht_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_GM_length_faslen_licht_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_GM_length_faslen_licht_mean, STR_POST_GM_length_faslen_licht_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_GM_length_faslen_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_GM_length_faslen_licht_mean, CON_PRE_GM_length_faslen_licht_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_GM_length_faslen_licht_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_GM_length_faslen_licht_mean, CON_POST_GM_length_faslen_licht_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_faslen)
                xlabel('Gonio angle (°)')
                ylabel('Fascicle length (mm)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('GM fascicle length vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,32),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,32))
                end
                axis(axis_faslen)
                xlabel('Gonio angle (°)')
                ylabel('Fascicle length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('GM fascicle length vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,32),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,32))
                end
                axis(axis_faslen)
                xlabel('Gonio angle (°)')
                ylabel('Fascicle length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM fascicle length vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,32),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,32),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,32),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,32),'b','LineStyle','-','LineWidth',1)
                axis(axis_faslen)
                xlabel('Gonio angle (°)')
                ylabel('Fascicle length (mm)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% GM pennation angle (Lichtwark) vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('GM pennation angle vs ankle angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,33),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,33),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,33),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,33),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_GM_pennation_licht_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_GM_pennation_licht_mean, STR_PRE_GM_pennation_licht_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_GM_pennation_licht_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_GM_pennation_licht_mean, STR_POST_GM_pennation_licht_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_GM_pennation_licht_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_GM_pennation_licht_mean, CON_PRE_GM_pennation_licht_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_GM_pennation_licht_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_GM_pennation_licht_mean, CON_POST_GM_pennation_licht_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_pennation)
                xlabel('Gonio angle (°)')
                ylabel('Pennation angle (°)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('GM pennation angle vs ankle angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,33),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,33))
                end
                axis(axis_pennation)
                xlabel('Gonio angle (°)')
                ylabel('Pennation angle (°)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('GM pennation angle vs ankle angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,33),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,33))
                end
                axis(axis_pennation)
                xlabel('Gonio angle (°)')
                ylabel('Pennation angle (°)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND GM pennation angle vs ankle angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,33),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,33),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,33),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,33),'b','LineStyle','-','LineWidth',1)
                axis(axis_pennation)
                xlabel('Gonio angle (°)')
                ylabel('Pennation angle (°)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% SOL MUSCLE portion: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle SOL length vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,20),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,20),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,20),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,20),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_L_msc_SOL_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_L_msc_SOL_mean, STR_PRE_L_msc_SOL_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_L_msc_SOL_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_L_msc_SOL_mean, STR_POST_L_msc_SOL_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_L_msc_SOL_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_L_msc_SOL_mean, CON_PRE_L_msc_SOL_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_L_msc_SOL_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_L_msc_SOL_mean, CON_POST_L_msc_SOL_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_len_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southeast')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('muscle SOL length vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,20),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,20))
                end
                axis(axis_len_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('muscle SOL length vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,20),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,20))
                end
                axis(axis_len_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle SOL length vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,20),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,20),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,20),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,20),'b','LineStyle','-','LineWidth',1)
                axis(axis_len_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Southeast')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% SOL MUSCLE portion: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle SOL elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,19),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,19),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,19),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,19),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_msc_SOL_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_msc_SOL_mean, STR_PRE_elong_msc_SOL_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_msc_SOL_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_msc_SOL_mean, STR_POST_elong_msc_SOL_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_msc_SOL_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_msc_SOL_mean, CON_PRE_elong_msc_SOL_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_msc_SOL_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_msc_SOL_mean, CON_POST_elong_msc_SOL_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_el_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('muscle SOL elongation vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,19),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,19))
                end
                axis(axis_el_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('muscle SOL elongation vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,19),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,19))
                end
                axis(axis_el_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle SOL elongation vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,19),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,19),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,19),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,19),'b','LineStyle','-','LineWidth',1)
                axis(axis_el_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% SOL MUSCLE portion: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('muscle SOL strain vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,27),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,27),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,27),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,27),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_strain_msc_SOL_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_strain_msc_SOL_mean, STR_PRE_strain_msc_SOL_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_strain_msc_SOL_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_strain_msc_SOL_mean, STR_POST_strain_msc_SOL_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_strain_msc_SOL_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_strain_msc_SOL_mean, CON_PRE_strain_msc_SOL_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_strain_msc_SOL_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_strain_msc_SOL_mean, CON_POST_strain_msc_SOL_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_str_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('muscle SOL strain vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,27),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,27))
                end
                axis(axis_str_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('muscle SOL strain vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,27),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,27))
                end
                axis(axis_str_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND muscle SOL strain vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,27),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,27),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,27),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,27),'b','LineStyle','-','LineWidth',1)
                axis(axis_str_SOL)
                xlabel('Gonio angle (°)')
                ylabel('Strain (% of initial length)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% Full MTU: Length vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('MTU length vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,14),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,14),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,14),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,14),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_L_MTU_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_L_MTU_mean, STR_PRE_L_MTU_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_L_MTU_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_L_MTU_mean, STR_POST_L_MTU_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_L_MTU_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_L_MTU_mean, CON_PRE_L_MTU_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_L_MTU_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_L_MTU_mean, CON_POST_L_MTU_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_len_MTU)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Southeast')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('MTU length vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,14),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,14))
                end
                axis(axis_len_MTU)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('MTU length vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,14),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,14))
                end
                axis(axis_len_MTU)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND MTU length vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,14),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,14),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,14),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,14),'b','LineStyle','-','LineWidth',1)
                axis(axis_len_MTU)
                xlabel('Gonio angle (°)')
                ylabel('Length (mm)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Southeast')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% Full MTU: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('MTU elongation vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,8),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,8),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,8),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,8),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_elong_MTU_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_elong_MTU_mean, STR_PRE_elong_MTU_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_elong_MTU_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_elong_MTU_mean, STR_POST_elong_MTU_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_elong_MTU_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_elong_MTU_mean, CON_PRE_elong_MTU_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_elong_MTU_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_elong_MTU_mean, CON_POST_elong_MTU_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_el_MTU)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('MTU elongation vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,8),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,8))
                end
                axis(axis_el_MTU)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('MTU elongation vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,8),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,8))
                end
                axis(axis_el_MTU)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND MTU elongation vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,8),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,8),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,8),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,8),'b','LineStyle','-','LineWidth',1)
                axis(axis_el_MTU)
                xlabel('Gonio angle (°)')
                ylabel('Elongation (mm)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            % commenting out MTU strain plots since they simply show the Grieve calculation factor
            
            %% Full MTU: Strain vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if plot_check
%                 plottitle = horzcat('MTU strain vs angle - 1');
%                 figure('Name',plottitle)
%                 hold on
%                 plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,23),'Color',col_lightred,'LineStyle','--','LineWidth',1)
%                 plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,23),'r','LineStyle','-','LineWidth',1)
%                 plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,23),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
%                 plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,23),'b','LineStyle','-','LineWidth',1)
%                 
%                 herrorbar(STR_PRE_ROM_mean, STR_PRE_strain_MTU_mean, STR_PRE_ROM_SD, '*m')
%                 errorbar(STR_PRE_ROM_mean, STR_PRE_strain_MTU_mean, STR_PRE_strain_MTU_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
%                 herrorbar(STR_POST_ROM_mean, STR_POST_strain_MTU_mean, STR_POST_ROM_SD, '*r')
%                 errorbar(STR_POST_ROM_mean, STR_POST_strain_MTU_mean, STR_POST_strain_MTU_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
%                 
%                 herrorbar(CON_PRE_ROM_mean, CON_PRE_strain_MTU_mean, CON_PRE_ROM_SD, '*c')
%                 errorbar(CON_PRE_ROM_mean, CON_PRE_strain_MTU_mean, CON_PRE_strain_MTU_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
%                 herrorbar(CON_POST_ROM_mean, CON_POST_strain_MTU_mean, CON_POST_ROM_SD, '*b')
%                 errorbar(CON_POST_ROM_mean, CON_POST_strain_MTU_mean, CON_POST_strain_MTU_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
%                 
%                 axis(axis_str_MTU)
%                 xlabel('Gonio angle (°)')
%                 ylabel('Strain (% of initial length)')
%                 title(plottitle)
%                 legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
%                 saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
%             end
%             
%             if plot_check
%                 plottitle = horzcat('MTU strain vs angle - 2 STRETCH PRE-POST');
%                 figure('Name',plottitle)
%                 hold on
%                 for i = 1:STR_PRE_count
%                     plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,23),'LineStyle','--')
%                 end
%                 set(gca,'ColorOrderIndex',1)
%                 for i = 1:STR_POST_count
%                     plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,23))
%                 end
%                 axis(axis_str_MTU)
%                 xlabel('Gonio angle (°)')
%                 ylabel('Strain (% of initial length)')
%                 title(plottitle)
%                 %legend
%                 saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
%             end
%             if plot_check
%                 plottitle = horzcat('MTU strain vs angle - 3 CONTROL PRE-POST');
%                 figure('Name',plottitle)
%                 hold on
%                 for i = 1:CON_PRE_count
%                     plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,23),'LineStyle','--')
%                 end
%                 set(gca,'ColorOrderIndex',1)
%                 for i = 1:CON_POST_count
%                     plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,23))
%                 end
%                 axis(axis_str_MTU)
%                 xlabel('Gonio angle (°)')
%                 ylabel('Strain (% of initial length)')
%                 title(plottitle)
%                 %legend
%                 saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
%             end
%             
%             % the following plots require that trials are analysed in
%             % systematic order, and all subjects have all 4 legs/timepoints
%             if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
%                 for i = 1:CON_PRE_count
%                     plottitle = horzcat('IND MTU strain vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
%                     figure('Name',plottitle)
%                     hold on
% 
%                     plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,23),'Color',col_lightred,'LineStyle','--','LineWidth',1)
%                     plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,23),'r','LineStyle','-','LineWidth',1)
%                     plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,23),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
%                     plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,23),'b','LineStyle','-','LineWidth',1)
%                 axis(axis_str_MTU)
%                 xlabel('Gonio angle (°)')
%                 ylabel('Strain (% of initial length)')
%                     title(plottitle)
%                     legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
%                     saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
%                 end
%             end
            
            %% GMFAS scans: Displacement vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('displacement GMFAS vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,9),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,9),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,9),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,9),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_displ_GMFAS_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_displ_GMFAS_mean, STR_PRE_displ_GMFAS_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_displ_GMFAS_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_displ_GMFAS_mean, STR_POST_displ_GMFAS_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_displ_GMFAS_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_displ_GMFAS_mean, CON_PRE_displ_GMFAS_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_displ_GMFAS_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_displ_GMFAS_mean, CON_POST_displ_GMFAS_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_displ_GMFAS)
                xlabel('Gonio angle (°)')
                ylabel('Displacement (mm)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('displacement GMFAS vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,9),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,9))
                end
                axis(axis_displ_GMFAS)
                xlabel('Gonio angle (°)')
                ylabel('Displacement (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('displacement GMFAS vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,9),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,9))
                end
                axis(axis_displ_GMFAS)
                xlabel('Gonio angle (°)')
                ylabel('Displacement (mm)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND displacement GMFAS vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,9),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,9),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,9),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,9),'b','LineStyle','-','LineWidth',1)
                axis(axis_displ_GMFAS)
                xlabel('Gonio angle (°)')
                ylabel('Displacement (mm)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% Length +- SD, 3 MTU components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % MMM LATER if needed?
            
            %% EMG vs angle GM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('EMG gas.med. vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,3),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,3),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,3),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,3),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_EMG_gm_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_EMG_gm_mean, STR_PRE_EMG_gm_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_EMG_gm_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_EMG_gm_mean, STR_POST_EMG_gm_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_EMG_gm_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_EMG_gm_mean, CON_PRE_EMG_gm_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_EMG_gm_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_EMG_gm_mean, CON_POST_EMG_gm_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('EMG gas.med. vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,3),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,3))
                end
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('EMG gas.med. vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,3),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,3))
                end
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND EMG gas.med. vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,3),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,3),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,3),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,3),'b','LineStyle','-','LineWidth',1)
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% EMG vs angle GL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('EMG gas.lat. vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,4),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,4),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,4),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,4),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_EMG_gl_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_EMG_gl_mean, STR_PRE_EMG_gl_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_EMG_gl_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_EMG_gl_mean, STR_POST_EMG_gl_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_EMG_gl_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_EMG_gl_mean, CON_PRE_EMG_gl_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_EMG_gl_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_EMG_gl_mean, CON_POST_EMG_gl_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('EMG gas.lat. vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,4),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,4))
                end
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('EMG gas.lat. vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,4),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,4))
                end
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND EMG gas.lat. vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,4),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,4),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,4),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,4),'b','LineStyle','-','LineWidth',1)
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
            
            %% EMG vs angle SOL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_check
                plottitle = horzcat('EMG soleus vs angle - 1');
                figure('Name',plottitle)
                hold on
                plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,5),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,5),'r','LineStyle','-','LineWidth',1)
                plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,5),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,5),'b','LineStyle','-','LineWidth',1)
                
                herrorbar(STR_PRE_ROM_mean, STR_PRE_EMG_sol_mean, STR_PRE_ROM_SD, '*m')
                errorbar(STR_PRE_ROM_mean, STR_PRE_EMG_sol_mean, STR_PRE_EMG_sol_SD, 'Color', col_lightred, 'Marker', '*', 'MarkerFaceColor', col_lightred)
                herrorbar(STR_POST_ROM_mean, STR_POST_EMG_sol_mean, STR_POST_ROM_SD, '*r')
                errorbar(STR_POST_ROM_mean, STR_POST_EMG_sol_mean, STR_POST_EMG_sol_SD, 'Color', 'r', 'Marker', '*', 'MarkerFaceColor', 'r')
                
                herrorbar(CON_PRE_ROM_mean, CON_PRE_EMG_sol_mean, CON_PRE_ROM_SD, '*c')
                errorbar(CON_PRE_ROM_mean, CON_PRE_EMG_sol_mean, CON_PRE_EMG_sol_SD, 'Color', col_lightblue, 'Marker', '*', 'MarkerFaceColor', col_lightblue)
                herrorbar(CON_POST_ROM_mean, CON_POST_EMG_sol_mean, CON_POST_ROM_SD, '*b')
                errorbar(CON_POST_ROM_mean, CON_POST_EMG_sol_mean, CON_POST_EMG_sol_SD, 'Color', 'b', 'Marker', '*', 'MarkerFaceColor', 'b')
                
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                legend('STR PRE','STR POST','CON PRE','CON POST','STR PRE ind max','Location','Northwest')
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            if plot_check
                plottitle = horzcat('EMG soleus vs angle - 2 STRETCH PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:STR_PRE_count
                    plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,5),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:STR_POST_count
                    plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,5))
                end
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            if plot_check
                plottitle = horzcat('EMG soleus vs angle - 3 CONTROL PRE-POST');
                figure('Name',plottitle)
                hold on
                for i = 1:CON_PRE_count
                    plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,5),'LineStyle','--')
                end
                set(gca,'ColorOrderIndex',1)
                for i = 1:CON_POST_count
                    plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,5))
                end
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                title(plottitle)
                %legend
                saveas(gcf, horzcat('data_plots/GRP_INT ',plottitle,'.jpg'))
            end
            
            % the following plots require that trials are analysed in
            % systematic order, and all subjects have all 4 legs/timepoints
            if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
                for i = 1:CON_PRE_count
                    plottitle = horzcat('IND EMG soleus vs angle - SUBJECT ', num2str(CON_PRE_subject_ID(i)) ,' PRE-POST');
                    figure('Name',plottitle)
                    hold on

                    plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,5),'Color',col_lightred,'LineStyle','--','LineWidth',1)
                    plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,5),'r','LineStyle','-','LineWidth',1)
                    plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,5),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
                    plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,5),'b','LineStyle','-','LineWidth',1)
                axis(axis_EMG)
                xlabel('Gonio angle (°)')
                ylabel('EMG (% of MVC)')
                    title(plottitle)
                    legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
                    saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
                end
            end
        
    end
    %% PLOT GROUP FIGURES FINISHED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
   

    clear all
end