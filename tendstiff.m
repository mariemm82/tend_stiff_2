%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main file for tendon stiffness
% Marie Moltubakk 17.5.2013
% 
% Computes tendon stiffness at indidividual + common (across all trials) 
% force level, creates force-elongation plots.
% Analyses two sets of data in parallell:
%    "SOL" = free achilles tendon
%    "GM"  = whole achilles tendon (up to GM MTJ)
% 
% input argument 1 = project selection (1 = BD, 2 = intervent)
% input argument 2 = plot selection (0 = none, 1 = group plots, 2 = ind plots)
% input argument 3 = resume running of loop (0 = start over, 1 = resume)
% 
% Stiffness cutoff rules:
%   1. 6 trials, find lowest produced force
%   2. 90% of the lowest produced force
%   3. IF set lower than #2 above - manual cutoff level, where force-elong curves back
%   4. Individual stiffness = 90-100% of above, and 80-100% of above
%   5. Common stiffness = 90% of the lowest produced force or WEAKEST subject
% 
% The scripts assume a certain structure for input files from Tracker and
% Noraxon. I.e. number of and order of EMG and other channels. If these
% scripts are to be used for other projects, some of the code must be
% modified. These lines are marked % PROJECTSPECIFIC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% MMM TODO 
% - catch error in solve_sec_poly with fzero when equation never reaches 0.0
%   --- under header:   GRP: plot FIT, force-elong per subject (groupwise)
    



function [] = tendstiff(input_project, input_plot, input_resumerun)
    close all
    warning('off','MATLAB:xlswrite:AddSheet')
    

    forceintervals = 50; %VAR - Average stiffness across X N
    


    %% PLOTS - determine which plots to display
    global plot_achilles plot_norm plot_emg plot_check plot_us plot_conversion subject_id

    if input_plot >= 1 
        plot_check = 1; % plot only group summaries
    else
        plot_check = 0;
    end
    if input_plot >= 2
        plot_achilles = 1; % plot stiffness fits / troubleshoot plots
    else
        plot_achilles = 0;
    end
    if input_plot >= 3
        plot_norm = 1; % show torque before and after initial lowpass filter / onset force / ankle rotation fit plots
    else
        plot_norm = 0;
    end
    plot_conversion = 0;
    plot_us = 0; % ultrasound before and after extmark correction
    plot_emg = 0; % RMS 3 EMG channels per trial


    %% Set constants and globals % PROJECTSPECIFIC

    % declare for later use:
    % variables for NORM conversion factors calculated from actual data
    global convert_norm_angle_a convert_norm_angle_b % convert_norm_torque_a convert_norm_torque_b convert_norm_velocity_a convert_norm_velocity_b convert_norm_direction_b

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

    % preallocate "line", for reloading if "input_resumerun" == 1
    line = 0;
    

    %% set AXES etc for plots
    col_lightblue = [0.6 0.8 1];
    col_lightred = [1 0.6 0.8];
    axis_stiff_free = [0 14 -100 Inf];
    axis_stiff_whole = [0 24 -100 Inf];
    


    %% Read datamaster file (connect corresponding data files) EVENTUALLY resume running of loop
    % Produces arrays with file names and variables per trial, to be retrieved later
    global dm_subjectno dm_timepoint dm_side dm_trial dm_group
    global dm_stiff1_NX dm_stiff1_US dm_stiff1_US_frame dm_stiff2_NX dm_stiff2_US dm_stiff2_US_frame dm_stiff3_NX dm_stiff3_US dm_stiff3_US_frame 
%    global dm_heel1_NX dm_heel1_US dm_heel1_US_frame dm_heel2_NX dm_heel2_US dm_heel2_US_frame dm_heel3_NX dm_heel3_US dm_heel3_US_frame
    global dm_MVC_PF dm_MVC_DF dm_CPM_calc_NX dm_CPM_sol_NX dm_leg_length % dm_CPM_calc_US dm_CPM_calc_US_frame 
    global dm_tendonlength dm_cutforce
    global dm_rot_const dm_CSA
    global filepath
    dm_filename = 'data/datamaster_stiff.tsv';
    if input_resumerun == 1 % resume running of loop, with new datamaster version (filenames may be edited, line order NOT!)
        load all_data_stiff_SOL %_inloop
        line_start = line+1; % all_data_stiff_inloop ended on a line - resume with next line
        input_resumerun = 1; % overwrite variable coming from all_data_stiff_inloop
    else
        line_start = 1;
    end
    linestotal = read_datamaster_stiff(dm_filename);

    
    %% preallocate
    % common arrays for all subjects:
    if input_resumerun == 0
        all_stiff_output_head = {'Subject', 'Time', 'Side', 'Trial', ...
            'AT moment arm, m', 'Rotation correction, mm per deg', ... % 1-2
            'StiffEQ coeff 1', 'StiffEQ coeff 2', 'StiffEQ coeff 3', 'StiffEQ R2', ... % 3-6
            'Force Ramp ind F cutoff, N', 'Force Ramp ind F max', 'Force Ramp F at max-elong', 'Force MVC F plantflex', ... % 7-10
            'Force Ramp common F cutoff', 'Force Ramp common F max',... % 11-12
            'Stiff ind 80-100, N per mm', 'Stiff ind 90-100', ... % 13-14
            'Stiff common cutoff 80-100', 'Stiff common cutoff 90-100', 'Stiff common max 80-100', 'Stiff common max 90-100', ... % 15-18
            'Elong at ind F cutoff, mm', 'Elong ind max-elong', 'Elong at common F cutoff', 'Elong at common F max', ... % 19-22
            'Tend-length 0deg, mm', ... % 23
            'Strain at ind F cutoff, percent', 'Strain ind max-elong', 'Strain at common F cutoff', 'Strain at common F max', ... % 24-27
            'Young at ind F cutoff', 'Young at ind max-elong', 'Young at common F cutoff', 'Young at common F max', ... % 28-31
            'Tend-CSA, mm^2',... % 32
            }; % PROJECTSPECIFIC
        loc_stiff_a = 3;
        loc_stiff_b = 4;
        loc_stiff_c = 5;
        loc_ind_force_cut = 7;
        loc_ind_force_max = 8;
        % loc_ind_force_elongmax = 9;
        loc_common_force_cut_array = 11;
        loc_common_force_max_array = 12;
        loc_stiff_common_cut_80 = 15;
        loc_stiff_common_cut_90 = 16;
        loc_stiff_common_max_80 = 17;
        loc_stiff_common_max_90 = 18;
        loc_elong_common_cut = 21;
        loc_elong_common_max = 22;
        loc_tendonlength_initial = 23;
        loc_strain_common_cut = 26;
        loc_strain_common_max = 27;
        loc_young_common_cut = 30;
        loc_young_common_max = 31;
        loc_tendonCSA = 32;
        
        all_stiff_output = zeros(ceil(linestotal),length(all_stiff_output_head)-4); 
        all_stiff_output_txt = cell(ceil(linestotal),4);
    
        force_elong_ALL = cell(1,ceil(linestotal));
        
        if input_project == 1 % BD study
            BD_SOL_count = 0; % # of ballet dancer subjects
            CON_SOL_count = 0; % # of controls = intervention study subjects
            BD_GM_count = 0;
            CON_GM_count = 0;
            
            BD_SOL_no(ceil(linestotal)) = zeros;
            CON_SOL_no(ceil(linestotal)) = zeros;
            BD_GM_no(ceil(linestotal)) = zeros;
            CON_GM_no(ceil(linestotal)) = zeros;
            
            force_elong_SOL_BD = cell(1,ceil(linestotal));
            force_elong_SOL_CON = cell(1,ceil(linestotal));
            force_elong_GM_BD = cell(1,ceil(linestotal));
            force_elong_GM_CON = cell(1,ceil(linestotal));
            
            stiffness_SOL_BD(1:ceil(linestotal),1:5) = NaN;
            stiffness_SOL_CON(1:ceil(linestotal),1:5) = NaN;
            stiffness_GM_BD(1:ceil(linestotal),1:5) = NaN;
            stiffness_GM_CON(1:ceil(linestotal),1:5) = NaN;
        else % intervention
            STR_PRE_SOL_count = 0;
            STR_POST_SOL_count = 0;
            CON_PRE_SOL_count = 0;
            CON_POST_SOL_count = 0;
            STR_PRE_GM_count = 0;
            STR_POST_GM_count = 0;
            CON_PRE_GM_count = 0;
            CON_POST_GM_count = 0;
            
            STR_PRE_SOL_ID{ceil(linestotal)} = [];
            STR_POST_SOL_ID{ceil(linestotal)} = [];
            CON_PRE_SOL_ID{ceil(linestotal)} = [];
            CON_POST_SOL_ID{ceil(linestotal)} = [];
            STR_PRE_GM_ID{ceil(linestotal)} = [];
            STR_POST_GM_ID{ceil(linestotal)} = [];
            CON_PRE_GM_ID{ceil(linestotal)} = [];
            CON_POST_GM_ID{ceil(linestotal)} = [];
            
            force_elong_STR_PRE_SOL = cell(1,ceil(linestotal));
            force_elong_STR_POST_SOL = cell(1,ceil(linestotal));
            force_elong_CON_PRE_SOL = cell(1,ceil(linestotal));
            force_elong_CON_POST_SOL = cell(1,ceil(linestotal));
            force_elong_STR_PRE_GM = cell(1,ceil(linestotal));
            force_elong_STR_POST_GM = cell(1,ceil(linestotal));
            force_elong_CON_PRE_GM = cell(1,ceil(linestotal));
            force_elong_CON_POST_GM = cell(1,ceil(linestotal));
            
            stiffness_STR_PRE_SOL(1:ceil(linestotal),1:5) = NaN;
            stiffness_STR_POST_SOL(1:ceil(linestotal),1:5) = NaN;
            stiffness_CON_PRE_SOL(1:ceil(linestotal),1:5) = NaN;
            stiffness_CON_POST_SOL(1:ceil(linestotal),1:5) = NaN;
            stiffness_STR_PRE_GM(1:ceil(linestotal),1:5) = NaN;
            stiffness_STR_POST_GM(1:ceil(linestotal),1:5) = NaN;
            stiffness_CON_PRE_GM(1:ceil(linestotal),1:5) = NaN;
            stiffness_CON_POST_GM(1:ceil(linestotal),1:5) = NaN;
        end
    end
    
    
    %% Loop through all lines in datamaster file (except header line)
    for line = line_start:linestotal


        %% subject/trial identifier
        
        if input_project == 1 % BD study
            trial_subjectno = str2double(dm_subjectno{line});
            if trial_subjectno > 100
                filepath = 'data\BD\';
                subject_id = horzcat('Dancer ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
                if strcmp(dm_trial{line},'SOL')
                    BD_SOL_count = BD_SOL_count + 1;
                    BD_SOL_no(BD_SOL_count) = trial_subjectno;
                else % GM
                    BD_GM_count = BD_GM_count + 1;
                    BD_GM_no(BD_GM_count) = trial_subjectno;
                end
            else
                filepath = 'data\';
                subject_id = horzcat('Control ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
                if strcmp(dm_trial{line},'SOL')
                    CON_SOL_count = CON_SOL_count + 1;
                    CON_SOL_no(CON_SOL_count) = trial_subjectno;
                else % GM
                    CON_GM_count = CON_GM_count + 1;
                    CON_GM_no(CON_GM_count) = trial_subjectno;
                end
            end

        elseif input_project == 2 % intervention
            trial_subjectno = str2double(dm_subjectno{line});
            trial_timepoint = strcmp(dm_timepoint{line},'POST'); % 0 = PRE, 1 = POST
            %trial_leg = strcmp(dm_side{line},'L'); % 0 = RIGHT, 1 = LEFT
            trial_location = strcmp(dm_trial{line},'GM'); % 0 = SOL, 1 = GM
            trial_group = strcmp(dm_group{line},'STR'); % 0 = CON, 1 = STR
            filepath = 'data\';
            subject_id = horzcat('INT_', dm_subjectno{line}, '_', dm_trial{line}, '_', dm_timepoint{line}, '_', dm_group{line}, '_', dm_side{line});
            
            if trial_location == 0 % SOL
                if trial_timepoint == 0 && trial_group == 1 % PRE, STR
                    STR_PRE_SOL_count = STR_PRE_SOL_count + 1;
                    STR_PRE_SOL_ID{STR_PRE_SOL_count} = subject_id;
                elseif trial_timepoint == 1 && trial_group == 1 % POST, STR
                    STR_POST_SOL_count = STR_POST_SOL_count + 1;
                    STR_POST_SOL_ID{STR_POST_SOL_count} = subject_id;
                elseif trial_timepoint == 0 && trial_group == 0 % PRE, CON
                    CON_PRE_SOL_count = CON_PRE_SOL_count + 1;
                    CON_PRE_SOL_ID{CON_PRE_SOL_count} = subject_id;
                elseif trial_timepoint == 1 && trial_group == 0 % POST, CON
                    CON_POST_SOL_count = CON_POST_SOL_count + 1;
                    CON_POST_SOL_ID{CON_POST_SOL_count} = subject_id;
                end
            else % GM
                if trial_timepoint == 0 && trial_group == 1 % PRE, STR
                    STR_PRE_GM_count = STR_PRE_GM_count + 1;
                    STR_PRE_GM_ID{STR_PRE_GM_count} = subject_id;
                elseif trial_timepoint == 1 && trial_group == 1 % POST, STR
                    STR_POST_GM_count = STR_POST_GM_count + 1;
                    STR_POST_GM_ID{STR_POST_GM_count} = subject_id;
                elseif trial_timepoint == 0 && trial_group == 0 % PRE, CON
                    CON_PRE_GM_count = CON_PRE_GM_count + 1;
                    CON_PRE_GM_ID{CON_PRE_GM_count} = subject_id;
                elseif trial_timepoint == 1 && trial_group == 0 % POST, CON
                    CON_POST_GM_count = CON_POST_GM_count + 1;
                    CON_POST_GM_ID{CON_POST_GM_count} = subject_id;
                end
            end
        end
        cprintf('*black', horzcat('----------------', subject_id, '------------------\n'))

        
        %% Calculate conversion factors for tibialis anterior CO-ACTIVATION
        % Produce individual conversion factors for ANGLE
        % using only "angle_constants" not "calculate_norm_constants_ACTIVE" because this is for NORM data only, and tendstiff uses only passive trials from NORM
        
        % trial 1
        [convert_ind_angle_a1, convert_ind_angle_b1] = calculate_angle_constants(angle_cutoff, horzcat(filepath, dm_CPM_calc_NX{line}), dm_side{line});
        % trial 2
        [convert_ind_angle_a2, convert_ind_angle_b2] = calculate_angle_constants(angle_cutoff, horzcat(filepath, dm_CPM_sol_NX{line}), dm_side{line});
        % average
        convert_norm_angle_a = mean([convert_ind_angle_a1 convert_ind_angle_a2]);
        convert_norm_angle_b = mean([convert_ind_angle_b1 convert_ind_angle_b2]);
        
        % Read co-activation noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
        % Produce a new noraxon data array
        noraxon_coact = read_noraxon_stiffness(strcat(filepath, dm_MVC_DF{line}), freq_default, dm_side{line}, 'MVC dorsi');

        % Calculate co-activation constants
        % Read complete, prepared noraxon array + number of frames to average (freq * time)
        % Produce max torque, max EMG constants
        [coact_max_torque,coact_max_EMG] = calculate_coactivation(noraxon_coact, freq_default*(mvc_window_ms/1000), dm_side{line});
        

        %% Calculate achilles tendon MOMENT ARM

        % Read moment arm trial US data file, determine time stamps, set trigger frame as time = zero
        % Produce US sample frequency, create new US array containing time and displacement
%        [usdata_CPM, usfreq_CPM] = read_us_file(strcat(filepath, dm_CPM_calc_US{line}, '.txt'), str2double(dm_CPM_calc_US_frame{line}), 'CPM calcaneus');

        % Read moment arm noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
        % Produce a new noraxon data array
%        noraxon_CPM = read_noraxon_stiffness(strcat(filepath, dm_CPM_calc_NX{line}), usfreq_CPM, dm_side{line}, 'CPM calcaneus');

        % above variables are no longer used for momentarm, but are (re-)used for ankle rotation below
        
        % Read complete, prepared noraon array + prepared us array
        % Produce AT moment arm constant
        at_momentarm = calculate_momentarm(0,0, dm_leg_length{line});


        %% ANKLE ROTATION from CALC CPM: calc displacement vs ankle angle
        % Read complete, prepared noraxon array + prepared us array
        % Produce ankle rotation constant, mm/degree
        % 
        % REUSING us file from moment arm
        % REUSING noraxon file from moment arm
        
        % no longer correcting for rotation at this point - instead, "3x OTJ" contains rotation correction
        at_rotation_const = 0; 

% 3 lines below = enable to retrieve rotation constants from CPM files (i.e. SOL CPM trials). 2017-12-12
%        [usdata_CPM, usfreq_CPM] = read_us_file(strcat(filepath, dm_CPM_calc_US{line}, '.txt'), str2double(dm_CPM_calc_US_frame{line}), 'CPM nowSOL');
%        noraxon_CPM = read_noraxon_stiffness(strcat(filepath, dm_CPM_calc_NX{line}), usfreq_CPM, dm_side{line}, 'CPM nowSOL');
%        at_rotation_const = calculate_rotation_correction(noraxon_CPM, usdata_CPM);
        
%         if strcmp(subject_id,'Control 4 R PRE SOL') || strcmp(subject_id,'INT_4_SOL_PRE_STR_R') || strcmp(subject_id,'INT_4_GM_PRE_STR_R')
%             % ROUGH coding, manually calculated for below subjects
%             at_rotation_const = -0.04496; % BD predefined = -0.14; % auto-calculated value = -0.61/-0.35
%         elseif strcmp(subject_id,'INT_10_SOL_PRE_CON_L') || strcmp(subject_id,'INT_10_GM_PRE_CON_L')
%             at_rotation_const = -0.091707; % auto-calculated value = positive
%         else
%             % NORMALLY
%             at_rotation_const = calculate_rotation_correction(noraxon_CPM, usdata_CPM);
%         end


        %% Calculate MVC for plantarflexion

        % Read MVC noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
        % Produce a new noraxon data array
        % sending in a length corresponding to 9 seconds (delete anything after 9 sec)
        noraxon_MVC = read_noraxon_stiffness(strcat(filepath, dm_MVC_PF{line}), freq_default, dm_side{line}, 'MVC plantar flex');

        % Calculate co-activation constants
        % Read complete, prepared noraxon array + number of frames to average (freq * time)
        % Produce max torque, max EMG constants
        [plantflex_max_torque,~] = calculate_MVC(noraxon_MVC, freq_default*(mvc_window_ms/1000), at_momentarm, dm_side{line});


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
        % 2017-12-08: using rotation correction for MTJ instead of OTJ displacement
%         if(strcmp(dm_heel1_NX{line}, 'null'))
%             time_force_displ_otj1 = zeros(0);
%         else
%             [time_force_displ_otj1,trial_force_max(4)] = extract_force_displ_singletrial(dm_heel1_NX{line}, dm_heel1_US{line}, dm_heel1_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, 0, dm_side{line}, 'OTJ1');
%         end
%         if(strcmp(dm_heel2_NX{line}, 'null'))
%             time_force_displ_otj2 = zeros(0);
%         else
%             [time_force_displ_otj2,trial_force_max(5)] = extract_force_displ_singletrial(dm_heel2_NX{line}, dm_heel2_US{line}, dm_heel2_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, 0, dm_side{line}, 'OTJ2');
%         end
%         if(strcmp(dm_heel3_NX{line}, 'null'))
%             time_force_displ_otj3 = zeros(0);
%         else
%             [time_force_displ_otj3,trial_force_max(6)] = extract_force_displ_singletrial(dm_heel3_NX{line}, dm_heel3_US{line}, dm_heel3_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, 0, dm_side{line}, 'OTJ3');
%         end
        
        
        %% ankle rotation correction per trial
        % create time-force-displacement for 3 trials
        % displacement is based on: dm_rot_const * gonio angle
        % 
        % reading 3 MTJ trials (can combine US displacement and ankle rotation from same trial)
        % saving into OTJ variables for simplicity
        % 
        if(strcmp(dm_stiff1_NX{line}, 'null'))
            time_force_displ_otj1 = zeros(0);
        else 
            [time_force_displ_otj1,~] = extract_force_displ_singletrial_rotation(dm_stiff1_NX{line}, dm_stiff1_US{line}, dm_stiff1_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, dm_rot_const{line}, dm_side{line}, 'MTJ1-rotation');
        end
        if(strcmp(dm_stiff2_NX{line}, 'null'))
            time_force_displ_otj2 = zeros(0);
        else
            [time_force_displ_otj2,~] = extract_force_displ_singletrial_rotation(dm_stiff2_NX{line}, dm_stiff2_US{line}, dm_stiff2_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, dm_rot_const{line}, dm_side{line}, 'MTJ2-rotation');
        end
        if(strcmp(dm_stiff3_NX{line}, 'null'))
            time_force_displ_otj3 = zeros(0);
        else
            [time_force_displ_otj3,~] = extract_force_displ_singletrial_rotation(dm_stiff3_NX{line}, dm_stiff3_US{line}, dm_stiff3_US_frame{line}, coact_max_torque, coact_max_EMG, at_momentarm, dm_rot_const{line}, dm_side{line}, 'MTJ3-rotation');
        end
        


        %% Final stiffness
        % Read time-force-displacement data
        %    sending 3+3 trials + max forces from trials + manually set cutoff force
        % Produce stiffness equation (based on CUT data set) +
        %    force-elong-array (up to defined force level of 90% of 6-trial-common-force or 90% of manual cutoff)
        [stiff_eq, stiff_gof, force_elong_array, loc_cutoff, elong_indmax, force_indelongmax, strain_indmax, young_indmax, elong_indcut, force_indcut, strain_indcut, young_indcut] = final_stiffness(time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3, time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3, forceintervals, dm_cutforce{line}, trial_force_max, dm_tendonlength{line}, dm_CSA{line});


        %% calculate stiffness for last 10 and 20% of ind max:
        force_cutoff_ind = force_elong_array(loc_cutoff,2); % defined cutoff point of 90% of 6-trial-common-force or 90% of manually set force
        stiff_ind_80 = calculate_stiffness(stiff_eq, force_cutoff_ind, 0.8, 1.0, 'ind max'); % last two variables are percent range, from 0.00 to 1.00
        stiff_ind_90 = calculate_stiffness(stiff_eq, force_cutoff_ind, 0.9, 1.0, 'ind max');
       
        
        %% save individual data to common array
        % add all individual variables to a common array for all subjects    

        % headers
        all_stiff_output_txt(line,:) = [dm_subjectno(line) dm_timepoint(line) dm_group(line) dm_trial(line)];
        
        % max force, 2 options:
        % -- min(trial_force_max) = absolute value, common of 6 trials, eg 4235 M
        % -- force_elong_array(end,2) = rounded value, common of 6 trials, eg 4200 N
             
        all_stiff_output(line,:) = [...
            at_momentarm at_rotation_const ... % 1-2
            coeffvalues(stiff_eq) stiff_gof.rsquare ... % 3-6 -- 3x coeffisients + R^2
            force_indcut min(trial_force_max) force_indelongmax plantflex_max_torque ... % 7-10 -- FORCE INDIVIDUAL 90%-of-6-trial-common/manual-force / 100%-6-trial-force / force at maxelong / MVC-force
            NaN NaN ... % 11-12 -- FORCE COMMON: common cutoff / common force max
            stiff_ind_80 stiff_ind_90 ... % 13-14
            NaN NaN NaN NaN ... % 15-18 --- NaN for 4x stiffness at force levels common to all subjects
            elong_indcut elong_indmax NaN NaN ... % 19-22: ELONG, indmax + COMMON: NaN for 2x elong at force levels common to all subjects
            str2double(dm_tendonlength{line}) ... % 23
            strain_indcut strain_indmax NaN NaN... % 24-27 --- STRAIN
            young_indcut young_indmax NaN NaN... % 28-31 --- YOUNG'S MODULUS
            str2double(dm_CSA{line}) ... % 32
            ];
        
        
        %% save stiffness fits and force-elongation arrays to groupwise cells, for future averaging
        % contents of stiffness_(((SOL_BD))): 
        %     2nd order equation coeffs / elongation max (= xlim for plots) / force max

        force_elong_ALL{line} = force_elong_array; % old: save only until cut force - (1:loc_cutoff);
        
         if input_project == 1 % BD study
            if trial_subjectno > 100 % BD subject
                if strcmp(dm_trial{line}, 'SOL')
                    force_elong_SOL_BD{BD_SOL_count} = force_elong_array(1:loc_cutoff,:);
                    stiffness_SOL_BD(BD_SOL_count,:) = [coeffvalues(stiff_eq), force_elong_array(loc_cutoff,1), force_elong_array(loc_cutoff,2)];
                else % GM
                    force_elong_GM_BD{BD_GM_count} = force_elong_array(1:loc_cutoff,:);
                    stiffness_GM_BD(BD_GM_count,:) = [coeffvalues(stiff_eq), force_elong_array(loc_cutoff,1), force_elong_array(loc_cutoff,2)];
                end
            else % CON subject
                if strcmp(dm_trial{line}, 'SOL') == 1
                    force_elong_SOL_CON{CON_SOL_count} = force_elong_array(1:loc_cutoff,:);
                    stiffness_SOL_CON(CON_SOL_count,:) = [coeffvalues(stiff_eq), force_elong_array(loc_cutoff,1), force_elong_array(loc_cutoff,2)];
                else % GM
                    force_elong_GM_CON{CON_GM_count} = force_elong_array(1:loc_cutoff,:);
                    stiffness_GM_CON(CON_GM_count,:) = [coeffvalues(stiff_eq), force_elong_array(loc_cutoff,1), force_elong_array(loc_cutoff,2)];
                end
            end
         else % intervention study
            if trial_location == 0 % SOL
                if trial_timepoint == 0 && trial_group == 1 % PRE, STR
                    force_elong_STR_PRE_SOL{STR_PRE_SOL_count} = force_elong_array(1:loc_cutoff,:);
                    stiffness_STR_PRE_SOL(STR_PRE_SOL_count,:) = [coeffvalues(stiff_eq), force_elong_array(loc_cutoff,1), force_elong_array(loc_cutoff,2)];
                elseif trial_timepoint == 1 && trial_group == 1 % POST, STR
                    force_elong_STR_POST_SOL{STR_POST_SOL_count} = force_elong_array(1:loc_cutoff,:);
                    stiffness_STR_POST_SOL(STR_POST_SOL_count,:) = [coeffvalues(stiff_eq), force_elong_array(loc_cutoff,1), force_elong_array(loc_cutoff,2)];
                elseif trial_timepoint == 0 && trial_group == 0 % PRE, CON
                    force_elong_CON_PRE_SOL{CON_PRE_SOL_count} = force_elong_array(1:loc_cutoff,:);
                    stiffness_CON_PRE_SOL(CON_PRE_SOL_count,:) = [coeffvalues(stiff_eq), force_elong_array(loc_cutoff,1), force_elong_array(loc_cutoff,2)];
                elseif trial_timepoint == 1 && trial_group == 0 % POST, CON
                    force_elong_CON_POST_SOL{CON_POST_SOL_count} = force_elong_array(1:loc_cutoff,:);
                    stiffness_CON_POST_SOL(CON_POST_SOL_count,:) = [coeffvalues(stiff_eq), force_elong_array(loc_cutoff,1), force_elong_array(loc_cutoff,2)];
                end
            else % GM
                if trial_timepoint == 0 && trial_group == 1 % PRE, STR
                    force_elong_STR_PRE_GM{STR_PRE_GM_count} = force_elong_array(1:loc_cutoff,:);
                    stiffness_STR_PRE_GM(STR_PRE_GM_count,:) = [coeffvalues(stiff_eq), force_elong_array(loc_cutoff,1), force_elong_array(loc_cutoff,2)];
                elseif trial_timepoint == 1 && trial_group == 1 % POST, STR
                    force_elong_STR_POST_GM{STR_POST_GM_count} = force_elong_array(1:loc_cutoff,:);
                    stiffness_STR_POST_GM(STR_POST_GM_count,:) = [coeffvalues(stiff_eq), force_elong_array(loc_cutoff,1), force_elong_array(loc_cutoff,2)];
                elseif trial_timepoint == 0 && trial_group == 0 % PRE, CON
                    force_elong_CON_PRE_GM{CON_PRE_GM_count} = force_elong_array(1:loc_cutoff,:);
                    stiffness_CON_PRE_GM(CON_PRE_GM_count,:) = [coeffvalues(stiff_eq), force_elong_array(loc_cutoff,1), force_elong_array(loc_cutoff,2)];
                elseif trial_timepoint == 1 && trial_group == 0 % POST, CON
                    force_elong_CON_POST_GM{CON_POST_GM_count} = force_elong_array(1:loc_cutoff,:);
                    stiffness_CON_POST_GM(CON_POST_GM_count,:) = [coeffvalues(stiff_eq), force_elong_array(loc_cutoff,1), force_elong_array(loc_cutoff,2)];
                end
            end
         end
         save all_data_stiff_inloop
         close all
    end
    %% LOOP finished --- loop end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cprintf('*black', horzcat('----------------', ' Loop finished ', '------------------\n'))
    
        
    %% truncate cells
    if input_project == 1 % BD study
        force_elong_SOL_BD = force_elong_SOL_BD(~cellfun('isempty',force_elong_SOL_BD));
        force_elong_GM_BD = force_elong_GM_BD(~cellfun('isempty',force_elong_GM_BD));
        force_elong_SOL_CON = force_elong_SOL_CON(~cellfun('isempty',force_elong_SOL_CON));
        force_elong_GM_CON = force_elong_GM_CON(~cellfun('isempty',force_elong_GM_CON));
        
        stiffness_SOL_BD(BD_SOL_count+1:end,:) = [];
        stiffness_GM_BD(BD_GM_count+1:end,:) = [];
        stiffness_SOL_CON(CON_SOL_count+1:end,:) = [];
        stiffness_GM_CON(CON_GM_count+1:end,:) = [];
        
    else % intervention
        force_elong_STR_PRE_SOL = force_elong_STR_PRE_SOL(~cellfun('isempty',force_elong_STR_PRE_SOL));
        force_elong_STR_POST_SOL = force_elong_STR_POST_SOL(~cellfun('isempty',force_elong_STR_POST_SOL));
        force_elong_CON_PRE_SOL = force_elong_CON_PRE_SOL(~cellfun('isempty',force_elong_CON_PRE_SOL));
        force_elong_CON_POST_SOL = force_elong_CON_POST_SOL(~cellfun('isempty',force_elong_CON_POST_SOL));
        force_elong_STR_PRE_GM = force_elong_STR_PRE_GM(~cellfun('isempty',force_elong_STR_PRE_GM));
        force_elong_STR_POST_GM = force_elong_STR_POST_GM(~cellfun('isempty',force_elong_STR_POST_GM));
        force_elong_CON_PRE_GM = force_elong_CON_PRE_GM(~cellfun('isempty',force_elong_CON_PRE_GM));
        force_elong_CON_POST_GM = force_elong_CON_POST_GM(~cellfun('isempty',force_elong_CON_POST_GM));
        
        stiffness_STR_PRE_SOL(STR_PRE_SOL_count+1:end,:) = [];
        stiffness_STR_POST_SOL(STR_POST_SOL_count+1:end,:) = [];
        stiffness_CON_PRE_SOL(CON_PRE_SOL_count+1:end,:) = [];
        stiffness_CON_POST_SOL(CON_POST_SOL_count+1:end,:) = [];
        stiffness_STR_PRE_GM(STR_PRE_GM_count+1:end,:) = [];
        stiffness_STR_POST_GM(STR_POST_GM_count+1:end,:) = [];
        stiffness_CON_PRE_GM(CON_PRE_GM_count+1:end,:) = [];
        stiffness_CON_POST_GM(CON_POST_GM_count+1:end,:) = [];
    end
    
    
    %% IND data: calculate stiffness at group COMMON FORCE
    % select common force (maximal force reached by all subjects)
    stiff_common_force = min(all_stiff_output(:,loc_ind_force_cut));
    stiff_common_force_max = min(all_stiff_output(:,loc_ind_force_max));
    
    % create stiffness equation as cfit
    f = fittype('a*x^2+b*x+c');
    a = 0; % temporary values, will be filled/replaced inside loop
    b = 0;
    c = 0;
    stiff_eq_group = cfit(f, a, b, c);

    for i = 1:size(all_stiff_output,1)
        stiff_eq_group.a = all_stiff_output(i,loc_stiff_a); % a 
        stiff_eq_group.b = all_stiff_output(i,loc_stiff_b); % b
        stiff_eq_group.c = all_stiff_output(i,loc_stiff_c); % c

        % calculate stiffness
        stiff_common_80 = calculate_stiffness(stiff_eq_group, stiff_common_force, 0.8, 1.0, horzcat(all_stiff_output_txt{i,1}, ' ', all_stiff_output_txt{i,4}, ' ', all_stiff_output_txt{i,2}, ' ', all_stiff_output_txt{i,3}, ' common cutoff force')); % last two = percent of submitted force - %VAR
        stiff_common_90 = calculate_stiffness(stiff_eq_group, stiff_common_force, 0.9, 1.0, horzcat( all_stiff_output_txt{i,1}, ' ', all_stiff_output_txt{i,4}, ' ', all_stiff_output_txt{i,2}, ' ', all_stiff_output_txt{i,3}, ' common cutoff force')); %VAR
        stiff_common_80_max = calculate_stiffness(stiff_eq_group, stiff_common_force_max, 0.8, 1.0, horzcat(all_stiff_output_txt{i,1}, ' ', all_stiff_output_txt{i,4}, ' ', all_stiff_output_txt{i,2}, ' ', all_stiff_output_txt{i,3}, ' common max force')); % last two = percent of submitted force - %VAR
        stiff_common_90_max = calculate_stiffness(stiff_eq_group, stiff_common_force_max, 0.9, 1.0, horzcat(all_stiff_output_txt{i,1}, ' ', all_stiff_output_txt{i,4}, ' ', all_stiff_output_txt{i,2}, ' ', all_stiff_output_txt{i,3}, ' common max force')); %VAR

        % add to array across subjects
        all_stiff_output(i,loc_stiff_common_cut_80) = stiff_common_80;
        all_stiff_output(i,loc_stiff_common_cut_90) = stiff_common_90;
        all_stiff_output(i,loc_stiff_common_max_80) = stiff_common_80_max;
        all_stiff_output(i,loc_stiff_common_max_90) = stiff_common_90_max;
        all_stiff_output(i,loc_common_force_cut_array) = stiff_common_force;
        all_stiff_output(i,loc_common_force_max_array) = stiff_common_force_max;
    end
    
    cprintf('blue*',horzcat('Stiffness: Common cutoff force = ', num2str(stiff_common_force), ' N, common max force = ', num2str(round(stiff_common_force_max,0)), ' N.\n'))

        
    %% IND data: calculate STRAIN - ELONGATION - YOUNG'S MODULUS at COMMON force levels
    % re-use force levels from stiffness:
    %    stiff_common_force - 90%/cutoff
    %    stiff_common_force_max - common max
    
    % young's modulus = stress / strain
    % stress = force per CSA (in m^2) between 80 and 100% force
    % strain = % elongation between 80 and 100% force
    young_percentage = 0.8; %VAR calculate between 80% and defined end (max/cutoff)
    
    for i = 1:size(all_stiff_output,1)
        % extract ELONGATION and STRAIN from force-elong-cell containing all subjects:
        
        ym_csa = (all_stiff_output(i,loc_tendonCSA)/1000/1000); % converting mm^2 to m^2
        ym_tendlength = all_stiff_output(i,loc_tendonlength_initial);
        
        % CUT: elongation and strain at common force which is cut to 90%:
        loc_common_force_cut = find(force_elong_ALL{i}(:,2)>=stiff_common_force,1,'first');
        if isempty(loc_common_force_cut)
            elong_common_cut = NaN;
            strain_common_cut = NaN;
            young_common_cut = NaN;
            cprintf('red',horzcat('WARNING: Line ', num2str(i), ' (subj ', all_stiff_output_txt{i,1}, '): Elongation @ 90%% common force not found.\n'))
        else
            % end values
            elong_common_cut = force_elong_ALL{i}(loc_common_force_cut,1);
            strain_common_cut = elong_common_cut / ym_tendlength * 100;
            force_common_cut = stiff_common_force;
            % 80% of end force - for Young
            loc_force_common_cut_80 = find(force_elong_ALL{i}(:,2) >= (force_common_cut * young_percentage), 1, 'first');
            elong_common_cut_80 = force_elong_ALL{i}(loc_force_common_cut_80,1);
            force_common_cut_80 = force_elong_ALL{i}(loc_force_common_cut_80,2);
            % Young data
            ym_stress_common_80cut = (force_common_cut - force_common_cut_80) / ym_csa; 
            ym_strain_common_80cut = (elong_common_cut - elong_common_cut_80) / ym_tendlength * 100;
            young_common_cut = ym_stress_common_80cut / ym_strain_common_80cut;
        end
        
        % MAX: elongation and strain at common force BEFORE cut to 90%:
        %    The below variable will not be extracted (will be NaN) for several subjects, since
        %    force_elong_ALL contains force-elongation only up to individual 90%/manual-cut force level.
        %    Not deleted from data analysis, just to save time.
        loc_common_force_max = find(force_elong_ALL{i}(:,2)>=stiff_common_force_max,1,'first');
        if isempty(loc_common_force_max)
            elong_common_max = NaN;
            strain_common_max = NaN;
            young_common_max = NaN;
        else
            % end values
            elong_common_max = force_elong_ALL{i}(loc_common_force_max,1);
            strain_common_max = elong_common_max / ym_tendlength * 100;
            force_common_max = force_elong_ALL{i}(loc_common_force_max,2); % cannot use exact max, do not have corresponding elong
            % 80% of end force - for Young
            loc_force_common_max_80 = find(force_elong_ALL{i}(:,2) >= (stiff_common_force_max * young_percentage), 1, 'first');
            elong_common_max_80 = force_elong_ALL{i}(loc_force_common_max_80,1);
            force_common_max_80 = force_elong_ALL{i}(loc_force_common_max_80,2);
            % Young data
            ym_stress_common_80max = (force_common_max - force_common_max_80) / ym_csa; 
            ym_strain_common_80max = (elong_common_max - elong_common_max_80) / ym_tendlength * 100;
            young_common_max = ym_stress_common_80max / ym_strain_common_80max;
        end
        
        % add to array across subjects
        all_stiff_output(i,loc_elong_common_cut) = elong_common_cut;
        all_stiff_output(i,loc_elong_common_max) = elong_common_max;
        all_stiff_output(i,loc_strain_common_cut) = strain_common_cut;
        all_stiff_output(i,loc_strain_common_max) = strain_common_max;
        all_stiff_output(i,loc_young_common_cut) = young_common_cut;
        all_stiff_output(i,loc_young_common_max) = young_common_max;
    end
    
    
    %% save matlab workspace to file
    save all_data_stiff_endloop

    
    %% IND: save array with individual variables to XLS

    if ispc
        filename_output = strcat('data_output/all_stiff_output_', datestr(now, 'yyyymmdd_HHMM'), '.xlsx');
        xlswrite(filename_output, all_stiff_output_head, 1, 'A1')
        xlswrite(filename_output, all_stiff_output_txt, 1, 'A2')
        xlswrite(filename_output, all_stiff_output, 1, 'E2')
    else
        filename_output = strcat('data_output/all_stiff_output_', datestr(now, 'yyyymmdd_HHMM'), '.csv');
        csvwrite(filename_output, all_stiff_output)
    end
        
    
    %% IND: write graphpad prism files
    % TODO? currently, the script must be run once for all GM trials + once
    % for all SOL trials, and filenames do not indicate what the contents
    % is
    
    if input_project == 2
        filename_basis = strcat('data_output/prism_stiff/all_stiff_', datestr(now, 'yyyymmdd_HHMM'), '_');

        % create output file header
        j = 1;
        prism_array_col1 = {'Subj'; all_stiff_output_txt{j,3}; all_stiff_output_txt{j+2,3} };

        % for each output variable (column):
        loc_relevant_var = 7; % first column that has data for prism. 7 = start of force
        no_subjects = size(all_stiff_output_txt,1)/4;
        for i = loc_relevant_var:size(all_stiff_output,2)

            % reset
            prism_array(1:3,1:(2*no_subjects)) = NaN;
            write_col = 1;

            % filename with output variable
            filename_variable = strcat(filename_basis, all_stiff_output_head{i+4}, '.xlsx');

            for j = 1:4:size(all_stiff_output_txt,1) % number of lines / trials

            % row 1 = subject no
            prism_array(1,write_col) = str2double(all_stiff_output_txt{j,1});
            prism_array(1,write_col+no_subjects) = str2double(all_stiff_output_txt{j,1});
            
            % rows 2-3 in columns 1 and no_subjects+1 = data
            prism_array(2,write_col) = all_stiff_output(j,i); % datamaster 1 = pre con
            prism_array(2,write_col+no_subjects) = all_stiff_output(j+1,i); % datamaster 2 = post con
            prism_array(3,write_col) = all_stiff_output(j+2,i); % datamaster 3 = pre str
            prism_array(3,write_col+no_subjects) = all_stiff_output(j+3,i); % datamaster 4 = post str
            write_col = write_col + 1;
            end

            % write header
            xlswrite(filename_variable, prism_array_col1, 1, 'A1')
            % write data
            xlswrite(filename_variable, prism_array, 1, 'B1')
        end
    end
    
    
    %% GRP: plot FITS: force-elong per group  //  mean figures
    % ind: plot the fit up until IND cutoff ELONGATION (x axis determines)
    % grp: mean fit curve, up until weakest subject FORCE
    %      reusing variable: stiff_common_force
    
    fit_elong_guess = 4; %VAR guessing approximate elongation (X value) for 2nd order stiffness fit Y values
    
    % create stiffness equation as cfit
    f = fittype('a*x^2+b*x+c');
    a = 0; % temporary values, will be filled/replaced inside for loop
    b = 0;
    c = 0;
    stiff_eq_group2 = cfit(f, a, b, c);
    
    if input_project == 1 
        %% PLOTS - BD study
        if BD_SOL_count > 0
            % prepare data
            BD_SOL_force_common = stiff_common_force; % min(stiffness_SOL_BD(:,5));
            BD_SOL_force_array = (0:forceintervals:BD_SOL_force_common)';
            BD_SOL_elong_fit(1:length(BD_SOL_force_array),1:BD_SOL_count) = NaN;
            BD_SOL_elongmax_mean = mean(stiffness_SOL_BD(:,4));
            BD_SOL_elongmax_SD = std(stiffness_SOL_BD(:,4));
            BD_SOL_forcemax_mean = mean(stiffness_SOL_BD(:,5));
            BD_SOL_forcemax_SD = std(stiffness_SOL_BD(:,5));
            % prepare plot
            plottitle = horzcat('Free AT, stiffness curve fit, dancers');
            figure('Name',plottitle)
            hold on
            % individual data
            for i = 1:BD_SOL_count
                stiff_eq_group2.a = stiffness_SOL_BD(i,1); % a 
                stiff_eq_group2.b = stiffness_SOL_BD(i,2); % b
                stiff_eq_group2.c = stiffness_SOL_BD(i,3); % c
                xlim([0,stiffness_SOL_BD(i,4)]);
                plot(stiff_eq_group2)
                BD_SOL_elong_fit(:,i) = solve_sec_poly(stiffness_SOL_BD(i,1), stiffness_SOL_BD(i,2), stiffness_SOL_BD(i,3), 0, BD_SOL_force_common, forceintervals, fit_elong_guess);
            end
            % mean data
            BD_SOL_elong_fit_mean = nanmean(BD_SOL_elong_fit,2); % nanmean if fits through zero, mean if firt not through zero (4 places in total)
            h1 = plot(BD_SOL_elong_fit_mean,BD_SOL_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(BD_SOL_elongmax_mean,BD_SOL_forcemax_mean,BD_SOL_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(BD_SOL_elongmax_mean,BD_SOL_forcemax_mean,BD_SOL_elongmax_SD, 'k.')
            % visual
            axis(axis_stiff_free)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            saveas(gcf, horzcat('data_plots_stiff/GRP_BD_freeAT.jpg'))
        end

        if BD_GM_count > 0
            % prepare data
            BD_GM_force_common = stiff_common_force;
            BD_GM_force_array = (0:forceintervals:BD_GM_force_common)';
            BD_GM_elong_fit(1:length(BD_GM_force_array),1:BD_GM_count) = NaN;
            BD_GM_elongmax_mean = mean(stiffness_GM_BD(:,4));
            BD_GM_elongmax_SD = std(stiffness_GM_BD(:,4));
            BD_GM_forcemax_mean = mean(stiffness_GM_BD(:,5));
            BD_GM_forcemax_SD = std(stiffness_GM_BD(:,5));
            % prepare plot
            plottitle = horzcat('Whole AT, stiffness curve fit, dancers');
            figure('Name',plottitle)
            hold on
            % individual data
            for i = 1:BD_GM_count
                stiff_eq_group2.a = stiffness_GM_BD(i,1); % a 
                stiff_eq_group2.b = stiffness_GM_BD(i,2); % b
                stiff_eq_group2.c = stiffness_GM_BD(i,3); % c
                xlim([0,stiffness_GM_BD(i,4)]);
                plot(stiff_eq_group2)
                BD_GM_elong_fit(:,i) = solve_sec_poly(stiffness_GM_BD(i,1), stiffness_GM_BD(i,2), stiffness_GM_BD(i,3), 0, BD_GM_force_common, forceintervals, fit_elong_guess);
            end
            % mean data
            BD_GM_elong_fit_mean = nanmean(BD_GM_elong_fit,2);
            h1 = plot(BD_GM_elong_fit_mean,BD_GM_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(BD_GM_elongmax_mean,BD_GM_forcemax_mean,BD_GM_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(BD_GM_elongmax_mean,BD_GM_forcemax_mean,BD_GM_elongmax_SD, 'k.')
            % visual
            axis(axis_stiff_whole)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            saveas(gcf, horzcat('data_plots_stiff/GRP_BD_wholeAT.jpg'))
        end

        if CON_SOL_count > 0
            % prepare data
            CON_SOL_force_common = stiff_common_force;
            CON_SOL_force_array = (0:forceintervals:CON_SOL_force_common)';
            CON_SOL_elong_fit(1:length(CON_SOL_force_array),1:CON_SOL_count) = NaN;
            CON_SOL_elongmax_mean = mean(stiffness_SOL_CON(:,4));
            CON_SOL_elongmax_SD = std(stiffness_SOL_CON(:,4));
            CON_SOL_forcemax_mean = mean(stiffness_SOL_CON(:,5));
            CON_SOL_forcemax_SD = std(stiffness_SOL_CON(:,5));
            % prepare plot
            plottitle = horzcat('Free AT, stiffness curve fit, controls');
            figure('Name',plottitle)
            hold on
            % individual data
            for i = 1:CON_SOL_count
                stiff_eq_group2.a = stiffness_SOL_CON(i,1); % a 
                stiff_eq_group2.b = stiffness_SOL_CON(i,2); % b
                stiff_eq_group2.c = stiffness_SOL_CON(i,3); % c
                xlim([0,stiffness_SOL_CON(i,4)]);
                plot(stiff_eq_group2)
                CON_SOL_elong_fit(:,i) = solve_sec_poly(stiffness_SOL_CON(i,1), stiffness_SOL_CON(i,2), stiffness_SOL_CON(i,3), 0, CON_SOL_force_common, forceintervals, fit_elong_guess);
            end
            % mean data
            CON_SOL_elong_fit_mean = nanmean(CON_SOL_elong_fit,2);
            h1 = plot(CON_SOL_elong_fit_mean,CON_SOL_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(CON_SOL_elongmax_mean,CON_SOL_forcemax_mean,CON_SOL_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(CON_SOL_elongmax_mean,CON_SOL_forcemax_mean,CON_SOL_elongmax_SD, 'k.')
            % visual
            axis(axis_stiff_free)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            saveas(gcf, horzcat('data_plots_stiff/GRP_CON_freeAT.jpg'))
        end

        if CON_GM_count > 0
            % prepare data
            CON_GM_force_common = stiff_common_force; 
            CON_GM_force_array = (0:forceintervals:CON_GM_force_common)';
            CON_GM_elong_fit(1:length(CON_GM_force_array),1:CON_GM_count) = NaN;
            CON_GM_elongmax_mean = mean(stiffness_GM_CON(:,4));
            CON_GM_elongmax_SD = std(stiffness_GM_CON(:,4));
            CON_GM_forcemax_mean = mean(stiffness_GM_CON(:,5));
            CON_GM_forcemax_SD = std(stiffness_GM_CON(:,5));
            % prepare plot
            plottitle = horzcat('Whole AT, stiffness curve fit, controls');
            figure('Name',plottitle)
            hold on
            % individual data
            for i = 1:CON_GM_count
                stiff_eq_group2.a = stiffness_GM_CON(i,1); % a 
                stiff_eq_group2.b = stiffness_GM_CON(i,2); % b
                stiff_eq_group2.c = stiffness_GM_CON(i,3); % c
                xlim([0,stiffness_GM_CON(i,4)]);
                plot(stiff_eq_group2)
                CON_GM_elong_fit(:,i) = solve_sec_poly(stiffness_GM_CON(i,1), stiffness_GM_CON(i,2), stiffness_GM_CON(i,3), 0, CON_GM_force_common, forceintervals, fit_elong_guess);
            end
            % mean data
            CON_GM_elong_fit_mean = nanmean(CON_GM_elong_fit,2);
            h1 = plot(CON_GM_elong_fit_mean,CON_GM_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(CON_GM_elongmax_mean,CON_GM_forcemax_mean,CON_GM_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(CON_GM_elongmax_mean,CON_GM_forcemax_mean,CON_GM_elongmax_SD, 'k.')
            % visual
            axis(axis_stiff_whole)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            saveas(gcf, horzcat('data_plots_stiff/GRP_CON_wholeAT.jpg'))
        end



        if BD_GM_count > 0 && BD_SOL_count > 0 && CON_GM_count > 0 && CON_SOL_count > 0 
            plottitle = horzcat('stiffness curves MEAN');
            figure('Name',plottitle)
            hold on
            % mean data
            h1 = plot(BD_SOL_elong_fit_mean,BD_SOL_force_array,'r','Linewidth',2);
            h2 = plot(BD_GM_elong_fit_mean,BD_GM_force_array,'r--','Linewidth',2);
            h3 = plot(CON_SOL_elong_fit_mean,CON_SOL_force_array,'b','Linewidth',2);
            h4 = plot(CON_GM_elong_fit_mean,CON_GM_force_array,'b--','Linewidth',2);
            % error bars
            errorbar(BD_SOL_elongmax_mean,BD_SOL_forcemax_mean,BD_SOL_forcemax_SD, 'ro', 'MarkerFaceColor', 'r', 'Markersize',4); % avg of ind max force/elong
            herrorbar(BD_SOL_elongmax_mean,BD_SOL_forcemax_mean,BD_SOL_elongmax_SD, 'r.')
            errorbar(BD_GM_elongmax_mean,BD_GM_forcemax_mean,BD_GM_forcemax_SD, 'ro', 'MarkerFaceColor', 'r', 'Markersize',4); % avg of ind max force/elong
            herrorbar(BD_GM_elongmax_mean,BD_GM_forcemax_mean,BD_GM_elongmax_SD, 'r.')
            errorbar(CON_SOL_elongmax_mean,CON_SOL_forcemax_mean,CON_SOL_forcemax_SD, 'bo', 'MarkerFaceColor', 'b', 'Markersize',4); % avg of ind max force/elong
            herrorbar(CON_SOL_elongmax_mean,CON_SOL_forcemax_mean,CON_SOL_elongmax_SD, 'b.')
            errorbar(CON_GM_elongmax_mean,CON_GM_forcemax_mean,CON_GM_forcemax_SD, 'bo', 'MarkerFaceColor', 'b', 'Markersize',4); % avg of ind max force/elong
            herrorbar(CON_GM_elongmax_mean,CON_GM_forcemax_mean,CON_GM_elongmax_SD, 'b.')
            % visual
            axis(axis_stiff_whole)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend([h1 h2 h3 h4], 'BD free AT', 'BD whole AT', 'CON free AT', 'CON whole AT', 'Location','Southeast')
            saveas(gcf, horzcat('data_plots_stiff/GRP_mean_all.jpg'))



            plottitle = horzcat('Free AT, stiffness curves MEAN');
            figure('Name',plottitle)
            hold on
            % mean data
            h1 = plot(BD_SOL_elong_fit_mean,BD_SOL_force_array,'r','Linewidth',2);
            h3 = plot(CON_SOL_elong_fit_mean,CON_SOL_force_array,'b','Linewidth',2);
            % error bars
            errorbar(BD_SOL_elongmax_mean,BD_SOL_forcemax_mean,BD_SOL_forcemax_SD, 'ro', 'MarkerFaceColor', 'r', 'Markersize',4); % avg of ind max force/elong
            herrorbar(BD_SOL_elongmax_mean,BD_SOL_forcemax_mean,BD_SOL_elongmax_SD, 'r.')
            errorbar(CON_SOL_elongmax_mean,CON_SOL_forcemax_mean,CON_SOL_forcemax_SD, 'bo', 'MarkerFaceColor', 'b', 'Markersize',4); % avg of ind max force/elong
            herrorbar(CON_SOL_elongmax_mean,CON_SOL_forcemax_mean,CON_SOL_elongmax_SD, 'b.')
            % visual
            axis(axis_stiff_free)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend([h1 h3], 'BD free AT', 'CON free AT', 'Location','Southeast')
            saveas(gcf, horzcat('data_plots_stiff/GRP_mean_freeAT.jpg'))


            plottitle = horzcat('Whole AT, stiffness curves MEAN');
            figure('Name',plottitle)
            hold on
            % mean data
            h2 = plot(BD_GM_elong_fit_mean,BD_GM_force_array,'r--','Linewidth',2);
            h4 = plot(CON_GM_elong_fit_mean,CON_GM_force_array,'b--','Linewidth',2);
            % error bars
            errorbar(BD_GM_elongmax_mean,BD_GM_forcemax_mean,BD_GM_forcemax_SD, 'ro', 'MarkerFaceColor', 'r', 'Markersize',4); % avg of ind max force/elong
            herrorbar(BD_GM_elongmax_mean,BD_GM_forcemax_mean,BD_GM_elongmax_SD, 'r.')
            errorbar(CON_GM_elongmax_mean,CON_GM_forcemax_mean,CON_GM_forcemax_SD, 'bo', 'MarkerFaceColor', 'b', 'Markersize',4); % avg of ind max force/elong
            herrorbar(CON_GM_elongmax_mean,CON_GM_forcemax_mean,CON_GM_elongmax_SD, 'b.')
            % visual
            axis(axis_stiff_whole)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend([h2 h4], 'BD whole AT', 'CON whole AT', 'Location','Southeast')
            saveas(gcf, horzcat('data_plots_stiff/GRP_mean_wholeAT.jpg'))
        end
        
    else
        %% plots - intervention %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if STR_PRE_SOL_count > 0
            % prepare data
            STR_PRE_SOL_force_common = stiff_common_force;
            STR_PRE_SOL_force_array = (0:forceintervals:STR_PRE_SOL_force_common)';
            STR_PRE_SOL_elong_fit(1:length(STR_PRE_SOL_force_array),1:STR_PRE_SOL_count) = NaN;
            STR_PRE_SOL_elongmax_mean = mean(stiffness_STR_PRE_SOL(:,4));
            STR_PRE_SOL_elongmax_SD = std(stiffness_STR_PRE_SOL(:,4));
            STR_PRE_SOL_forcemax_mean = mean(stiffness_STR_PRE_SOL(:,5));
            STR_PRE_SOL_forcemax_SD = std(stiffness_STR_PRE_SOL(:,5));
            % prepare plot
            plottitle = horzcat('Free AT, stiffness curve fit, STR PRE');
            figure('Name',plottitle)
            hold on
            % individual data
            for i = 1:STR_PRE_SOL_count
                stiff_eq_group2.a = stiffness_STR_PRE_SOL(i,1); % a 
                stiff_eq_group2.b = stiffness_STR_PRE_SOL(i,2); % b
                stiff_eq_group2.c = stiffness_STR_PRE_SOL(i,3); % c
                xlim([0,stiffness_STR_PRE_SOL(i,4)]);
                plot(stiff_eq_group2)
                STR_PRE_SOL_elong_fit(:,i) = solve_sec_poly(stiffness_STR_PRE_SOL(i,1), stiffness_STR_PRE_SOL(i,2), stiffness_STR_PRE_SOL(i,3), 0, STR_PRE_SOL_force_common, forceintervals, fit_elong_guess);
                if isnan(STR_PRE_SOL_elong_fit(3,i)) % checking 3rd value of force-elongation (normally @ 100N)
                    cprintf('red', horzcat('WARNING: ', STR_PRE_SOL_ID{i}, ' force-elong curve fit gives NaN after first value.\n Average across subjects will be corrupted.\n'))
                end
            % mean data
            STR_PRE_SOL_elong_fit_mean = nanmean(STR_PRE_SOL_elong_fit,2); % nanmean if fits through zero, mean if fit not through zero (4 places in total)
            end
            h1 = plot(STR_PRE_SOL_elong_fit_mean,STR_PRE_SOL_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(STR_PRE_SOL_elongmax_mean,STR_PRE_SOL_forcemax_mean,STR_PRE_SOL_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(STR_PRE_SOL_elongmax_mean,STR_PRE_SOL_forcemax_mean,STR_PRE_SOL_elongmax_SD, 'k.')
            % visual
            axis(axis_stiff_free)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            print(horzcat('data_plots_stiff/GRP_',plottitle),'-dpng')
        end

        if STR_POST_SOL_count > 0
            % prepare data
            STR_POST_SOL_force_common = stiff_common_force;
            STR_POST_SOL_force_array = (0:forceintervals:STR_POST_SOL_force_common)';
            STR_POST_SOL_elong_fit(1:length(STR_POST_SOL_force_array),1:STR_POST_SOL_count) = NaN;
            STR_POST_SOL_elongmax_mean = mean(stiffness_STR_POST_SOL(:,4));
            STR_POST_SOL_elongmax_SD = std(stiffness_STR_POST_SOL(:,4));
            STR_POST_SOL_forcemax_mean = mean(stiffness_STR_POST_SOL(:,5));
            STR_POST_SOL_forcemax_SD = std(stiffness_STR_POST_SOL(:,5));
            % prepare plot
            plottitle = horzcat('Free AT, stiffness curve fit, STR POST');
            figure('Name',plottitle)
            hold on
            % individual data
            for i = 1:STR_POST_SOL_count
                stiff_eq_group2.a = stiffness_STR_POST_SOL(i,1); % a 
                stiff_eq_group2.b = stiffness_STR_POST_SOL(i,2); % b
                stiff_eq_group2.c = stiffness_STR_POST_SOL(i,3); % c
                xlim([0,stiffness_STR_POST_SOL(i,4)]);
                plot(stiff_eq_group2)
                STR_POST_SOL_elong_fit(:,i) = solve_sec_poly(stiffness_STR_POST_SOL(i,1), stiffness_STR_POST_SOL(i,2), stiffness_STR_POST_SOL(i,3), 0, STR_POST_SOL_force_common, forceintervals, fit_elong_guess);
                if isnan(STR_POST_SOL_elong_fit(3,i)) % checking 3rd value of force-elongation (normally @ 100N)
                    cprintf('red', horzcat('WARNING: ', STR_POST_SOL_ID{i}, ' force-elong curve fit gives NaN after first value.\n Average across subjects will be corrupted.\n'))
                end

            end
            % mean data
            STR_POST_SOL_elong_fit_mean = nanmean(STR_POST_SOL_elong_fit,2);
            h1 = plot(STR_POST_SOL_elong_fit_mean,STR_POST_SOL_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(STR_POST_SOL_elongmax_mean,STR_POST_SOL_forcemax_mean,STR_POST_SOL_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(STR_POST_SOL_elongmax_mean,STR_POST_SOL_forcemax_mean,STR_POST_SOL_elongmax_SD, 'k.')
            % visual
            axis(axis_stiff_free)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            print(horzcat('data_plots_stiff/GRP_',plottitle),'-dpng')
        end
        
        if CON_PRE_SOL_count > 0
            % prepare data
            CON_PRE_SOL_force_common = stiff_common_force;
            CON_PRE_SOL_force_array = (0:forceintervals:CON_PRE_SOL_force_common)';
            CON_PRE_SOL_elong_fit(1:length(CON_PRE_SOL_force_array),1:CON_PRE_SOL_count) = NaN;
            CON_PRE_SOL_elongmax_mean = mean(stiffness_CON_PRE_SOL(:,4));
            CON_PRE_SOL_elongmax_SD = std(stiffness_CON_PRE_SOL(:,4));
            CON_PRE_SOL_forcemax_mean = mean(stiffness_CON_PRE_SOL(:,5));
            CON_PRE_SOL_forcemax_SD = std(stiffness_CON_PRE_SOL(:,5));
            % prepare plot
            plottitle = horzcat('Free AT, stiffness curve fit, CON PRE');
            figure('Name',plottitle)
            hold on
            % individual data
            for i = 1:CON_PRE_SOL_count
                stiff_eq_group2.a = stiffness_CON_PRE_SOL(i,1); % a 
                stiff_eq_group2.b = stiffness_CON_PRE_SOL(i,2); % b
                stiff_eq_group2.c = stiffness_CON_PRE_SOL(i,3); % c
                xlim([0,stiffness_CON_PRE_SOL(i,4)]);
                plot(stiff_eq_group2)
                CON_PRE_SOL_elong_fit(:,i) = solve_sec_poly(stiffness_CON_PRE_SOL(i,1), stiffness_CON_PRE_SOL(i,2), stiffness_CON_PRE_SOL(i,3), 0, CON_PRE_SOL_force_common, forceintervals, fit_elong_guess);
                if isnan(CON_PRE_SOL_elong_fit(3,i)) % checking 3rd value of force-elongation (normally @ 100N)
                    cprintf('red', horzcat('WARNING: ', CON_PRE_SOL_ID{i}, ' force-elong curve fit gives NaN after first value.\n Average across subjects will be corrupted.\n'))
                end
            end
            % mean data
            CON_PRE_SOL_elong_fit_mean = nanmean(CON_PRE_SOL_elong_fit,2);
            h1 = plot(CON_PRE_SOL_elong_fit_mean,CON_PRE_SOL_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(CON_PRE_SOL_elongmax_mean,CON_PRE_SOL_forcemax_mean,CON_PRE_SOL_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(CON_PRE_SOL_elongmax_mean,CON_PRE_SOL_forcemax_mean,CON_PRE_SOL_elongmax_SD, 'k.')
            % visual
            axis(axis_stiff_free)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            print(horzcat('data_plots_stiff/GRP_',plottitle),'-dpng')
        end

        if CON_POST_SOL_count > 0
            % prepare data
            CON_POST_SOL_force_common = stiff_common_force;
            CON_POST_SOL_force_array = (0:forceintervals:CON_POST_SOL_force_common)';
            CON_POST_SOL_elong_fit(1:length(CON_POST_SOL_force_array),1:CON_POST_SOL_count) = NaN;
            CON_POST_SOL_elongmax_mean = mean(stiffness_CON_POST_SOL(:,4));
            CON_POST_SOL_elongmax_SD = std(stiffness_CON_POST_SOL(:,4));
            CON_POST_SOL_forcemax_mean = mean(stiffness_CON_POST_SOL(:,5));
            CON_POST_SOL_forcemax_SD = std(stiffness_CON_POST_SOL(:,5));
            % prepare plot
            plottitle = horzcat('Free AT, stiffness curve fit, CON POST');
            figure('Name',plottitle)
            hold on
            % individual data
            for i = 1:CON_POST_SOL_count
                stiff_eq_group2.a = stiffness_CON_POST_SOL(i,1); % a 
                stiff_eq_group2.b = stiffness_CON_POST_SOL(i,2); % b
                stiff_eq_group2.c = stiffness_CON_POST_SOL(i,3); % c
                xlim([0,stiffness_CON_POST_SOL(i,4)]);
                plot(stiff_eq_group2)
                CON_POST_SOL_elong_fit(:,i) = solve_sec_poly(stiffness_CON_POST_SOL(i,1), stiffness_CON_POST_SOL(i,2), stiffness_CON_POST_SOL(i,3), 0, CON_POST_SOL_force_common, forceintervals, fit_elong_guess);
                if isnan(CON_POST_SOL_elong_fit(3,i)) % checking 3rd value of force-elongation (normally @ 100N)
                    cprintf('red', horzcat('WARNING: ', CON_POST_SOL_ID{i}, ' force-elong curve fit gives NaN after first value.\n Average across subjects will be corrupted.\n'))
                end
            end
            % mean data
            CON_POST_SOL_elong_fit_mean = nanmean(CON_POST_SOL_elong_fit,2);
            h1 = plot(CON_POST_SOL_elong_fit_mean,CON_POST_SOL_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(CON_POST_SOL_elongmax_mean,CON_POST_SOL_forcemax_mean,CON_POST_SOL_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(CON_POST_SOL_elongmax_mean,CON_POST_SOL_forcemax_mean,CON_POST_SOL_elongmax_SD, 'k.')
            % visual
            axis(axis_stiff_free)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            print(horzcat('data_plots_stiff/GRP_',plottitle),'-dpng')
        end
        

        if STR_PRE_GM_count > 0
            % prepare data
            STR_PRE_GM_force_common = stiff_common_force;
            STR_PRE_GM_force_array = (0:forceintervals:STR_PRE_GM_force_common)';
            STR_PRE_GM_elong_fit(1:length(STR_PRE_GM_force_array),1:STR_PRE_GM_count) = NaN;
            STR_PRE_GM_elongmax_mean = mean(stiffness_STR_PRE_GM(:,4));
            STR_PRE_GM_elongmax_SD = std(stiffness_STR_PRE_GM(:,4));
            STR_PRE_GM_forcemax_mean = mean(stiffness_STR_PRE_GM(:,5));
            STR_PRE_GM_forcemax_SD = std(stiffness_STR_PRE_GM(:,5));
            % prepare plot
            plottitle = horzcat('Whole AT, stiffness curve fit, STR PRE');
            figure('Name',plottitle)
            hold on
            % individual data
            for i = 1:STR_PRE_GM_count
                stiff_eq_group2.a = stiffness_STR_PRE_GM(i,1); % a 
                stiff_eq_group2.b = stiffness_STR_PRE_GM(i,2); % b
                stiff_eq_group2.c = stiffness_STR_PRE_GM(i,3); % c
                xlim([0,stiffness_STR_PRE_GM(i,4)]);
                plot(stiff_eq_group2)
                STR_PRE_GM_elong_fit(:,i) = solve_sec_poly(stiffness_STR_PRE_GM(i,1), stiffness_STR_PRE_GM(i,2), stiffness_STR_PRE_GM(i,3), 0, STR_PRE_GM_force_common, forceintervals, fit_elong_guess);
                if isnan(STR_PRE_GM_elong_fit(3,i)) % checking 3rd value of force-elongation (normally @ 100N)
                    cprintf('red', horzcat('WARNING: ', STR_PRE_GM_ID{i}, ' force-elong curve fit gives NaN after first value.\n Average across subjects will be corrupted.\n'))
                end
            end
            % mean data
            STR_PRE_GM_elong_fit_mean = nanmean(STR_PRE_GM_elong_fit,2);
            h1 = plot(STR_PRE_GM_elong_fit_mean,STR_PRE_GM_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(STR_PRE_GM_elongmax_mean,STR_PRE_GM_forcemax_mean,STR_PRE_GM_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(STR_PRE_GM_elongmax_mean,STR_PRE_GM_forcemax_mean,STR_PRE_GM_elongmax_SD, 'k.')
            % visual
            axis(axis_stiff_whole)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            print(horzcat('data_plots_stiff/GRP_',plottitle),'-dpng')
        end

        if STR_POST_GM_count > 0
            % prepare data
            STR_POST_GM_force_common = stiff_common_force;
            STR_POST_GM_force_array = (0:forceintervals:STR_POST_GM_force_common)';
            STR_POST_GM_elong_fit(1:length(STR_POST_GM_force_array),1:STR_POST_GM_count) = NaN;
            STR_POST_GM_elongmax_mean = mean(stiffness_STR_POST_GM(:,4));
            STR_POST_GM_elongmax_SD = std(stiffness_STR_POST_GM(:,4));
            STR_POST_GM_forcemax_mean = mean(stiffness_STR_POST_GM(:,5));
            STR_POST_GM_forcemax_SD = std(stiffness_STR_POST_GM(:,5));
            % prepare plot
            plottitle = horzcat('Whole AT, stiffness curve fit, STR POST');
            figure('Name',plottitle)
            hold on
            % individual data
            for i = 1:STR_POST_GM_count
                stiff_eq_group2.a = stiffness_STR_POST_GM(i,1); % a 
                stiff_eq_group2.b = stiffness_STR_POST_GM(i,2); % b
                stiff_eq_group2.c = stiffness_STR_POST_GM(i,3); % c
                xlim([0,stiffness_STR_POST_GM(i,4)]);
                plot(stiff_eq_group2)
                STR_POST_GM_elong_fit(:,i) = solve_sec_poly(stiffness_STR_POST_GM(i,1), stiffness_STR_POST_GM(i,2), stiffness_STR_POST_GM(i,3), 0, STR_POST_GM_force_common, forceintervals, fit_elong_guess);
                if isnan(STR_POST_GM_elong_fit(3,i)) % checking 3rd value of force-elongation (normally @ 100N)
                    cprintf('red', horzcat('WARNING: ', STR_POST_GM_ID{i}, ' force-elong curve fit gives NaN after first value.\n Average across subjects will be corrupted.\n'))
                end
            end
            % mean data
            STR_POST_GM_elong_fit_mean = nanmean(STR_POST_GM_elong_fit,2);
            h1 = plot(STR_POST_GM_elong_fit_mean,STR_POST_GM_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(STR_POST_GM_elongmax_mean,STR_POST_GM_forcemax_mean,STR_POST_GM_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(STR_POST_GM_elongmax_mean,STR_POST_GM_forcemax_mean,STR_POST_GM_elongmax_SD, 'k.')
            % visual
            axis(axis_stiff_whole)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            print(horzcat('data_plots_stiff/GRP_',plottitle),'-dpng')
        end
        
        if CON_PRE_GM_count > 0
            % prepare data
            CON_PRE_GM_force_common = stiff_common_force;
            CON_PRE_GM_force_array = (0:forceintervals:CON_PRE_GM_force_common)';
            CON_PRE_GM_elong_fit(1:length(CON_PRE_GM_force_array),1:CON_PRE_GM_count) = NaN;
            CON_PRE_GM_elongmax_mean = mean(stiffness_CON_PRE_GM(:,4));
            CON_PRE_GM_elongmax_SD = std(stiffness_CON_PRE_GM(:,4));
            CON_PRE_GM_forcemax_mean = mean(stiffness_CON_PRE_GM(:,5));
            CON_PRE_GM_forcemax_SD = std(stiffness_CON_PRE_GM(:,5));
            % prepare plot
            plottitle = horzcat('Whole AT, stiffness curve fit, CON PRE');
            figure('Name',plottitle)
            hold on
            % individual data
            for i = 1:CON_PRE_GM_count
                stiff_eq_group2.a = stiffness_CON_PRE_GM(i,1); % a 
                stiff_eq_group2.b = stiffness_CON_PRE_GM(i,2); % b
                stiff_eq_group2.c = stiffness_CON_PRE_GM(i,3); % c
                xlim([0,stiffness_CON_PRE_GM(i,4)]);
                plot(stiff_eq_group2)
                CON_PRE_GM_elong_fit(:,i) = solve_sec_poly(stiffness_CON_PRE_GM(i,1), stiffness_CON_PRE_GM(i,2), stiffness_CON_PRE_GM(i,3), 0, CON_PRE_GM_force_common, forceintervals, fit_elong_guess);
                if isnan(CON_PRE_GM_elong_fit(3,i)) % checking 3rd value of force-elongation (normally @ 100N)
                    cprintf('red', horzcat('WARNING: ', CON_PRE_GM_ID{i}, ' force-elong curve fit gives NaN after first value.\n Average across subjects will be corrupted.\n'))
                end
            end
            % mean data
            CON_PRE_GM_elong_fit_mean = nanmean(CON_PRE_GM_elong_fit,2);
            h1 = plot(CON_PRE_GM_elong_fit_mean,CON_PRE_GM_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(CON_PRE_GM_elongmax_mean,CON_PRE_GM_forcemax_mean,CON_PRE_GM_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(CON_PRE_GM_elongmax_mean,CON_PRE_GM_forcemax_mean,CON_PRE_GM_elongmax_SD, 'k.')
            % visual
            axis(axis_stiff_whole)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            print(horzcat('data_plots_stiff/GRP_',plottitle),'-dpng')
        end

        if CON_POST_GM_count > 0
            % prepare data
            CON_POST_GM_force_common = stiff_common_force;
            CON_POST_GM_force_array = (0:forceintervals:CON_POST_GM_force_common)';
            CON_POST_GM_elong_fit(1:length(CON_POST_GM_force_array),1:CON_POST_GM_count) = NaN;
            CON_POST_GM_elongmax_mean = mean(stiffness_CON_POST_GM(:,4));
            CON_POST_GM_elongmax_SD = std(stiffness_CON_POST_GM(:,4));
            CON_POST_GM_forcemax_mean = mean(stiffness_CON_POST_GM(:,5));
            CON_POST_GM_forcemax_SD = std(stiffness_CON_POST_GM(:,5));
            % prepare plot
            plottitle = horzcat('Whole AT, stiffness curve fit, CON POST');
            figure('Name',plottitle)
            hold on
            % individual data
            for i = 1:CON_POST_GM_count
                stiff_eq_group2.a = stiffness_CON_POST_GM(i,1); % a 
                stiff_eq_group2.b = stiffness_CON_POST_GM(i,2); % b
                stiff_eq_group2.c = stiffness_CON_POST_GM(i,3); % c
                xlim([0,stiffness_CON_POST_GM(i,4)]);
                plot(stiff_eq_group2)
                CON_POST_GM_elong_fit(:,i) = solve_sec_poly(stiffness_CON_POST_GM(i,1), stiffness_CON_POST_GM(i,2), stiffness_CON_POST_GM(i,3), 0, CON_POST_GM_force_common, forceintervals, fit_elong_guess);
                if isnan(CON_POST_GM_elong_fit(3,i)) % checking 3rd value of force-elongation (normally @ 100N)
                    cprintf('red', horzcat('WARNING: ', CON_POST_GM_ID{i}, ' force-elong curve fit gives NaN after first value.\n Average across subjects will be corrupted.\n'))
                end
            end
            % mean data
            CON_POST_GM_elong_fit_mean = nanmean(CON_POST_GM_elong_fit,2);
            h1 = plot(CON_POST_GM_elong_fit_mean,CON_POST_GM_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(CON_POST_GM_elongmax_mean,CON_POST_GM_forcemax_mean,CON_POST_GM_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(CON_POST_GM_elongmax_mean,CON_POST_GM_forcemax_mean,CON_POST_GM_elongmax_SD, 'k.')
            % visual
            axis(axis_stiff_whole)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            print(horzcat('data_plots_stiff/GRP_',plottitle),'-dpng')
        end
        
        if STR_PRE_SOL_count > 0 && STR_POST_SOL_count > 0 && CON_PRE_SOL_count > 0 && CON_POST_SOL_count > 0
            plottitle = horzcat('Free AT, stiffness curves MEAN');
            figure('Name',plottitle)
            hold on
            % mean data
            h1 = plot(STR_PRE_SOL_elong_fit_mean,STR_PRE_SOL_force_array,'Color',col_lightred,'LineStyle','--','LineWidth',2);
            h2 = plot(STR_POST_SOL_elong_fit_mean,STR_POST_SOL_force_array,'r','LineStyle','-','LineWidth',2);
            h3 = plot(CON_PRE_SOL_elong_fit_mean,CON_PRE_SOL_force_array,'Color',col_lightblue,'LineStyle','--','LineWidth',2);
            h4 = plot(CON_POST_SOL_elong_fit_mean,CON_POST_SOL_force_array,'b','LineStyle','-','LineWidth',2);
            % error bars
            herrorbar(STR_PRE_SOL_elongmax_mean,STR_PRE_SOL_forcemax_mean,STR_PRE_SOL_elongmax_SD, '*m')
            errorbar(STR_PRE_SOL_elongmax_mean,STR_PRE_SOL_forcemax_mean,STR_PRE_SOL_forcemax_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred) % avg of ind max force/elong
            herrorbar(STR_POST_SOL_elongmax_mean,STR_POST_SOL_forcemax_mean,STR_POST_SOL_elongmax_SD, 'r.')
            errorbar(STR_POST_SOL_elongmax_mean,STR_POST_SOL_forcemax_mean,STR_POST_SOL_forcemax_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
            herrorbar(CON_PRE_SOL_elongmax_mean,CON_PRE_SOL_forcemax_mean,CON_PRE_SOL_elongmax_SD,  '*c')
            errorbar(CON_PRE_SOL_elongmax_mean,CON_PRE_SOL_forcemax_mean,CON_PRE_SOL_forcemax_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
            herrorbar(CON_POST_SOL_elongmax_mean,CON_POST_SOL_forcemax_mean,CON_POST_SOL_elongmax_SD, 'b.')
            errorbar(CON_POST_SOL_elongmax_mean,CON_POST_SOL_forcemax_mean,CON_POST_SOL_forcemax_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
            % visual
            axis(axis_stiff_free)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend([h1 h2 h3 h4], 'STR PRE', 'STR POST', 'CON PRE', 'CON POST', 'Location','Southeast')
            print(horzcat('data_plots_stiff/GRP_',plottitle),'-dpng')
        end
        
        if STR_PRE_GM_count > 0 && STR_POST_GM_count > 0 && CON_PRE_GM_count > 0 && CON_POST_GM_count > 0
            plottitle = horzcat('Whole AT, stiffness curves MEAN');
            figure('Name',plottitle)
            hold on
            % mean data
            h1 = plot(STR_PRE_GM_elong_fit_mean,STR_PRE_GM_force_array,'Color',col_lightred,'LineStyle','--','LineWidth',2);
            h2 = plot(STR_POST_GM_elong_fit_mean,STR_POST_GM_force_array,'r','LineStyle','-','LineWidth',2);
            h3 = plot(CON_PRE_GM_elong_fit_mean,CON_PRE_GM_force_array,'Color',col_lightblue,'LineStyle','--','LineWidth',2);
            h4 = plot(CON_POST_GM_elong_fit_mean,CON_POST_GM_force_array,'b','LineStyle','-','LineWidth',2);
            % error bars
            herrorbar(STR_PRE_GM_elongmax_mean,STR_PRE_GM_forcemax_mean,STR_PRE_GM_elongmax_SD, '*m')
            errorbar(STR_PRE_GM_elongmax_mean,STR_PRE_GM_forcemax_mean,STR_PRE_GM_forcemax_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred) % avg of ind max force/elong
            herrorbar(STR_POST_GM_elongmax_mean,STR_POST_GM_forcemax_mean,STR_POST_GM_elongmax_SD, 'r.')
            errorbar(STR_POST_GM_elongmax_mean,STR_POST_GM_forcemax_mean,STR_POST_GM_forcemax_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
            herrorbar(CON_PRE_GM_elongmax_mean,CON_PRE_GM_forcemax_mean,CON_PRE_GM_elongmax_SD,  '*c')
            errorbar(CON_PRE_GM_elongmax_mean,CON_PRE_GM_forcemax_mean,CON_PRE_GM_forcemax_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
            herrorbar(CON_POST_GM_elongmax_mean,CON_POST_GM_forcemax_mean,CON_POST_GM_elongmax_SD, 'b.')
            errorbar(CON_POST_GM_elongmax_mean,CON_POST_GM_forcemax_mean,CON_POST_GM_forcemax_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
            % visual
            axis(axis_stiff_whole)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend([h1 h2 h3 h4], 'STR PRE', 'STR POST', 'CON PRE', 'CON POST', 'Location','Southeast')
            print(horzcat('data_plots_stiff/GRP_',plottitle),'-dpng')
        end
                
    end
      
        
    %% IND + GRP: extract common force-elong, calculate group means
    % output = 
    %   group variables within common force region:
    %   force
    %   elong mean
    %   elong SD
    
    % find lowest force (shortest cell array among all groups and all SOL/GM)
    len = 10000;
    
    if input_project == 1 % BD study
        for i = 1:BD_SOL_count
            tmp_size = size(force_elong_SOL_BD{1,i});
            if tmp_size(1) < len
                len = tmp_size(1);
            end
        end
        for i = 1:BD_GM_count
            tmp_size = size(force_elong_GM_BD{1,i});
            if tmp_size(1) < len
                len = tmp_size(1);
            end
        end
        for i = 1:CON_SOL_count
            tmp_size = size(force_elong_SOL_CON{1,i});
            if tmp_size(1) < len
                len = tmp_size(1);
            end
        end
        for i = 1:CON_GM_count
            tmp_size = size(force_elong_GM_CON{1,i});
            if tmp_size(1) < len
                len = tmp_size(1);
            end
        end

        % collect elong up to common force
        if BD_SOL_count > 0
            BD_SOL_force = force_elong_SOL_BD{1,1}(1:len,2);
            BD_SOL_elong_data(1:len,1:BD_SOL_count) = NaN;
            for i = 1:length(force_elong_SOL_BD)
                BD_SOL_elong_data(:,i) = force_elong_SOL_BD{1,i}(1:len,1);
            end
            BD_SOL_elong_data_mean = nanmean(BD_SOL_elong_data,2);
            BD_SOL_elong_data_SD = nanstd(BD_SOL_elong_data,0,2);
        end
        if BD_GM_count > 0
            BD_GM_force = force_elong_GM_BD{1,1}(1:len,2);
            BD_GM_elong_data(1:len,1:BD_GM_count) = NaN;
            for i = 1:length(force_elong_GM_BD)
                BD_GM_elong_data(:,i) = force_elong_GM_BD{1,i}(1:len,1);
            end
            BD_GM_elong_data_mean = nanmean(BD_GM_elong_data,2);
            BD_GM_elong_data_SD = nanstd(BD_GM_elong_data,0,2);
        end
        if CON_SOL_count > 0
            CON_SOL_force = force_elong_SOL_CON{1,1}(1:len,2);
            CON_SOL_elong_data(1:len,1:CON_SOL_count) = NaN;
            for i = 1:length(force_elong_SOL_CON)
                CON_SOL_elong_data(:,i) = force_elong_SOL_CON{1,i}(1:len,1);
            end
            CON_SOL_elong_data_mean = nanmean(CON_SOL_elong_data,2);
            CON_SOL_elong_data_SD = nanstd(CON_SOL_elong_data,0,2);
        end
        if CON_GM_count > 0
            CON_GM_force = force_elong_GM_CON{1,1}(1:len,2);
            CON_GM_elong_data(1:len,1:CON_GM_count) = NaN;
            for i = 1:length(force_elong_GM_CON)
                CON_GM_elong_data(:,i) = force_elong_GM_CON{1,i}(1:len,1);
            end
            CON_GM_elong_data_mean = nanmean(CON_GM_elong_data,2);
            CON_GM_elong_data_SD = nanstd(CON_GM_elong_data,0,2);
        end
        
    else % intervention
        for i = 1:STR_PRE_SOL_count
            tmp_size = size(force_elong_STR_PRE_SOL{1,i});
            if tmp_size(1) < len
                len = tmp_size(1);
            end
        end
        for i = 1:STR_POST_SOL_count
            tmp_size = size(force_elong_STR_POST_SOL{1,i});
            if tmp_size(1) < len
                len = tmp_size(1);
            end
        end
        for i = 1:CON_PRE_SOL_count
            tmp_size = size(force_elong_CON_PRE_SOL{1,i});
            if tmp_size(1) < len
                len = tmp_size(1);
            end
        end
        for i = 1:CON_POST_SOL_count
            tmp_size = size(force_elong_CON_POST_SOL{1,i});
            if tmp_size(1) < len
                len = tmp_size(1);
            end
        end
        for i = 1:STR_PRE_GM_count
            tmp_size = size(force_elong_STR_PRE_GM{1,i});
            if tmp_size(1) < len
                len = tmp_size(1);
            end
        end
        for i = 1:STR_POST_GM_count
            tmp_size = size(force_elong_STR_POST_GM{1,i});
            if tmp_size(1) < len
                len = tmp_size(1);
            end
        end
        for i = 1:CON_PRE_GM_count
            tmp_size = size(force_elong_CON_PRE_GM{1,i});
            if tmp_size(1) < len
                len = tmp_size(1);
            end
        end
        for i = 1:CON_POST_GM_count
            tmp_size = size(force_elong_CON_POST_GM{1,i});
            if tmp_size(1) < len
                len = tmp_size(1);
            end
        end

        % collect elong up to common force
        if STR_PRE_SOL_count > 0
            STR_PRE_SOL_force = force_elong_STR_PRE_SOL{1,1}(1:len,2);
            STR_PRE_SOL_elong_data(1:len,1:STR_PRE_SOL_count) = NaN;
            for i = 1:length(force_elong_STR_PRE_SOL)
                STR_PRE_SOL_elong_data(:,i) = force_elong_STR_PRE_SOL{1,i}(1:len,1);
            end
            STR_PRE_SOL_elong_data_mean = nanmean(STR_PRE_SOL_elong_data,2);
            STR_PRE_SOL_elong_data_SD = nanstd(STR_PRE_SOL_elong_data,0,2);
        end
        if STR_POST_SOL_count > 0
            STR_POST_SOL_force = force_elong_STR_POST_SOL{1,1}(1:len,2);
            STR_POST_SOL_elong_data(1:len,1:STR_POST_SOL_count) = NaN;
            for i = 1:length(force_elong_STR_POST_SOL)
                STR_POST_SOL_elong_data(:,i) = force_elong_STR_POST_SOL{1,i}(1:len,1);
            end
            STR_POST_SOL_elong_data_mean = nanmean(STR_POST_SOL_elong_data,2);
            STR_POST_SOL_elong_data_SD = nanstd(STR_POST_SOL_elong_data,0,2);
        end

        if CON_PRE_SOL_count > 0
            CON_PRE_SOL_force = force_elong_CON_PRE_SOL{1,1}(1:len,2);
            CON_PRE_SOL_elong_data(1:len,1:CON_PRE_SOL_count) = NaN;
            for i = 1:length(force_elong_CON_PRE_SOL)
                CON_PRE_SOL_elong_data(:,i) = force_elong_CON_PRE_SOL{1,i}(1:len,1);
            end
            CON_PRE_SOL_elong_data_mean = nanmean(CON_PRE_SOL_elong_data,2);
            CON_PRE_SOL_elong_data_SD = nanstd(CON_PRE_SOL_elong_data,0,2);
        end
        if CON_POST_SOL_count > 0
            CON_POST_SOL_force = force_elong_CON_POST_SOL{1,1}(1:len,2);
            CON_POST_SOL_elong_data(1:len,1:CON_POST_SOL_count) = NaN;
            for i = 1:length(force_elong_CON_POST_SOL)
                CON_POST_SOL_elong_data(:,i) = force_elong_CON_POST_SOL{1,i}(1:len,1);
            end
            CON_POST_SOL_elong_data_mean = nanmean(CON_POST_SOL_elong_data,2);
            CON_POST_SOL_elong_data_SD = nanstd(CON_POST_SOL_elong_data,0,2);
        end

        if STR_PRE_GM_count > 0
            STR_PRE_GM_force = force_elong_STR_PRE_GM{1,1}(1:len,2);
            STR_PRE_GM_elong_data(1:len,1:STR_PRE_GM_count) = NaN;
            for i = 1:length(force_elong_STR_PRE_GM)
                STR_PRE_GM_elong_data(:,i) = force_elong_STR_PRE_GM{1,i}(1:len,1);
            end
            STR_PRE_GM_elong_data_mean = nanmean(STR_PRE_GM_elong_data,2);
            STR_PRE_GM_elong_data_SD = nanstd(STR_PRE_GM_elong_data,0,2);
        end
        if STR_POST_GM_count > 0
            STR_POST_GM_force = force_elong_STR_POST_GM{1,1}(1:len,2);
            STR_POST_GM_elong_data(1:len,1:STR_POST_GM_count) = NaN;
            for i = 1:length(force_elong_STR_POST_GM)
                STR_POST_GM_elong_data(:,i) = force_elong_STR_POST_GM{1,i}(1:len,1);
            end
            STR_POST_GM_elong_data_mean = nanmean(STR_POST_GM_elong_data,2);
            STR_POST_GM_elong_data_SD = nanstd(STR_POST_GM_elong_data,0,2);
        end

        if CON_PRE_GM_count > 0
            CON_PRE_GM_force = force_elong_CON_PRE_GM{1,1}(1:len,2);
            CON_PRE_GM_elong_data(1:len,1:CON_PRE_GM_count) = NaN;
            for i = 1:length(force_elong_CON_PRE_GM)
                CON_PRE_GM_elong_data(:,i) = force_elong_CON_PRE_GM{1,i}(1:len,1);
            end
            CON_PRE_GM_elong_data_mean = nanmean(CON_PRE_GM_elong_data,2);
            CON_PRE_GM_elong_data_SD = nanstd(CON_PRE_GM_elong_data,0,2);
        end
        if CON_POST_GM_count > 0
            CON_POST_GM_force = force_elong_CON_POST_GM{1,1}(1:len,2);
            CON_POST_GM_elong_data(1:len,1:CON_POST_GM_count) = NaN;
            for i = 1:length(force_elong_CON_POST_GM)
                CON_POST_GM_elong_data(:,i) = force_elong_CON_POST_GM{1,i}(1:len,1);
            end
            CON_POST_GM_elong_data_mean = nanmean(CON_POST_GM_elong_data,2);
            CON_POST_GM_elong_data_SD = nanstd(CON_POST_GM_elong_data,0,2);
        end
    end
    
        
    %% IND: save individual force-elong arrays (to common force) to one file common to all subjects
    if ispc
        filename_output = strcat('data_output/prism_stiff/arrays_stiff_across_force-elong_', datestr(now, 'yyyymmdd_HHMM'), '.xls');
        if input_project == 1 % bd study
            if BD_SOL_count > 0
                filename_subj_BD_SOL = cellstr(num2str(BD_SOL_no'))';
                xlswrite(filename_output, ['Force (N)' filename_subj_BD_SOL], 'Force-el data BD SOL', 'A1')
                xlswrite(filename_output, BD_SOL_force, 'Force-el data BD SOL', 'A2')
                xlswrite(filename_output, BD_SOL_elong_data, 'Force-el data BD SOL', 'B2')
                xlswrite(filename_output, ['Force (N)' filename_subj_BD_SOL], 'Force-el fit BD SOL', 'A1')
                xlswrite(filename_output, BD_SOL_force, 'Force-el fit BD SOL', 'A2')
                xlswrite(filename_output, BD_SOL_elong_fit, 'Force-el fit BD SOL', 'B2')
            end
            if CON_SOL_count > 0
                filename_subj_CON_SOL = cellstr(num2str(CON_SOL_no'))';
                xlswrite(filename_output, ['Force (N)' filename_subj_CON_SOL], 'Force-el data CON SOL', 'A1')
                xlswrite(filename_output, CON_SOL_force, 'Force-el data CON SOL', 'A2')
                xlswrite(filename_output, CON_SOL_elong_data, 'Force-el data CON SOL', 'B2')
                xlswrite(filename_output, ['Force (N)' filename_subj_CON_SOL], 'Force-el fit CON SOL', 'A1')
                xlswrite(filename_output, CON_SOL_force, 'Force-el fit CON SOL', 'A2')
                xlswrite(filename_output, CON_SOL_elong_fit, 'Force-el fit CON SOL', 'B2')
            end
            if BD_GM_count > 0
                filename_subj_BD_GM = cellstr(num2str(BD_GM_no'))';
                xlswrite(filename_output, ['Force (N)' filename_subj_BD_GM], 'Force-el data BD GM', 'A1')
                xlswrite(filename_output, BD_GM_force, 'Force-el data BD GM', 'A2')
                xlswrite(filename_output, BD_GM_elong_data, 'Force-el data BD GM', 'B2')
                xlswrite(filename_output, ['Force (N)' filename_subj_BD_GM], 'Force-el fit BD GM', 'A1')
                xlswrite(filename_output, BD_GM_force, 'Force-el fit BD GM', 'A2')
                xlswrite(filename_output, BD_GM_elong_fit, 'Force-el fit BD GM', 'B2')
            end
            if CON_GM_count > 0
                filename_subj_CON_GM = cellstr(num2str(CON_GM_no'))';
                xlswrite(filename_output, ['Force (N)' filename_subj_CON_GM], 'Force-el data CON GM', 'A1')
                xlswrite(filename_output, CON_GM_force, 'Force-el data CON GM', 'A2')
                xlswrite(filename_output, CON_GM_elong_data, 'Force-el data CON GM', 'B2')
                xlswrite(filename_output, ['Force (N)' filename_subj_CON_GM], 'Force-el fit CON GM', 'A1')
                xlswrite(filename_output, CON_GM_force, 'Force-el fit CON GM', 'A2')
                xlswrite(filename_output, CON_GM_elong_fit, 'Force-el fit CON GM', 'B2')
            end
            
        else % intervention
            if STR_PRE_SOL_count > 0
                filename_subj_STR_PRE_SOL = STR_PRE_SOL_ID;
                xlswrite(filename_output, ['Force (N)' filename_subj_STR_PRE_SOL], 'Force-el data STR PRE SOL', 'A1')
                xlswrite(filename_output, STR_PRE_SOL_force, 'Force-el data STR PRE SOL', 'A2')
                xlswrite(filename_output, STR_PRE_SOL_elong_data, 'Force-el data STR PRE SOL', 'B2')
                xlswrite(filename_output, ['Force (N)' filename_subj_STR_PRE_SOL], 'Force-el fit STR PRE SOL', 'A1')
                xlswrite(filename_output, STR_PRE_SOL_force, 'Force-el fit STR PRE SOL', 'A2')
                xlswrite(filename_output, STR_PRE_SOL_elong_fit, 'Force-el fit STR PRE SOL', 'B2')
            end
            if STR_POST_SOL_count > 0
                filename_subj_STR_POST_SOL = STR_POST_SOL_ID;
                xlswrite(filename_output, ['Force (N)' filename_subj_STR_POST_SOL], 'Force-el data STR POST SOL', 'A1')
                xlswrite(filename_output, STR_POST_SOL_force, 'Force-el data STR POST SOL', 'A2')
                xlswrite(filename_output, STR_POST_SOL_elong_data, 'Force-el data STR POST SOL', 'B2')
                xlswrite(filename_output, ['Force (N)' filename_subj_STR_POST_SOL], 'Force-el fit STR POST SOL', 'A1')
                xlswrite(filename_output, STR_POST_SOL_force, 'Force-el fit STR POST SOL', 'A2')
                xlswrite(filename_output, STR_POST_SOL_elong_fit, 'Force-el fit STR POST SOL', 'B2')
            end
            if CON_PRE_SOL_count > 0
                filename_subj_CON_PRE_SOL = CON_PRE_SOL_ID;
                xlswrite(filename_output, ['Force (N)' filename_subj_CON_PRE_SOL], 'Force-el data CON PRE SOL', 'A1')
                xlswrite(filename_output, CON_PRE_SOL_force, 'Force-el data CON PRE SOL', 'A2')
                xlswrite(filename_output, CON_PRE_SOL_elong_data, 'Force-el data CON PRE SOL', 'B2')
                xlswrite(filename_output, ['Force (N)' filename_subj_CON_PRE_SOL], 'Force-el fit CON PRE SOL', 'A1')
                xlswrite(filename_output, CON_PRE_SOL_force, 'Force-el fit CON PRE SOL', 'A2')
                xlswrite(filename_output, CON_PRE_SOL_elong_fit, 'Force-el fit CON PRE SOL', 'B2')
            end
            if CON_POST_SOL_count > 0
                filename_subj_CON_POST_SOL = CON_POST_SOL_ID;
                xlswrite(filename_output, ['Force (N)' filename_subj_CON_POST_SOL], 'Force-el data CON POST SOL', 'A1')
                xlswrite(filename_output, CON_POST_SOL_force, 'Force-el data CON POST SOL', 'A2')
                xlswrite(filename_output, CON_POST_SOL_elong_data, 'Force-el data CON POST SOL', 'B2')
                xlswrite(filename_output, ['Force (N)' filename_subj_CON_POST_SOL], 'Force-el fit CON POST SOL', 'A1')
                xlswrite(filename_output, CON_POST_SOL_force, 'Force-el fit CON POST SOL', 'A2')
                xlswrite(filename_output, CON_POST_SOL_elong_fit, 'Force-el fit CON POST SOL', 'B2')
            end

            if STR_PRE_GM_count > 0
                filename_subj_STR_PRE_GM = STR_PRE_GM_ID;
                xlswrite(filename_output, ['Force (N)' filename_subj_STR_PRE_GM], 'Force-el data STR PRE GM', 'A1')
                xlswrite(filename_output, STR_PRE_GM_force, 'Force-el data STR PRE GM', 'A2')
                xlswrite(filename_output, STR_PRE_GM_elong_data, 'Force-el data STR PRE GM', 'B2')
                xlswrite(filename_output, ['Force (N)' filename_subj_STR_PRE_GM], 'Force-el fit STR PRE GM', 'A1')
                xlswrite(filename_output, STR_PRE_GM_force, 'Force-el fit STR PRE GM', 'A2')
                xlswrite(filename_output, STR_PRE_GM_elong_fit, 'Force-el fit STR PRE GM', 'B2')
            end
            if STR_POST_GM_count > 0
                filename_subj_STR_POST_GM = STR_POST_GM_ID;
                xlswrite(filename_output, ['Force (N)' filename_subj_STR_POST_GM], 'Force-el data STR POST GM', 'A1')
                xlswrite(filename_output, STR_POST_GM_force, 'Force-el data STR POST GM', 'A2')
                xlswrite(filename_output, STR_POST_GM_elong_data, 'Force-el data STR POST GM', 'B2')
                xlswrite(filename_output, ['Force (N)' filename_subj_STR_POST_GM], 'Force-el fit STR POST GM', 'A1')
                xlswrite(filename_output, STR_POST_GM_force, 'Force-el fit STR POST GM', 'A2')
                xlswrite(filename_output, STR_POST_GM_elong_fit, 'Force-el fit STR POST GM', 'B2')
            end
            if CON_PRE_GM_count > 0
                filename_subj_CON_PRE_GM = CON_PRE_GM_ID;
                xlswrite(filename_output, ['Force (N)' filename_subj_CON_PRE_GM], 'Force-el data CON PRE GM', 'A1')
                xlswrite(filename_output, CON_PRE_GM_force, 'Force-el data CON PRE GM', 'A2')
                xlswrite(filename_output, CON_PRE_GM_elong_data, 'Force-el data CON PRE GM', 'B2')
                xlswrite(filename_output, ['Force (N)' filename_subj_CON_PRE_GM], 'Force-el fit CON PRE GM', 'A1')
                xlswrite(filename_output, CON_PRE_GM_force, 'Force-el fit CON PRE GM', 'A2')
                xlswrite(filename_output, CON_PRE_GM_elong_fit, 'Force-el fit CON PRE GM', 'B2')
            end
            if CON_POST_GM_count > 0
                filename_subj_CON_POST_GM = CON_POST_GM_ID;
                xlswrite(filename_output, ['Force (N)' filename_subj_CON_POST_GM], 'Force-el data CON POST GM', 'A1')
                xlswrite(filename_output, CON_POST_GM_force, 'Force-el data CON POST GM', 'A2')
                xlswrite(filename_output, CON_POST_GM_elong_data, 'Force-el data CON POST GM', 'B2')
                xlswrite(filename_output, ['Force (N)' filename_subj_CON_POST_GM], 'Force-el fit CON POST GM', 'A1')
                xlswrite(filename_output, CON_POST_GM_force, 'Force-el fit CON POST GM', 'A2')
                xlswrite(filename_output, CON_POST_GM_elong_fit, 'Force-el fit CON POST GM', 'B2')
            end
        end
    else
        % csvwrite
    end
    
    
    %% GRP: plot DATA, force-elong avg + mean (GM and SOL separately)
    
    if input_project == 1
        % Free AT (SOL):
        if (BD_SOL_count > 0 || CON_SOL_count > 0)
            fig_f_e_legend = [];
            plottitle = horzcat('Free AT, force-elongation (data, no fit)');
            figure('Name',plottitle)
            hold on
            if BD_SOL_count > 0
                bd1 = plot(BD_SOL_elong_data_mean, BD_SOL_force, 'r', 'LineWidth',2);
                fig_f_e_legend{end+1} = 'BD avg';
            end
            if CON_SOL_count > 0
                con1 = plot(CON_SOL_elong_data_mean, CON_SOL_force, 'b', 'LineWidth',2);
                fig_f_e_legend{end+1} = 'CON avg';
            end
            if BD_SOL_count > 0
                bd2 = plot(BD_SOL_elong_data_mean+BD_SOL_elong_data_SD, BD_SOL_force, 'r--', 'LineWidth',0.5);
                plot(BD_SOL_elong_data_mean-BD_SOL_elong_data_SD, BD_SOL_force, 'r--', 'LineWidth',0.5);
                fig_f_e_legend{end+1} = 'BD SD';
            end
            if CON_SOL_count > 0
                con2 = plot(CON_SOL_elong_data_mean+CON_SOL_elong_data_SD, CON_SOL_force, 'b--', 'LineWidth',0.5);
                plot(CON_SOL_elong_data_mean-CON_SOL_elong_data_SD, CON_SOL_force, 'b--', 'LineWidth',0.5);
                fig_f_e_legend{end+1} = 'CON SD';
            end
            %axis(axis_force)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            if CON_SOL_count > 0 && BD_SOL_count > 0
                legend([bd1 con1 bd2 con2],fig_f_e_legend, 'Location','Northwest')
            elseif BD_SOL_count > 0
                legend([bd1 bd2],fig_f_e_legend, 'Location','Northwest')
            else
                legend([con1 con2],fig_f_e_legend, 'Location','Northwest')
            end
            saveas(gcf, horzcat('data_plots_stiff/GRP_',plottitle,'.jpg'))
        end



        % WholeAT (GM):
        if (BD_GM_count > 0 || CON_GM_count > 0)
            fig_f_e_legend = [];
            plottitle = horzcat('Whole AT, force-elongation (data, no fit)');
            figure('Name',plottitle)
            hold on
            if BD_GM_count > 0
                bd1 = plot(BD_GM_elong_data_mean, BD_GM_force, 'r', 'LineWidth',2);
                fig_f_e_legend{end+1} = 'BD avg';
            end
            if CON_GM_count > 0
                con1 = plot(CON_GM_elong_data_mean, CON_GM_force, 'b', 'LineWidth',2);
                fig_f_e_legend{end+1} = 'CON avg';
            end
            if BD_GM_count > 0
                bd2 = plot(BD_GM_elong_data_mean+BD_GM_elong_data_SD, BD_GM_force, 'r--', 'LineWidth',0.5);
                plot(BD_GM_elong_data_mean-BD_GM_elong_data_SD, BD_GM_force, 'r--', 'LineWidth',0.5);
                fig_f_e_legend{end+1} = 'BD SD';
            end
            if CON_GM_count > 0
                con2 = plot(CON_GM_elong_data_mean+CON_GM_elong_data_SD, CON_GM_force, 'b--', 'LineWidth',0.5);
                plot(CON_GM_elong_data_mean-CON_GM_elong_data_SD, CON_GM_force, 'b--', 'LineWidth',0.5);
                fig_f_e_legend{end+1} = 'CON SD';
            end
            %axis(axis_force)
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            if CON_GM_count > 0 && BD_GM_count > 0
                legend([bd1 con1 bd2 con2],fig_f_e_legend, 'Location','Northwest')
            elseif BD_GM_count > 0
                legend([bd1 bd2],fig_f_e_legend, 'Location','Northwest')
            elseif CON_GM_count
                legend([con1 con2],fig_f_e_legend, 'Location','Northwest')
            end
            saveas(gcf, horzcat('data_plots_stiff/GRP_',plottitle,'.jpg'))
        end
        
    else % intervention 
        % Free AT (SOL), STR PRE-POST:
        if (STR_PRE_SOL_count > 0 || STR_POST_SOL_count > 0)
            fig_f_e_legend = [];
            plottitle = horzcat('Free AT, force-elongation (data, no fit) STR');
            figure('Name',plottitle)
            hold on
            if STR_PRE_SOL_count > 0
                fig_pre_1 = plot(STR_PRE_SOL_elong_data_mean, STR_PRE_SOL_force, 'r', 'LineWidth',2);
                fig_f_e_legend{end+1} = 'STR PRE avg';
            end
            if STR_POST_SOL_count > 0
                fig_post_1 = plot(STR_POST_SOL_elong_data_mean, STR_POST_SOL_force, 'b', 'LineWidth',2);
                fig_f_e_legend{end+1} = 'STR POST avg';
            end
            if STR_PRE_SOL_count > 0
                fig_pre_2 = plot(STR_PRE_SOL_elong_data_mean+STR_PRE_SOL_elong_data_SD, STR_PRE_SOL_force, 'r--', 'LineWidth',0.5);
                plot(STR_PRE_SOL_elong_data_mean-STR_PRE_SOL_elong_data_SD, STR_PRE_SOL_force, 'r--', 'LineWidth',0.5);
                fig_f_e_legend{end+1} = 'STR PRE SD';
            end
            if STR_POST_SOL_count > 0
                fig_post_2 = plot(STR_POST_SOL_elong_data_mean+STR_POST_SOL_elong_data_SD, STR_POST_SOL_force, 'b--', 'LineWidth',0.5);
                plot(STR_POST_SOL_elong_data_mean-STR_POST_SOL_elong_data_SD, STR_POST_SOL_force, 'b--', 'LineWidth',0.5);
                fig_f_e_legend{end+1} = 'STR POST SD';
            end
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            if STR_POST_SOL_count > 0 && STR_PRE_SOL_count > 0
                legend([fig_pre_1 fig_post_1 fig_pre_2 fig_post_2],fig_f_e_legend, 'Location','Northwest')
            elseif STR_PRE_SOL_count > 0
                legend([fig_pre_1 fig_pre_2],fig_f_e_legend, 'Location','Northwest')
            else
                legend([fig_post_1 fig_post_2],fig_f_e_legend, 'Location','Northwest')
            end
            print(horzcat('data_plots_stiff/GRP_',plottitle),'-dpng')
        end

        % Free AT (SOL), CON PRE-POST:
        if (CON_PRE_SOL_count > 0 || CON_POST_SOL_count > 0)
            fig_f_e_legend = [];
            plottitle = horzcat('Free AT, force-elongation (data, no fit) CON');
            figure('Name',plottitle)
            hold on
            if CON_PRE_SOL_count > 0
                fig_pre_1 = plot(CON_PRE_SOL_elong_data_mean, CON_PRE_SOL_force, 'r', 'LineWidth',2);
                fig_f_e_legend{end+1} = 'CON PRE avg';
            end
            if CON_POST_SOL_count > 0
                fig_post_1 = plot(CON_POST_SOL_elong_data_mean, CON_POST_SOL_force, 'b', 'LineWidth',2);
                fig_f_e_legend{end+1} = 'CON POST avg';
            end
            if CON_PRE_SOL_count > 0
                fig_pre_2 = plot(CON_PRE_SOL_elong_data_mean+CON_PRE_SOL_elong_data_SD, CON_PRE_SOL_force, 'r--', 'LineWidth',0.5);
                plot(CON_PRE_SOL_elong_data_mean-CON_PRE_SOL_elong_data_SD, CON_PRE_SOL_force, 'r--', 'LineWidth',0.5);
                fig_f_e_legend{end+1} = 'CON PRE SD';
            end
            if CON_POST_SOL_count > 0
                fig_post_2 = plot(CON_POST_SOL_elong_data_mean+CON_POST_SOL_elong_data_SD, CON_POST_SOL_force, 'b--', 'LineWidth',0.5);
                plot(CON_POST_SOL_elong_data_mean-CON_POST_SOL_elong_data_SD, CON_POST_SOL_force, 'b--', 'LineWidth',0.5);
                fig_f_e_legend{end+1} = 'CON POST SD';
            end
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            if CON_POST_SOL_count > 0 && CON_PRE_SOL_count > 0
                legend([fig_pre_1 fig_post_1 fig_pre_2 fig_post_2],fig_f_e_legend, 'Location','Northwest')
            elseif CON_PRE_SOL_count > 0
                legend([fig_pre_1 fig_pre_2],fig_f_e_legend, 'Location','Northwest')
            else
                legend([fig_post_1 fig_post_2],fig_f_e_legend, 'Location','Northwest')
            end
            print(horzcat('data_plots_stiff/GRP_',plottitle),'-dpng')
        end


        % Whole AT (GM), STR PRE-POST:
        if (STR_PRE_GM_count > 0 || STR_POST_GM_count > 0)
            fig_f_e_legend = [];
            plottitle = horzcat('Whole AT, force-elongation (data, no fit) STR');
            figure('Name',plottitle)
            hold on
            if STR_PRE_GM_count > 0
                fig_pre_1 = plot(STR_PRE_GM_elong_data_mean, STR_PRE_GM_force, 'r', 'LineWidth',2);
                fig_f_e_legend{end+1} = 'STR PRE avg';
            end
            if STR_POST_GM_count > 0
                fig_post_1 = plot(STR_POST_GM_elong_data_mean, STR_POST_GM_force, 'b', 'LineWidth',2);
                fig_f_e_legend{end+1} = 'STR POST avg';
            end
            if STR_PRE_GM_count > 0
                fig_pre_2 = plot(STR_PRE_GM_elong_data_mean+STR_PRE_GM_elong_data_SD, STR_PRE_GM_force, 'r--', 'LineWidth',0.5);
                plot(STR_PRE_GM_elong_data_mean-STR_PRE_GM_elong_data_SD, STR_PRE_GM_force, 'r--', 'LineWidth',0.5);
                fig_f_e_legend{end+1} = 'STR PRE SD';
            end
            if STR_POST_GM_count > 0
                fig_post_2 = plot(STR_POST_GM_elong_data_mean+STR_POST_GM_elong_data_SD, STR_POST_GM_force, 'b--', 'LineWidth',0.5);
                plot(STR_POST_GM_elong_data_mean-STR_POST_GM_elong_data_SD, STR_POST_GM_force, 'b--', 'LineWidth',0.5);
                fig_f_e_legend{end+1} = 'STR POST SD';
            end
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            if STR_POST_GM_count > 0 && STR_PRE_GM_count > 0
                legend([fig_pre_1 fig_post_1 fig_pre_2 fig_post_2],fig_f_e_legend, 'Location','Northwest')
            elseif STR_PRE_GM_count > 0
                legend([fig_pre_1 fig_pre_2],fig_f_e_legend, 'Location','Northwest')
            elseif STR_POST_GM_count
                legend([fig_post_1 fig_post_2],fig_f_e_legend, 'Location','Northwest')
            end
            print(horzcat('data_plots_stiff/GRP_',plottitle),'-dpng')
        end

        % Whole AT (GM), CON PRE-POST:
        if (CON_PRE_GM_count > 0 || CON_POST_GM_count > 0)
            fig_f_e_legend = [];
            plottitle = horzcat('Whole AT, force-elongation (data, no fit) CON');
            figure('Name',plottitle)
            hold on
            if CON_PRE_GM_count > 0
                fig_pre_1 = plot(CON_PRE_GM_elong_data_mean, CON_PRE_GM_force, 'r', 'LineWidth',2);
                fig_f_e_legend{end+1} = 'CON PRE avg';
            end
            if CON_POST_GM_count > 0
                fig_post_1 = plot(CON_POST_GM_elong_data_mean, CON_POST_GM_force, 'b', 'LineWidth',2);
                fig_f_e_legend{end+1} = 'CON POST avg';
            end
            if CON_PRE_GM_count > 0
                fig_pre_2 = plot(CON_PRE_GM_elong_data_mean+CON_PRE_GM_elong_data_SD, CON_PRE_GM_force, 'r--', 'LineWidth',0.5);
                plot(CON_PRE_GM_elong_data_mean-CON_PRE_GM_elong_data_SD, CON_PRE_GM_force, 'r--', 'LineWidth',0.5);
                fig_f_e_legend{end+1} = 'CON PRE SD';
            end
            if CON_POST_GM_count > 0
                fig_post_2 = plot(CON_POST_GM_elong_data_mean+CON_POST_GM_elong_data_SD, CON_POST_GM_force, 'b--', 'LineWidth',0.5);
                plot(CON_POST_GM_elong_data_mean-CON_POST_GM_elong_data_SD, CON_POST_GM_force, 'b--', 'LineWidth',0.5);
                fig_f_e_legend{end+1} = 'CON POST SD';
            end
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            if CON_POST_GM_count > 0 && CON_PRE_GM_count > 0
                legend([fig_pre_1 fig_post_1 fig_pre_2 fig_post_2],fig_f_e_legend, 'Location','Northwest')
            elseif CON_PRE_GM_count > 0
                legend([fig_pre_1 fig_pre_2],fig_f_e_legend, 'Location','Northwest')
            elseif CON_POST_GM_count
                legend([fig_post_1 fig_post_2],fig_f_e_legend, 'Location','Northwest')
            end
            print(horzcat('data_plots_stiff/GRP_',plottitle),'-dpng')
        end
        
    end
end