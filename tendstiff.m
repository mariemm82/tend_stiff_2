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

    
function [] = tendstiff(input_project, input_plot)
    close all

    forceintervals = 50; %VAR - Average stiffness across X N
    


    %% PLOTS - determine which plots to display
    global plot_achilles plot_norm plot_emg plot_check plot_us plot_conversion subject_id

    if input_plot >= 1 
        plot_check = 1; % turn on/off checkpoint plots (leave only stiffness)
    else
        plot_check = 0;
    end
    if input_plot >= 2
        plot_achilles = 1; % turn on/off all troubleshoot plots
    else
        plot_achilles = 0;
    end
    if input_plot >= 3
        plot_conversion = 1;
        plot_norm = 1; % show torque before and after initial lowpass filter / ankle rotation fit plots
    else
        plot_conversion = 0;
        plot_norm = 0; % show torque before and after initial lowpass filter / ankle rotation fit plots
    end
    plot_us = 0;
    plot_emg = 0;  % RMS 3 EMG channels per trial

    

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

    
    %% set AXES etc for plots
    col_lightblue = [0.6 0.8 1];
    col_lightred = [1 0.6 0.8];

    
    %% Read datamaster file, to connect corresponding data files
    % Produces arrays with file names and variables per trial, to be retrieved later

    global dm_subjectno dm_timepoint dm_side dm_trial dm_group
    global dm_stiff1_NX dm_stiff1_US dm_stiff1_US_frame dm_stiff2_NX dm_stiff2_US dm_stiff2_US_frame dm_stiff3_NX dm_stiff3_US dm_stiff3_US_frame 
    global dm_heel1_NX dm_heel1_US dm_heel1_US_frame dm_heel2_NX dm_heel2_US dm_heel2_US_frame dm_heel3_NX dm_heel3_US dm_heel3_US_frame
    global dm_MVC_PF dm_MVC_DF dm_CPM_calc_NX dm_CPM_sol_NX dm_CPM_calc_US dm_CPM_calc_US_frame dm_leg_length
    global dm_cutforce %new2014-04-14
    global filepath
    dm_filename = 'data/datamaster_stiff.tsv';
    linestotal = read_datamaster_stiff(dm_filename);

    
    %% preallocate
    % common arrays for all subjects:
    all_stiff_output_head = {'Subject', 'Time', 'Side', 'Trial', ...
        'Stiff coeff 1', 'Stiff coeff 2', 'Stiff coeff 3', 'Stiff R2', ...
        'Ramp force cutoff (N)', 'Ramp force max (N)', 'PF MVC (N)', ...
        'Stiffness ind 80-100 (N/mm)', 'Stiff ind 90-100', 'Stiff common cutoff 80-100', 'Stiff common cutoff 90-100', 'Stiff common max 80-100', 'Stiff common max 90-100', ...
         'AT moment arm (m)', 'Rotation correction (mm/deg)', 'Common force cutoff', 'Common force max'}; % PROJECTSPECIFIC
    all_stiff_output = zeros(ceil(linestotal),length(all_stiff_output_head)-4); 
    all_stiff_output_txt = cell(ceil(linestotal),4);

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

    
    %% Loop through all lines in datamaster file (except header line)
    for line = 1:linestotal


        %% subject/trial identifier
        trial_subjectno = str2double(dm_subjectno{line});
        
        if input_project == 1 % BD study
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
        [usdata_CPM, usfreq_CPM] = read_us_file(strcat(filepath, dm_CPM_calc_US{line}, '.txt'), str2double(dm_CPM_calc_US_frame{line}), 'CPM calcaneus');

        % Read moment arm noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
        % Produce a new noraxon data array
        noraxon_CPM = read_noraxon_stiffness(strcat(filepath, dm_CPM_calc_NX{line}), usfreq_CPM, dm_side{line}, 'CPM calcaneus');

        % above variables are no longer used for momentarm, but are (re-)used for ankle rotation below
        
        % Read complete, prepared noraon array + prepared us array
        % Produce AT moment arm constant
        at_momentarm = calculate_momentarm(0,0, dm_leg_length{line});


        %% Calculate relationship between calc displacement and ANKLE ROTATION
        % Read complete, prepared noraxon array + prepared us array
        % Produce ankle rotation constant, mm/degree
        
        % REUSING us file from moment arm
        % REUSING noraxon file from moment arm
        if strcmp(subject_id,'Control 4 R PRE SOL') || strcmp(subject_id,'INT_4_SOL_PRE_STR_R')
            % rough coding, manually calculated for this trial
            at_rotation_const = -0.14;
        else %normally
            at_rotation_const = calculate_rotation_correction(noraxon_CPM, usdata_CPM);
            % MMM TODO - same rot const pre and post??
        end


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
        %    sending 3+3 trials + max forces from trials + manually set cutoff force
        % Produce stiffness equation (upon full data set) +
        %    force-elong-array (up to defined force level of 90% of 6-trial-common-force or 90% of manual cutoff)
        [stiff_eq, stiff_gof, force_elong_array] = final_stiffness(time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3, time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3, forceintervals, dm_cutforce{line}, trial_force_max);
        % MMM TODO - is it right that stiff equation is made upon FULL
        % dataset (up to 100% of common force or cutoff force), while
        % stiffness is calculated in defined regions up to 90% of the
        % common/cutoff force? E.g. stiffness at 90-100% of 90% of manual
        % cutoff force

        %% calculate stiffness for last 10 and 20% of ind max:
        force_cutoff_ind = force_elong_array(end,2); % defined cutoff point of 90% of 6-trial-common-force or 90% of manually set force
        stiff_ind_80 = calculate_stiffness(stiff_eq, force_cutoff_ind, 0.8, 1.0, 'ind max'); % last two variables are percent range, from 0.00 to 1.00
        stiff_ind_90 = calculate_stiffness(stiff_eq, force_cutoff_ind, 0.9, 1.0, 'ind max');


        %% Other calculations? %LATER
        % strain_rate = calculate_strain_rate(time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3, time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3);
        % tendon elongation (mm)?
        % young's modulus?
        
        
        %% save individual data to common array
        % add all individual variables to a common array for all subjects    
        all_stiff_output_txt(line,:) = [dm_subjectno(line) dm_timepoint(line) dm_side(line) dm_trial(line)];
        % add NaN at locations for stiffness at force levels common to all subjects
        all_stiff_output(line,:) = [coeffvalues(stiff_eq) stiff_gof.rsquare ...
            force_elong_array(end,2) min(trial_force_max) plantflex_max_torque ... % 90%-of-common/manual-force / 100%-common-force / MVC-force
            stiff_ind_80 stiff_ind_90 NaN NaN NaN NaN at_momentarm at_rotation_const NaN NaN];
        
        
        %% add force-elongation arrays to cell for future averaging
        % add stiffness fits to array for future averaging, containing:
        %     2nd order equation coeffs / elongation max (= xlim for plots) / force max


         if input_project == 1 % BD study
            if trial_subjectno > 100 % BD subject
                if strcmp(dm_trial{line}, 'SOL')
                    force_elong_SOL_BD{BD_SOL_count} = force_elong_array;
                    stiffness_SOL_BD(BD_SOL_count,:) = [coeffvalues(stiff_eq), force_elong_array(end,1), force_elong_array(end,2)];
                else % GM
                    force_elong_GM_BD{BD_GM_count} = force_elong_array;
                    stiffness_GM_BD(BD_GM_count,:) = [coeffvalues(stiff_eq), force_elong_array(end,1), force_elong_array(end,2)];
                end
            else % CON subject
                if strcmp(dm_trial{line}, 'SOL') == 1
                    force_elong_SOL_CON{CON_SOL_count} = force_elong_array;
                    stiffness_SOL_CON(CON_SOL_count,:) = [coeffvalues(stiff_eq), force_elong_array(end,1), force_elong_array(end,2)];
                else % GM
                    force_elong_GM_CON{CON_GM_count} = force_elong_array;
                    stiffness_GM_CON(CON_GM_count,:) = [coeffvalues(stiff_eq), force_elong_array(end,1), force_elong_array(end,2)];
                end
            end
         else % intervention study
            if trial_location == 0 % SOL
                if trial_timepoint == 0 && trial_group == 1 % PRE, STR
                    force_elong_STR_PRE_SOL{STR_PRE_SOL_count} = force_elong_array;
                    stiffness_STR_PRE_SOL(STR_PRE_SOL_count,:) = [coeffvalues(stiff_eq), force_elong_array(end,1), force_elong_array(end,2)];
                elseif trial_timepoint == 1 && trial_group == 1 % POST, STR
                    force_elong_STR_POST_SOL{STR_POST_SOL_count} = force_elong_array;
                    stiffness_STR_POST_SOL(STR_POST_SOL_count,:) = [coeffvalues(stiff_eq), force_elong_array(end,1), force_elong_array(end,2)];
                elseif trial_timepoint == 0 && trial_group == 0 % PRE, CON
                    force_elong_CON_PRE_SOL{CON_PRE_SOL_count} = force_elong_array;
                    stiffness_CON_PRE_SOL(CON_PRE_SOL_count,:) = [coeffvalues(stiff_eq), force_elong_array(end,1), force_elong_array(end,2)];
                elseif trial_timepoint == 1 && trial_group == 0 % POST, CON
                    force_elong_CON_POST_SOL{CON_POST_SOL_count} = force_elong_array;
                    stiffness_CON_POST_SOL(CON_POST_SOL_count,:) = [coeffvalues(stiff_eq), force_elong_array(end,1), force_elong_array(end,2)];
                end
            else % GM
                if trial_timepoint == 0 && trial_group == 1 % PRE, STR
                    force_elong_STR_PRE_GM{STR_PRE_GM_count} = force_elong_array;
                    stiffness_STR_PRE_GM(STR_PRE_GM_count,:) = [coeffvalues(stiff_eq), force_elong_array(end,1), force_elong_array(end,2)];
                elseif trial_timepoint == 1 && trial_group == 1 % POST, STR
                    force_elong_STR_POST_GM{STR_POST_GM_count} = force_elong_array;
                    stiffness_STR_POST_GM(STR_POST_GM_count,:) = [coeffvalues(stiff_eq), force_elong_array(end,1), force_elong_array(end,2)];
                elseif trial_timepoint == 0 && trial_group == 0 % PRE, CON
                    force_elong_CON_PRE_GM{CON_PRE_GM_count} = force_elong_array;
                    stiffness_CON_PRE_GM(CON_PRE_GM_count,:) = [coeffvalues(stiff_eq), force_elong_array(end,1), force_elong_array(end,2)];
                elseif trial_timepoint == 1 && trial_group == 0 % POST, CON
                    force_elong_CON_POST_GM{CON_POST_GM_count} = force_elong_array;
                    stiffness_CON_POST_GM(CON_POST_GM_count,:) = [coeffvalues(stiff_eq), force_elong_array(end,1), force_elong_array(end,2)];
                end
            end
         end
         save all_data_stiff_inloop
    end
    %% LOOP finished %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    save all_data_stiff
    
    
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
    
    
    %% IND: calculate stiffness at group common force levels
    
    % select common force (maximal force reached by all subjects)
    all_stiff_col = 5; % 5 for cutoff force level
    stiff_common_force = min(all_stiff_output(:,all_stiff_col));
    all_stiff_col = 6; % 6 for max force level
    stiff_common_force_max = min(all_stiff_output(:,all_stiff_col));
    
    % create stiffness equation as cfit
    f = fittype('a*x^2+b*x+c');
    a = 0; % temporary values, will be filled/replaced inside for loop
    b = 0;
    c = 0;
    stiff_eq_group = cfit(f, a, b, c);

    for i = 1:size(all_stiff_output,1)
        stiff_eq_group.a = all_stiff_output(i,1); % a 
        stiff_eq_group.b = all_stiff_output(i,2); % b
        stiff_eq_group.c = all_stiff_output(i,3); % c

        % calculate stiffness
        stiff_common_80 = calculate_stiffness(stiff_eq_group, stiff_common_force, 0.8, 1.0, horzcat('FP', all_stiff_output_txt{i,1}, ' ', all_stiff_output_txt{i,4}, ' common cutoff force')); % last two = percent of submitted force - %VAR
        stiff_common_90 = calculate_stiffness(stiff_eq_group, stiff_common_force, 0.9, 1.0, horzcat('FP', all_stiff_output_txt{i,1}, ' ', all_stiff_output_txt{i,4}, ' common cutoff force')); %VAR
        stiff_common_80_max = calculate_stiffness(stiff_eq_group, stiff_common_force_max, 0.8, 1.0, horzcat('FP', all_stiff_output_txt{i,1}, ' ', all_stiff_output_txt{i,4}, ' common max force')); % last two = percent of submitted force - %VAR
        stiff_common_90_max = calculate_stiffness(stiff_eq_group, stiff_common_force_max, 0.9, 1.0, horzcat('FP', all_stiff_output_txt{i,1}, ' ', all_stiff_output_txt{i,4}, ' common max force')); %VAR

        % add to array across subjects
        all_stiff_output(i,10) = stiff_common_80;
        all_stiff_output(i,11) = stiff_common_90;
        all_stiff_output(i,12) = stiff_common_80_max;
        all_stiff_output(i,13) = stiff_common_90_max;
        all_stiff_output(i,16) = stiff_common_force;
        all_stiff_output(i,17) = stiff_common_force_max;
    end
    
    cprintf('blue*',horzcat('Stiffness: Common cutoff force = ', num2str(stiff_common_force), ' N, common max force = ', num2str(round(stiff_common_force_max,0)), ' N.\n'))

 
    %% IND: save array with individual variables to file

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
      
    
    %% GRP: plot stiff FIT lines per subject (groupwise)
    % ind: plot the fit up until ind cutoff force/cutoff elong
    % grp: mean fit curve, up until weakest subject
    
    fit_elong_guess = 4; %VAR guessing approximate elongation (X value) for 2nd order stiffness fit Y values
    
    if input_project == 1 % BD study
        % 5 for cutoff force level
        allgroups_force_common = min( [min(stiffness_SOL_BD(:,5)) min(stiffness_GM_BD(:,5)) min(stiffness_SOL_CON(:,5)) min(stiffness_GM_CON(:,5))] );
    else % intervention
        allgroups_force_common = min( [...
            min(stiffness_STR_PRE_SOL(:,5)) min(stiffness_STR_PRE_GM(:,5)) min(stiffness_CON_PRE_SOL(:,5)) min(stiffness_CON_PRE_GM(:,5)) ...
            min(stiffness_STR_POST_SOL(:,5)) min(stiffness_STR_POST_GM(:,5)) min(stiffness_CON_POST_SOL(:,5)) min(stiffness_CON_POST_GM(:,5)) ...
            ] );
        tmp_forcelevels = [min(stiffness_STR_PRE_SOL(:,5)) min(stiffness_STR_PRE_GM(:,5)) min(stiffness_CON_PRE_SOL(:,5)) min(stiffness_CON_PRE_GM(:,5)) ...
            min(stiffness_STR_POST_SOL(:,5)) min(stiffness_STR_POST_GM(:,5)) min(stiffness_CON_POST_SOL(:,5)) min(stiffness_CON_POST_GM(:,5))]% TMP MMM;
    end
    
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
            BD_SOL_force_common = allgroups_force_common; % min(stiffness_SOL_BD(:,5));
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
            axis([0 14 -100 3600])
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle)
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            saveas(gcf, horzcat('data_plots_stiff/GRP_stiff_fit_BD_freeAT.jpg'))
        end

        if BD_GM_count > 0
            % prepare data
            BD_GM_force_common = allgroups_force_common;
            BD_GM_force_array = (0:forceintervals:BD_GM_force_common)';
            BD_GM_elong_fit(1:length(BD_GM_force_array),1:BD_GM_count) = NaN;
            BD_GM_elongmax_mean = mean(stiffness_GM_BD(:,4));
            BD_GM_elongmax_SD = std(stiffness_GM_BD(:,4));
            BD_GM_forcemax_mean = mean(stiffness_GM_BD(:,5));
            BD_GM_forcemax_SD = std(stiffness_GM_BD(:,5));
            % prepare plot
            plottitle = horzcat('Entire AT, stiffness curve fit, dancers');
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
            axis([0 24 -100 3600])
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle)
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            saveas(gcf, horzcat('data_plots_stiff/GRP_stiff_fit_BD_entireAT.jpg'))
        end

        if CON_SOL_count > 0
            % prepare data
            CON_SOL_force_common = allgroups_force_common;
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
            axis([0 14 -100 3600])
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle)
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            saveas(gcf, horzcat('data_plots_stiff/GRP_stiff_fit_CON_freeAT.jpg'))
        end

        if CON_GM_count > 0
            % prepare data
            CON_GM_force_common = allgroups_force_common; 
            CON_GM_force_array = (0:forceintervals:CON_GM_force_common)';
            CON_GM_elong_fit(1:length(CON_GM_force_array),1:CON_GM_count) = NaN;
            CON_GM_elongmax_mean = mean(stiffness_GM_CON(:,4));
            CON_GM_elongmax_SD = std(stiffness_GM_CON(:,4));
            CON_GM_forcemax_mean = mean(stiffness_GM_CON(:,5));
            CON_GM_forcemax_SD = std(stiffness_GM_CON(:,5));
            % prepare plot
            plottitle = horzcat('Entire AT, stiffness curve fit, controls');
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
            axis([0 24 -100 3600])
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle)
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            saveas(gcf, horzcat('data_plots_stiff/GRP_stiff_fit_CON_entireAT.jpg'))
        end



        if BD_GM_count > 0 && BD_SOL_count > 0 && CON_GM_count > 0 && CON_SOL_count > 0 
            plottitle = horzcat('Mean stiffness curves');
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
            axis([0 22 0 3600])
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle)
            legend([h1 h2 h3 h4], 'BD free AT', 'BD entire AT', 'CON free AT', 'CON entire AT', 'Location','Southeast')
            saveas(gcf, horzcat('data_plots_stiff/GRP_stiff_fit_mean_all.jpg'))



            plottitle = horzcat('Free AT, mean stiffness curves');
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
            axis([0 12 0 3600])
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle)
            legend([h1 h3], 'BD free AT', 'CON free AT', 'Location','Southeast')
            saveas(gcf, horzcat('data_plots_stiff/GRP_stiff_fit_mean_freeAT.jpg'))


            plottitle = horzcat('Entire AT, mean stiffness curves');
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
            axis([0 22 0 3600])
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle)
            legend([h2 h4], 'BD entire AT', 'CON entire AT', 'Location','Southeast')
            saveas(gcf, horzcat('data_plots_stiff/GRP_stiff_fit_mean_entireAT.jpg'))
        end
        
    else
        %% plots - intervention %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if STR_PRE_SOL_count > 0
            % prepare data
            STR_PRE_SOL_force_common = allgroups_force_common;
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
            end
            % mean data
            STR_PRE_SOL_elong_fit_mean = nanmean(STR_PRE_SOL_elong_fit,2); % nanmean if fits through zero, mean if firt not through zero (4 places in total)
            h1 = plot(STR_PRE_SOL_elong_fit_mean,STR_PRE_SOL_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(STR_PRE_SOL_elongmax_mean,STR_PRE_SOL_forcemax_mean,STR_PRE_SOL_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(STR_PRE_SOL_elongmax_mean,STR_PRE_SOL_forcemax_mean,STR_PRE_SOL_elongmax_SD, 'k.')
            % visual
            axis([0 14 0 4000])
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle)
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            print(horzcat('data_plots_stiff/GRP_stiff_',plottitle),'-dpng')
        end

        if STR_POST_SOL_count > 0
            % prepare data
            STR_POST_SOL_force_common = allgroups_force_common;
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
            end
            % mean data
            STR_POST_SOL_elong_fit_mean = nanmean(STR_POST_SOL_elong_fit,2);
            h1 = plot(STR_POST_SOL_elong_fit_mean,STR_POST_SOL_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(STR_POST_SOL_elongmax_mean,STR_POST_SOL_forcemax_mean,STR_POST_SOL_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(STR_POST_SOL_elongmax_mean,STR_POST_SOL_forcemax_mean,STR_POST_SOL_elongmax_SD, 'k.')
            % visual
            axis([0 14 0 4000])
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle)
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            print(horzcat('data_plots_stiff/GRP_stiff_',plottitle),'-dpng')
        end
        
        if CON_PRE_SOL_count > 0
            % prepare data
            CON_PRE_SOL_force_common = allgroups_force_common;
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
            end
            % mean data
            CON_PRE_SOL_elong_fit_mean = nanmean(CON_PRE_SOL_elong_fit,2);
            h1 = plot(CON_PRE_SOL_elong_fit_mean,CON_PRE_SOL_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(CON_PRE_SOL_elongmax_mean,CON_PRE_SOL_forcemax_mean,CON_PRE_SOL_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(CON_PRE_SOL_elongmax_mean,CON_PRE_SOL_forcemax_mean,CON_PRE_SOL_elongmax_SD, 'k.')
            % visual
            axis([0 14 0 4000])
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle)
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            print(horzcat('data_plots_stiff/GRP_stiff_',plottitle),'-dpng')
        end

        if CON_POST_SOL_count > 0
            % prepare data
            CON_POST_SOL_force_common = allgroups_force_common;
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
            end
            % mean data
            CON_POST_SOL_elong_fit_mean = nanmean(CON_POST_SOL_elong_fit,2);
            h1 = plot(CON_POST_SOL_elong_fit_mean,CON_POST_SOL_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(CON_POST_SOL_elongmax_mean,CON_POST_SOL_forcemax_mean,CON_POST_SOL_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(CON_POST_SOL_elongmax_mean,CON_POST_SOL_forcemax_mean,CON_POST_SOL_elongmax_SD, 'k.')
            % visual
            axis([0 14 0 4000])
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle)
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            print(horzcat('data_plots_stiff/GRP_stiff_',plottitle),'-dpng')
        end
        

        if STR_PRE_GM_count > 0
            % prepare data
            STR_PRE_GM_force_common = allgroups_force_common;
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
            end
            % mean data
            STR_PRE_GM_elong_fit_mean = nanmean(STR_PRE_GM_elong_fit,2);
            h1 = plot(STR_PRE_GM_elong_fit_mean,STR_PRE_GM_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(STR_PRE_GM_elongmax_mean,STR_PRE_GM_forcemax_mean,STR_PRE_GM_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(STR_PRE_GM_elongmax_mean,STR_PRE_GM_forcemax_mean,STR_PRE_GM_elongmax_SD, 'k.')
            % visual
            axis([0 24 0 4000])
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle)
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            print(horzcat('data_plots_stiff/GRP_stiff_',plottitle),'-dpng')
        end

        if STR_POST_GM_count > 0
            % prepare data
            STR_POST_GM_force_common = allgroups_force_common;
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
            end
            % mean data
            STR_POST_GM_elong_fit_mean = nanmean(STR_POST_GM_elong_fit,2);
            h1 = plot(STR_POST_GM_elong_fit_mean,STR_POST_GM_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(STR_POST_GM_elongmax_mean,STR_POST_GM_forcemax_mean,STR_POST_GM_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(STR_POST_GM_elongmax_mean,STR_POST_GM_forcemax_mean,STR_POST_GM_elongmax_SD, 'k.')
            % visual
            axis([0 24 0 4000])
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle)
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            print(horzcat('data_plots_stiff/GRP_stiff_',plottitle),'-dpng')
        end
        
        if CON_PRE_GM_count > 0
            % prepare data
            CON_PRE_GM_force_common = allgroups_force_common;
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
            end
            % mean data
            CON_PRE_GM_elong_fit_mean = nanmean(CON_PRE_GM_elong_fit,2);
            h1 = plot(CON_PRE_GM_elong_fit_mean,CON_PRE_GM_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(CON_PRE_GM_elongmax_mean,CON_PRE_GM_forcemax_mean,CON_PRE_GM_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(CON_PRE_GM_elongmax_mean,CON_PRE_GM_forcemax_mean,CON_PRE_GM_elongmax_SD, 'k.')
            % visual
            axis([0 24 0 4000])
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle)
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            print(horzcat('data_plots_stiff/GRP_stiff_',plottitle),'-dpng')
        end

        if CON_POST_GM_count > 0
            % prepare data
            CON_POST_GM_force_common = allgroups_force_common;
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
            end
            % mean data
            CON_POST_GM_elong_fit_mean = nanmean(CON_POST_GM_elong_fit,2);
            h1 = plot(CON_POST_GM_elong_fit_mean,CON_POST_GM_force_array,'k','Linewidth',2); % average curve
            h2 = errorbar(CON_POST_GM_elongmax_mean,CON_POST_GM_forcemax_mean,CON_POST_GM_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',4); % avg of ind max force/elong
            herrorbar(CON_POST_GM_elongmax_mean,CON_POST_GM_forcemax_mean,CON_POST_GM_elongmax_SD, 'k.')
            % visual
            axis([0 24 0 4000])
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle)
            legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
            print(horzcat('data_plots_stiff/GRP_stiff_',plottitle),'-dpng')
        end
        
        if STR_PRE_SOL_count > 0 && STR_POST_SOL_count > 0 && CON_PRE_SOL_count > 0 && CON_POST_SOL_count > 0
            plottitle = horzcat('Free AT, mean stiffness curves');
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
            axis([0 14 0 4000])
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle)
            legend([h1 h2 h3 h4], 'STR PRE', 'STR POST', 'CON PRE', 'CON POST', 'Location','Southeast')
            print(horzcat('data_plots_stiff/GRP_stiff_',plottitle),'-dpng')
        end
        
        if STR_PRE_GM_count > 0 && STR_POST_GM_count > 0 && CON_PRE_GM_count > 0 && CON_POST_GM_count > 0
            plottitle = horzcat('Whole AT, mean stiffness curves');
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
            axis([0 24 0 4000])
            xlabel('Tendon elongation (mm)')
            ylabel('Force (N)')
            title(plottitle)
            legend([h1 h2 h3 h4], 'STR PRE', 'STR POST', 'CON PRE', 'CON POST', 'Location','Southeast')
            print(horzcat('data_plots_stiff/GRP_stiff_',plottitle),'-dpng')
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
        filename_output = strcat('data_output/stiff_GRP_force-elong_', datestr(now, 'yyyy-mm-dd HH-MM'), '.xls');
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
    
    
    %% GRP: plot force-elong per subject + mean (GM and SOL separately)
    
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
            title(plottitle)
            if CON_SOL_count > 0 && BD_SOL_count > 0
                legend([bd1 con1 bd2 con2],fig_f_e_legend, 'Location','Northwest')
            elseif BD_SOL_count > 0
                legend([bd1 bd2],fig_f_e_legend, 'Location','Northwest')
            else
                legend([con1 con2],fig_f_e_legend, 'Location','Northwest')
            end
            saveas(gcf, horzcat('data_plots_stiff/GRP_',plottitle,'.jpg'))
        end



        % Entire AT (GM):
        if (BD_GM_count > 0 || CON_GM_count > 0)
            fig_f_e_legend = [];
            plottitle = horzcat('Entire AT, force-elongation (data, no fit)');
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
            title(plottitle)
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
            plottitle = horzcat('Free AT, force-elongation (data, no fit)');
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
            title(plottitle)
            if STR_POST_SOL_count > 0 && STR_PRE_SOL_count > 0
                legend([fig_pre_1 fig_post_1 fig_pre_2 fig_post_2],fig_f_e_legend, 'Location','Northwest')
            elseif STR_PRE_SOL_count > 0
                legend([fig_pre_1 fig_pre_2],fig_f_e_legend, 'Location','Northwest')
            else
                legend([fig_post_1 fig_post_2],fig_f_e_legend, 'Location','Northwest')
            end
            saveas(gcf, horzcat('data_plots_stiff/GRP_',plottitle,'.jpg'))
        end

        % Free AT (SOL), CON PRE-POST:
        if (CON_PRE_SOL_count > 0 || CON_POST_SOL_count > 0)
            fig_f_e_legend = [];
            plottitle = horzcat('Free AT, force-elongation (data, no fit)');
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
            title(plottitle)
            if CON_POST_SOL_count > 0 && CON_PRE_SOL_count > 0
                legend([fig_pre_1 fig_post_1 fig_pre_2 fig_post_2],fig_f_e_legend, 'Location','Northwest')
            elseif CON_PRE_SOL_count > 0
                legend([fig_pre_1 fig_pre_2],fig_f_e_legend, 'Location','Northwest')
            else
                legend([fig_post_1 fig_post_2],fig_f_e_legend, 'Location','Northwest')
            end
            saveas(gcf, horzcat('data_plots_stiff/GRP_',plottitle,'.jpg'))
        end


        % Entire AT (GM), STR PRE-POST:
        if (STR_PRE_GM_count > 0 || STR_POST_GM_count > 0)
            fig_f_e_legend = [];
            plottitle = horzcat('Entire AT, force-elongation (data, no fit)');
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
            title(plottitle)
            if STR_POST_GM_count > 0 && STR_PRE_GM_count > 0
                legend([fig_pre_1 fig_post_1 fig_pre_2 fig_post_2],fig_f_e_legend, 'Location','Northwest')
            elseif STR_PRE_GM_count > 0
                legend([fig_pre_1 fig_pre_2],fig_f_e_legend, 'Location','Northwest')
            elseif STR_POST_GM_count
                legend([fig_post_1 fig_post_2],fig_f_e_legend, 'Location','Northwest')
            end
            saveas(gcf, horzcat('data_plots_stiff/GRP_',plottitle,'.jpg'))
        end

        % Entire AT (GM), CON PRE-POST:
        if (CON_PRE_GM_count > 0 || CON_POST_GM_count > 0)
            fig_f_e_legend = [];
            plottitle = horzcat('Entire AT, force-elongation (data, no fit)');
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
            title(plottitle)
            if CON_POST_GM_count > 0 && CON_PRE_GM_count > 0
                legend([fig_pre_1 fig_post_1 fig_pre_2 fig_post_2],fig_f_e_legend, 'Location','Northwest')
            elseif CON_PRE_GM_count > 0
                legend([fig_pre_1 fig_pre_2],fig_f_e_legend, 'Location','Northwest')
            elseif CON_POST_GM_count
                legend([fig_post_1 fig_post_2],fig_f_e_legend, 'Location','Northwest')
            end
            saveas(gcf, horzcat('data_plots_stiff/GRP_',plottitle,'.jpg'))
        end
        
    end
end