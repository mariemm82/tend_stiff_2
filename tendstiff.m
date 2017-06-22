%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main file for tendon stiffness
% Marie Moltubakk 17.5.2013
% 
% Run for all trials for which stiffness should be computed at common
% level + force-elongation plot, together.
% Separates SOL from GM trials.
% 
% Note 1:
% The scripts assume a certain structure for input files from Tracker and
% Noraxon. I.e. number of and order of EMG and other channels. If these
% scripts are to be used for other projects, some of the code must be
% modified. These lines are marked % PROJECTSPECIFIC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
function [] = tendstiff(input_project, input_plot)
    close all



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

    % Average stiffness across X N
    forceintervals = 50; %VAR
    
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
    angle_cutoff = 20/(noraxonfreq/2); % 9/(noraxonfreq/2); % cutoff freq, Norm angle PASSIVE - ref Winter 1990 = 15 hz. Kongsgaard = 8hz
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


    %% Read datamaster file, to connect corresponding data files
    % Produces arrays with file names and variables per trial, to be retrieved later

    global dm_subjectno dm_timepoint dm_side dm_trial 
    global dm_stiff1_NX dm_stiff1_US dm_stiff1_US_frame dm_stiff2_NX dm_stiff2_US dm_stiff2_US_frame dm_stiff3_NX dm_stiff3_US dm_stiff3_US_frame 
    global dm_heel1_NX dm_heel1_US dm_heel1_US_frame dm_heel2_NX dm_heel2_US dm_heel2_US_frame dm_heel3_NX dm_heel3_US dm_heel3_US_frame
    global dm_MVC_PF dm_MVC_DF dm_CPM_calc_NX dm_CPM_calc_US dm_CPM_calc_US_frame dm_leg_length
    global dm_cutforce %new2014-04-14
    global filepath
    dm_filename = 'data/datamaster_stiff.tsv';
    dm_columns = 29; % number of data columns entered per subject % PROJECTSPECIFIC %new2014-04-14 - was 28
    linestotal = read_datamaster(dm_filename,dm_columns);


    %% preallocate
    %preallocate output arrays
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
        BD_GM_count = 0; % # of ballet dancer subjects
        CON_GM_count = 0; % # of controls = intervention study subjects
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
    else
        %LATER
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
                    BD_SOL_no(BD_SOL_count) = str2double(dm_subjectno{line});
                else % GM
                    BD_GM_count = BD_GM_count + 1;
                    BD_GM_no(BD_GM_count) = str2double(dm_subjectno{line});
                end
            else
                filepath = 'data\';
                subject_id = horzcat('Control ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
                if strcmp(dm_trial{line},'SOL')
                    CON_SOL_count = CON_SOL_count + 1;
                    CON_SOL_no(CON_SOL_count) = str2double(dm_subjectno{line});
                else % GM
                    CON_GM_count = CON_GM_count + 1;
                    CON_GM_no(CON_GM_count) = str2double(dm_subjectno{line});
                end
            end
        elseif input_project == 2 % intervention
            filepath = 'data\';
            subject_id = horzcat('Intervent ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
            %LATER
%             if trial_timepoint == 0 && trial_leg == 0 % PRE, CON
%                 CON_PRE_count = CON_PRE_count + 1;
%                 CON_PRE_no(CON_PRE_count) = str2double(dm_subjectno{line});
%                 CON_PRE_subject_ID(CON_PRE_count) = trial_subjectno;
%             elseif trial_timepoint == 0 && trial_leg == 1 % PRE, STR
%                 STR_PRE_count = STR_PRE_count + 1;
%                 STR_PRE_no(STR_PRE_count) = str2double(dm_subjectno{line});
%             elseif trial_timepoint == 1 && trial_leg == 0 % POST, CON
%                 CON_POST_count = CON_POST_count + 1;
%                 CON_POST_no(CON_POST_count) = str2double(dm_subjectno{line});
%             elseif trial_timepoint == 1 && trial_leg == 1 % POST, STR
%                 STR_POST_count = STR_POST_count + 1;
%                 STR_POST_no(STR_POST_count) = str2double(dm_subjectno{line});
%             end
        end
        
        cprintf('*black', horzcat('----------------', subject_id, '------------------\n'))

        

        %% Calculate conversion factors for tibialis anterior CO-ACTIVATION

        % Produce individual conversion factors for angle
        [convert_norm_angle_a, convert_norm_angle_b] = calculate_angle_constants(angle_cutoff, horzcat(filepath, dm_CPM_calc_NX{line}), dm_side{line});
        % using only "angle_constants" not "constants_ACTIVE" because this is for NORM data only, and tendstiff uses only passive trials from NORM

        % Read co-activation noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
        % Produce a new noraxon data array
        % sending in a length corresponding to 9 seconds (delete anything after 9 sec)
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

        % above variables are no longer used for momentarm, but are (re-)used for ankle rot
        
        % Read complete, prepared noraon array + prepared us array
        % Produce AT moment arm constant
        at_momentarm = calculate_momentarm(0,0, dm_leg_length{line});


        %% Calculate relationship between calc displacement and ANKLE ROTATION

        % REUSING us file from moment arm
        % REUSING noraxon file from moment arm

        % Read complete, prepared noraxon array + prepared us array
        % Produce ankle rotation constant, mm/degree
        if strcmp(subject_id,'Control 4 R PRE SOL') %very special case for this trial...
            at_rotation_const = -0.14;
        else %normally
            at_rotation_const = calculate_rotation_correction(noraxon_CPM, usdata_CPM);
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
        % Produce stiffness equation
        [stiff_eq, stiff_gof, force_elong_array] = final_stiffness(time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3, time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3, forceintervals, dm_cutforce{line}, trial_force_max, at_momentarm);


        %% calculate stiffness for last 10 and 20% of ind max:
        stiff_ind_80 = calculate_stiffness(stiff_eq, force_elong_array(end,2), 0.8, 1.0, 'ind max'); % last two variables are percent range, from 0.00 to 1.00
        stiff_ind_90 = calculate_stiffness(stiff_eq, force_elong_array(end,2), 0.9, 1.0, 'ind max'); %force_elong... = cutoff force level


        %% Strain rate
        %LATER
        % strain_rate = calculate_strain_rate(time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3, time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3);
        
        
        %% other calculations from Melina datasheet?
        %LATER
        % output tendon elongation / strain / young's modulus?
        
    
        %% save individual data to common array + force-elong-file

        % add all individual variables to a common array for all subjects    
        all_stiff_output_txt(line,:) = [dm_subjectno(line) dm_timepoint(line) dm_side(line) dm_trial(line)];
        % NaN at locations for stiffness at force levels common to all subjects
        all_stiff_output(line,:) = [coeffvalues(stiff_eq) stiff_gof.rsquare force_elong_array(end,2) min(trial_force_max) plantflex_max_torque stiff_ind_80 stiff_ind_90 NaN NaN NaN NaN at_momentarm at_rotation_const NaN NaN];

%         % add force-elongation arrays to xls file per subject
%         if ispc
%             filename_output = strcat('data_output/stiff_forceoutput_', subject_id, '.xls');
%             xlswrite(filename_output, {'Elongation (mm)','Force (N)'}, 'Final stiff', 'A1')
%             xlswrite(filename_output, force_elong_array, 'Final stiff', 'A2')
%         else
%             filename_output = strcat('data_output/stiff_forceoutput_', subject_id, '_averaged.csv');
%             csvwrite(filename_output, force_elong_array)
%         end
        
        
        
        %% add force-elongation arrays to cell for future averaging
        % add stiffness fits to array for future averaging

         if input_project == 1 % BD study
            if trial_subjectno > 100 % BD subject
                if strcmp(dm_trial{line}, 'SOL')
                    force_elong_SOL_BD{BD_SOL_count} = force_elong_array;
                    stiffness_SOL_BD(BD_SOL_count,:) = [coeffvalues(stiff_eq), force_elong_array(end,1), force_elong_array(end,2)]; % 2nd order equation coeffs, elongation max (= xlim), force max
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
         else
             %LATER
         end
         
    end
    %% LOOP finished %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    save loop_end
    
    
    %% truncate cells
    force_elong_SOL_BD = force_elong_SOL_BD(~cellfun('isempty',force_elong_SOL_BD));
    force_elong_GM_BD = force_elong_GM_BD(~cellfun('isempty',force_elong_GM_BD));
    force_elong_SOL_CON = force_elong_SOL_CON(~cellfun('isempty',force_elong_SOL_CON));
    force_elong_GM_CON = force_elong_GM_CON(~cellfun('isempty',force_elong_GM_CON));
    stiffness_SOL_BD(BD_SOL_count+1:end,:) = [];
    stiffness_GM_BD(BD_GM_count+1:end,:) = [];
    stiffness_SOL_CON(CON_SOL_count+1:end,:) = [];
    stiffness_GM_CON(CON_GM_count+1:end,:) = [];
    
    
    %% IND: calculate stiffness at group common force levels
    
    % select common force (maximal force reached by all subjects)
    all_stiff_col = 5; % 5 for cutoff force level, 6 for max force level
    stiff_common_force = min(all_stiff_output(:,all_stiff_col));
    all_stiff_col = 6; % 5 for cutoff force level, 6 for max force level
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
    
    allgroups_force_common = min( [min(stiffness_SOL_BD(:,5)) min(stiffness_GM_BD(:,5)) min(stiffness_SOL_CON(:,5)) min(stiffness_GM_CON(:,5))] );
    
    % create stiffness equation as cfit
    f = fittype('a*x^2+b*x+c');
    a = 0; % temporary values, will be filled/replaced inside for loop
    b = 0;
    c = 0;
    stiff_eq_group2 = cfit(f, a, b, c);
    
    if plot_check && BD_SOL_count > 0
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
        h2 = errorbar(BD_SOL_elongmax_mean,BD_SOL_forcemax_mean,BD_SOL_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',6); % avg of ind max force/elong
        herrorbar(BD_SOL_elongmax_mean,BD_SOL_forcemax_mean,BD_SOL_elongmax_SD, 'ko')
        % visual
        axis([0 14 -100 3600])
        xlabel('Tendon elongation (mm)')
        ylabel('Force (N)')
        title(plottitle)
        legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
        saveas(gcf, horzcat('data_plots_stiff/GRP_stiff_fit_BD_freeAT.jpg'))
    end
    
    if plot_check && BD_GM_count > 0
        % prepare data
        BD_GM_force_common = allgroups_force_common; % min(stiffness_GM_BD(:,5));
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
        h2 = errorbar(BD_GM_elongmax_mean,BD_GM_forcemax_mean,BD_GM_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',6); % avg of ind max force/elong
        herrorbar(BD_GM_elongmax_mean,BD_GM_forcemax_mean,BD_GM_elongmax_SD, 'ko')
        % visual
        axis([0 24 -100 3600])
        xlabel('Tendon elongation (mm)')
        ylabel('Force (N)')
        title(plottitle)
        legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
        saveas(gcf, horzcat('data_plots_stiff/GRP_stiff_fit_BD_entireAT.jpg'))
    end
    
    if plot_check && CON_SOL_count > 0
        % prepare data
        CON_SOL_force_common = allgroups_force_common; % min(stiffness_SOL_CON(:,5));
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
        h2 = errorbar(CON_SOL_elongmax_mean,CON_SOL_forcemax_mean,CON_SOL_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',6); % avg of ind max force/elong
        herrorbar(CON_SOL_elongmax_mean,CON_SOL_forcemax_mean,CON_SOL_elongmax_SD, 'ko')
        % visual
        axis([0 14 -100 3600])
        xlabel('Tendon elongation (mm)')
        ylabel('Force (N)')
        title(plottitle)
        legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
        saveas(gcf, horzcat('data_plots_stiff/GRP_stiff_fit_CON_freeAT.jpg'))
    end
    
    if plot_check && CON_GM_count > 0
        % prepare data
        CON_GM_force_common = allgroups_force_common; % min(stiffness_GM_CON(:,5));
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
        h2 = errorbar(CON_GM_elongmax_mean,CON_GM_forcemax_mean,CON_GM_forcemax_SD, 'ko', 'MarkerFaceColor', 'k', 'Markersize',6); % avg of ind max force/elong
        herrorbar(CON_GM_elongmax_mean,CON_GM_forcemax_mean,CON_GM_elongmax_SD, 'ko')
        % visual
        axis([0 24 -100 3600])
        xlabel('Tendon elongation (mm)')
        ylabel('Force (N)')
        title(plottitle)
        legend([h1 h2], 'Mean curve', 'Ind max', 'Location','Northwest')
        saveas(gcf, horzcat('data_plots_stiff/GRP_stiff_fit_CON_entireAT.jpg'))
    end
    
    if plot_check && BD_GM_count > 0 && BD_SOL_count > 0 && CON_GM_count > 0 && CON_SOL_count > 0 
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
        herrorbar(BD_SOL_elongmax_mean,BD_SOL_forcemax_mean,BD_SOL_elongmax_SD, 'ro')
        errorbar(BD_GM_elongmax_mean,BD_GM_forcemax_mean,BD_GM_forcemax_SD, 'ro', 'MarkerFaceColor', 'r', 'Markersize',4); % avg of ind max force/elong
        herrorbar(BD_GM_elongmax_mean,BD_GM_forcemax_mean,BD_GM_elongmax_SD, 'ro')
        errorbar(CON_SOL_elongmax_mean,CON_SOL_forcemax_mean,CON_SOL_forcemax_SD, 'bo', 'MarkerFaceColor', 'b', 'Markersize',4); % avg of ind max force/elong
        herrorbar(CON_SOL_elongmax_mean,CON_SOL_forcemax_mean,CON_SOL_elongmax_SD, 'bo')
        errorbar(CON_GM_elongmax_mean,CON_GM_forcemax_mean,CON_GM_forcemax_SD, 'bo', 'MarkerFaceColor', 'b', 'Markersize',4); % avg of ind max force/elong
        herrorbar(CON_GM_elongmax_mean,CON_GM_forcemax_mean,CON_GM_elongmax_SD, 'bo')
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
        herrorbar(BD_SOL_elongmax_mean,BD_SOL_forcemax_mean,BD_SOL_elongmax_SD, 'ro')
        errorbar(CON_SOL_elongmax_mean,CON_SOL_forcemax_mean,CON_SOL_forcemax_SD, 'bo', 'MarkerFaceColor', 'b', 'Markersize',4); % avg of ind max force/elong
        herrorbar(CON_SOL_elongmax_mean,CON_SOL_forcemax_mean,CON_SOL_elongmax_SD, 'bo')
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
        herrorbar(BD_GM_elongmax_mean,BD_GM_forcemax_mean,BD_GM_elongmax_SD, 'ro')
        errorbar(CON_GM_elongmax_mean,CON_GM_forcemax_mean,CON_GM_forcemax_SD, 'bo', 'MarkerFaceColor', 'b', 'Markersize',4); % avg of ind max force/elong
        herrorbar(CON_GM_elongmax_mean,CON_GM_forcemax_mean,CON_GM_elongmax_SD, 'bo')
        % visual
        axis([0 22 0 3600])
        xlabel('Tendon elongation (mm)')
        ylabel('Force (N)')
        title(plottitle)
        legend([h2 h4], 'BD entire AT', 'CON entire AT', 'Location','Southeast')
        saveas(gcf, horzcat('data_plots_stiff/GRP_stiff_fit_mean_entireAT.jpg'))
        
    end
    
      
    
    %% IND + GRP: extract common force-elong, calculate group means
    % output = 
    %   group variables within common force region:
    %   force
    %   elong mean
    %   elong SD
    
    % find lowest force (shortest cell array among all groups and all SOL/GM)
    len = 10000;
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
    
    
    
    %% IND: save individual force-elong arrays (to common force) to one file common to all subjects
    if ispc
        filename_output = strcat('data_output/stiff_GRP_force-elong_', datestr(now, 'yyyy-mm-dd HH-MM'), '.xls');
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
    else
        % csvwrite
    end
    
    
    %% GRP: plot force-elong per subject + mean (GM and SOL separately)
    
    % Free AT (SOL):
    if plot_check && (BD_SOL_count > 0 || CON_SOL_count > 0)
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
    if plot_check && (BD_GM_count > 0 || CON_GM_count > 0)
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
            
end