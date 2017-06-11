%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main file for tendon stiffness
% Marie Moltubakk 17.5.2013
% 
% Run for all trials for which stiffness should be computed at common
% level + force-elongation plot, together.
% E.g. one run for BD+CON SOL trials, another run for BD+CON GM trials.
% 
% Note 1:
% The scripts assume a certain structure for input files from Tracker and
% Noraxon. I.e. number of and order of EMG and other channels. If these
% scripts are to be used for other projects, some of the code must be
% modified. These lines are marked % PROJECTSPECIFIC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % TODO MMM: 
    % CON 5 SOL - stiff coeff negative
    % BD files existing???
    % 50 or 100N force intervals?
    % check predetermined cutoff for each subject (3 trials sent O%J)
    % common stiffness = greatest common force level = 2300 N ?
    % GM vs SOL trials - same force levels for stiffness, different F/E figures?
    
    
    
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
        plot_norm = 0; % show torque before and after initial lowpass filter / ankle rotation fit plots
        plot_conversion = 1;
    else
        plot_norm = 0; % show torque before and after initial lowpass filter / ankle rotation fit plots
        plot_conversion = 0;
    end
    plot_us = 0;
    plot_emg = 0;  % RMS 3 EMG channels per trial


    %% Set constants and globals % PROJECTSPECIFIC

    % Average stiffness across X N
    forceintervals = 50; %TMP
    
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
        BD_count = 0; % # of ballet dancer subjects
        CON_count = 0; % # of controls = intervention study subjects
        BD_no(ceil(linestotal)) = zeros;
        CON_no(ceil(linestotal)) = zeros;
        BD_stiff{ceil(linestotal)} = zeros;
        CON_stiff{ceil(linestotal)} = zeros;
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
                BD_count = BD_count + 1;
                BD_no(BD_count) = str2double(dm_subjectno{line});
            else
                filepath = 'data\';
                subject_id = horzcat('Control ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
                CON_count = CON_count + 1;
                CON_no(CON_count) = str2double(dm_subjectno{line});
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
        [usdata_momentarm, usfreq_momentarm] = read_us_file(strcat(filepath, dm_CPM_calc_US{line}, '.txt'), str2double(dm_CPM_calc_US_frame{line}), 'CPM calcaneus');

        % Read moment arm noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
        % Produce a new noraxon data array
        noraxon_momentarm = read_noraxon_stiffness(strcat(filepath, dm_CPM_calc_NX{line}), usfreq_momentarm, dm_side{line}, 'CPM calcaneus');

        % Read complete, prepared noraon array + prepared us array
        % Produce AT moment arm constant
        at_momentarm = calculate_momentarm(0,0, dm_leg_length{line});


        %% Calculate relationship between calc displacement and ANKLE ROTATION

        % REUSING us file from moment arm
        % REUSING noraxon file from moment arm

        % Read complete, prepared noraxon array + prepared us array
        % Produce ankle rotation constant, mm/degree
        if strcmp(subject_id,'subject 4 R PRE SOL') %very special case for this trial...
            at_rotation_const = -0.14;
        else %normally
            at_rotation_const = calculate_rotation_correction(noraxon_momentarm, usdata_momentarm);
        end


        %% Calculate MVC for plantarflexion

        % Read MVC noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
        % Produce a new noraxon data array
        % sending in a length corresponding to 9 seconds (delete anything after 9 sec)
        noraxon_MVC = read_noraxon_stiffness(strcat(filepath, dm_MVC_PF{line}), freq_default, dm_side{line}, 'MVC dorsi');

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
        [stiff_eq, stiff_gof, force_elong_array, stiff_force_cutoff] = final_stiffness(time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3, time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3, forceintervals, dm_cutforce{line}, trial_force_max, at_momentarm);


        %% calculate stiffness for last 10 and 20% of ind max:
        stiff_ind_80 = calculate_stiffness(stiff_eq, stiff_force_cutoff, 0.8, 1.0, 'ind max'); % last two variables are percent range, from 0.00 to 1.00
        stiff_ind_90 = calculate_stiffness(stiff_eq, stiff_force_cutoff, 0.9, 1.0, 'ind max');


        %% Strain rate
        %LATER
        % strain_rate = calculate_strain_rate(time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3, time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3);
        
        
        %% other calculations from Melina datasheet?
        %LATER
        % AT/GM tendon elongation / strain / young's modulus
        
    
        %% save individual data to common array + force-elong-file

        % add all individual variables to a common array for all subjects    
        all_stiff_output_txt(line,:) = [dm_subjectno(line) dm_timepoint(line) dm_side(line) dm_trial(line)];
        % NaN at locations for stiffness at force levels common to all subjects
        all_stiff_output(line,:) = [coeffvalues(stiff_eq) stiff_gof.rsquare stiff_force_cutoff min(trial_force_max) plantflex_max_torque stiff_ind_80 stiff_ind_90 NaN NaN NaN NaN at_momentarm at_rotation_const NaN NaN];

        % add force-elongation arrays to xls file per subject
        if ispc
            filename_output = strcat('data_output/stiff_forceoutput_', subject_id, '.xls');
            xlswrite(filename_output, {'Elongation (mm)','Force (N)'}, 'Final stiff', 'A1')
            xlswrite(filename_output, force_elong_array, 'Final stiff', 'A2')
        else
            filename_output = strcat('data_output/stiff_forceoutput_', subject_id, '_averaged.csv');
            csvwrite(filename_output, force_elong_array)
        end
        
        % add force-elongation arrays to cell for future averaging 
        % MMM TODO - separate GM and SOL
         if input_project == 1 % BD study
            if trial_subjectno > 100 % BD subject
                BD_stiff{BD_count} = force_elong_array;
            else
                CON_stiff{CON_count} = force_elong_array;
            end
         else
             %LATER
         end
        
    end
    %% LOOP finished %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %% truncate cells
    BD_stiff(BD_count+1:end) = [];
    CON_stiff(CON_count+1:end) = [];
    
    
    %% calculate stiffness at group common force levels

    % select common force (maximal force reached by all subjects)
    all_stiff_col = 5; % 5 for cutoff force level, 6 for max force level
    stiff_common_force = min(all_stiff_output(:,all_stiff_col));
    all_stiff_col = 6; % 5 for cutoff force level, 6 for max force level
    stiff_common_force_max = min(all_stiff_output(:,all_stiff_col));

    for i = 1:size(all_stiff_output,1)
        % reuse last stiff_eq (cfit object)
        stiff_eq.p1 = all_stiff_output(i,1);
        stiff_eq.p2 = all_stiff_output(i,2);

        % calculate stiffness
        stiff_common_80 = calculate_stiffness(stiff_eq, stiff_common_force, 0.8, 1.0, horzcat('FP', all_stiff_output_txt{i,1}, ' ', all_stiff_output_txt{i,4}, ' common cutoff force')); % last two = percent of submitted force - %VAR
        stiff_common_90 = calculate_stiffness(stiff_eq, stiff_common_force, 0.9, 1.0, horzcat('FP', all_stiff_output_txt{i,1}, ' ', all_stiff_output_txt{i,4}, ' common cutoff force')); %VAR
        stiff_common_80_max = calculate_stiffness(stiff_eq, stiff_common_force_max, 0.8, 1.0, horzcat('FP', all_stiff_output_txt{i,1}, ' ', all_stiff_output_txt{i,4}, ' common max force')); % last two = percent of submitted force - %VAR
        stiff_common_90_max = calculate_stiffness(stiff_eq, stiff_common_force_max, 0.9, 1.0, horzcat('FP', all_stiff_output_txt{i,1}, ' ', all_stiff_output_txt{i,4}, ' common max force')); %VAR

        % add to array across subjects
        all_stiff_output(i,10) = stiff_common_80;
        all_stiff_output(i,11) = stiff_common_90;
        all_stiff_output(i,12) = stiff_common_80_max;
        all_stiff_output(i,13) = stiff_common_90_max;
        all_stiff_output(i,14) = stiff_common_force;
        all_stiff_output(i,15) = stiff_common_force_max;
    end
    
    cprintf('blue*',horzcat('Stiffness: Common cutoff force = ', num2str(stiff_common_force), ' N, common max force = ', num2str(round(stiff_common_force_max,0)), ' N.\n'))



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
    
    
    %% calculate group force-elongation
    if BD_count > 0
        % find lowest force (shortest cell array)
        len = 10000;
        for i = 1:BD_count
            tmp_size = size(BD_stiff{1,i});
            if tmp_size(1) < len
                len = tmp_size(1);
            end
        end

        % collect elong up to common force
        BD_force = BD_stiff{1,1}(1:len,2);
        BD_elong(1:len,1:BD_count) = NaN;
        for i = 1:length(BD_stiff)
            BD_elong(:,i) = BD_stiff{1,i}(1:len,1);
        end
        BD_elong_mean = nanmean(BD_elong,2);
        BD_elong_SD = nanstd(BD_elong,0,2);
    end
    
    if CON_count > 0
        len = 10000;
        for i = 1:CON_count
            tmp_size = size(CON_stiff{1,i});
            if tmp_size(1) < len
                len = tmp_size(1);
            end
        end

        % collect elong up to common force
        CON_force = CON_stiff{1,1}(1:len,2);
        CON_elong(1:len,1:CON_count) = NaN;
        for i = 1:length(CON_stiff)
            CON_elong(:,i) = CON_stiff{1,i}(1:len,1);
        end
        CON_elong_mean = nanmean(CON_elong,2);
        CON_elong_SD = nanstd(CON_elong,0,2);
    end
    
    
    %% save force-elong arrays to one file common to all subjects
    if ispc
        filename_output = strcat('data_output/stiff_GRP_force-elong_', datestr(now, 'yyyy-mm-dd HH-MM'), '.xls');

        if BD_count > 0
            filename_subj_BD = cellstr(num2str(BD_no'))';
            xlswrite(filename_output, ['Force (N)' filename_subj_BD], 'Force-elong BD', 'A1')
            xlswrite(filename_output, BD_force, 'Force-elong BD', 'A2')
            xlswrite(filename_output, BD_elong, 'Force-elong BD', 'B2')
        end
        if CON_count > 0
            filename_subj_CON = cellstr(num2str(CON_no'))';
            xlswrite(filename_output, ['Force (N)' filename_subj_CON], 'Force-elong CON', 'A1')
            xlswrite(filename_output, CON_force, 'Force-elong CON', 'A2')
            xlswrite(filename_output, CON_elong, 'Force-elong CON', 'B2')
        end
    else
        
        % csvwrite
    end
    
    
    %% plot group force-elong figure
    if plot_check
        fig_f_e_legend = [];

        plottitle = horzcat('Force-elongation, common angles');
        figure('Name',plottitle)
        hold on
        if BD_count > 0
            h1 = plot(BD_elong_mean, BD_force, 'r', 'LineWidth',2);
            fig_f_e_legend{end+1} = 'BD avg';
        end
        if CON_count > 0
            h2 = plot(CON_elong_mean, CON_force, 'b', 'LineWidth',2);
            fig_f_e_legend{end+1} = 'CON avg';
        end
        if BD_count > 0
            h3 = plot(BD_elong_mean+BD_elong_SD, BD_force, 'r--', 'LineWidth',0.5);
            plot(BD_elong_mean-BD_elong_SD, BD_force, 'r--', 'LineWidth',0.5);
            fig_f_e_legend{end+1} = 'BD SD';
        end
        if CON_count > 0
            h5 = plot(CON_elong_mean+CON_elong_SD, CON_force, 'b--', 'LineWidth',0.5);
            plot(CON_elong_mean-CON_elong_SD, CON_force, 'b--', 'LineWidth',0.5);
            fig_f_e_legend{end+1} = 'CON SD';
        end
        %axis(axis_force)
        xlabel('Tendon elongation (mm)')
        ylabel('Force (N)')
        title(plottitle)
        legend([h1 h2 h3 h5],fig_f_e_legend, 'Location','Northwest')
        saveas(gcf, horzcat('data_plots_stiff/GRP_BD ',plottitle,'.jpg'))
    end
            
end