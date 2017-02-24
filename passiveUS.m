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
%
% 03.08.16: INPUT ARGUMENT: pass 0/1/2 for amount of plots (see below)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [] = passiveUS(input_plot)
    close all
    
    
    
    %%% Determine which plots to produce during script running
    global plot_achilles plot_norm plot_emg plot_check plot_us subject_id plot_licht plot_individual plot_conversion

    if input_plot >= 1 
        plot_check = 1; % toggle main / summarizing plots (LEVEL 1)
    else
        plot_check = 0;
    end
    if input_plot >= 2
        plot_individual = 1; % plot all trials per subject (force-angle etc)
    else
        plot_individual = 0;
    end
    plot_norm = 0; % turn on/off Norm/checkpoint data plots (LEVEL 2)

    plot_conversion = 0; % turn on/off plots for data conversion Norm
    plot_us = 0; % tracked feature vs external marker 
    plot_emg = 0; % RMS 3 EMG channels per trial
    plot_achilles = 0; % turn on/off Achilles machine plots
    plot_licht = 0; % plot averaging of trials from Lichtwark US fascicle tracking





    %%% Set constants and globals % PROJECTSPECIFIC

    % sampling frequencies
    global us_zerodispframes noraxonfreq freq_default
    us_zerodispframes = 1; % No of US frames to average as zero displacement
    noraxonfreq = 1500; % sampling frequency of noraxon data
    freq_default = 100; % output frequency for noraxon data without US video (MVC, etc)

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

    % Average active tendon stiffness across X N
    % forceintervals = 100; 

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
    

    
    
    
    
    
    
    
    %%% Read files from "create_angles_passive.m", to extract max angles and forces per trial/subject/common
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
    




    %%% Read datamaster file, to connect corresponding data files
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
    
    
    
    
    %%% preallocate output arrays
    % common arrays for all subjects:
    all_passive_output = zeros(ceil(linestotal),50); 
    all_passive_output_txt = cell(ceil(linestotal),4);

    % BD-SPECIFIC
    BD_count = 0;
    CON_count = 0;
    
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



    
    
    %%%%%%%%%%%%%%%% LOOP through all lines in datamaster file (except header line)
    for line = 1:linestotal
        clear emg_all;
        emg_all{3,6} = zeros; % 6 EMG channels for 1 subject, to be reused across subjects




        %%% subject/trial identifier
        subjectno = str2double(dm_subjectno{line});
        if subjectno > 100
            filepath = 'data\BD\';
            subject_id = horzcat('dancer ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
            BD_count = BD_count + 1; %BD-SPECIFIC
        else
            filepath = 'data\';
            subject_id = horzcat('control ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
            CON_count = CON_count + 1; %BD-SPECIFIC
        end
        cprintf('*black', horzcat('----------------', subject_id, '------------------\n'))




        %%% Calculate data for muscle activation (EMG)

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



        
        %%% calculate NORM conversion factors
        % retrieve conversion constants for Norm data
        [convert_norm_angle_a, convert_norm_angle_b, convert_norm_torque_a, convert_norm_torque_b, convert_norm_velocity_a, convert_norm_velocity_b, convert_norm_direction_b] = calculate_norm_constants_complete(dm_subjectno{line}, dm_side{line}, dm_timepoint{line}, line, 'passive');




        %%% calculate ACHILLES TENDON MOMENT ARM
        at_momentarm = calculate_momentarm(0, 0, dm_leg_length{line});

        


        

        %%% Calculations for 2x SOL trials
        
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
            axis([-2 35 -3 11]) %VAR
            ylabel('Displacement (mm)')
            xlabel('Gonio angle (°)')
            title(plottitle)
            legend('Trial 1','Trial 2','Location','Southeast')
            saveas(gcf, horzcat('data_plots/', plottitle,'.jpg'))
        end




        %%% Calculations for 2x GM MTJ trials

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
            axis([-2 35 -3 15]) %VAR
            ylabel('Displacement (mm)')
            xlabel('Gonio angle (°)')
            title(plottitle)
            legend('Trial 1','Trial 2','Location','Southeast')
            saveas(gcf, horzcat('data_plots/', plottitle,'.jpg'))
        end


        

        %%% Calculations for 2x GM fascicle trials - ORIGINAL calculations, same as for GMmtj and for SOL

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
            axis([-2 35 -1 2.5]) %VAR
            ylabel('Displacement (mm)')
            xlabel('Gonio angle (°)')
            title(plottitle)
            legend('Trial 1','Trial 2','Location','Southeast')
            saveas(gcf, horzcat('data_plots/', plottitle,'.jpg'))
        end
        
        
        
      
        
        
        
        
        
        
        
        
        %%% Calculations for 2x GM fascicle trials - LICHTWARK analysis
        
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

        
        
        
        
        %%% GM length change ---- MMM TODO as below:
        
        % Fukunaga 2001: In vivo behavor of human ...

        %   data_GMFAS_licht_GM
        % containing:
        %   averaged angle (currently calculated from gonio)
        %   averaged fasicle length
        %   averaged pennation angle
        % OR containing zeros (if nonexistent)
        
        % muscle longitudinal length change = faslen * cos penn_ang, at each joint angle
        % no data about initial length, so only elongation possible, not total muscle length or strain
        % current "muscle" length change can be split into GM and proximal tendon
        
        
        
        
        
        
        

        %%% Average SOL + GMMTJ + GMFAS trials for force, gonio, angle, EMG
        
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
            axis([-2 35 0 1400]) %VAR
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
            
            axis([-2 35 -3 30]) %VAR
            ylabel('Muscle activation (% of MVC)')
            xlabel('Gonio angle (°)')
            title(plottitle)
            legend('Gastr.med.','Gastr.lat.','Soleus','Location','Northwest')
            saveas(gcf, horzcat('data_plots/', plottitle,'.jpg'))
        end
        
        
        
        
        
        %%% check conformation of goniometer to norm angle
        if plot_check && plot_norm  
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
        
        
        
        
        
        
                
        %%% Convert displacements to lengths - plot MTU length across angles
        
        % MTU_length_array =                      MTU_elong_array =
        %1       angle
        %2       calc to SOL ins = free tendon                 
        %3       calc go GM ins = GM tendon
        %4       calc to knee = entire calf
        %5          ~                           GMFAS_displ   = displacement of GMFAS point
        %6       SOL to GM = GM apo
        %7       GM to knee = GM msc
        %7       SOL length (from Hawkins & Hull) minus free AT = SOL msc
        
        % send max gonio angle
        if strcmpi(dm_timepoint{line}, 'PRE') == 1
            if strcmpi(dm_side{line}, 'R') == 1
                out_ROM_trial_max = str2double(input_gon_pre_r{subjectno}) - 0.00000001; % tweak for computer storage of numbers, where 8.3000 is stored as 8.2999999999999999
            else % L
                out_ROM_trial_max = str2double(input_gon_pre_l{subjectno}) - 0.00000001;
            end
        else % POST
            if strcmpi(dm_side{line}, 'R') == 1
                out_ROM_trial_max = str2double(input_gon_post_r{subjectno}) - 0.00000001;
            else % L
                out_ROM_trial_max = str2double(input_gon_post_l{subjectno}) - 0.00000001;
            end
        end
        
        % calculate MTU lengths
        [MTU_length_array, MTU_elong_array, MTU_strain_array] = calculate_mtu_length(data_SOL(:,2:3), data_GMMTJ(:,2:3), data_GMFAS(:,2:3), dm_at_SOL_length{line}, dm_at_GM_length{line}, dm_leg_length{line}, out_ROM_trial_max);

        if plot_check && plot_individual
            % raw lengths (mm)
            plottitle = horzcat('IND MTU length vs angle, ', subject_id);
            figure('Name',plottitle);
            hold on
            plot(MTU_length_array(:,1),MTU_length_array(:,4),'LineWidth',1.5,'Color','blue') % AT + GM apo + GM msc
            plot(MTU_length_array(:,1),MTU_length_array(:,8)+MTU_length_array(:,2),'LineWidth',1.5,'Color','black') % AT + SOL msc
            plot(MTU_length_array(:,1),MTU_length_array(:,3),'LineWidth',1.5,'Color','red') % AT + GM apo
            plot(MTU_length_array(:,1),MTU_length_array(:,2),'LineWidth',1.5,'Color','green') % AT
            % plot vertical lines to show color coding of subsequent figs
            plot([1,1], [MTU_length_array(1,4), MTU_length_array(1,3)],'Color','cyan')
            plot([1,1], [MTU_length_array(1,3), MTU_length_array(1,2)],'Color',[1 0.75 0])
            plot([1,1], [MTU_length_array(1,2), 0],'Color','green')
            axis([-2 35 0 550]) %VAR
            ylabel('Length (mm)')
            xlabel('Gonio angle (°)')
            title(plottitle)
            legend('GM MTU', 'SOL MTU', 'AT up to GM insert.', 'Free AT', 'Location','East');
            saveas(gcf, horzcat('data_plots/', plottitle,'.jpg'))
            
            % raw elongation (mm)
            plottitle = horzcat('IND MTU elongation vs angle, ', subject_id);
            figure('Name',plottitle);
            hold on
            plot(MTU_elong_array(:,1),MTU_elong_array(:,4),'LineWidth',1.5,'Color','blue') % MTU = AT + GM apo + msc
            plot(MTU_elong_array(:,1),MTU_elong_array(:,7),'LineWidth',1.5,'Color','cyan') % GM msc
            plot(MTU_elong_array(:,1),MTU_elong_array(:,8),'LineWidth',1.5,'Color','black') % SOL msc
            plot(MTU_elong_array(:,1),MTU_elong_array(:,6),'LineWidth',1.5,'Color',[1 0.75 0]) % GM apo
            plot(MTU_elong_array(:,1),MTU_elong_array(:,2),'LineWidth',1.5,'Color','green') % AT
            plot(MTU_elong_array(:,1),MTU_elong_array(:,3),'LineWidth',1.5,'Color','red') % GM tendon = AT + GM apo
            plot(MTU_elong_array(:,1),MTU_elong_array(:,5),'LineWidth',1.5,'Color','yellow') % GM FAS displacement
            axis([-2 35 -3 30]) %VAR
            ylabel('Elongation (mm)')
            xlabel('Gonio angle (°)')
            title(plottitle)
            legend('GM MTU', 'GM msc', 'SOL msc', 'GM aponeur.', 'Free AT', 'GM tendon', 'GM fascicle DISPL.', 'Location','Northwest');
            saveas(gcf, horzcat('data_plots/', plottitle,'.jpg'))
            
            % strain (percent of initial length)
            plottitle = horzcat('IND MTU strain vs angle, ', subject_id);
            figure('Name',plottitle);
            hold on
            plot(MTU_strain_array(:,1),MTU_strain_array(:,4),'LineWidth',1.5,'Color','blue') % MTU = AT + GM apo + msc
            plot(MTU_strain_array(:,1),MTU_strain_array(:,7),'LineWidth',1.5,'Color','cyan') % GM msc
            plot(MTU_strain_array(:,1),MTU_strain_array(:,8),'LineWidth',1.5,'Color','black') % SOL msc
            plot(MTU_strain_array(:,1),MTU_strain_array(:,6),'LineWidth',1.5,'Color',[1 0.75 0]) % GM apo
            plot(MTU_strain_array(:,1),MTU_strain_array(:,2),'LineWidth',1.5,'Color','green') % AT
            plot(MTU_strain_array(:,1),MTU_strain_array(:,3),'LineWidth',1.5,'Color','red') % GM tendon = AT + GM apo
            axis([-2 35 -1 24]) %VAR
            ylabel('Strain (% of initial length)')
            xlabel('Gonio angle (°)')
            title(plottitle)
            legend('GM MTU', 'GM msc', 'SOL msc', 'GM aponeur.', 'Free AT', 'GM tendon', 'Location','Northwest');
            saveas(gcf, horzcat('data_plots/', plottitle,'.jpg'))
        end
        
        
        
        
        
        
        
        
        
        %%% Extract FORCE, ANGLE, DISPLACEMENT, EMG at various sets of max joint angles
        %   force, angle, EMG from all 3 scan locations / 6 trials, averaged
        %   displacement of SOL, GMMTJ, GMFAS from 2 trials per scan location
        %   elongation (abs and normalized) of AT, GM apo, msc portion, from displacement data
        %   joint angles = 
        %       out_ROM_trial_max = trial max (different PRE, POST, L, R)
        %       out_ROM_ind_max = subject ind max (lowest PRE, POST, L, R)
        %       out_ROM_common_max = common max (lowest of all subjects)
        %       out_ROM_submax_1 = additional predetermined angle, e.g. 1/3 and 2/3 of trial max ROM
        %       out_ROM_submax_2 = additional predetermined angle, e.g. 1/3 and 2/3 of trial max ROM
                
        % column choices
        col_force = 1;
        col_angle = 2;
        col_displ = 3;
        col_lichtfas = 2;
        col_lichtpenn = 3;
        
        % print error if max angles do not exist --> create_angles_passive needs to be run
        if str2double(input_gon_ind_max{subjectno}) == 100
            cprintf('*red', 'ERROR: Max ROM values are not calculated for subject. Run create_angles_passive first.\n')
        end

        % goniometer angles 
        if strcmpi(dm_timepoint{line}, 'PRE') == 1
            if strcmpi(dm_side{line}, 'R') == 1
                out_ROM_trial_max = str2double(input_gon_pre_r{subjectno}) - 0.00000001; % tweak for computer storage of numbers, where 8.3000 is stored as 8.2999999999999999
            else % L
                out_ROM_trial_max = str2double(input_gon_pre_l{subjectno}) - 0.00000001;
            end
        else % POST
            if strcmpi(dm_side{line}, 'R') == 1
                out_ROM_trial_max = str2double(input_gon_post_r{subjectno}) - 0.00000001;
            else % L
                out_ROM_trial_max = str2double(input_gon_post_l{subjectno}) - 0.00000001;
            end
        end
        out_ROM_ind_max = str2double(input_gon_ind_max{subjectno}) - 0.00000001;
        out_ROM_common_max = str2double(input_gon_common_max{subjectno}) - 0.00000001;
        out_ROM_submax_1 = out_ROM_trial_max * 1/3; %VAR
        out_ROM_submax_2 =  out_ROM_trial_max * 2/3; %VAR
        
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
        out_F_submax_1 = data_force_gonio(loc_frame,col_force);
        loc_frame = find(data_force_gonio(:,col_angle)>=out_ROM_submax_2,1,'first'); 
        out_F_submax_2 = data_force_gonio(loc_frame,col_force);
        
        % torques
        out_T_trial_max_ROM = out_F_trial_max_ROM *at_momentarm;
        out_T_trial_max_F = out_F_trial_max_F *at_momentarm;
        out_T_ind_max = out_F_ind_max *at_momentarm;
        out_T_common_max = out_F_common_max *at_momentarm;
        out_T_zero = out_F_zero *at_momentarm;
        out_T_submax_1 = out_F_submax_1 *at_momentarm;
        out_T_submax_2 = out_F_submax_2 *at_momentarm;

        % displacements SOL
        loc_frame = find(data_SOL(:,col_angle)>=out_ROM_trial_max,1,'first'); 
        out_displ_SOL_trial_max = data_SOL(loc_frame,col_displ);
        loc_frame = find(data_SOL(:,col_angle)>=out_ROM_ind_max,1,'first'); 
        out_displ_SOL_ind_max = data_SOL(loc_frame,col_displ); 
        loc_frame = find(data_SOL(:,col_angle)>=out_ROM_common_max,1,'first'); 
        out_displ_SOL_common_max = data_SOL(loc_frame,col_displ); 
        loc_frame = find(data_SOL(:,col_angle)>=out_ROM_submax_1,1,'first'); 
        out_displ_SOL_submax_1 = data_SOL(loc_frame,col_displ); 
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
        out_displ_GMMTJ_submax_1 = data_GMMTJ(loc_frame,col_displ); 
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
        out_displ_GMFAS_submax_1 = data_GMFAS(loc_frame,col_displ); 
        loc_frame = find(data_GMFAS(:,col_angle)>=out_ROM_submax_2,1,'first'); 
        out_displ_GMFAS_submax_2 = data_GMFAS(loc_frame,col_displ); 

        % column choices ELONGATION / LENGTH / STRAIN
        col_angle_elong = 1;
        col_AT = 2;
        col_GM_tend = 3;
        col_leg = 4;
        col_GM_apo = 6;
        col_msc_GM = 7;
        col_msc_SOL = 8;
        
        % absolute elong AT - should be same as displ SOL
        % absolute elong GM apo
        % absolute elong msc part (prox of GM insertion)
        % relative elong SOL
        % relative elong GM apo
        % relative elong msc part
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
        loc_frame = find(MTU_elong_array(:,col_angle_elong)>=out_ROM_submax_1,1,'first'); 
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
        
        
        

        
        
        
        
        
        
        % EMG (from all 3 scan locations / 6 trials, averaged)
        emg_step = 9; %VAR - number of EMG values BEFORE relevant angle, to include in average. 9 values = a span of 0.5 degrees.
        % cannot average values AROUND relevant angle, because for a few (least flexible) trials, we want the data AT the last available angle.

        % identify locations of the relevant goniometer angles, in the array with all 2+2+2 trials averaged:
        loc_frame_trial_max = find(data_force_gonio(:,col_angle)>=out_ROM_trial_max,1,'first'); % when averaging angles across trials, intervals of 0.05 degrees are used. This means that any required angle will exist in all data series
        loc_frame_ind_max = find(data_force_gonio(:,col_angle)>=out_ROM_ind_max,1,'first'); 
        loc_frame_common_max = find(data_force_gonio(:,col_angle)>=out_ROM_common_max,1,'first'); 
        loc_frame_submax_1 = find(data_force_gonio(:,col_angle)>=out_ROM_submax_1,1,'first'); 
        loc_frame_submax_2 = find(data_force_gonio(:,col_angle)>=out_ROM_submax_2,1,'first'); 

        % break if frame not found
        if isempty(loc_frame_ind_max)
            % cprintf('error', horzcat('ERROR: Computed trial max ROM (', num2str(out_ROM_trial_max), ') does not exist in the current data series (max = ', num2str(max(data_force_gonio(:,col_angle))), '). Check "create_angles_passive" vs current data.\n' ));
            error(horzcat('ERROR: Computed trial max ROM (', num2str(out_ROM_trial_max), ') does not exist in the current data series (max = ', num2str(max(data_force_gonio(:,col_angle))), '). Check "create_angles_passive" vs current data.\n' ))
        end
        
        % EMG gm
        out_emg_gm_trial_max = mean(data_force_gonio(loc_frame_trial_max-emg_step:loc_frame_trial_max,3)); % EMG gm = column 3
        out_emg_gm_ind_max = mean(data_force_gonio(loc_frame_ind_max-emg_step:loc_frame_ind_max,3));
        out_emg_gm_common_max = mean(data_force_gonio(loc_frame_common_max-emg_step:loc_frame_common_max,3));
        out_emg_gm_submax_1 = mean(data_force_gonio(loc_frame_submax_1-emg_step:loc_frame_submax_1,3));
        out_emg_gm_submax_2 = mean(data_force_gonio(loc_frame_submax_2-emg_step:loc_frame_submax_2,3));

        % EMG gl
        out_emg_gl_trial_max = mean(data_force_gonio(loc_frame_trial_max-emg_step:loc_frame_trial_max,4)); % EMG gl = column 4
        out_emg_gl_ind_max = mean(data_force_gonio(loc_frame_ind_max-emg_step:loc_frame_ind_max,4)); 
        out_emg_gl_common_max = mean(data_force_gonio(loc_frame_common_max-emg_step:loc_frame_common_max,4));
        out_emg_gl_submax_1 = mean(data_force_gonio(loc_frame_submax_1-emg_step:loc_frame_submax_1,4));
        out_emg_gl_submax_2 = mean(data_force_gonio(loc_frame_submax_2-emg_step:loc_frame_submax_2,4));

        % EMG sol
        out_emg_sol_trial_max = mean(data_force_gonio(loc_frame_trial_max-emg_step:loc_frame_trial_max,5)); % EMG sol = column 5
        out_emg_sol_ind_max = mean(data_force_gonio(loc_frame_ind_max-emg_step:loc_frame_ind_max,5));
        out_emg_sol_common_max = mean(data_force_gonio(loc_frame_common_max-emg_step:loc_frame_common_max,5));
        out_emg_sol_submax_1 = mean(data_force_gonio(loc_frame_submax_1-emg_step:loc_frame_submax_1,5));
        out_emg_sol_submax_2 = mean(data_force_gonio(loc_frame_submax_2-emg_step:loc_frame_submax_2,5));


        
        

        %%% extract FASCICLE LENGTH, PENNATION ANGLE at 3 sets of max joint angles and at zero angle
        %   data_GMFAS_licht_GM
        %   data_GMFAS_licht_SOL
        % containing:
        %   averaged angle (currently calculated from gonio)
        %   averaged fasicle length
        %   averaged pennation angle
        % OR containing zeros (if nonexistent)

        if data_GMFAS_licht_GM == 0
            % no licht data existing
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
        else
            % at trial max angle:
            loc_frame = find(data_GMFAS(:,col_angle)>=out_ROM_trial_max,1,'first'); 
            out_licht_faslen_GM_trial_max = data_GMFAS_licht_GM(loc_frame,col_lichtfas);
            out_licht_pennation_GM_trial_max = data_GMFAS_licht_GM(loc_frame,col_lichtpenn);
            if data_GMFAS_licht_SOL == 0 % SOL data don't exist
                out_licht_faslen_SOL_trial_max = 100;
                out_licht_pennation_SOL_trial_max = 100;
            else
                % SOL exists
                out_licht_faslen_SOL_trial_max = data_GMFAS_licht_SOL(loc_frame,col_lichtfas); 
                out_licht_pennation_SOL_trial_max = data_GMFAS_licht_SOL(loc_frame,col_lichtpenn); 
            end

            % at individual max angle (across sides/timepoints):
            loc_frame = find(data_GMFAS(:,col_angle)>=out_ROM_ind_max,1,'first'); 
            out_licht_faslen_GM_ind_max = data_GMFAS_licht_GM(loc_frame,col_lichtfas); 
            out_licht_pennation_GM_ind_max = data_GMFAS_licht_GM(loc_frame,col_lichtpenn); 
            if data_GMFAS_licht_SOL == 0 % SOL data don't exist
                out_licht_faslen_SOL_ind_max = 100;
                out_licht_pennation_SOL_ind_max = 100;
            else
                % SOL exists
                out_licht_faslen_SOL_ind_max = data_GMFAS_licht_SOL(loc_frame,col_lichtfas); 
                out_licht_pennation_SOL_ind_max = data_GMFAS_licht_SOL(loc_frame,col_lichtpenn); 
            end

            % at subject common max angle:
            loc_frame = find(data_GMFAS(:,col_angle)>=out_ROM_common_max,1,'first'); 
            out_licht_faslen_GM_common_max = data_GMFAS_licht_GM(loc_frame,col_lichtfas); 
            out_licht_pennation_GM_common_max = data_GMFAS_licht_GM(loc_frame,col_lichtpenn); 
            if data_GMFAS_licht_SOL == 0 % SOL data don't exist
                out_licht_faslen_SOL_common_max = 100;
                out_licht_pennation_SOL_common_max = 100;
            else
                % SOL exists
                out_licht_faslen_SOL_common_max = data_GMFAS_licht_SOL(loc_frame,col_lichtfas); 
                out_licht_pennation_SOL_common_max = data_GMFAS_licht_SOL(loc_frame,col_lichtpenn); 
            end

            % at zero angle:
            loc_frame = find(data_GMFAS(:,2)>=0,1,'first'); 
            out_licht_faslen_GM_zero = data_GMFAS_licht_GM(loc_frame,col_lichtfas); 
            out_licht_pennation_GM_zero = data_GMFAS_licht_GM(loc_frame,col_lichtpenn); 
            if data_GMFAS_licht_SOL == 0 % SOL data don't exist
                out_licht_faslen_SOL_zero = 100;
                out_licht_pennation_SOL_zero = 100;
            else
                % SOL exists
                out_licht_faslen_SOL_zero = data_GMFAS_licht_SOL(loc_frame,col_lichtfas); 
                out_licht_pennation_SOL_zero = data_GMFAS_licht_SOL(loc_frame,col_lichtpenn); 
            end

            % at submax_1 angle:
            loc_frame = find(data_GMFAS(:,col_angle)>=out_ROM_submax_1,1,'first'); 
            out_licht_faslen_GM_submax_1 = data_GMFAS_licht_GM(loc_frame,col_lichtfas); 
            out_licht_pennation_GM_submax_1 = data_GMFAS_licht_GM(loc_frame,col_lichtpenn); 
            if data_GMFAS_licht_SOL == 0 % SOL data don't exist
                out_licht_faslen_SOL_submax_1 = 100;
                out_licht_pennation_SOL_submax_1 = 100;
            else
                % SOL exists
                out_licht_faslen_SOL_submax_1 = data_GMFAS_licht_SOL(loc_frame,col_lichtfas); 
                out_licht_pennation_SOL_submax_1 = data_GMFAS_licht_SOL(loc_frame,col_lichtpenn); 
            end

            % at submax_2 angle:
            loc_frame = find(data_GMFAS(:,col_angle)>=out_ROM_submax_2,1,'first'); 
            out_licht_faslen_GM_submax_2 = data_GMFAS_licht_GM(loc_frame,col_lichtfas); 
            out_licht_pennation_GM_submax_2 = data_GMFAS_licht_GM(loc_frame,col_lichtpenn); 
            if data_GMFAS_licht_SOL == 0 % SOL data don't exist
                out_licht_faslen_SOL_submax_2 = 100;
                out_licht_pennation_SOL_submax_2 = 100;
            else
                % SOL exists
                out_licht_faslen_SOL_submax_2 = data_GMFAS_licht_SOL(loc_frame,col_lichtfas); 
                out_licht_pennation_SOL_submax_2 = data_GMFAS_licht_SOL(loc_frame,col_lichtpenn); 
            end
        end







        %%%%%%%%%%%% PASSIVE STIFFNESS and STIFFNESS INDEX (Nordez 2006) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % gonio angle = data_force_gonio(:,col_angle)
        % force = data_force_gonio(:,col_force)
        % multiplying force with at_momentarm to convert to torque
        
        
        
        %%% PASSIVE STIFFNESS:
        
        % passive stiffness = delta torque / delta angle, at various angles
        % fit 4th order polynomial to averaged torque-angle curve, using data from start to:
        
        %    - subject's individual max ROM (across both legs)
        fit_ind_max = polyfit(data_force_gonio(1:loc_frame_ind_max,col_angle), at_momentarm*data_force_gonio(1:loc_frame_ind_max,col_force), 4);
        %   currently using torque-angle data up to the INDIVIDUAL max, to calculate the fit. Then calculating stiffness at common max angle, and at 15 degrees.
        %    - trial max ROM
        % fit_trial_max = polyfit(data_force_gonio(1:loc_angle_trial_max,col_angle), at_momentarm*data_force_gonio(1:loc_angle_trial_max,col_force), 4);
        %    - common max ROM
        % fit_common_max = polyfit(data_force_gonio(1:loc_angle_common_max,col_angle), at_momentarm*data_force_gonio(1:loc_angle_common_max,col_force), 4);

        % extract passive stiffness (derivate of 4th order poly) at:
        out_pstiff_trial_max = (4 * fit_ind_max(1) * out_ROM_trial_max^3) + (3 * fit_ind_max(2) * out_ROM_trial_max^2) + (2 * fit_ind_max(3) * out_ROM_trial_max) + fit_ind_max(4);
        out_pstiff_ind_max = (4 * fit_ind_max(1) * out_ROM_ind_max^3) + (3 * fit_ind_max(2) * out_ROM_ind_max^2) + (2 * fit_ind_max(3) * out_ROM_ind_max) + fit_ind_max(4);
        out_pstiff_common_max = (4 * fit_ind_max(1) * out_ROM_common_max^3) + (3 * fit_ind_max(2) * out_ROM_common_max^2) + (2 * fit_ind_max(3) * out_ROM_common_max) + fit_ind_max(4);
        out_pstiff_submax_1 = (4 * fit_ind_max(1) * out_ROM_submax_1^3) + (3 * fit_ind_max(2) * out_ROM_submax_1^2) + (2 * fit_ind_max(3) * out_ROM_submax_1) + fit_ind_max(4);
        out_pstiff_submax_2 = (4 * fit_ind_max(1) * out_ROM_submax_2^3) + (3 * fit_ind_max(2) * out_ROM_submax_2^2) + (2 * fit_ind_max(3) * out_ROM_submax_2) + fit_ind_max(4);
        out_pstiff_angle = 15; %VAR
        out_pstiff_15 = (4 * fit_ind_max(1) * out_pstiff_angle^3) + (3 * fit_ind_max(2) * out_pstiff_angle^2) + (2 * fit_ind_max(3) * out_pstiff_angle) + fit_ind_max(4);

        
        
        %%% STIFFNESS INDEX:
        
        % fit 2nd order polynomial to averaged torque-angle curve, using data from start to:
        
        %    - ind max ROM
        fit_ind_max = polyfit(data_force_gonio(1:loc_frame_ind_max,col_angle), data_force_gonio(1:loc_frame_ind_max,col_force), 2);
        %   currently using torque-angle data up to the INDIVIDUAL max, to calculate the fit
        %    - trial max ROM
        % fit_trial_max = polyfit(data_force_gonio(1:loc_angle_trial_max,col_angle), data_force_gonio(1:loc_angle_trial_max,col_force), 2);
        %    - common max ROM
        % fit_common_max = polyfit(data_force_gonio(1:loc_angle_common_max,col_angle), data_force_gonio(1:loc_angle_common_max,col_force), 2);

        % extract stiffness index as 2 * a
        out_pstiff_index = 2 * fit_ind_max(1);
        
        
        
        %%% AT STIFFNESS:     % MMM TODO GOON
        
        % stiffness = delta force / delta length for AT ? check literature
        
        
        
        
        
        
        %%% Extract ANGLES at specific FORCE levels
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
        if str2double(input_for_ind_max{subjectno}) == 10000
            cprintf('*red', 'ERROR: Max ROM values are not calculated for current subject. Run create_angles_passive first.\n')
        end

        % load the predetermined forces to be searched for
        if strcmpi(dm_timepoint{line}, 'PRE') == 1
            if strcmpi(dm_side{line}, 'R') == 1
                loc_F_trial_max = str2double(input_for_pre_r{subjectno}) - 0.0001; % tweak for computer storage of numbers, where 8.3000 is stored as 8.2999999999999999
            else % L
                loc_F_trial_max = str2double(input_for_pre_l{subjectno}) - 0.0001;
            end
        else % POST
            if strcmpi(dm_side{line}, 'R') == 1
                loc_F_trial_max = str2double(input_for_post_r{subjectno}) - 0.0001;
            else % L
                loc_F_trial_max = str2double(input_for_post_l{subjectno}) - 0.0001;
            end
        end
        loc_F_ind_max = str2double(input_for_ind_max{subjectno}) - 0.0001;
        loc_F_common_max = str2double(input_for_common_max{subjectno}) - 0.0001;
        loc_F_ind_rmax = str2double(input_for_ind_rmax{subjectno}) - 0.0001;
        loc_F_ind_lmax = str2double(input_for_ind_lmax{subjectno}) - 0.0001;

        % find goniometer angles (from all 3 scan locations / 6 trials, averaged)
        loc_frame = find(data_force_gonio(:,col_force)>=loc_F_trial_max,1,'first'); % when averaging angles across trials, intervals of 0.05 degrees are used. This means that any required angle will exist in all data series
        if isempty(loc_frame) % catch error of angle not found:
            loc_frame = find(data_force_gonio(:,col_force)>=max(data_force_gonio(:,col_force)),1,'first'); 
            cprintf('*red', 'ERROR: Recorded ind max angle not found in data. Run create_angles_passive?\n')
        end
        out_angle_trial_max = data_force_gonio(loc_frame,col_angle); 
        loc_frame = find(data_force_gonio(:,col_force)>=loc_F_ind_max,1,'first'); 
        if isempty(loc_frame) % catch error of angle not found:
            loc_frame = find(data_force_gonio(:,col_force)>=max(data_force_gonio(:,col_force)),1,'first'); 
            cprintf('*red', 'ERROR: Recorded ind max angle not found in data. Run create_angles_passive?\n')
        end
        out_angle_ind_max = data_force_gonio(loc_frame,col_angle);
        loc_frame = find(data_force_gonio(:,col_force)>=loc_F_common_max,1,'first');
        out_angle_common_max = data_force_gonio(loc_frame,col_angle);
        if str2double(input_for_ind_rmax(subjectno)) > 9000
            % data do not exist for the right side (array preloaded with "empty" values of 10000)
            out_angle_ind_rmax = 100;
        else
            loc_frame = find(data_force_gonio(:,col_force)>=loc_F_ind_rmax,1,'first'); 
                    if isempty(loc_frame) % catch error of angle not found:
                        loc_frame = find(data_force_gonio(:,col_force)>=max(data_force_gonio(:,col_force)),1,'first'); 
                        cprintf('*red', 'ERROR: Recorded ind max angle not found in data. Run create_angles_passive?\n')
                    end
            out_angle_ind_rmax = data_force_gonio(loc_frame,col_angle);
        end
        if str2double(input_for_ind_lmax(subjectno)) > 9000
            % data do not exist for the left side (array preloaded with "empty" values of 10000)
            out_angle_ind_lmax = 100;
        else
            loc_frame = find(data_force_gonio(:,col_force)>=loc_F_ind_lmax,1,'first'); 
                    if isempty(loc_frame) % catch error of angle not found:
                        loc_frame = find(data_force_gonio(:,col_force)>=max(data_force_gonio(:,col_force)),1,'first'); 
                        cprintf('*red', 'ERROR: Recorded ind max angle not found in data. Run create_angles_passive?\n')
                    end
            out_angle_ind_lmax = data_force_gonio(loc_frame,col_angle);
        end
        
        
        
        
        
        
        
        
        
        
        %%% OUTPUT data for plot across subjects

        % data_force_gonio
        % from M-file average_passive_forces_EMG
        % average of 6 trials
        % output array "data_force_gonio" contains: 
        %   average_force_gonio
        %   average_angle_array
        %   average_emg_gm_gonio
        %   average_emg_gl_gonio
        %	average_emg_sol_gonio
        
        % data_SOL data_GMMTJ data_GMFAS
        % from M-file average_passive_trials_EMG
        % average of 2 trials
        % output array "data_SOL" contains:
        %   average_force_gonio
        %   average_angle_array
        %   average_displ_gonio
        %   average_emg_gm_gonio
        %   average_emg_gl_gonio
        %	average_emg_sol_gonio
        
        % extract angle range common to all trials for current subject
        angle_start = 0 - 0.00000001; % tweak for computer storage of numbers, where 8.3000 is stored as 8.2999999999999999 %VAR
        angle_stop = out_ROM_trial_max;

        % identify locations of start/stop angles in above mentioned arrays
        loc_angle_start = find(data_force_gonio(:,col_angle)>=angle_start,1,'first');
        loc_angle_stop = find(data_force_gonio(:,col_angle)>=angle_stop,1,'first');
%         loc_angle_start_SOL = find(data_SOL(:,col_angle)>=angle_start,1,'first');
%         loc_angle_stop_SOL = find(data_SOL(:,col_angle)>=angle_stop,1,'first');
%         loc_angle_start_GMMTJ = find(data_GMMTJ(:,col_angle)>=angle_start,1,'first');
%         loc_angle_stop_GMMTJ = find(data_GMMTJ(:,col_angle)>=angle_stop,1,'first');
%         loc_angle_start_GMFAS = find(data_GMFAS(:,col_angle)>=angle_start,1,'first');
%         loc_angle_stop_GMFAS = find(data_GMFAS(:,col_angle)>=angle_stop,1,'first');
            
        if subjectno > 100 % BD subject
            
            % angle_vars contains: 
            %   1 angle
            %   2 F 
            %   3 EMG_gm 
            %   4 EMG_gl 
            %   5 EMG_sol 
            
            %   6 elong AT / displ SOL     
            %   7 elong GMtend / displ GMMTJ
            % * 8 elong leg        
            % * 9 ---- / displ GMFAS
            % *10 elong GMapo       
            % *11 elong msc         
            
            % *12 L_AT              
            % *13 L_GMtend          
            % *14 L_leg             
            % *15 ----- (GMFAS)
            % *16 L_GMapo           
            % *17 L_msc             
            
            % *18 Torque            
            
            % *19 elong msc SOL
            % *20 L msc SOL
            
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
                ];
            
            % all data in ONE cell, NORMALIZED data:
            %    - EMG only normalized to %MVC, not to max EMG @ trial
            %    - elongations (6-11) not normalized to anything
            %    NB: in #12-17, normalizing to initial length in array, instead of using constants: str2double(dm_at_SOL_length{line}), dm_at_GM_length, (cm to mm) 10*dm_leg_length. Tested and values are the same except rounding at 2nd or 3rd decimal.
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
                ];

            % reshape
            BD_angle_vars_norm{BD_count} = [ (0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,2), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,3), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,4), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,5), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,6), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,7), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,8), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,9), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,10), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,11), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,12), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,13), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,14), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,15), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,16), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,17), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,18), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,19), 0:0.05:100)', ...
                spline(BD_angle_vars_norm_indlength{1,BD_count}(:,1), BD_angle_vars_norm_indlength{1,BD_count}(:,20), 0:0.05:100)'];
            
            
            
            
            
                        
            
        else % CON subject
            
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
                ];

            % reshape
            CON_angle_vars_norm{CON_count} = [ (0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,2), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,3), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,4), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,5), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,6), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,7), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,8), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,9), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,10), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,11), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,12), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,13), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,14), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,15), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,16), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,17), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,18), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,19), 0:0.05:100)', ...
                spline(CON_angle_vars_norm_indlength{1,CON_count}(:,1), CON_angle_vars_norm_indlength{1,CON_count}(:,20), 0:0.05:100)'];
            
        end
        
        
        
        
        
        
        
        
        %%% OUTPUT final individual data to file

        % add data to a common array for all subjects    
        i = 1;
        
        % txt trial ID
        all_passive_output_txt(line,1) = dm_subjectno(line);
        all_passive_output_txt(line,2) = dm_timepoint(line);
        all_passive_output_txt(line,3) = dm_side(line);
        all_passive_output_txt(line,4) = dm_trial(line);
        
        % ROM
        all_passive_output(line,i) = out_ROM_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_ROM_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_ROM_common_max;
        i = i+1;
        all_passive_output(line,i) = out_ROM_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_ROM_submax_2;
        i = i+1;
        
        % Force
        all_passive_output(line,i) = out_F_trial_max_ROM;
        i = i+1;
        all_passive_output(line,i) = out_F_trial_max_F;
        i = i+1;
        all_passive_output(line,i) = out_F_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_F_common_max;
        i = i+1;
        all_passive_output(line,i) = out_F_zero; % special
        i = i+1;
        all_passive_output(line,i) = out_F_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_F_submax_2;
        i = i+1;
        
        % Force
        all_passive_output(line,i) = out_T_trial_max_ROM;
        i = i+1;
        all_passive_output(line,i) = out_T_trial_max_F;
        i = i+1;
        all_passive_output(line,i) = out_T_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_T_common_max;
        i = i+1;
        all_passive_output(line,i) = out_T_zero; % special
        i = i+1;
        all_passive_output(line,i) = out_T_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_T_submax_2;
        i = i+1;
        
        % angle @ F
        all_passive_output(line,i) = out_angle_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_angle_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_angle_common_max;
        i = i+1;
        all_passive_output(line,i) = out_angle_ind_rmax; % special
        i = i+1;
        all_passive_output(line,i) = out_angle_ind_lmax; % special
        i = i+1;
        
        % stiff
        all_passive_output(line,i) = out_pstiff_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_pstiff_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_pstiff_common_max;
        i = i+1;
        all_passive_output(line,i) = out_pstiff_15; % special
        i = i+1;
        all_passive_output(line,i) = out_pstiff_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_pstiff_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_pstiff_index; % special
        i = i+1;
        
        % US displ -- should be same as elongation
        all_passive_output(line,i) = out_displ_SOL_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_displ_SOL_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_displ_SOL_common_max;
        i = i+1;
        all_passive_output(line,i) = out_displ_SOL_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_displ_SOL_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_displ_GMMTJ_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_displ_GMMTJ_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_displ_GMMTJ_common_max;
        i = i+1;
        all_passive_output(line,i) = out_displ_GMMTJ_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_displ_GMMTJ_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_displ_GMFAS_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_displ_GMFAS_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_displ_GMFAS_common_max;
        i = i+1;
        all_passive_output(line,i) = out_displ_GMFAS_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_displ_GMFAS_submax_2;
        i = i+1;
        
        % elongation AT, GM tend, leg, GM apo, msc part
        all_passive_output(line,i) = out_elong_AT_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_GMapo_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_GMtend_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_leg_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_msc_GM_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_msc_SOL_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_AT_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_GMtend_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_leg_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_GMapo_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_msc_GM_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_msc_SOL_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_AT_common_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_GMtend_common_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_leg_common_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_GMapo_common_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_msc_GM_common_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_msc_SOL_common_max;
        i = i+1;
        all_passive_output(line,i) = out_elong_AT_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_elong_GMtend_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_elong_leg_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_elong_GMapo_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_elong_msc_GM_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_elong_msc_SOL_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_elong_AT_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_elong_GMtend_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_elong_leg_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_elong_GMapo_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_elong_msc_GM_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_elong_msc_SOL_submax_2;
        i = i+1;
        
        % strain
        all_passive_output(line,i) = out_strain_AT_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_GMtend_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_leg_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_GMapo_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_msc_GM_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_msc_SOL_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_AT_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_GMtend_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_leg_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_GMapo_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_msc_GM_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_msc_SOL_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_AT_common_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_GMtend_common_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_leg_common_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_GMapo_common_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_msc_GM_common_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_msc_SOL_common_max;
        i = i+1;
        all_passive_output(line,i) = out_strain_AT_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_strain_GMtend_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_strain_leg_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_strain_GMapo_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_strain_msc_GM_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_strain_msc_SOL_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_strain_AT_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_strain_GMtend_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_strain_leg_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_strain_GMapo_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_strain_msc_GM_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_strain_msc_SOL_submax_2;
        i = i+1;
        
        
        
        % length
        all_passive_output(line,i) = out_length_AT_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_length_GMtend_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_length_leg_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_length_GMapo_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_length_msc_GM_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_length_msc_SOL_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_length_AT_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_length_GMtend_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_length_leg_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_length_GMapo_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_length_msc_GM_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_length_msc_SOL_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_length_AT_common_max;
        i = i+1;
        all_passive_output(line,i) = out_length_GMtend_common_max;
        i = i+1;
        all_passive_output(line,i) = out_length_leg_common_max;
        i = i+1;
        all_passive_output(line,i) = out_length_GMapo_common_max;
        i = i+1;
        all_passive_output(line,i) = out_length_msc_GM_common_max;
        i = i+1;
        all_passive_output(line,i) = out_length_msc_SOL_common_max;
        i = i+1;
        all_passive_output(line,i) = out_length_AT_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_length_GMtend_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_length_leg_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_length_GMapo_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_length_msc_GM_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_length_msc_SOL_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_length_AT_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_length_GMtend_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_length_leg_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_length_GMapo_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_length_msc_GM_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_length_msc_SOL_submax_2;
        i = i+1;
        
        % US lichtwark GM
        all_passive_output(line,i) = out_licht_faslen_GM_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_licht_faslen_GM_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_licht_faslen_GM_common_max;
        i = i+1;
        all_passive_output(line,i) = out_licht_faslen_GM_zero; % special
        i = i+1;
        all_passive_output(line,i) = out_licht_faslen_GM_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_licht_faslen_GM_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_licht_pennation_GM_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_licht_pennation_GM_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_licht_pennation_GM_common_max;
        i = i+1;
        all_passive_output(line,i) = out_licht_pennation_GM_zero; % special
        i = i+1;
        all_passive_output(line,i) = out_licht_pennation_GM_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_licht_pennation_GM_submax_2;
        i = i+1;
        % US lichtwark SOL
        all_passive_output(line,i) = out_licht_faslen_SOL_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_licht_faslen_SOL_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_licht_faslen_SOL_common_max;
        i = i+1;
        all_passive_output(line,i) = out_licht_faslen_SOL_zero; % special
        i = i+1;
        all_passive_output(line,i) = out_licht_faslen_SOL_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_licht_faslen_SOL_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_licht_pennation_SOL_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_licht_pennation_SOL_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_licht_pennation_SOL_common_max;
        i = i+1;
        all_passive_output(line,i) = out_licht_pennation_SOL_zero; % special
        i = i+1;
        all_passive_output(line,i) = out_licht_pennation_SOL_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_licht_pennation_SOL_submax_2;
        i = i+1;

        % EMG
        all_passive_output(line,i) = out_emg_gm_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_emg_gm_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_emg_gm_common_max;
        i = i+1;
        all_passive_output(line,i) = out_emg_gm_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_emg_gm_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_emg_gl_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_emg_gl_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_emg_gl_common_max;
        i = i+1;
        all_passive_output(line,i) = out_emg_gl_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_emg_gl_submax_2;
        i = i+1;
        all_passive_output(line,i) = out_emg_sol_trial_max;
        i = i+1;
        all_passive_output(line,i) = out_emg_sol_ind_max;
        i = i+1;
        all_passive_output(line,i) = out_emg_sol_common_max;
        i = i+1;
        all_passive_output(line,i) = out_emg_sol_submax_1;
        i = i+1;
        all_passive_output(line,i) = out_emg_sol_submax_2;
        

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%% LOOP FINISHED
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATIONS ACROSS SUBJECTS - NUMBERS

    %%%  mean and stdav of each subject's INDIVIDUAL MAX ROM, force, elong, EMG
    n_o_array_elements = 20; %VAR - number of elements in BD_angle_xxx arrays. OLD version = 13
    
    if BD_count > 0
        % preallocate array
        BD_max(BD_count,n_o_array_elements) = zeros;
        BD_max_norm(BD_count,n_o_array_elements) = zeros;
        % collect variables per subject
        for i = 1:BD_count % per subject
            for j = 1:n_o_array_elements % per element in arrays
                BD_max(i,j) = max(BD_angle_vars{1,i}(:,j));
                BD_max_norm(i,j) = max(BD_angle_vars_norm{1,i}(:,j));
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
        BD_displ_SOL_mean = mean(BD_max(:,6)); % / elong AT
        BD_displ_SOL_SD = std(BD_max(:,6));
        BD_displ_GMMTJ_mean = mean(BD_max(:,7)); % / elong GM tend
        BD_displ_GMMTJ_SD = std(BD_max(:,7));
        BD_displ_MTU_mean = mean(BD_max(:,8)); % / elong leg
        BD_displ_MTU_SD = std(BD_max(:,8));
        BD_displ_GMFAS_mean = mean(BD_max(:,9)); 
        BD_displ_GMFAS_SD = std(BD_max(:,9));
        BD_elong_GMapo_mean = mean(BD_max(:,10)); 
        BD_elong_GMapo_SD = std(BD_max(:,10));
        BD_elong_msc_mean = mean(BD_max(:,11)); % GM msc
        BD_elong_msc_SD = std(BD_max(:,11));
        BD_L_at_SOL_mean = mean(BD_max(:,12)); % L AT
        BD_L_at_SOL_SD = std(BD_max(:,12));
        BD_L_at_GM_mean = mean(BD_max(:,13)); % L GM tend
        BD_L_at_GM_SD = std(BD_max(:,13));
        BD_L_MTU_mean = mean(BD_max(:,14)); % L leg
        BD_L_MTU_SD = std(BD_max(:,14));
        BD_L_GMapo_mean = mean(BD_max(:,16)); 
        BD_L_GMapo_SD = std(BD_max(:,16));
        BD_L_msc_mean = mean(BD_max(:,17)); 
        BD_L_msc_SD = std(BD_max(:,17));
        BD_torque_mean = mean(BD_max(:,18));
        BD_torque_SD = std(BD_max(:,18));
        BD_elong_msc_SOL_mean = mean(BD_max(:,19)); 
        BD_elong_msc_SOL_SD = std(BD_max(:,19));
        BD_L_msc_SOL_mean = mean(BD_max(:,20)); 
        BD_L_msc_SOL_SD = std(BD_max(:,20));
        % determine common angle range
        BD_common_ROM = min(BD_max(:,1));
    end
    
    if CON_count > 0
        % preallocate array
        CON_max(CON_count,n_o_array_elements) = zeros;
        CON_max_norm(CON_count,n_o_array_elements) = zeros;
        % collect variables per subject
        for i = 1:CON_count % per subject
            for j = 1:n_o_array_elements % per element in arrays
                CON_max(i,j) = max(CON_angle_vars{1,i}(:,j));
                CON_max_norm(i,j) = max(CON_angle_vars_norm{1,i}(:,j));
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
        CON_displ_SOL_mean = mean(CON_max(:,6)); % / elong AT
        CON_displ_SOL_SD = std(CON_max(:,6));
        CON_displ_GMMTJ_mean = mean(CON_max(:,7)); % / elong GM tend
        CON_displ_GMMTJ_SD = std(CON_max(:,7));
        CON_displ_MTU_mean = mean(CON_max(:,8)); % / elong leg
        CON_displ_MTU_SD = std(CON_max(:,8));
        CON_displ_GMFAS_mean = mean(CON_max(:,9)); 
        CON_displ_GMFAS_SD = std(CON_max(:,9));
        CON_elong_GMapo_mean = mean(CON_max(:,10)); 
        CON_elong_GMapo_SD = std(CON_max(:,10));
        CON_elong_msc_mean = mean(CON_max(:,11)); % GM msc
        CON_elong_msc_SD = std(CON_max(:,11));
        CON_L_at_SOL_mean = mean(CON_max(:,12)); % L AT
        CON_L_at_SOL_SD = std(CON_max(:,12));
        CON_L_at_GM_mean = mean(CON_max(:,13)); % L GM tend
        CON_L_at_GM_SD = std(CON_max(:,13));
        CON_L_MTU_mean = mean(CON_max(:,14)); % L leg
        CON_L_MTU_SD = std(CON_max(:,14));
        CON_L_GMapo_mean = mean(CON_max(:,16)); 
        CON_L_GMapo_SD = std(CON_max(:,16));
        CON_L_msc_mean = mean(CON_max(:,17)); 
        CON_L_msc_SD = std(CON_max(:,17));
        CON_torque_mean = mean(CON_max(:,18));
        CON_torque_SD = std(CON_max(:,18));
        CON_elong_msc_SOL_mean = mean(CON_max(:,19)); 
        CON_elong_msc_SOL_SD = std(CON_max(:,19));
        CON_L_msc_SOL_mean = mean(CON_max(:,20)); 
        CON_L_msc_SOL_SD = std(CON_max(:,20));
        % determine common angle range
        CON_common_ROM = min(CON_max(:,1));
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%% OUTPUT KEY VARIABLES FOR ALL SUBJECTS TO FILE
    % write xls
    if ispc
        filename_output = strcat('data_output/all_passive_output_', datestr(now, 'yyyy-mm-dd HH-MM'));
        
        all_passive_output_head = {'Subject', 'Time', 'Side', 'Trial', ...
            'ROM trial (°)', 'ROM subject (L/R/PRE/POST)', 'ROM common (all subjects)', '33% ROM', '67% ROM', ...
            'force @ trial ROM (N)', 'force @ subject max ROM', 'force @ subject max force', 'force @ common max ROM', 'force @ 0 deg', 'force @ 33% ROM', 'force @ 67% ROM', ...
            'torque @ trial ROM (N)', 'torque @ subject max ROM', 'torque @ subject max torque', 'torque @ common max ROM', 'torque @ 0 deg', 'torque @ 33% ROM', 'torque @ 67% ROM', ...
            'angle @ trial max force (°)', 'angle @ subject max force', 'angle @ common max force', 'angle @ ind R max force', 'angle @ ind L max force', ...            
            'Stiffness @ trial ROM (Nm/°)', 'Stiffness @ subject max ROM', 'Stiffness @ common max ROM', 'Stiffness @ 15 deg', 'Stiffness @ 33% ROM', 'Stiffness @ 67% ROM',...
            'Stiffness index (Nm/°^2)', ...
            'displ SOL @ trial ROM (mm)', 'displ SOL @ subject max ROM', 'displ SOL @ common max ROM', 'displ SOL @ 33% ROM', 'displ SOL @ 67% ROM', ...
            'displ GMMTJ @ trial ROM', 'displ GMMTJ @ subject max ROM', 'displ GMMTJ @ common max ROM', 'displ GMMTJ @ 33% ROM', 'displ GMMTJ @ 67% ROM', ...
            'displ GMfas @ trial ROM', 'displ GMfas @ subject max ROM', 'displ GMfas @ common max ROM', 'displ GMFAS @ 33% ROM', 'displ GMFAS @ 67% ROM', ...
            'elong_AT_trial_max', 'elong_GMtend_trial_max', 'elong_msc_trial_max', 'elong_GMapo_trial_max', 'elong_msc_GM_trial_max', 'elong_msc_SOL_trial_max',...
            'elong_AT_ind_max', 'elong_GMtend_ind_max', 'elong_leg_ind_max', 'elong_GMapo_ind_max', 'elong_msc_GM_ind_max', 'elong_msc_SOL_ind_max', ...
            'elong_AT_common_max', 'elong_GMtend_common_max', 'elong_leg_common_max', 'elong_GMapo_common_max', 'elong_msc_GM_common_max', 'elong_msc_SOL_common_max', ...
            'elong_AT_submax_1', 'elong_GMtend_submax_1', 'elong_leg_submax_1', 'elong_GMapo_submax_1', 'elong_msc_GM_submax_1', 'elong_msc_SOL_submax_1', ...
            'elong_AT_submax_2', 'elong_GMtend_submax_2', 'elong_leg_submax_2', 'elong_GMapo_submax_2', 'elong_msc_GM_submax_2', 'elong_msc_SOL_submax_2', ...
            'strain_AT_trial_max', 'strain_GMtend_trial_max', 'strain_leg_trial_max', 'strain_GMapo_trial_max', 'strain_msc_GM_trial_max', 'strain_msc_SOL_trial_max', ...
            'strain_AT_ind_max', 'strain_GMtend_ind_max', 'strain_leg_ind_max', 'strain_GMapo_ind_max', 'strain_msc_GM_ind_max', 'strain_msc_SOL_ind_max', ...
            'strain_AT_common_max', 'strain_GMtend_common_max', 'strain_leg_common_max', 'strain_GMapo_common_max', 'strain_msc_GM_common_max', 'strain_msc_SOL_common_max', ...
            'strain_AT_submax_1', 'strain_GMtend_submax_1', 'strain_leg_submax_1', 'strain_GMapo_submax_1', 'strain_msc_GM_submax_1', 'strain_msc_SOL_submax_1', ...
            'strain_AT_submax_2', 'strain_GMtend_submax_2', 'strain_leg_submax_2', 'strain_GMapo_submax_2', 'strain_msc_GM_submax_2', 'strain_msc_SOL_submax_2', ...
            'length_AT_trial_max', 'length_GMtend_trial_max', 'length_leg_trial_max', 'length_GMapo_trial_max', 'length_msc_GM_trial_max', 'length_msc_SOL_trial_max', ...
            'length_AT_ind_max', 'length_GMtend_ind_max', 'length_leg_ind_max', 'length_GMapo_ind_max', 'length_msc_GM_ind_max', 'length_msc_SOL_ind_max', ...
            'length_AT_common_max', 'length_GMtend_common_max', 'length_leg_common_max', 'length_GMapo_common_max', 'length_msc_GM_common_max', 'length_msc_SOL_common_max', ...
            'length_AT_submax_1', 'length_GMtend_submax_1', 'length_leg_submax_1', 'length_GMapo_submax_1', 'length_msc_GM_submax_1', 'length_msc_SOL_submax_1', ...
            'length_AT_submax_2', 'length_GMtend_submax_2', 'length_leg_submax_2', 'length_GMapo_submax_2', 'length_msc_GM_submax_2', 'length_msc_SOL_submax_2', ...
            'faslen GM @ trial ROM (mm)', 'faslen GM @ subject max ROM', 'faslen GM @ common max ROM', 'faslen GM @ 0 deg', 'faslen GM @ 33% ROM', 'faslen GM @ 67% ROM', ...
            'pennation GM @ trial ROM (°)', 'pennation GM @ subject max ROM', 'pennation GM @ common max ROM', 'pennation GM @ 0 deg', 'pennation GM @ 33% ROM', 'pennation GM @ 67% ROM', ...
            'faslen SOL @ trial ROM', 'faslen SOL @ subject max ROM', 'faslen SOL @ common max ROM', 'faslen SOL @ 0 deg', 'faslen SOL @ 33% ROM', 'faslen SOL @ 67% ROM', ...
            'pennation SOL @ trial ROM', 'pennation SOL @ subject max ROM', 'pennation SOL @ common max ROM', 'pennation SOL @ 0 deg', 'pennation SOL @ 33% ROM', 'pennation SOL @ 67% ROM', ...
            'EMG GM @ trial ROM (%)', 'EMG GM @ subject max ROM', 'EMG GM @ common max ROM', 'EMG GM @ 33% ROM', 'EMG GM @ 67% ROM', ...
            'EMG GL @ trial ROM', 'EMG GL @ subject max ROM', 'EMG GL @ common max ROM', 'EMG GL @ 33% ROM', 'EMG GL @ 67% ROM', ...
            'EMG SOL @ trial ROM', 'EMG SOL @ subject max ROM', 'EMG SOL @ common max ROM', 'EMG SOL @ 33% ROM', 'EMG SOL @ 67% ROM', ...
            }; % PROJECTSPECIFIC
        
        xlswrite(filename_output, all_passive_output_head, 1, 'A1')
        xlswrite(filename_output, all_passive_output_txt, 1, 'A2')
        xlswrite(filename_output, all_passive_output, 1, 'E2')
    else
        filename_output = strcat('data_output/all_passive_output_', datestr(now, 'yyyy-mm-dd HH-MM'), '.csv');
        csvwrite(filename_output, all_passive_output)
    end

    
    
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATIONS ACROSS SUBJECTS - ARRAYS
    
    
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%% PLOT GROUP FIGURES
    % 
    % group mean absolute               common ROM + mean + SD
    % group mean normalized             common ROM
    % group BD absolute separately      ind ROM
    % group BD normalized separately    ind ROM
    % group CON absolute separately     ind ROM
    % group CON normalized separately   ind ROM
    
    
    
    % force-angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP force vs angle - 1 ABSOLUTE');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,2),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,2),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_F_mean, BD_F_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_F_mean, CON_F_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_F_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_F_mean, CON_ROM_SD, '*b')
            axis([-2 35 0 1400]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Force (N)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP force vs angle - 2 ABSOLUTE IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,2))
            end
            axis([-2 35 0 1400]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Force (N)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP force vs angle - 3 ABSOLUTE IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,2))
            end
            axis([-2 35 0 1400]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Force (N)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && CON_count > 1 && plot_check            
            plottitle = horzcat('GRP force vs angle - 4 NORMALIZED');
            figure('Name',plottitle)
            plot(BD_angle_vars_norm_mean(:,1), BD_angle_vars_norm_mean(:,2),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_norm_mean(:,1), CON_angle_vars_norm_mean(:,2),'b','LineWidth',2)
            axis([-1 100 0 110]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Force (% of ind max)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP force vs angle - 5 NORMALIZED IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars_norm{1,i}(:,1),BD_angle_vars_norm{1,i}(:,2))
            end
            axis([-1 100 0 110]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Force (% of ind max)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP force vs angle - 6 NORMALIZED IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars_norm{1,i}(:,1),CON_angle_vars_norm{1,i}(:,2))
            end
            axis([-1 100 0 110]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Force (% of ind max)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    
    
    
    % torque-angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP torque vs angle - 1 ABSOLUTE');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,18),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,18),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_torque_mean, BD_torque_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_torque_mean, CON_torque_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_torque_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_torque_mean, CON_ROM_SD, '*b')
            axis([-2 35 0 80]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Torque (Nm)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP torque vs angle - 2 ABSOLUTE IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,18))
            end
            axis([-2 35 0 80]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Torque (Nm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP torque vs angle - 3 ABSOLUTE IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,18))
            end
            axis([-2 35 0 80]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Torque (Nm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP torque vs angle - 4 NORMALIZED');
            figure('Name',plottitle)
            plot(BD_angle_vars_norm_mean(:,1), BD_angle_vars_norm_mean(:,18),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_norm_mean(:,1), CON_angle_vars_norm_mean(:,18),'b','LineWidth',2)
            axis([-1 100 0 110]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Torque (% of ind max)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP torque vs angle - 5 NORMALIZED IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars_norm{1,i}(:,1),BD_angle_vars_norm{1,i}(:,18))
            end
            axis([-1 100 0 110]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Torque (% of ind max)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP torque vs angle - 6 NORMALIZED IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars_norm{1,i}(:,1),CON_angle_vars_norm{1,i}(:,18))
            end
            axis([-1 100 0 110]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Torque (% of ind max)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    
    
    
    
    
    
    % SOL SCANS - FREE AT: LENGTH vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP free AT length vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,12),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,12),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_L_at_SOL_mean, BD_L_at_SOL_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_L_at_SOL_mean, CON_L_at_SOL_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_L_at_SOL_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_L_at_SOL_mean, CON_ROM_SD, '*b')
            axis([-2 35 0 120]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Southeast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP free AT length vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,12))
            end
            axis([-2 35 0 120]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP free AT length vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,12))
            end
              axis([-2 35 0 120]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP free AT strain vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_norm_mean(:,1), BD_angle_vars_norm_mean(:,12),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_norm_mean(:,1), CON_angle_vars_norm_mean(:,12),'b','LineWidth',2)
            axis([-10 100 -1 24]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP free AT strain vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars_norm{1,i}(:,1),BD_angle_vars_norm{1,i}(:,12))
            end
            axis([-10 100 -1 24]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP free AT strain vs angle - 2B IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count 
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars_norm_indlength{1,i}(:,12))
            end
            axis([-2 35 -1 24]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Strain (% of initial length)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP free AT strain vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars_norm{1,i}(:,1),CON_angle_vars_norm{1,i}(:,12))
            end
            axis([-10 100 -1 24]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP free AT strain vs angle - 3B IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars_norm_indlength{1,i}(:,12))
            end
            axis([-2 35 -1 24]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Strain (% of initial length)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    
    
    
    
    
    
    
    % SOL SCANS: Displacement vs angle = elongation AT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP free AT elongation vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,6),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,6),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_displ_SOL_mean, BD_displ_SOL_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_displ_SOL_mean, CON_displ_SOL_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_displ_SOL_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_displ_SOL_mean, CON_ROM_SD, '*b')
            axis([-2 35 -1 10]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP free AT elongation vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,6))
            end
            axis([-2 35 -1 10]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP free AT elongation vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,6))
            end
            axis([-2 35 -1 10]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end

    
    
    
    
    

    
    
    
    % GM TENDON: LENGTH vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP GM tendon (from calc) length vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,13),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,13),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_L_at_GM_mean, BD_L_at_GM_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_L_at_GM_mean, CON_L_at_GM_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_L_at_GM_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_L_at_GM_mean, CON_ROM_SD, '*b')
            axis([-2 35 0 300]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Southeast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP GM tendon (from calc) length vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,13))
            end
            axis([-2 35 0 300]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP GM tendon (from calc) length vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,13))
            end
            axis([-2 35 0 300]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP GM tendon (from calc) strain vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_norm_mean(:,1), BD_angle_vars_norm_mean(:,13),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_norm_mean(:,1), CON_angle_vars_norm_mean(:,13),'b','LineWidth',2)
            axis([-10 100 -1 7]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP GM tendon (from calc) strain vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars_norm{1,i}(:,1),BD_angle_vars_norm{1,i}(:,13))
            end
            axis([-10 100 -1 7]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP GM tendon (from calc) strain vs angle - 2B IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars_norm_indlength{1,i}(:,13))
            end
            axis([-2 35 -1 7]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Strain (% of initial length)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP GM tendon (from calc) strain vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars_norm{1,i}(:,1),CON_angle_vars_norm{1,i}(:,13))
            end
            axis([-10 100 -1 7]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP GM tendon (from calc) strain vs angle - 3B IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars_norm_indlength{1,i}(:,13))
            end
            axis([-2 35 -1 7]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Strain (% of initial length)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    
    
    
    
    
    % GM TENDON: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP GM tendon (from calc) elongation vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,7),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,7),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_displ_GMMTJ_mean, BD_displ_GMMTJ_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_displ_GMMTJ_mean, CON_displ_GMMTJ_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_displ_GMMTJ_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_displ_GMMTJ_mean, CON_ROM_SD, '*b')
            axis([-2 35 -1 13]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP GM tendon (from calc) elongation vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,7))
            end
            axis([-2 35 -1 13]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP GM tendon (from calc) elongation vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,7))
            end
            axis([-2 35 -1 13]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end

    
    
    
    
    

    
    % GM APO: LENGTH vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP GM apo (SOL ins-GM ins) length vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,16),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,16),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_L_GMapo_mean, BD_L_GMapo_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_L_GMapo_mean, CON_L_GMapo_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_L_GMapo_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_L_GMapo_mean, CON_ROM_SD, '*b')
            axis([-2 35 0 190]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Southeast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP GM apo (SOL ins-GM ins) length vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,16))
            end
            axis([-2 35 0 190]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP GM apo (SOL ins-GM ins) length vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,16))
            end
            axis([-2 35 0 190]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP GM apo (SOL ins-GM ins) strain vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_norm_mean(:,1), BD_angle_vars_norm_mean(:,16),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_norm_mean(:,1), CON_angle_vars_norm_mean(:,16),'b','LineWidth',2)
            axis([-10 100 -2 7]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)') 
            title(plottitle)
            legend('Dancer avg', 'Control avg','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP GM apo (SOL ins-GM ins) strain vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars_norm{1,i}(:,1),BD_angle_vars_norm{1,i}(:,16))
            end
            axis([-10 100 -2 7]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)') 
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP GM apo (SOL ins-GM ins) strain vs angle - 2B IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars_norm_indlength{1,i}(:,16))
            end
            axis([-2 35 -2 7]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Strain (% of initial length)') 
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP GM apo (SOL ins-GM ins) strain vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars_norm{1,i}(:,1),CON_angle_vars_norm{1,i}(:,16))
            end
            axis([-10 100 -2 7]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP GM apo (SOL ins-GM ins) strain vs angle - 3B IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars_norm_indlength{1,i}(:,16))
            end
            axis([-2 35 -2 7]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Strain (% of initial length)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    
    
    
    
    
    
    % GM APO: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP GM apo (SOL ins-GM ins) elongation vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,10),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,10),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_elong_GMapo_mean, BD_elong_GMapo_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_elong_GMapo_mean, CON_elong_GMapo_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_elong_GMapo_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_elong_GMapo_mean, CON_ROM_SD, '*b')
            axis([-2 35 -3 9]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP GM apo (SOL ins-GM ins) elongation vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,10))
            end
            axis([-2 35 -3 9]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP GM apo (SOL ins-GM ins) elongation vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,10))
            end
            axis([-2 35 -3 9]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end

    
    
    
    
    
    
    
    

    % GM MUSCLE portion: LENGTH vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP muscle GM (ins-knee) length vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,17),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,17),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_L_msc_mean, BD_L_msc_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_L_msc_mean, CON_L_msc_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_L_msc_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_L_msc_mean, CON_ROM_SD, '*b')
            axis([-2 35 0 350]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Southeast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP muscle GM (ins-knee) length vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,17))
            end
            axis([-2 35 0 350]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP muscle GM (ins-knee) length vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,17))
            end
            axis([-2 35 0 350]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP muscle GM (ins-knee) strain vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_norm_mean(:,1), BD_angle_vars_norm_mean(:,17),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_norm_mean(:,1), CON_angle_vars_norm_mean(:,17),'b','LineWidth',2)
            axis([-10 100 -1 10]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)') 
            title(plottitle)
            legend('Dancer avg', 'Control avg','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP muscle GM (ins-knee) strain vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars_norm{1,i}(:,1),BD_angle_vars_norm{1,i}(:,17))
            end
            axis([-10 100 -1 10]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)') 
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP muscle GM (ins-knee) strain vs angle - 2B IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars_norm_indlength{1,i}(:,17))
            end
            axis([-2 35 -1 10]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Strain (% of initial length)') 
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP muscle GM (ins-knee) strain vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars_norm{1,i}(:,1),CON_angle_vars_norm{1,i}(:,17))
            end
            axis([-10 100 -1 10]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)') 
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP muscle GM (ins-knee) strain vs angle - 3B IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars_norm_indlength{1,i}(:,17))
            end
            axis([-2 35 -1 10]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Strain (% of initial length)') 
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    
    
    
    
    
    
    % GM MUSCLE portion: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP muscle GM (ins-knee) elongation vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,11),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,11),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_elong_msc_mean, BD_elong_msc_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_elong_msc_mean, CON_elong_msc_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_elong_msc_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_elong_msc_mean, CON_ROM_SD, '*b')
            axis([-2 35 -2 24]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP muscle GM (ins-knee) elongation vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,11))
            end
            axis([-2 35 -2 24]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP muscle GM (ins-knee) elongation vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,11))
            end
            axis([-2 35 -2 24]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end

    
    
   
    
    
    
    
    
    

    % SOL MUSCLE portion: LENGTH vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP muscle SOL length vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,20),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,20),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_L_msc_SOL_mean, BD_L_msc_SOL_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_L_msc_SOL_mean, CON_L_msc_SOL_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_L_msc_SOL_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_L_msc_SOL_mean, CON_ROM_SD, '*b')
            axis([-2 35 0 350]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Southeast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP muscle SOL length vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,20))
            end
            axis([-2 35 0 350]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP muscle SOL length vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,20))
            end
            axis([-2 35 0 350]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP muscle SOL strain vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_norm_mean(:,1), BD_angle_vars_norm_mean(:,20),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_norm_mean(:,1), CON_angle_vars_norm_mean(:,20),'b','LineWidth',2)
            axis([-10 100 -1 11]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)') 
            title(plottitle)
            legend('Dancer avg', 'Control avg','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP muscle SOL strain vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars_norm{1,i}(:,1),BD_angle_vars_norm{1,i}(:,20))
            end
            axis([-10 100 -1 11]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)') 
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP muscle SOL strain vs angle - 2B IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars_norm_indlength{1,i}(:,20))
            end
            axis([-2 35 -1 11]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Strain (% of initial length)') 
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP muscle SOL strain vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars_norm{1,i}(:,1),CON_angle_vars_norm{1,i}(:,20))
            end
            axis([-10 100 -1 11]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)') 
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP muscle SOL strain vs angle - 3B IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars_norm_indlength{1,i}(:,20))
            end
            axis([-2 35 -1 11]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Strain (% of initial length)') 
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    
    
    
    
    
    
    % SOL MUSCLE portion: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP muscle SOL elongation vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,19),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,19),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_elong_msc_SOL_mean, BD_elong_msc_SOL_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_elong_msc_SOL_mean, CON_elong_msc_SOL_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_elong_msc_SOL_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_elong_msc_SOL_mean, CON_ROM_SD, '*b')
            axis([-2 35 -2 28]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP muscle SOL elongation vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,19))
            end
            axis([-2 35 -2 28]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP muscle SOL elongation vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,19))
            end
            axis([-2 35 -2 28]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end

    
    
   
    
    
    
    
    
    
    
    
    
    % GMFAS scans: Displacement vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP GMFAS displacement vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,9),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,9),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_displ_GMFAS_mean, BD_displ_GMFAS_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_displ_GMFAS_mean, CON_displ_GMFAS_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_displ_GMFAS_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_displ_GMFAS_mean, CON_ROM_SD, '*b')
            axis([-2 35 -1 2.5]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Displacement (mm)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP GMFAS displacement vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,9))
            end
            axis([-2 35 -1 2.5]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Displacement (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP GMFAS displacement vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,9))
            end
            axis([-2 35 -1 2.5]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Displacement (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    
    
    
    
    
    
    

    % Total MTU: LENGTH vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP MTU length vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,14),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,14),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_L_MTU_mean, BD_L_MTU_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_L_MTU_mean, CON_L_MTU_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_L_MTU_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_L_MTU_mean, CON_ROM_SD, '*b')
            axis([-2 35 0 550]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Southeast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP MTU length vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,14))
            end
            axis([-2 35 0 550]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP MTU length vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,14))
            end
            axis([-2 35 0 550]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP MTU strain vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_norm_mean(:,1), BD_angle_vars_norm_mean(:,14),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_norm_mean(:,1), CON_angle_vars_norm_mean(:,14),'b','LineWidth',2)
            axis([-10 100 0 7]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP MTU strain vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars_norm{1,i}(:,1),BD_angle_vars_norm{1,i}(:,14))
            end
            axis([-10 100 0 7]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP MTU strain vs angle - 2B IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars_norm_indlength{1,i}(:,14))
            end
            axis([-2 35 0 7]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Strain (% of initial length)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP MTU strain vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars_norm{1,i}(:,1),CON_angle_vars_norm{1,i}(:,14))
            end
            axis([-10 100 0 7]) %VAR
            xlabel('Gonio angle (% of ind max)')
            ylabel('Strain (% of initial length)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP MTU strain vs angle - 3B IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars_norm_indlength{1,i}(:,14))
            end
            axis([-2 35 0 7]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Strain (% of initial length)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    
    
    
    
    
    
    
    % Total MTU: Elongation vs angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP MTU elongation vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,8),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,8),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_displ_MTU_mean, BD_displ_MTU_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_displ_MTU_mean, CON_displ_MTU_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_displ_MTU_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_displ_MTU_mean, CON_ROM_SD, '*b')
            axis([-2 35 -1 30]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP MTU elongation vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,8))
            end
            axis([-2 35 -1 30]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP MTU elongation vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,8))
            end
            axis([-2 35 -1 30]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Elongation (mm)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    

    
    
    
    
    
    
    
    
    % lengths 3 MTU components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP length ALL MTU components vs angle');
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
            
            axis([-2 22 -0 500]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('Length (mm)')
            title(plottitle)
            legend('Dancer MTU', 'Dancer GM tend','Dancer free AT','Control MTU','Control GM tend','Control free AT','Location','Southeast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    
    
   
    
    
    
        
    
        
    % EMG vs angle GM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP EMG gas.med. vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,3),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,3),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_EMG_gm_mean, BD_EMG_gm_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_EMG_gm_mean, CON_EMG_gm_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_EMG_gm_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_EMG_gm_mean, CON_ROM_SD, '*b')
            axis([-2 35 0 20]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('EMG (% of MVC)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP EMG gas.med. vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,3))
            end
            axis([-2 35 0 20]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('EMG (% of MVC)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP EMG gas.med. vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,3))
            end
            axis([-2 35 0 20]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('EMG (% of MVC)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    
    
    
    
    % EMG vs angle GL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP EMG gas.lat. vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,4),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,4),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_EMG_gl_mean, BD_EMG_gl_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_EMG_gl_mean, CON_EMG_gl_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_EMG_gl_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_EMG_gl_mean, CON_ROM_SD, '*b')
            axis([-2 35 0 10]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('EMG (% of MVC)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP EMG gas.lat. vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,4))
            end
            axis([-2 35 0 10]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('EMG (% of MVC)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP EMG gas.lat. vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,4))
            end
            axis([-2 35 0 10]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('EMG (% of MVC)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    
    
    
    
    % EMG vs angle SOL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP EMG soleus vs angle - 1');
            figure('Name',plottitle)
            plot(BD_angle_vars_mean(:,1), BD_angle_vars_mean(:,5),'r','LineWidth',2)
            hold on
            plot(CON_angle_vars_mean(:,1), CON_angle_vars_mean(:,5),'b','LineWidth',2)
            errorbar(BD_ROM_mean, BD_EMG_sol_mean, BD_EMG_sol_SD, '*r', 'MarkerFaceColor', 'r')
            errorbar(CON_ROM_mean, CON_EMG_sol_mean, CON_EMG_sol_SD, '*b', 'MarkerFaceColor', 'b')
            herrorbar(BD_ROM_mean, BD_EMG_sol_mean, BD_ROM_SD, '*r')
            herrorbar(CON_ROM_mean, CON_EMG_sol_mean, CON_ROM_SD, '*b')
            axis([-2 35 0 30]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('EMG (% of MVC)')
            title(plottitle)
            legend('Dancer avg', 'Control avg','Dancer ind max','Control ind max','Location','Northwest')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && plot_check
            plottitle = horzcat('GRP EMG soleus vs angle - 2 IND, dancers');
            figure('Name',plottitle)
            hold on
            for i = 1:BD_count
                plot(BD_angle_vars{1,i}(:,1),BD_angle_vars{1,i}(:,5))
            end
            axis([-2 35 0 30]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('EMG (% of MVC)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if CON_count > 1 && plot_check
            plottitle = horzcat('GRP EMG soleus vs angle - 3 IND, controls');
            figure('Name',plottitle)
            hold on
            for i = 1:CON_count
                plot(CON_angle_vars{1,i}(:,1),CON_angle_vars{1,i}(:,5))
            end
            axis([-2 35 0 30]) %VAR
            xlabel('Gonio angle (°)')
            ylabel('EMG (% of MVC)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    
    
    
    
    
    
    
    
    


end