%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main file for analysis of passive dorsiflexion with US
% Marie Moltubakk 4.2.2015
% 
% input argument 1 = project selection (1 = BD, 2 = intervent)
% input argument 2 = plot selection (0 = none, 1 = group plots, 2 = ind plots)
% input argument 3 = but set in beginning of function: toggle_normalization (0 = absolute, 1 = normalized)
%
% The scripts assume a certain structure for input files from Tracker and
% Noraxon. I.e. number of and order of EMG and other channels. If these
% scripts are to be used for other projects, some of the code must be
% modified. These lines are marked % PROJECTSPECIFIC
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MMM TODO!!! change tendon length calculations to be based on prone length, not zero angle length?
% --- length, pennation and angle angle in input data = prone position
% --- for LICHT:
%            % version 2 - length PRONE is zero elong, zero strain
%            data_GMFAS_licht_GM(:,4) = data_GMFAS_licht_GM(:,2) - prone_GMfas_length; %fascicle elong 
%            data_GMFAS_licht_GM(:,5) = (data_GMFAS_licht_GM(:,2) - prone_GMfas_length) / prone_GMfas_length * 100; %fascicle strain



function [] = passiveUS_rotcorr(input_project, input_plot)
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
    
    % LEVEL 2: main checkpoint plots
    plot_norm = 0;
    
    % LEVEL 3:
    plot_conversion = 0; % turn on/off plots for data conversion Norm
    plot_us = 0; % tracked feature vs external marker 
    plot_emg = 0; % RMS 3 EMG channels per trial
    plot_achilles = 0; % turn on/off Achilles machine plots
    plot_licht = 0; % plot averaging of trials from Lichtwark US fascicle tracking



    %% Set constants and globals % PROJECTSPECIFIC

    % sampling frequencies
    global us_zerodispframes noraxonfreq freq_default
    us_zerodispframes = 1; % No of US frames to average as zero displacement
    noraxonfreq = 1500; % sampling frequency of noraxon data
    freq_default = 100; % output frequency for noraxon data without US video (MVC, etc)
%     angle_step_plots = 0.05; % resampled, averaged data extracted every x degrees for PLOTS
%     angle_step_stats_abs = 0.5; % every x degrees
%     angle_step_stats_norm = 10;  % every x percent
    

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
    %% Set constants and globals 
    
           
    
    %% set AXES etc for plots
    txt_displ = 'Displacement (mm)';
    txt_gonio = 'Gonio angle (°)';
    
    
    
    
    %% Read datamaster file, to connect corresponding data files
    % Produces arrays with file names and variables per trial

    global dm_subjectno dm_timepoint dm_side dm_trial 
    global dm_ROM_sol1_NX dm_ROM_sol1_US dm_ROM_sol1_US_frame dm_ROM_sol2_NX dm_ROM_sol2_US dm_ROM_sol2_US_frame
    global dm_ROM_gmmtj1_NX dm_ROM_gmmtj1_US dm_ROM_gmmtj1_US_frame dm_ROM_gmmtj2_NX dm_ROM_gmmtj2_US dm_ROM_gmmtj2_US_frame
%    global dm_ROM_gmfas1_NX dm_ROM_gmfas1_US dm_ROM_gmfas1_US_frame dm_ROM_gmfas2_NX dm_ROM_gmfas2_US dm_ROM_gmfas2_US_frame  dm_ROM_gmfas1_licht dm_ROM_gmfas2_licht
    global dm_CPM_calc_NX dm_CPM_calc_US dm_CPM_calc_US_frame %dm_MVC_PF dm_MVC_DF dm_CPM_sol_NX dm_CPM_sol_US dm_CPM_sol_US_frame 
%    global dm_leg_length dm_GMmsc_penn dm_GMmsc_faslen dm_ankle_angle_rest % dm_at_SOL_length dm_at_GM_length 
%    global at_momentarm 
    global filepath
    dm_filename = 'data/datamaster_passive.tsv';
    linestotal = read_datamaster_passive(dm_filename);
    %% Read datamaster file, to connect corresponding data files
        
    
    
    %% preallocate output arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % common arrays for all subjects:
        all_passive_output_head = {'Subject', 'Time', 'Side', 'Trial', ...
            'Gonio start SOL1', 'Gonio start SOL2', 'Gonio start GMMTJ1', 'Gonio start GMMTJ2', ...
            'Norm start SOL1', 'Norm start SOL2', 'Norm start GMMTJ1', 'Norm start GMMTJ2', ...
            'Rot const SOL1', 'Rot const SOL2', 'Rot const GMMTJ1', 'Rot const GMMTJ2', ...
            'Rsquare SOL1', 'Rsquare SOL2', 'Rsquare GMMTJ1', 'Rsquare GMMTJ2', ...
            'Rot const CALC', ...
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
        % old CON_PRE_subject_ID(ceil(linestotal/4)) = zeros; % holds subject numbers for intervention study

        STR_PRE_no(ceil(linestotal)) = zeros;
        STR_POST_no(ceil(linestotal)) = zeros;
        CON_PRE_no(ceil(linestotal)) = zeros;
        CON_POST_no(ceil(linestotal)) = zeros;

    
    
    
    
    %% LOOP through all lines in datamaster file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for line = 1:linestotal 
        
        
        
        %% subject/trial identifier
%        trial_subjectno = str2double(dm_subjectno{line});
        trial_timepoint = strcmp(dm_timepoint{line},'POST'); % 0 = PRE, 1 = POST
        trial_leg = strcmp(dm_trial{line},'STR'); % 0 = CON, 1 = STR
%         trial_calf_length = str2double(dm_leg_length{line}) * 10;  % convert from input in cm to calf length in mm
%         prone_GMfas_length = str2double(dm_GMmsc_faslen{line}); % resting (prone) length of GM fascicle 
%         prone_GMfas_penn = str2double(dm_GMmsc_penn{line}); % resting (prone) pennation angle of GM 
%         prone_GMfas_ankle = -str2double(dm_ankle_angle_rest{line}); % plantarflexed angle - input is positive value, but must be neg for plots
        

            if input_project == 1
                filepath = 'data\BD\';
                subject_id = horzcat('BD', dm_subjectno{line}, '_', dm_trial{line}, '_', dm_timepoint{line}, '_', dm_side{line});
            else
                filepath = 'data\';
                subject_id = horzcat('INT', dm_subjectno{line}, '_', dm_trial{line}, '_', dm_timepoint{line}, '_', dm_side{line});
            end
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
        
        
        
        
        %% calculate NORM conversion factors
        % retrieve conversion constants for Norm data
        [convert_norm_angle_a, convert_norm_angle_b, convert_norm_torque_a, convert_norm_torque_b, convert_norm_velocity_a, convert_norm_velocity_b, convert_norm_direction_b] = calculate_norm_constants_complete(dm_subjectno{line}, dm_side{line}, dm_timepoint{line}, line, 'passive');
        
        
        %% Calculate relationship between calc displacement and ANKLE ROTATION

        if strcmp(dm_CPM_calc_US{line},'null')
            CALC_at_rotation_const = NaN;
        else
            % Read moment arm trial US data file, determine time stamps, set trigger frame as time = zero
            % Produce US sample frequency, create new US array containing time and displacement
            [usdata_CPM, usfreq_CPM] = read_us_file(strcat(filepath, dm_CPM_calc_US{line}, '.txt'), str2double(dm_CPM_calc_US_frame{line}), 'CPM calcaneus');

            % Read moment arm noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
            % Produce a new noraxon data array
            noraxon_CPM = read_noraxon_stiffness(strcat(filepath, dm_CPM_calc_NX{line}), usfreq_CPM, dm_side{line}, 'CPM calcaneus');

            CALC_at_rotation_const = calculate_rotation_correction(noraxon_CPM, usdata_CPM);
        end
        
        
        %% preparations
        % create displ/angle equation as cfit
        f = fittype('a*x+b');
        a = 0; % temporary values, will be filled/replaced
        b = 0;
        fit_ankle_rot = cfit(f, a, b);
        
                
        %% calculations for 2x SOL trials
        
        % extract force, gonio, angle, displacement for EACH TRIAL
        % NB: extract_force_displ_singletrial_passive_EMG is where torque is converted to force
        if(strcmpi(dm_ROM_sol1_NX{line}, 'null'))
            % allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case in final_stiffness
            SOL_gonio_1 = NaN;
            SOL_angle_1 = NaN;
            SOL_displ_1 = NaN;
            SOL_rsquare_1 = NaN;
            SOL_at_rotation_const_1 = [NaN NaN];
        else 
            [SOL_at_rotation_const_1, SOL_rsquare_1, SOL_gonio_1, SOL_angle_1, SOL_displ_1] = extract_rot_corr(dm_ROM_sol1_NX{line}, dm_ROM_sol1_US{line}, dm_ROM_sol1_US_frame{line}, dm_side{line}, 'SOL1');
        end
        if(strcmpi(dm_ROM_sol2_NX{line}, 'null'))
            SOL_gonio_2 = NaN;
            SOL_angle_2 = NaN;
            SOL_displ_2 = NaN;
            SOL_rsquare_2 = NaN;
            SOL_at_rotation_const_2 = [NaN NaN];
        else 
            [SOL_at_rotation_const_2, SOL_rsquare_2, SOL_gonio_2, SOL_angle_2, SOL_displ_2] = extract_rot_corr(dm_ROM_sol2_NX{line}, dm_ROM_sol2_US{line}, dm_ROM_sol2_US_frame{line}, dm_side{line}, 'SOL2');
        end

        % plot summary of 2 trials: Displacement-angle
        if plot_check
            plottitle = horzcat('IND ankle rot NEW SOL, ', subject_id);
            figure('Name',plottitle);
            hold on
            plot(SOL_gonio_1,SOL_displ_1,'b.')
            plot(SOL_gonio_2,SOL_displ_2,'r.')

            fit_ankle_rot.a = SOL_at_rotation_const_1(1); % a
            fit_ankle_rot.b = SOL_at_rotation_const_1(2); % b
            plot(fit_ankle_rot,'b-')
            fit_ankle_rot.a = SOL_at_rotation_const_2(1); % a
            fit_ankle_rot.b = SOL_at_rotation_const_2(2); % b
            plot(fit_ankle_rot,'r-')
            text(-3,0.5,horzcat('Coeffs: ', num2str(SOL_at_rotation_const_1(1)), ' / ', num2str(SOL_at_rotation_const_2(1))))
            
            ylabel(txt_displ)
            xlabel(txt_gonio)
            title(plottitle,'Interpreter', 'none')
            legend('Trial 1','Trial 2','Location','Northeast')
            print(horzcat('data_plots/',plottitle),'-dpng')
        end
        

        %% calculations for 2x GM MTJ trials
        % extract force, gonio, angle, displacement for EACH TRIAL
        % NB: extract_force_displ_singletrial_passive_EMG is where torque is converted to force
        if(strcmpi(dm_ROM_gmmtj1_NX{line}, 'null'))
            % allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case in final_stiffness
            GMMTJ_gonio_1 = NaN;
            GMMTJ_angle_1 = NaN;
            GMMTJ_displ_1 = NaN;
            GMMTJ_rsquare_1 = NaN;
            GMMTJ_at_rotation_const_1 = [NaN NaN];
        else 
            [GMMTJ_at_rotation_const_1, GMMTJ_rsquare_1, GMMTJ_gonio_1, GMMTJ_angle_1, GMMTJ_displ_1] = extract_rot_corr(dm_ROM_gmmtj1_NX{line}, dm_ROM_gmmtj1_US{line}, dm_ROM_gmmtj1_US_frame{line}, dm_side{line}, 'GMMTJ1');
        end
        if(strcmpi(dm_ROM_gmmtj2_NX{line}, 'null'))
            GMMTJ_gonio_2 = NaN;
            GMMTJ_angle_2 = NaN;
            GMMTJ_displ_2 = NaN;
            GMMTJ_rsquare_2 = NaN;
            GMMTJ_at_rotation_const_2 = [NaN NaN];
        else 
            [GMMTJ_at_rotation_const_2, GMMTJ_rsquare_2, GMMTJ_gonio_2, GMMTJ_angle_2, GMMTJ_displ_2] = extract_rot_corr(dm_ROM_gmmtj2_NX{line}, dm_ROM_gmmtj2_US{line}, dm_ROM_gmmtj2_US_frame{line}, dm_side{line}, 'GMMTJ2');
        end
        
        % plot summary of 2 trials: Displacement-angle
        if plot_check
            plottitle = horzcat('IND ankle rot NEW GMMTJ, ', subject_id);
            figure('Name',plottitle);
            hold on
            plot(GMMTJ_gonio_1,GMMTJ_displ_1,'b.')
            plot(GMMTJ_gonio_2,GMMTJ_displ_2,'r.')

            fit_ankle_rot.a = GMMTJ_at_rotation_const_1(1); % a
            fit_ankle_rot.b = GMMTJ_at_rotation_const_1(2); % b
            plot(fit_ankle_rot,'b-')
            fit_ankle_rot.a = GMMTJ_at_rotation_const_2(1); % a
            fit_ankle_rot.b = GMMTJ_at_rotation_const_2(2); % b
            plot(fit_ankle_rot,'r-')
            text(-3,0.5,horzcat('Coeffs: ', num2str(GMMTJ_at_rotation_const_1(1)), ' / ', num2str(GMMTJ_at_rotation_const_2(1))))
            
            %axis(axis_displ_GMMTJ)
            ylabel(txt_displ)
            xlabel(txt_gonio)
            title(plottitle,'Interpreter', 'none')
            legend('Trial 1','Trial 2','Location','Northeast')
            print(horzcat('data_plots/',plottitle),'-dpng')
        end
        
        
        %% prepare prone arrays for group plots & stats
        % TODO - extract arrays?
%          %% intervention study
%             if trial_timepoint == 0 && trial_leg == 1 % PRE, STR
%                 STR_PRE_prone(STR_PRE_count,:) = MTU_prone_vars;
%             elseif trial_timepoint == 1 && trial_leg == 1 % POST, STR
%                 STR_POST_prone(STR_POST_count,:) = MTU_prone_vars;
%             elseif trial_timepoint == 0 && trial_leg == 0 % PRE, CON
%                 CON_PRE_prone(CON_PRE_count,:) = MTU_prone_vars;
%             elseif trial_timepoint == 1 && trial_leg == 0 % POST, CON
%                 CON_POST_prone(CON_POST_count,:) = MTU_prone_vars;
%             end

        
        %% prepare arrays for individual trial data to file

        % txt trial ID
        all_passive_output_txt(line,:) = [dm_subjectno(line) dm_timepoint(line) dm_side(line) dm_trial(line)];
        
        % TODO: choose only one of two constants per site?
        
        % add data to a common array for all subjects    
        all_passive_output(line,:) = [ SOL_gonio_1(1) SOL_gonio_2(1) GMMTJ_gonio_1(1) GMMTJ_gonio_2(1) ...
            SOL_angle_1(1) SOL_angle_2(1) GMMTJ_angle_1(1) GMMTJ_angle_2(1) ...
            SOL_at_rotation_const_1(1) SOL_at_rotation_const_2(1) GMMTJ_at_rotation_const_1(1) GMMTJ_at_rotation_const_2(1) ...
            SOL_rsquare_1 SOL_rsquare_2 GMMTJ_rsquare_1 GMMTJ_rsquare_2 ...
            CALC_at_rotation_const ...
            ];
        
        save all_data_ankle_inloop
        close all
        
    end
    %% LOOP FINISHED  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save all_data_ankle
    
    
    %% OUTPUT individual trial data TO FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % write xls
    if ispc
        filename_output = strcat('data_output/all_anklerot_output_', datestr(now, 'yyyymmdd_HHMM'), '.xlsx');
        
        xlswrite(filename_output, all_passive_output_head, 1, 'A1')
        xlswrite(filename_output, all_passive_output_txt, 1, 'A2')
        xlswrite(filename_output, all_passive_output, 1, 'E2')
    else
        filename_output = strcat('data_output/all_anklerot_output_', datestr(now, 'yyyymmdd_HHMM'), '.csv');
        csvwrite(filename_output, all_passive_output)
    end
    
end