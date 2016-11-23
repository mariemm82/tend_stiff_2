%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main file for analysis of strength trials (isometric and isokinetic)
% Marie Moltubakk 26.10.2016
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




function [] = strength_analysis(input_plot)
    close all
    
    
    
    %%% Determine which plots to produce during script running
    global plot_check subject_id plot_individual plot_conversion

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
    plot_conversion = 1;
    



    %%% Set constants and globals % PROJECTSPECIFIC

    % sampling frequencies
    global noraxonfreq freq_default
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
    norm_volt_per_nm_a = 0.07490071; % from file "MMM conversion volt-torque 2014-DES NEW DATA.xlsx"
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
    

    
    
    
    
    
    
    
    %%% Read datamaster file, to connect corresponding data files
    %%% Produces arrays with file names and variables per trial, to be retrieved later

    global dm_subjectno dm_timepoint dm_side dm_trial 
    global dm_isomet_P10_1 dm_isomet_P10_2 dm_isomet_D00_1 dm_isomet_D00_2 dm_isomet_D05_1 dm_isomet_D05_2 dm_isomet_D10_1 dm_isomet_D10_2 dm_isomet_D15_1 dm_isomet_D15_2
    global dm_isokinD30 dm_isokinP30 dm_isokinP45 dm_isokinP60 dm_isokinP90
    global dm_CPM_calc_NX dm_CPM_sol_NX
    global dm_leg_length
    global filepath
    dm_filename = 'data_strength/datamaster_strength.tsv';
    dm_columns = 22; % number of data columns entered per subject % PROJECTSPECIFIC
    linestotal = read_datamaster_strength(dm_filename,dm_columns);
    
    
    
    
    % MMM TODO
    %%% preallocate output arrays
    % common arrays for all subjects:
    all_strength_output = zeros(ceil(linestotal),50); 
    all_strength_output_txt = cell(ceil(linestotal),4);

    % BD-SPECIFIC
    BD_count = 0;
    CON_count = 0;
    
%    BD_angle_vars{ceil(linestotal)} = zeros;
%    CON_angle_vars{ceil(linestotal)} = zeros;


    



    
    %%%%%%%%%%%%%%%% LOOP through all lines in datamaster file (except header line)
    for line = 1:linestotal;




        %%% subject/trial identifier
        subjectno = str2double(dm_subjectno{line});
        if subjectno > 100
            filepath = 'data\BD\'; % MMM TODO
            subject_id = horzcat('dancer ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
            BD_count = BD_count + 1; %BD-SPECIFIC
        else
            filepath = 'data_strength\';
            subject_id = horzcat('control ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
            CON_count = CON_count + 1; %BD-SPECIFIC
        end
        cprintf('*black', horzcat('----------------', subject_id, '------------------\n'))




        %%% Calculate data for muscle activation (EMG)
        % script below is from passiveUS.m, not modified for strength data,
        % and disabled because for now, we do not analyse EMG with strength
        % data

%         % prepare column placement
%         if strcmpi(dm_side{line},'R') == 1
%             column_tibant = column_r_tibant;
%             column_gm = column_r_gm;
%             column_gl = column_r_gl;
%             column_sol = column_r_sol;
%         else % left
%             column_tibant = column_l_tibant;
%             column_gm = column_l_gm;
%             column_gl = column_l_gl;
%             column_sol = column_l_sol;
%         end
% 
%         % Read noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample - DORSIFLEXION
%         % Produce a new noraxon data array
%         noraxon_mvc_dorsi = read_noraxon_stiffness(strcat(filepath, dm_MVC_DF{line}), freq_default, dm_side{line}, 'MVC dorsi');
% 
%         % Calculate co-activation constants
%         % Read complete, prepared noraxon array + number of frames to average (freq * time)
%         % Produce max torque, max EMG constants
%         [~,EMG_max_TA] = calculate_EMG_max(noraxon_mvc_dorsi, freq_default*(mvc_window_ms/1000), column_tibant, 1); % 1 = invert torque for dorsiflexion
% 
%         % Read noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample - PLANTAR FLEXION
%         % Produce a new noraxon data array
%         noraxon_mvc_plantar = read_noraxon_stiffness(strcat(filepath, dm_MVC_PF{line}), freq_default, dm_side{line}, 'MVC plantar');
% 
%         % Calculate co-activation constants
%         % Read complete, prepared noraxon array + number of frames to average (freq * time)
%         % Produce max torque, max EMG constants
%         [~,EMG_max_gm] = calculate_EMG_max(noraxon_mvc_plantar, freq_default*(mvc_window_ms/1000), column_gm, 0);
%         [~,EMG_max_gl] = calculate_EMG_max(noraxon_mvc_plantar, freq_default*(mvc_window_ms/1000), column_gl, 0);
%         [~,EMG_max_sol] = calculate_EMG_max(noraxon_mvc_plantar, freq_default*(mvc_window_ms/1000), column_sol, 0);



        
        %%% calculate NORM conversion factors
        % retrieve conversion constants for Norm data
        [convert_norm_angle_a, convert_norm_angle_b, convert_norm_torque_a, convert_norm_torque_b, convert_norm_velocity_a, convert_norm_velocity_b, convert_norm_direction_b] = calculate_norm_constants_complete(dm_subjectno{line}, dm_side{line}, dm_timepoint{line}, line, 'active');
        

        
        


        
        
        

        
        %%% calculations for ISOKINETIC trials
        
        % identify files
        isokinetic_data = {dm_isokinD30{line} dm_isokinP30{line} dm_isokinP45{line} dm_isokinP60{line} dm_isokinP90{line}};
        isokinetic_labels = {'isokin DF 30' 'isokin PF 30' 'isokin PF 45' 'isokin PF 60' 'isokin PF 90'};

        % preallocate % MMM TODO move out of loop per subject?
        clear isokinetic_torque_angle_work isokinetic_arrays
        isokinetic_torque_angle_work(length(isokinetic_data),4) = zeros();
        isokinetic_arrays{length(isokinetic_data)} = zeros();
        
        % extract data from Noraxon files
        for i = 1:length(isokinetic_data)
            if(strcmpi(isokinetic_data{i}, 'null'))
                % allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case later
                isokinetic_torque_angle_work(i,1) = 0;
            else 
                [isokinetic_torque_angle_work(i,1), isokinetic_torque_angle_work(i,2), isokinetic_torque_angle_work(i,3), isokinetic_torque_angle_work(i,4), isokinetic_arrays{i}] = extract_isokinetic(isokinetic_data{i}, dm_side{line}, isokinetic_labels{i});
                % torque, angle, velocity, work, data arrays
            end
        end
        
        % plot curves
        if plot_check
            plottitle = horzcat('Isokinetic torque-angle trials ', subject_id);
            figure('Name',plottitle)
            hold on
            for i = 1:length(isokinetic_data)
                plot(isokinetic_arrays{i}(:,2),isokinetic_arrays{i}(:,1))
            end
            % axis([-2 35 -1 2.5]) %VAR
            % set(ax, 'xdir','reverse')
            xlabel('Norm angle (deg)')
            ylabel('Isokinetic torque (Nm)')
            title(plottitle)
            legend(isokinetic_labels)
            % saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
        end
        
        
        % plot summary of trials
        if plot_check
            plottitle = horzcat('Isokinetic PLANTAR FLEXION torque/velocity ', subject_id);
            figure('Name',plottitle)
            plot(isokinetic_torque_angle_work(1,3), isokinetic_torque_angle_work(1,1), 'LineStyle','none','Marker','o','MarkerSize',6) % dorsiflexion - tweak to get colors of remaining lines to match other plots
            hold on
            for i = 2:length(isokinetic_torque_angle_work) % not including first entry (dorsiflexion)
                plot(isokinetic_torque_angle_work(i,3), isokinetic_torque_angle_work(i,1), 'LineStyle','none','Marker','o','MarkerSize',6)
            end
             axis([-100 20 0 max(isokinetic_torque_angle_work(2:end,1))*1.1]) %VAR
            ax = gca;
            set(ax, 'xdir','reverse')
            xlabel('Velocity (deg/s)')
            ylabel('Isokinetic peak torque (Nm)')
            title(plottitle)
            legend(isokinetic_labels(1:end),'location','SouthWest')
            % saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
        end
        
        
        % MMM TODO GOON
        % further calculations and averages for isokinetic output
        
        
        
        
        
        

        %%% calculations for ISOMETRIC trials
        % dm_isomet_P10_1 dm_isomet_P10_2 dm_isomet_D00_1 dm_isomet_D00_2 dm_isomet_D05_1 dm_isomet_D05_2 dm_isomet_D10_1 dm_isomet_D10_2 dm_isomet_D15_1 dm_isomet_D15_2
        
        % preallocate
        isometric_torque_angle(10,2) = zeros();
        
        isometric_data = {dm_isomet_P10_1{line} dm_isomet_P10_2{line} dm_isomet_D00_1{line} dm_isomet_D00_2{line} dm_isomet_D05_1{line} dm_isomet_D05_2{line} dm_isomet_D10_1{line} dm_isomet_D10_2{line} dm_isomet_D15_1{line} dm_isomet_D15_2{line}};
        isometric_labels = {'P10 trial 1' 'P10 trial 2' 'D00 trial 1' 'D00 trial 2' 'D05 trial 1' 'D05 trial 2' 'D10 trial 1' 'D10 trial 2' 'D15 trial 1' 'D15 trial 2'};
        
        for i = 1:length(isometric_data)
            if(strcmpi(isometric_data{i}, 'null'))
                % allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case later
            else 
                [isometric_torque_angle(i,1), isometric_torque_angle(i,2)] = extract_isometric(isometric_data{i}, dm_side{line}, isometric_labels{i});
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%% MMM TODO GOON
        % average two trials

        % plot summary of trials
        if plot_check
            plottitle = horzcat('Isometric torque-angle ', subject_id);
            figure('Name',plottitle)
            hold on
            for i = 1:length(isometric_torque_angle(:,2))
                if isometric_torque_angle(i,1) ~= 0
                    plot(isometric_torque_angle(i,2), isometric_torque_angle(i,1), 'LineStyle','none','Marker','o','MarkerSize',6)
                end
            end
            axis([-17 12 min(isometric_torque_angle(:,1))*.9 max(isometric_torque_angle(:,1))*1.1]) %VAR
            ax = gca;
            set(ax, 'xdir','reverse')
            xlabel('Norm angle (deg)')
            ylabel('Isometric peak torque (Nm)')
            title(plottitle)
            % legend(isometric_labels)
            % saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
        end
        

        
        
        
        break
        % MMM TODO GOON
        % final output data
        
            
            % all data in ONE cell, common angles, RAW data:
            BD_angle_vars{BD_count} = [ ...
                data_force_gonio(loc_angle_start:loc_angle_stop,col_angle) ...  1
                data_force_gonio(loc_angle_start:loc_angle_stop,col_force) ...  2
                data_force_gonio(loc_angle_start:loc_angle_stop,3) ...          3
                data_force_gonio(loc_angle_start:loc_angle_stop,4)...           4
                data_force_gonio(loc_angle_start:loc_angle_stop,5) ...          5
                data_SOL(loc_angle_start_SOL:loc_angle_stop_SOL,3) ...          6
                data_GMMTJ(loc_angle_start_GMMTJ:loc_angle_stop_GMMTJ,3) ...    7
                data_GMFAS(loc_angle_start_GMFAS:loc_angle_stop_GMFAS,3) ...    8
                MTU_length_array(:,4)-MTU_length_array(1,4) ...                 9
                MTU_length_array(:,2) ...                                   10
                MTU_length_array(:,3) ...                                   11
                MTU_length_array(:,4) ...                                   12
                data_force_gonio(loc_angle_start:loc_angle_stop,col_force)*at_momentarm ... % 13
                ];
                        
            
        
        
        
        
        %%% OUTPUT final individual data to file
        % MMM TODO GOON - below are just examples from PASSIVE
        
        % add data to a common array for all subjects    
        i = 1;
        
        % txt trial ID
        all_strength_output_txt(line,1) = dm_subjectno(line);
        all_strength_output_txt(line,2) = dm_timepoint(line);
        all_strength_output_txt(line,3) = dm_side(line);
        all_strength_output_txt(line,4) = dm_trial(line);
        
        % ROM
        all_strength_output(line,i) = out_ROM_trial_max;
        i = i+1;
        all_strength_output(line,i) = out_emg_sol_submax_1;
        i = i+1;
        all_strength_output(line,i) = out_emg_sol_submax_2;
        

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%% LOOP FINISHED
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATIONS ACROSS SUBJECTS - NUMBERS
    % MMM TODO  - below are just examples from PASSIVE
    
    %%%  mean and stdav of each subject's INDIVIDUAL MAX ROM, force, elong, EMG
    n_o_array_elements = 13; %VAR - number of elements in BD_angle_xxx arrays
    
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
    end
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%% OUTPUT KEY VARIABLES FOR ALL SUBJECTS TO FILE
    % write xls
    if ispc
        filename_output = strcat('data_output/all_strength_output_', datestr(now, 'yyyy-mm-dd HH-MM'));
        
        % MMM TODO  - below are just examples from PASSIVE
        all_sterngth_output_head = {'Subject', 'Time', 'Side', 'Trial', ...
            'ROM trial (°)', 'ROM subject (L/R/PRE/POST)', 'ROM common (all subjects)', '33% ROM', '67% ROM', ...
            }; % PROJECTSPECIFIC
        
        xlswrite(filename_output, all_sterngth_output_head, 1, 'A1')
        xlswrite(filename_output, all_strength_output_txt, 1, 'A2')
        xlswrite(filename_output, all_strength_output, 1, 'E2')
    else
        filename_output = strcat('data_output/all_strength_output_', datestr(now, 'yyyy-mm-dd HH-MM'), '.csv');
        csvwrite(filename_output, all_strength_output)
    end

    
    
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATIONS ACROSS SUBJECTS - ARRAYS
    
    
    
    %%% average ABSOLUTE arrays, up to all subjects' COMMON MAX ROM, for force, elong, EMG
    % MMM TODO  - below are just examples from PASSIVE
    
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
    
    
    clear BD_angle_vars_mean_tmp BD_angle_vars_norm_mean_tmp CON_angle_vars_mean_tmp CON_angle_vars_norm_mean_tmp
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%% PLOT GROUP FIGURES
    % MMM TODO  - below are just examples from PASSIVE
    
    
    
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
            xlabel('Gonio angle (deg)')
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
            xlabel('Gonio angle (deg)')
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
            xlabel('Gonio angle (deg)')
            ylabel('Force (N)')
            title(plottitle)
            %legend
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    if BD_count > 1 && CON_count > 1 && plot_check            
            % NORMALIZED
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
    
    
    
end