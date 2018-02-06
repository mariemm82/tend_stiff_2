%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate max ROM angles and torques across all subjects, for loading into passiveUS script
%    basis = main file for analysis of passive dorsiflexion with US
% Marie Moltubakk 4.2.2015
% 
% The scripts assume a certain structure for input files from Tracker and
% Noraxon. I.e. number of and order of EMG and other channels. If these
% scripts are to be used for other projects, some of the code must be
% modified. These lines are marked % PROJECTSPECIFIC
%
% INPUT ARGUMENT: 1 for checkup plots or 0 for no plots. e.g. create_angles_passive(1). 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function []=create_angles_passive(input_plot)
    close all
    
    load forces_output 
    load angles_output
    
    
    
    %% Determine which plots to produce during script running
    global plot_achilles plot_norm plot_emg plot_check plot_us subject_id plot_conversion

    plot_check = input_plot; % turn on/off main checkpoint plots
    plot_achilles = 0; % turn on/off troubleshoot plots
    plot_norm = 0;
    plot_emg = 0;  % RMS 3 EMG channels per trial
    plot_us = 0;
    plot_conversion = 0;
    


    %% Set constants % PROJECTSPECIFIC

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
    norm_volt_per_nm_a = 0.07490071; % from file "MMM conversion volt-torque 2014-DES NEW DATA.xlsx"
    norm_volt_per_nm_b = 0.69283161;

    % cutoff frequencies for filters
    global emg_bandpass emg_rms_ms mvc_window_ms 
    emg_bandpass = [10/(noraxonfreq/2) 500/(noraxonfreq/2)]; % cutoff frequencies for EMG butterworth bandpass filter
    emg_rms_ms = 100; % milliseconds RMS window, must be divisible by 4 - ref Basmajian 1985 = 50 ms, ref Aagaard paper Passive tensile stress = 200 ms
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

    

    %% Read datamaster file, to connect corresponding data files
    % Produces arrays with file names and variables per trial, to be retrieved later

    global dm_subjectno dm_timepoint dm_side dm_trial 
    global dm_ROM_sol1_NX dm_ROM_sol1_US dm_ROM_sol1_US_frame dm_ROM_sol2_NX dm_ROM_sol2_US dm_ROM_sol2_US_frame
    global dm_ROM_gmmtj1_NX dm_ROM_gmmtj1_US dm_ROM_gmmtj1_US_frame dm_ROM_gmmtj2_NX dm_ROM_gmmtj2_US dm_ROM_gmmtj2_US_frame
    global dm_ROM_gmfas1_NX dm_ROM_gmfas1_US dm_ROM_gmfas1_US_frame dm_ROM_gmfas2_NX dm_ROM_gmfas2_US dm_ROM_gmfas2_US_frame
    % global dm_MVC_PF dm_MVC_DF % dm_CPM_calc_NX dm_CPM_calc_US dm_CPM_calc_US_frame dm_CPM_sol_NX dm_CPM_sol_US dm_CPM_sol_US_frame 
    global dm_leg_length
    global at_momentarm
    global filepath
    dm_filename = 'data/datamaster_passive.tsv';
    linestotal = read_datamaster_passive(dm_filename);
    subject_last = max(str2double(dm_subjectno)); % highest subject number in table

    
    
    %% predefine output arrays
    % prepare output array for angles with values of 100, to avoid reporting zero angles when data are missing
    angles_output(1:subject_last, 1:13) = 100;
    % output array for forces
    forces_output(1:subject_last, 1:9) = 100000;

    
    %% Loop through all lines in datamaster file (except header line)
    for line = 1:linestotal



        %% subject/trial identifier
        subjectno = str2double(dm_subjectno{line});
        if subjectno > 100
            filepath = 'data\BD\';
            subject_id = horzcat('dancer ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
        else
            filepath = 'data\';
            subject_id = horzcat('intervention ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
        end
        cprintf('*black', horzcat('----------------', subject_id, '------------------\n'))



        %% Calculate data for muscle activation (EMG)
% EMG not needed to create max angles and forces. Only used for output to
% screen, and by setting EMG_max to 0, output will report "inf". Does not
% affect data.
        EMG_max_TA = 0;
        EMG_max_GM = 0;
        EMG_max_GL = 0;
        EMG_max_SOL = 0;
% 
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
%     %    [coact_max_torque,coact_max_EMG] = calculate_coactivation(noraxon_coact, freq_default*(mvc_window_ms/1000), dm_side{line});
%         [~,EMG_max_TA] = calculate_EMG_max(noraxon_mvc_dorsi, freq_default*(mvc_window_ms/1000), column_tibant, 1); % 1 = invert torque for dorsiflexion
% 
%         % Read noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample - PLANTAR FLEXION
%         % Produce a new noraxon data array
%         noraxon_mvc_plantar = read_noraxon_stiffness(strcat(filepath, dm_MVC_PF{line}), freq_default, dm_side{line}, 'MVC plantarflex'); %TMP BD
% 
%         % Calculate co-activation constants
%         % Read complete, prepared noraxon array + number of frames to average (freq * time)
%         % Produce max torque, max EMG constants
%     %    [coact_max_torque,coact_max_EMG] = calculate_coactivation(noraxon_coact, freq_default*(mvc_window_ms/1000), dm_side{line});
%         [~,EMG_max_GM] = calculate_EMG_max(noraxon_mvc_plantar, freq_default*(mvc_window_ms/1000), column_gm, 0);
%         [~,EMG_max_GL] = calculate_EMG_max(noraxon_mvc_plantar, freq_default*(mvc_window_ms/1000), column_gl, 0);
%         [~,EMG_max_SOL] = calculate_EMG_max(noraxon_mvc_plantar, freq_default*(mvc_window_ms/1000), column_sol, 0);



        %% calculate norm conversion factors
        % retrieve conversion constants for Norm data
        [convert_norm_angle_a, convert_norm_angle_b, convert_norm_torque_a, convert_norm_torque_b, convert_norm_velocity_a, convert_norm_velocity_b, convert_norm_direction_b] = calculate_norm_constants_complete(dm_subjectno{line}, dm_side{line}, dm_timepoint{line}, line, 'passive');



        %% calculate achilles tendon moment arm
        at_momentarm = calculate_momentarm(0, 0, dm_leg_length{line});


        %% Calculations for 2x GM MTJ trials

        % extract force, gonio, angle, displacement
        if(strcmp(dm_ROM_gmmtj1_NX{line}, 'null'))
            % allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case in final_stiffness
            force_gmmtj_1 = zeros(0);
        else 
            [force_gmmtj_1, gonio_gmmtj_1, angle_gmmtj_1, displacement_gmmtj_1, ~, ~, ~, time_gmmtj_1, torque_gmmtj_1] = extract_force_displ_singletrial_passive_EMG(dm_ROM_gmmtj1_NX{line}, dm_ROM_gmmtj1_US{line}, dm_ROM_gmmtj1_US_frame{line}, EMG_max_TA, EMG_max_GM, EMG_max_GL, EMG_max_SOL, dm_leg_length{line}, dm_side{line}, line, 'GMMTJ1');
        end
        if(strcmp(dm_ROM_gmmtj2_NX{line}, 'null'))
            force_gmmtj_2 = zeros(0);
        else 
            [force_gmmtj_2, gonio_gmmtj_2, angle_gmmtj_2, displacement_gmmtj_2,  ~, ~, ~, time_gmmtj_2, torque_gmmtj_2] = extract_force_displ_singletrial_passive_EMG(dm_ROM_gmmtj2_NX{line}, dm_ROM_gmmtj2_US{line}, dm_ROM_gmmtj2_US_frame{line}, EMG_max_TA, EMG_max_GM, EMG_max_GL, EMG_max_SOL, dm_leg_length{line}, dm_side{line}, line, 'GMMTJ2');    
        end

        % average two passive trials
        if(strcmp(dm_ROM_gmmtj1_NX{line}, 'null')==0 && strcmp(dm_ROM_gmmtj2_NX{line}, 'null')==0)
            % if 2 trials exist (none of variables are 'null')

            % plot summary of 2 trials
            if plot_check && plot_norm
                max_force = max([max(force_gmmtj_1) max(force_gmmtj_2)]);
                max_displ = max([max(displacement_gmmtj_1) max(displacement_gmmtj_2)]);
                min_displ = min([min(displacement_gmmtj_1) min(displacement_gmmtj_2)]);
       %         min_angle = min([min(gonio_gmmtj_1) min(gonio_gmmtj_2)]);
                plottitle = horzcat('T-A-D, 2 trials, GM MTJ, ', subject_id);
                figure('Name',plottitle);
                % top panel = force
                AXa = subplot(2,1,1);
                plot(gonio_gmmtj_1,force_gmmtj_1)
                hold on
                plot(gonio_gmmtj_2,force_gmmtj_2)
                set(get(AXa,'Ylabel'),'String','Force (N)')
         %       set(AXa,'XLim',[min_angle max_angle+1])
                set(AXa,'YLim',[0 max_force],'YTick',(0:round(max_force/10,-2):max_force))
                title(plottitle)
                % bottom panel = displacement
                AXc = subplot(2,1,2);
                plot(gonio_gmmtj_1,displacement_gmmtj_1)
                hold on
                plot(gonio_gmmtj_2,displacement_gmmtj_2)
                set(get(AXc,'Ylabel'),'String','Displacement (mm)')
       %         set(AXc,'XLim',[min_angle max_angle+1])
                set(AXc,'YLim',[min_displ max_displ],'YTick',(0:5:1.1*max_displ))
                xlabel('Gonio angle (deg)');
                legend('Trial 1','Trial 2','Location','Southeast');
            end

            data_gmmtj = average_passive_trials_EMG(force_gmmtj_1, gonio_gmmtj_1, angle_gmmtj_1, displacement_gmmtj_1, 0,0,0, time_gmmtj_1, torque_gmmtj_1, force_gmmtj_2, gonio_gmmtj_2, angle_gmmtj_2, displacement_gmmtj_2, 0,0,0, time_gmmtj_2, torque_gmmtj_2);

        elseif(strcmp(dm_ROM_gmmtj1_NX{line}, 'null')==0)
            % keep first trial
            data_gmmtj = average_passive_trials_EMG(force_gmmtj_1, gonio_gmmtj_1, angle_gmmtj_1, displacement_gmmtj_1, 0,0,0, time_gmmtj_1, torque_gmmtj_1);
            angle_gmmtj_2 = 10000;
        else
            % keep second trial
            data_gmmtj = average_passive_trials_EMG(force_gmmtj_2, gonio_gmmtj_2, angle_gmmtj_2, displacement_gmmtj_2, 0,0,0, time_gmmtj_2, torque_gmmtj_2);
            angle_gmmtj_1 = 10000;
        end



        %% Calculations for 2x GM fascicle trials

        % extract force, gonio, angle, displacement
        if(strcmp(dm_ROM_gmfas1_NX{line}, 'null'))
            % allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case in final_stiffness
            force_gmfas_1 = zeros(0);
        else 
            [force_gmfas_1, gonio_gmfas_1, angle_gmfas_1, displacement_gmfas_1,  ~, ~, ~, time_gmfas_1, torque_gmfas_1] = extract_force_displ_singletrial_passive_EMG(dm_ROM_gmfas1_NX{line}, dm_ROM_gmfas1_US{line}, dm_ROM_gmfas1_US_frame{line}, EMG_max_TA, EMG_max_GM, EMG_max_GL, EMG_max_SOL, dm_leg_length{line}, dm_side{line}, line, 'GMfas1');
        end
        if(strcmp(dm_ROM_gmfas2_NX{line}, 'null'))
            force_gmfas_2 = zeros(0);
        else 
            [force_gmfas_2, gonio_gmfas_2, angle_gmfas_2, displacement_gmfas_2,  ~, ~, ~, time_gmfas_2, torque_gmfas_2] = extract_force_displ_singletrial_passive_EMG(dm_ROM_gmfas2_NX{line}, dm_ROM_gmfas2_US{line}, dm_ROM_gmfas2_US_frame{line}, EMG_max_TA, EMG_max_GM, EMG_max_GL, EMG_max_SOL, dm_leg_length{line}, dm_side{line}, line, 'GMfas2');    
        end

        % average two passive trials
        if(strcmp(dm_ROM_gmfas1_NX{line}, 'null')==0 && strcmp(dm_ROM_gmfas2_NX{line}, 'null')==0)
            % if 2 trials exist (none of variables are 'null')

            % plot summary of 2 trials
            if plot_check && plot_norm
                max_force = max([max(force_gmfas_1) max(force_gmfas_2)]);
                max_displ = max([max(displacement_gmfas_1) max(displacement_gmfas_2)]);
                min_displ = min([min(displacement_gmfas_1) min(displacement_gmfas_2)]);
                plottitle = horzcat('T-A-D, 2 trials, GM FAS, ', subject_id);
                figure('Name',plottitle);
                % top panel = force
                AXa = subplot(2,1,1);
                plot(gonio_gmfas_1,force_gmfas_1)
                hold on
                plot(gonio_gmfas_2,force_gmfas_2)
                set(get(AXa,'Ylabel'),'String','Force (N)')
                set(AXa,'YLim',[0 max_force],'YTick',(0:round(max_force/10,-2):max_force))
                title(plottitle)
                % bottom panel = displacement
                AXc = subplot(2,1,2);
                plot(gonio_gmfas_1,displacement_gmfas_1)
                hold on
                plot(gonio_gmfas_2,displacement_gmfas_2)
                set(get(AXc,'Ylabel'),'String','Displacement (mm)')
                set(AXc,'YLim',[min_displ max_displ],'YTick',(0:5:1.1*max_displ))
                xlabel('Gonio angle (deg)');
                legend('Trial 1','Trial 2','Location','Southeast');

            end

            data_gmfas = average_passive_trials_EMG(force_gmfas_1, gonio_gmfas_1, angle_gmfas_1, displacement_gmfas_1, 0,0,0, time_gmfas_1, torque_gmfas_1, force_gmfas_2, gonio_gmfas_2, angle_gmfas_2, displacement_gmfas_2, 0,0,0, time_gmfas_2, torque_gmfas_2);
        elseif(strcmp(dm_ROM_gmfas1_NX{line}, 'null')==0) % trial 1 not null
            % keep first trial
            data_gmfas = average_passive_trials_EMG(force_gmfas_1, gonio_gmfas_1, angle_gmfas_1, displacement_gmfas_1, 0,0,0, time_gmfas_1, torque_gmfas_1);
            angle_gmfas_2 = 10000;
        else
            % keep second trial
            data_gmfas = average_passive_trials_EMG(force_gmfas_2, gonio_gmfas_2, angle_gmfas_2, displacement_gmfas_2, 0,0,0, time_gmfas_2, torque_gmfas_2);
            angle_gmfas_1 = 10000;
        end


        %% Calculations for 2x SOL trials

        % extract force, gonio, angle, displacement
        if(strcmp(dm_ROM_sol1_NX{line}, 'null'))
            % allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case in final_stiffness
            force_sol_1 = zeros(0);
        else 
            [force_sol_1, gonio_sol_1, angle_sol_1, displacement_sol_1, ~, ~, ~, time_sol_1, torque_sol_1] = extract_force_displ_singletrial_passive_EMG(dm_ROM_sol1_NX{line}, dm_ROM_sol1_US{line}, dm_ROM_sol1_US_frame{line}, EMG_max_TA, EMG_max_GM, EMG_max_GL, EMG_max_SOL, dm_leg_length{line}, dm_side{line}, line, 'sol1');
        end
        if(strcmp(dm_ROM_sol2_NX{line}, 'null'))
            force_sol_2 = zeros(0);
        else 
            [force_sol_2, gonio_sol_2, angle_sol_2, displacement_sol_2, ~, ~, ~, time_sol_2, torque_sol_2] = extract_force_displ_singletrial_passive_EMG(dm_ROM_sol2_NX{line}, dm_ROM_sol2_US{line}, dm_ROM_sol2_US_frame{line}, EMG_max_TA, EMG_max_GM, EMG_max_GL, EMG_max_SOL, dm_leg_length{line}, dm_side{line}, line, 'sol2');    
        end

        % average two passive trials
        if(strcmp(dm_ROM_sol1_NX{line}, 'null')==0 && strcmp(dm_ROM_sol2_NX{line}, 'null')==0)
            % if 2 trials exist (none of variables are 'null')

            % plot summary of 2 trials
            if plot_check && plot_norm
                max_force = max([max(force_sol_1) max(force_sol_2)]);
                max_displ = max([max(displacement_sol_1) max(displacement_sol_2)]);
                min_displ = min([min(displacement_sol_1) min(displacement_sol_2)]);
                plottitle = horzcat('T-A-D, 2 trials, SOL, ', subject_id);
                figure('Name',plottitle);
                % top panel = force
                AXa = subplot(2,1,1);
                plot(gonio_sol_1,force_sol_1)
                hold on
                plot(gonio_sol_2,force_sol_2)
                set(get(AXa,'Ylabel'),'String','Force (N)')
                set(AXa,'YLim',[0 max_force],'YTick',(0:round(max_force/10,-2):max_force))
                title(plottitle)
                % bottom panel = displacement
                AXc = subplot(2,1,2);
                plot(gonio_sol_1,displacement_sol_1)
                hold on
                plot(gonio_sol_2,displacement_sol_2)
                set(get(AXc,'Ylabel'),'String','Displacement (mm)')
                set(AXc,'YLim',[min_displ max_displ],'YTick',(0:5:1.1*max_displ))
                xlabel('Gonio angle (deg)');
                legend('Trial 1','Trial 2','Location','Southeast');

            end

            data_sol = average_passive_trials_EMG(force_sol_1, gonio_sol_1, angle_sol_1, displacement_sol_1, 0,0,0, time_sol_1, torque_sol_1, force_sol_2, gonio_sol_2, angle_sol_2, displacement_sol_2, 0,0,0, time_sol_2, torque_sol_2);
        elseif(strcmp(dm_ROM_sol1_NX{line}, 'null')==0)
            % keep first trial
            data_sol = average_passive_trials_EMG(force_sol_1, gonio_sol_1, angle_sol_1, displacement_sol_1, 0,0,0, time_sol_1, torque_sol_1);
            angle_sol_2 = 10000;
        else
            % keep second trial
            data_sol = average_passive_trials_EMG(force_sol_2, gonio_sol_2, angle_sol_2, displacement_sol_2, 0,0,0, time_sol_2, torque_sol_2);
            angle_sol_1 = 10000;
        end



        %% Combine SOL, GMMTJ, GMFAS

         % all three displacements separately 
        if plot_check
            plottitle = horzcat('IND 3x force displacement, ', subject_id);
            figure('Name',plottitle);
            max_force = max([max(data_sol(:,1)) max(data_gmmtj(:,1)) max(data_gmfas(:,1))]);
            max_displ = max([max(data_sol(:,3)) max(data_gmmtj(:,3)) max(data_gmfas(:,3))]);
            min_displ = min([min(data_sol(:,3)) min(data_gmmtj(:,3)) min(data_gmfas(:,3))]);
 %           min_angle = min([min(data_sol(:,2)) min(data_gmmtj(:,2)) min(data_gmfas(:,2))]);
            % top panel = force
            AXa = subplot(2,1,1);
            plot(data_sol(:,2),data_sol(:,1))
            hold on
            plot(data_gmmtj(:,2),data_gmmtj(:,1))
            plot(data_gmfas(:,2),data_gmfas(:,1))
            set(get(AXa,'Ylabel'),'String','Force (N)')
 %           set(AXa,'XLim',[min_angle max_angle+1])
            set(AXa,'YLim',[0 max_force],'YTick',(0:round(max_force/10,-2):max_force))
            title(plottitle)
            %legend('SOL','GMMTJ','GMFAS','Location','Southeast');
            % bottom panel = displacement
            AXc = subplot(2,1,2);
            plot(data_sol(:,2),data_sol(:,3))
            hold on
            plot(data_gmmtj(:,2),data_gmmtj(:,3))
            plot(data_gmfas(:,2),data_gmfas(:,3))
            set(get(AXc,'Ylabel'),'String','Displacement (mm)')
   %         set(AXc,'XLim',[min_angle max_angle+1])
            set(AXc,'YLim',[min_displ max_displ],'YTick',(0:5:1.1*max_displ))
            xlabel('Gonio angle (deg)');
            legend('SOL','GMMTJ','GMFAS','Location','Northwest');
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
        end
        
        % average data from the three scan sites
        data_force_gonio = average_passive_forces(data_sol, data_gmmtj, data_gmfas);
        
        
        
        %% extract max ROM angles per trial, ind max and common max, write to array all_output

        % determine output column (L/R, PRE/POST)
        if strcmp(dm_timepoint{line}, 'PRE') == 1
            if strcmp(dm_side{line}, 'R') == 1
                output_column = 2;
            else % L
                output_column = 3;
            end
        else % POST
            if strcmp(dm_side{line}, 'R') == 1
                output_column = 4;
            else % L
                output_column = 5;
            end
        end
        
        % select lowest max ROM ANGLE after averaging 6 trials
        angle_gonio = max(data_force_gonio(:,2));
        % alternative method for Norm, because Norm angles are not currently contained in the "data" arrays.
        angle_norm = min([max(angle_sol_1) max(angle_sol_1) max(angle_gmmtj_1) max(angle_gmmtj_2) max(angle_gmfas_1) max(angle_gmfas_2)]);   % NORM angles probably not to be used

        % write max angle to array
        angles_output(subjectno, 1) = subjectno;
        angles_output(subjectno, output_column) = angle_norm;        % 'ANG_PRE_R', 'ANG_PRE_L', 'ANG_POST_R', 'ANG_POST_L'                   % NORM angles probably not to be used
        angles_output(subjectno, output_column + 6) = angle_gonio;   % 'GON_PRE_R', 'GON_PRE_L', 'GON_POST_R', 'GON_POST_L'
        
        save angles_output angles_output
        
        
        
        %% extract max FORCES per trial, ind max and common max, write to array all_output
        
        % write max averaged FORCE to array
        forces_output(subjectno, 1) = subjectno;
        forces_output(subjectno, output_column) = max(data_force_gonio(:,1));             % 'FORCE_PRE_R', 'FORCE_PRE_L', 'FORCE_POST_R', 'FORCE_POST_L', 
        
        save forces_output forces_output
        
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% loop finished %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    save all_data_createangles_endloop
    
    

    %% for all trials/subjects in datamaster, extract and store common and ind max

    % common max (across all subjects)
    angle_norm_common_max = min(min(angles_output(:, 2:5)));
    angle_gonio_common_max = min(min(angles_output(:, 8:11)));
    force_common_max = min(min(forces_output(:, 2:5)));

    % ind max (per subject)
    for ind = 1:length(angles_output(:, 1))
       angles_output(ind, 6) = min(angles_output(ind, 2:5));    % ANG_IND_MAX
       angles_output(ind, 12) = min(angles_output(ind, 8:11));  % GON_IND_MAX
       angles_output(ind, 7) = angle_norm_common_max;           % ANG_COMMON_MAX
       angles_output(ind, 13) = angle_gonio_common_max;         % GON_COMMON_MAX
       
       forces_output(ind, 6) = min(forces_output(ind, 2:5));    % FORCE_IND_MAX
       forces_output(ind, 7) = force_common_max;                % FORCE_COMMON_MAX
       forces_output(ind, 8) = min([forces_output(ind, 2) forces_output(ind, 4)]);    % FORCE_IND_R_MAX
       forces_output(ind, 9) = min([forces_output(ind, 3) forces_output(ind, 5)]);    % FORCE_IND_L_MAX
    end

    % output TSV file to be read by matlab
    dlmwrite('angles_output.tsv', angles_output, 'delimiter','\t', 'precision',10)
    dlmwrite('forces_output.tsv', forces_output, 'delimiter','\t', 'precision',10)

    % output XLS file with headers
    if ispc 
        filename_output = 'angles_output.xls';
        angles_output_head = {'SUBJECT', 'ANG_PRE_R', 'ANG_PRE_L', 'ANG_POST_R', 'ANG_POST_L', 'ANG_IND_MAX', 'ANG_COMMON_MAX', ...
             'GON_PRE_R', 'GON_PRE_L', 'GON_POST_R', 'GON_POST_L', 'GON_IND_MAX', 'GON_COMMON_MAX'};
        xlswrite(filename_output, angles_output_head, 1, 'A1')
        xlswrite(filename_output, angles_output, 1, 'A2')
        
        filename_output = 'forces_output.xls';
        forces_output_head = {'SUBJECT', 'FORCE_PRE_R', 'FORCE_PRE_L', 'FORCE_POST_R', 'FORCE_POST_L', 'FORCE_IND_MAX', 'FORCE_COMMON_MAX', ...
             'FORCE_IND_R_MAX', 'FORCE_IND_L_MAX'};
        xlswrite(filename_output, forces_output_head, 1, 'A1')
        xlswrite(filename_output, forces_output, 1, 'A2')
    end 
    
    % output MAT files
    save angles_output angles_output
    save forces_output forces_output

end