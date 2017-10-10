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

function [] = strength_analysis(input_project, input_plot)
    close all
 
    
    
    
    
    %% PLOTS: Determine which plots to output 
    global plot_check subject_id plot_individual plot_conversion plot_norm

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
    plot_norm = 0;



    %% Set CONSTANTS and globals % PROJECTSPECIFIC

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
    norm_volt_per_nm_a = 0.07490071; % from file "conversion volt-torque 2014-DES NEW DATA.xlsx"
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

    
    %% set AXES etc for plots
    global mat_version
    col_lightblue = [0.6 0.8 1];
    col_lightred = [1 0.6 0.8];
    %col_orange = [1 0.75 0];
    %col_grey = [0.3 0.3 0.3];
	mat_version = version('-release');


    
    %% Read datamaster file, to connect corresponding data files
    %%% Produces arrays with file names and variables per trial, to be retrieved later

    global dm_subjectno dm_timepoint dm_side dm_trial 
    global dm_isomet_P10_1 dm_isomet_P10_2 dm_isomet_D00_1 dm_isomet_D00_2 dm_isomet_D05_1 dm_isomet_D05_2 dm_isomet_D10_1 dm_isomet_D10_2 dm_isomet_D15_1 dm_isomet_D15_2
    global dm_isokinD30 dm_isokinP30 dm_isokinP45 dm_isokinP60 dm_isokinP90
%    global dm_CPM_calc_NX dm_CPM_sol_NX global dm_leg_length %  ---- not used in main script
    global filepath
    dm_filename = 'data/datamaster_strength.tsv';
    dm_columns = 22; % number of data columns entered per subject % PROJECTSPECIFIC
    linestotal = read_datamaster_strength(dm_filename,dm_columns);
    
        

    %% preallocate output arrays
    
    % common arrays for numbers across all subjects:
    all_strength_output = nan(ceil(linestotal),20); 
    all_strength_output_txt = cell(ceil(linestotal),4);
    all_strength_isokin_angles = nan(ceil(linestotal),15*5); 
    all_strength_isokin_angles_head = cell(1,75);

    if input_project == 1 % BD study
        BD_count = 0;
        CON_count = 0;
        
        BD_angle_vars{ceil(linestotal)} = [];
        CON_angle_vars{ceil(linestotal)} = [];
        BD_angle_intervals_vars{ceil(linestotal)} = [];
        CON_angle_intervals_vars{ceil(linestotal)} = [];
        BD_data = nan(ceil(linestotal),20);
        CON_data = nan(ceil(linestotal),20);
        BD_no(ceil(linestotal)) = zeros;
        CON_no(ceil(linestotal)) = zeros;
    else % intervention study
        STR_PRE_count = 0;
        CON_PRE_count = 0;
        STR_POST_count = 0;
        CON_POST_count = 0;

        STR_PRE_angle_vars{ceil(linestotal)} = [];
        CON_PRE_angle_vars{ceil(linestotal)} = [];
        STR_PRE_angle_intervals_vars{ceil(linestotal)} = [];
        CON_PRE_angle_intervals_vars{ceil(linestotal)} = [];
        STR_PRE_data = nan(ceil(linestotal),20);
        CON_PRE_data = nan(ceil(linestotal),20);
        STR_PRE_id{ceil(linestotal)} = [];
        CON_PRE_id{ceil(linestotal)} = [];

        STR_POST_angle_vars{ceil(linestotal)} = [];
        CON_POST_angle_vars{ceil(linestotal)} = [];
        STR_POST_angle_intervals_vars{ceil(linestotal)} = [];
        CON_POST_angle_intervals_vars{ceil(linestotal)} = [];
        STR_POST_data = nan(ceil(linestotal),20);
        CON_POST_data = nan(ceil(linestotal),20);
        STR_POST_id{ceil(linestotal)} = [];
        CON_POST_id{ceil(linestotal)} = [];
    end
        
    % below variables will be deleted and again preallocated within the loop
    isokinetic_torque_angle_work(1:5,1:4) = nan;
    isokinetic_arrays{5} = [];
    isokinetic_intervals{5} = [];
    


    
    %% LOOP through all lines in datamaster file (except header line)
    for line = 1:linestotal

        

        %% subject/trial identifier
        subjectno = str2double(dm_subjectno{line});
        
        if input_project == 1 % BD study
            if subjectno > 100
                filepath = 'data\BD\';
                subject_id = horzcat('dancer ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
                BD_count = BD_count + 1;
                BD_no(BD_count) = str2double(dm_subjectno{line});
                isokinetic_labels = {'isokin DF 30' 'isokin PF 45' 'isokin PF 60' 'isokin PF 90' 'isokin PF 120'};
                %isokinetic_speeds = [30 45 60 90 120];
            else
                filepath = 'data/stretcher_strength/';
                subject_id = horzcat('control ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
                CON_count = CON_count + 1;
                CON_no(CON_count) = str2double(dm_subjectno{line});
                isokinetic_labels = {'isokin DF 30' 'isokin PF 30' 'isokin PF 45' 'isokin PF 60' 'isokin PF 90'};
                %isokinetic_speeds = [30 30 45 60 90];
            end
            
        else
            isokinetic_labels = {'isokin DF 30' 'isokin PF 30' 'isokin PF 45' 'isokin PF 60' 'isokin PF 90'};
            %isokinetic_speeds = [30 30 45 60 90];
            filepath = 'data/stretcher_strength/';
            subject_id = horzcat('INT_', dm_subjectno{line}, '_', dm_trial{line}, '_', dm_timepoint{line}, '_', dm_side{line});
            if strcmp(dm_trial{line},'STR') %  == 1 --> STR
                if strcmp(dm_timepoint{line},'PRE') %  == 1 --> PRE
                    STR_PRE_count = STR_PRE_count + 1;
                    STR_PRE_id{STR_PRE_count} = subject_id;
                else % == 0 --> POST
                    STR_POST_count = STR_POST_count + 1;
                    STR_POST_id{STR_POST_count} = subject_id;
                end
            else %  == 0 --> CON
                if strcmp(dm_timepoint{line},'PRE') %  == 1 --> PRE
                    CON_PRE_count = CON_PRE_count + 1;
                    CON_PRE_id{CON_PRE_count} = subject_id;
                else % == 0 --> POST
                    CON_POST_count = CON_POST_count + 1;
                    CON_POST_id{CON_POST_count} = subject_id;
                end
            end
        end
        
        cprintf('*black', horzcat('----------------', subject_id, '------------------\n'))
        
        
        
        %% Calculate data for muscle activation (EMG)
        % script below is from passiveUS.m, not modified for strength data,
        % and disabled because for now, we do not analyse EMG with strength
        % data. Not upgraded to deal with intervention data, only BD study

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


        
        %% calculate NORM conversion factors
        % retrieve conversion constants for Norm data
        [convert_norm_angle_a, convert_norm_angle_b, convert_norm_torque_a, convert_norm_torque_b, convert_norm_velocity_a, convert_norm_velocity_b, convert_norm_direction_b] = calculate_norm_constants_complete(subjectno, dm_side{line}, dm_timepoint{line}, line, 'active');
        
        
        
        %% calculations for ISOKINETIC trials
        
        % identify files
        isokinetic_data = {dm_isokinD30{line} dm_isokinP30{line} dm_isokinP45{line} dm_isokinP60{line} dm_isokinP90{line}};
        
        % preallocate
        clear isokinetic_torque_angle_work isokinetic_arrays isokinetic_intervals
        isokinetic_torque_angle_work(1:length(isokinetic_data),1:4) = nan;
        isokinetic_arrays{length(isokinetic_data)} = [];
        isokinetic_intervals{length(isokinetic_data)} = [];
        
        % extract data from Noraxon files
        for i = 1:length(isokinetic_data)
            if(strcmpi(isokinetic_data{i}, 'null'))
                % allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case later
                isokinetic_torque_angle_work(i,1:4) = NaN;
            else 
                [isokinetic_torque_angle_work(i,1), isokinetic_torque_angle_work(i,2), isokinetic_torque_angle_work(i,3), isokinetic_torque_angle_work(i,4), isokinetic_arrays{i}, isokinetic_intervals{i}] = extract_isokinetic(isokinetic_data{i}, dm_side{line}, isokinetic_labels{i}, subject_id);
                % torque, angle, velocity, work, data arrays
            end
        end
        
        % plot torque-angle curves
        if plot_individual
            plottitle = horzcat('ISOKINETIC torque-angle trials, ', subject_id);
            figure('Name',plottitle)
            hold on
            for i = 1:length(isokinetic_data)
                if ~isempty(isokinetic_arrays{i})
                    plot(isokinetic_arrays{i}(:,2),isokinetic_arrays{i}(:,1))
                else
                    % no data for particular velocity
                    isokinetic_labels(i) = [];
                end
            end
            ax = gca;
            ax.ColorOrderIndex = 1;
            for i = 1:length(isokinetic_data)
                if ~isempty(isokinetic_arrays{i})
                    plot(isokinetic_torque_angle_work(i,2), isokinetic_torque_angle_work(i,1), '*') %torque/angle of peak torque
                end
            end
            % axis([-2 35 -1 2.5]) %VAR
            % set(ax, 'xdir','reverse')
            xlabel('Ankle angle (°)')
            ylabel('Isokinetic torque (Nm)')
            title(plottitle,'Interpreter', 'none')
            legend(isokinetic_labels)
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
        end
        
        % plot torque-velocity summary
        if plot_individual
            plottitle = horzcat('ISOKINETIC torque-velocity, ', subject_id);
            figure('Name',plottitle)
            if isokinetic_torque_angle_work(1,3) > 0 % trial exists
                plot(isokinetic_torque_angle_work(1,3), isokinetic_torque_angle_work(1,1), 'LineStyle','none','Marker','o','MarkerSize',6,'MarkerFaceColor','auto') % dorsiflexion - tweak to get colors of remaining lines to match other plots
            end
            hold on
            for i = 2:length(isokinetic_torque_angle_work) % not including first entry (dorsiflexion)
                plot(isokinetic_torque_angle_work(i,3), isokinetic_torque_angle_work(i,1), 'LineStyle','none','Marker','o','MarkerSize',6,'MarkerFaceColor','auto')
            end
             axis([20 130 0 max(isokinetic_torque_angle_work(2:end,1))*1.1]) %VAR
%            ax = gca;
%            set(ax, 'xdir','reverse')
            xlabel('Velocity (°/s)')
            ylabel('Isokinetic peak torque (Nm)')
            title(plottitle,'Interpreter', 'none')
            legend(isokinetic_labels(1:end),'location','SouthEast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
        end
        
        
        
        %% calculations for ISOMETRIC trials
        
        % identify files
        isometric_data = {dm_isomet_P10_1{line} dm_isomet_P10_2{line} dm_isomet_D00_1{line} dm_isomet_D00_2{line} dm_isomet_D05_1{line} dm_isomet_D05_2{line} dm_isomet_D10_1{line} dm_isomet_D10_2{line} dm_isomet_D15_1{line} dm_isomet_D15_2{line}};
        isometric_labels = {'P10 trial 1' 'P10 trial 2' 'D00 trial 1' 'D00 trial 2' 'D05 trial 1' 'D05 trial 2' 'D10 trial 1' 'D10 trial 2' 'D15 trial 1' 'D15 trial 2'};

        % preallocate
        isometric_torque_angle(1:10,1:2) = nan;
        isometric_arrays = []; % clear previous content
        isometric_arrays{length(isometric_data)} = [];

        % extract data from Noraxon files
        for i = 1:length(isometric_data)
            if(strcmpi(isometric_data{i}, 'null'))
                % allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case later
                isometric_torque_angle(i,1:2) = NaN;
            else 
                [isometric_torque_angle(i,1), isometric_torque_angle(i,2), isometric_arrays{i}] = extract_isometric(isometric_data{i}, dm_side{line}, isometric_labels{i});
            end
        end
        
        % plot torque-time curves
        if plot_individual
            plottitle = horzcat('ISOMETRIC angle checkup, ', subject_id);
            figure('Name',plottitle)
            hold on
            for i = 1:length(isometric_data)
                if ~isempty(isometric_arrays{i})
                    plot(isometric_arrays{i}(:,2))
                else
                    % no data for particular trial
                end
            end
            ax = gca;
            ax.ColorOrderIndex = 1;
            for i = 1:length(isometric_data)
                if ~isempty(isometric_arrays{i})
                    plot(length(isometric_arrays{i}(:,2)), isometric_arrays{i}(end,2), '*') %endpoint
                end
            end
            % axis([-2 35 -1 2.5]) %VAR
            % set(ax, 'xdir','reverse')
            xlabel('Frame')
            ylabel('Ankle angle (°)')
            title(plottitle,'Interpreter', 'none')
            legend(isometric_labels)
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
        end
        
        
        % plot summary of trials
        if plot_individual
            plottitle = horzcat('ISOMETRIC torque-angle trials, ', subject_id);
            figure('Name',plottitle)
            hold on
            for i = 1:length(isometric_torque_angle(:,2))
                if isometric_torque_angle(i,1) ~= 0
                    plot(isometric_torque_angle(i,2), isometric_torque_angle(i,1), 'LineStyle','none','Marker','o','MarkerSize',6,'MarkerFaceColor','auto')
                end
            end
            axis([-17 12 min(isometric_torque_angle(:,1))*.9 max(isometric_torque_angle(:,1))*1.1]) %VAR
            %ax = gca;
            %set(ax, 'xdir','reverse')
            xlabel('Ankle angle (°)')
            ylabel('Isometric peak torque (Nm)')
            title(plottitle,'Interpreter', 'none')
            % legend(isometric_labels)
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
        end
        
        % final data for output
        % isometric_torque_angle(i,1), isometric_torque_angle(i,2)
        % comes from extract_isometric: torque_max,angle_at_torque_max
        isometric_10 = max(isometric_torque_angle(1,1), isometric_torque_angle(2,1));
        isometric_0 = max(isometric_torque_angle(3,1), isometric_torque_angle(4,1));
        isometric_n5 = max(isometric_torque_angle(5,1), isometric_torque_angle(6,1));
        isometric_n10 = max(isometric_torque_angle(7,1), isometric_torque_angle(8,1));
        isometric_n15 = max(isometric_torque_angle(9,1), isometric_torque_angle(10,1));

        
        cprintf('blue', horzcat('Isometric angle check: ', num2str(round(isometric_torque_angle(:,2)',0)), '.\n'))
        
        

        %% collect OUTPUT DATA FOR AVERAGE PLOTS
        
        % arrays for averaging isokinetic curves across subjects
        if input_project == 1 % BD study
            if subjectno > 100
                BD_angle_vars{BD_count} = isokinetic_arrays;
                BD_angle_intervals_vars{BD_count} = isokinetic_intervals;
                BD_data(BD_count,1) = isokinetic_torque_angle_work(1,1);
                BD_data(BD_count,2) = isokinetic_torque_angle_work(2,1);
                BD_data(BD_count,3) = isokinetic_torque_angle_work(3,1);
                BD_data(BD_count,4) = isokinetic_torque_angle_work(4,1);
                BD_data(BD_count,5) = isokinetic_torque_angle_work(5,1);
                BD_data(BD_count,6) = isokinetic_torque_angle_work(1,2);
                BD_data(BD_count,7) = isokinetic_torque_angle_work(2,2);
                BD_data(BD_count,8) = isokinetic_torque_angle_work(3,2);
                BD_data(BD_count,9) = isokinetic_torque_angle_work(4,2);
                BD_data(BD_count,10) = isokinetic_torque_angle_work(5,2);
                BD_data(BD_count,11) = isokinetic_torque_angle_work(1,4);
                BD_data(BD_count,12) = isokinetic_torque_angle_work(2,4);
                BD_data(BD_count,13) = isokinetic_torque_angle_work(3,4);
                BD_data(BD_count,14) = isokinetic_torque_angle_work(4,4);
                BD_data(BD_count,15) = isokinetic_torque_angle_work(5,4);
                BD_data(BD_count,16) = isometric_10;
                BD_data(BD_count,17) = isometric_0;
                BD_data(BD_count,18) = isometric_n5;
                BD_data(BD_count,19) = isometric_n10;
                BD_data(BD_count,20) = isometric_n15;
            else
                CON_angle_vars{CON_count} = isokinetic_arrays;
                CON_angle_intervals_vars{CON_count} = isokinetic_intervals;
                CON_data(CON_count,1) = isokinetic_torque_angle_work(1,1);
                CON_data(CON_count,2) = isokinetic_torque_angle_work(2,1);
                CON_data(CON_count,3) = isokinetic_torque_angle_work(3,1);
                CON_data(CON_count,4) = isokinetic_torque_angle_work(4,1);
                CON_data(CON_count,5) = isokinetic_torque_angle_work(5,1);
                CON_data(CON_count,6) = isokinetic_torque_angle_work(1,2);
                CON_data(CON_count,7) = isokinetic_torque_angle_work(2,2);
                CON_data(CON_count,8) = isokinetic_torque_angle_work(3,2);
                CON_data(CON_count,9) = isokinetic_torque_angle_work(4,2);
                CON_data(CON_count,10) = isokinetic_torque_angle_work(5,2);
                CON_data(CON_count,11) = isokinetic_torque_angle_work(1,4);
                CON_data(CON_count,12) = isokinetic_torque_angle_work(2,4);
                CON_data(CON_count,13) = isokinetic_torque_angle_work(3,4);
                CON_data(CON_count,14) = isokinetic_torque_angle_work(4,4);
                CON_data(CON_count,15) = isokinetic_torque_angle_work(5,4);
                CON_data(CON_count,16) = isometric_10;
                CON_data(CON_count,17) = isometric_0;
                CON_data(CON_count,18) = isometric_n5;
                CON_data(CON_count,19) = isometric_n10;
                CON_data(CON_count,20) = isometric_n15;
            end
            
        else % intervention study
            if strcmp(dm_trial{line},'STR') %  == 1 --> STR
                if strcmp(dm_timepoint{line},'PRE') %  == 1 --> PRE
                    STR_PRE_angle_vars{STR_PRE_count} = isokinetic_arrays;
                    STR_PRE_angle_intervals_vars{STR_PRE_count} = isokinetic_intervals;
                    STR_PRE_data(STR_PRE_count,1) = isokinetic_torque_angle_work(1,1);
                    STR_PRE_data(STR_PRE_count,2) = isokinetic_torque_angle_work(2,1);
                    STR_PRE_data(STR_PRE_count,3) = isokinetic_torque_angle_work(3,1);
                    STR_PRE_data(STR_PRE_count,4) = isokinetic_torque_angle_work(4,1);
                    STR_PRE_data(STR_PRE_count,5) = isokinetic_torque_angle_work(5,1);
                    STR_PRE_data(STR_PRE_count,6) = isokinetic_torque_angle_work(1,2);
                    STR_PRE_data(STR_PRE_count,7) = isokinetic_torque_angle_work(2,2);
                    STR_PRE_data(STR_PRE_count,8) = isokinetic_torque_angle_work(3,2);
                    STR_PRE_data(STR_PRE_count,9) = isokinetic_torque_angle_work(4,2);
                    STR_PRE_data(STR_PRE_count,10) = isokinetic_torque_angle_work(5,2);
                    STR_PRE_data(STR_PRE_count,11) = isokinetic_torque_angle_work(1,4);
                    STR_PRE_data(STR_PRE_count,12) = isokinetic_torque_angle_work(2,4);
                    STR_PRE_data(STR_PRE_count,13) = isokinetic_torque_angle_work(3,4);
                    STR_PRE_data(STR_PRE_count,14) = isokinetic_torque_angle_work(4,4);
                    STR_PRE_data(STR_PRE_count,15) = isokinetic_torque_angle_work(5,4);
                    STR_PRE_data(STR_PRE_count,16) = isometric_10;
                    STR_PRE_data(STR_PRE_count,17) = isometric_0;
                    STR_PRE_data(STR_PRE_count,18) = isometric_n5;
                    STR_PRE_data(STR_PRE_count,19) = isometric_n10;
                    STR_PRE_data(STR_PRE_count,20) = isometric_n15;
                else % == 0 --> POST
                    STR_POST_angle_vars{STR_POST_count} = isokinetic_arrays;
                    STR_POST_angle_intervals_vars{STR_POST_count} = isokinetic_intervals;
                    STR_POST_data(STR_POST_count,1) = isokinetic_torque_angle_work(1,1);
                    STR_POST_data(STR_POST_count,2) = isokinetic_torque_angle_work(2,1);
                    STR_POST_data(STR_POST_count,3) = isokinetic_torque_angle_work(3,1);
                    STR_POST_data(STR_POST_count,4) = isokinetic_torque_angle_work(4,1);
                    STR_POST_data(STR_POST_count,5) = isokinetic_torque_angle_work(5,1);
                    STR_POST_data(STR_POST_count,6) = isokinetic_torque_angle_work(1,2);
                    STR_POST_data(STR_POST_count,7) = isokinetic_torque_angle_work(2,2);
                    STR_POST_data(STR_POST_count,8) = isokinetic_torque_angle_work(3,2);
                    STR_POST_data(STR_POST_count,9) = isokinetic_torque_angle_work(4,2);
                    STR_POST_data(STR_POST_count,10) = isokinetic_torque_angle_work(5,2);
                    STR_POST_data(STR_POST_count,11) = isokinetic_torque_angle_work(1,4);
                    STR_POST_data(STR_POST_count,12) = isokinetic_torque_angle_work(2,4);
                    STR_POST_data(STR_POST_count,13) = isokinetic_torque_angle_work(3,4);
                    STR_POST_data(STR_POST_count,14) = isokinetic_torque_angle_work(4,4);
                    STR_POST_data(STR_POST_count,15) = isokinetic_torque_angle_work(5,4);
                    STR_POST_data(STR_POST_count,16) = isometric_10;
                    STR_POST_data(STR_POST_count,17) = isometric_0;
                    STR_POST_data(STR_POST_count,18) = isometric_n5;
                    STR_POST_data(STR_POST_count,19) = isometric_n10;
                    STR_POST_data(STR_POST_count,20) = isometric_n15;
                end
            else %  == 0 --> CON
                if strcmp(dm_timepoint{line},'PRE') %  == 1 --> PRE
                    CON_PRE_angle_vars{CON_PRE_count} = isokinetic_arrays;
                    CON_PRE_angle_intervals_vars{CON_PRE_count} = isokinetic_intervals;
                    CON_PRE_data(CON_PRE_count,1) = isokinetic_torque_angle_work(1,1);
                    CON_PRE_data(CON_PRE_count,2) = isokinetic_torque_angle_work(2,1);
                    CON_PRE_data(CON_PRE_count,3) = isokinetic_torque_angle_work(3,1);
                    CON_PRE_data(CON_PRE_count,4) = isokinetic_torque_angle_work(4,1);
                    CON_PRE_data(CON_PRE_count,5) = isokinetic_torque_angle_work(5,1);
                    CON_PRE_data(CON_PRE_count,6) = isokinetic_torque_angle_work(1,2);
                    CON_PRE_data(CON_PRE_count,7) = isokinetic_torque_angle_work(2,2);
                    CON_PRE_data(CON_PRE_count,8) = isokinetic_torque_angle_work(3,2);
                    CON_PRE_data(CON_PRE_count,9) = isokinetic_torque_angle_work(4,2);
                    CON_PRE_data(CON_PRE_count,10) = isokinetic_torque_angle_work(5,2);
                    CON_PRE_data(CON_PRE_count,11) = isokinetic_torque_angle_work(1,4);
                    CON_PRE_data(CON_PRE_count,12) = isokinetic_torque_angle_work(2,4);
                    CON_PRE_data(CON_PRE_count,13) = isokinetic_torque_angle_work(3,4);
                    CON_PRE_data(CON_PRE_count,14) = isokinetic_torque_angle_work(4,4);
                    CON_PRE_data(CON_PRE_count,15) = isokinetic_torque_angle_work(5,4);
                    CON_PRE_data(CON_PRE_count,16) = isometric_10;
                    CON_PRE_data(CON_PRE_count,17) = isometric_0;
                    CON_PRE_data(CON_PRE_count,18) = isometric_n5;
                    CON_PRE_data(CON_PRE_count,19) = isometric_n10;
                    CON_PRE_data(CON_PRE_count,20) = isometric_n15;
                else % == 0 --> POST
                    CON_POST_angle_vars{CON_POST_count} = isokinetic_arrays;
                    CON_POST_angle_intervals_vars{CON_POST_count} = isokinetic_intervals;
                    CON_POST_data(CON_POST_count,1) = isokinetic_torque_angle_work(1,1);
                    CON_POST_data(CON_POST_count,2) = isokinetic_torque_angle_work(2,1);
                    CON_POST_data(CON_POST_count,3) = isokinetic_torque_angle_work(3,1);
                    CON_POST_data(CON_POST_count,4) = isokinetic_torque_angle_work(4,1);
                    CON_POST_data(CON_POST_count,5) = isokinetic_torque_angle_work(5,1);
                    CON_POST_data(CON_POST_count,6) = isokinetic_torque_angle_work(1,2);
                    CON_POST_data(CON_POST_count,7) = isokinetic_torque_angle_work(2,2);
                    CON_POST_data(CON_POST_count,8) = isokinetic_torque_angle_work(3,2);
                    CON_POST_data(CON_POST_count,9) = isokinetic_torque_angle_work(4,2);
                    CON_POST_data(CON_POST_count,10) = isokinetic_torque_angle_work(5,2);
                    CON_POST_data(CON_POST_count,11) = isokinetic_torque_angle_work(1,4);
                    CON_POST_data(CON_POST_count,12) = isokinetic_torque_angle_work(2,4);
                    CON_POST_data(CON_POST_count,13) = isokinetic_torque_angle_work(3,4);
                    CON_POST_data(CON_POST_count,14) = isokinetic_torque_angle_work(4,4);
                    CON_POST_data(CON_POST_count,15) = isokinetic_torque_angle_work(5,4);
                    CON_POST_data(CON_POST_count,16) = isometric_10;
                    CON_POST_data(CON_POST_count,17) = isometric_0;
                    CON_POST_data(CON_POST_count,18) = isometric_n5;
                    CON_POST_data(CON_POST_count,19) = isometric_n10;
                    CON_POST_data(CON_POST_count,20) = isometric_n15;
                end
            end
            
        end
        



        %% prepare arrays for individual trial data to file
        % save following arrays:
        %     all_strength_output_txt
        %     all_strength_output
        %     all_strength_isokin_angles
        
        % txt trial ID
        all_strength_output_txt(line,:) = [dm_subjectno(line) dm_timepoint(line) dm_side(line) dm_trial(line)];

        % isokinetic:
        % - comes from extract_isokinetic: torque_max, torque_max_angle, torque_max_velocity, work_max, array_output

        % add data to a common array for all subjects    
        all_strength_output(line,:) = [isokinetic_torque_angle_work(1,1) isokinetic_torque_angle_work(2,1) isokinetic_torque_angle_work(3,1) isokinetic_torque_angle_work(4,1) isokinetic_torque_angle_work(5,1) ... % peak torque
            isokinetic_torque_angle_work(1,2) isokinetic_torque_angle_work(2,2) isokinetic_torque_angle_work(3,2) isokinetic_torque_angle_work(4,2) isokinetic_torque_angle_work(5,2) ... % angle of peak torque
            isokinetic_torque_angle_work(1,4) isokinetic_torque_angle_work(2,4) isokinetic_torque_angle_work(3,4) isokinetic_torque_angle_work(4,4) isokinetic_torque_angle_work(5,4) ... % work
            isometric_10 isometric_0 isometric_n5 isometric_n10 isometric_n15 ... % isometric
            ];

        % isokinetic @ common angles:
        angles_common2 = (-7.5:2.5:27.5)'; % must match selection in extr_isokin
        n_o_angles = length(angles_common2);

        % data
        for j=1:length(isokinetic_intervals)
            i_start = n_o_angles*(j-1) + 1;
            i_end = n_o_angles + n_o_angles*(j-1);
            if isempty(isokinetic_intervals{1,j})
                all_strength_isokin_angles(line,i_start:i_end) = NaN;
            else % data exist
                all_strength_isokin_angles(line,i_start:i_end) = isokinetic_intervals{1,j};
            end
        end
        
        % save after each trial, in case of crash during next trial
        save all_data_strength_inloop
    end
    %% LOOP FINISHED ==============================================
    save all_data_strength
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %% HEADERS FOR DATA SAVED PER SUBJECT

    % isokinetic @ common angles
    for j=1:length(isokinetic_intervals)
        for k=1:(length(angles_common2))
            loc = k + (j-1)*length(angles_common2);
            all_strength_isokin_angles_head{loc} = horzcat('Trial ', num2str(j), ', ', num2str(angles_common2(k)), '°');
        end
    end
    
    % other strength data
    if input_project == 1 % BD study
        all_strength_output_head = {'Subject', 'Time', 'Side', 'Trial', ...
            'Peak torque DF 30 (Nm)', 'Peak torque PF 30/45 (Nm)', 'Peak torque PF 45/60 (Nm)', 'Peak torque PF 60/90 (Nm)', 'Peak torque PF 90/120 (Nm)' ...
            'Angle of PT DF 30 (°)', 'Angle of PT PF 30/ (°)', 'Angle of PT PF 45/ (°)', 'Angle of PT PF 60/ (°)', 'Angle of PT PF 90/ (°)' ...
            'Work DF 30 (J)', 'Work PF 30/ (J)', 'Work PF 45/ (J)', 'Work PF 60/ (J)', 'Work PF 90/ (J)' ...
            'Isometric PF torque +10° (Nm)', 'Isometric PF torque 0° (Nm)', 'Isometric PF torque -5° (Nm)', 'Isometric PF torque -10° (Nm)', 'Isometric PF torque -15° (Nm)' ...
            }; % PROJECTSPECIFIC
    else
        all_strength_output_head = {'Subject', 'Time', 'Side', 'Group', ...
            'Peak torque DF 30 (Nm)', 'Peak torque PF 30 (Nm)', 'Peak torque PF 45 (Nm)', 'Peak torque PF 60 (Nm)', 'Peak torque PF 90 (Nm)' ...
            'Angle of PT DF 30 (°)', 'Angle of PT PF 30 (°)', 'Angle of PT PF 45 (°)', 'Angle of PT PF 60 (°)', 'Angle of PT PF 90 (°)' ...
            'Work DF 30 (J)', 'Work PF 30 (J)', 'Work PF 45 (J)', 'Work PF 60 (J)', 'Work PF 90 (J)' ...
            'Isometric PF torque +10° (Nm)', 'Isometric PF torque 0° (Nm)', 'Isometric PF torque -5° (Nm)', 'Isometric PF torque -10° (Nm)', 'Isometric PF torque -15° (Nm)' ...
            }; % PROJECTSPECIFIC
    end

    
        
    %% GROUP CALCULATIONS - mean and stdav for plots
    
    %%%  mean and stdav of each subject's isometric torque, isokinetic torque/angle/work
    if input_project == 1 % BD study
        BD_data_mean = mean(BD_data,'omitnan');
        CON_data_mean = mean(CON_data,'omitnan');
        BD_data_SD = std(BD_data,'omitnan');
        CON_data_SD = std(CON_data,'omitnan');
    else % intervention
        STR_PRE_data_mean = mean(STR_PRE_data,'omitnan');
        CON_PRE_data_mean = mean(CON_PRE_data,'omitnan');
        STR_PRE_data_SD = std(STR_PRE_data,'omitnan');
        CON_PRE_data_SD = std(CON_PRE_data,'omitnan');
        STR_POST_data_mean = mean(STR_POST_data,'omitnan');
        CON_POST_data_mean = mean(CON_POST_data,'omitnan');
        STR_POST_data_SD = std(STR_POST_data,'omitnan');
        CON_POST_data_SD = std(CON_POST_data,'omitnan');
    end
    
    
    %% GROUP CALCULATIONS - average arrays for plots
    
    
    
    %%% average isokinetic arrays
    
    
    
    if input_project == 1 % BD study
        
        % find common angles across all subjects
        if BD_count > 0
            BD_angle_max(1:5,1:BD_count) = nan;
            BD_angle_min(1:5,1:BD_count) = nan;
            for i = 1:BD_count
                for n = 1:5 %5 trials
                    if ~isempty(BD_angle_vars{i}{n})
                        BD_angle_max(n,i) = max(BD_angle_vars{i}{n}(:,2));
                        BD_angle_min(n,i) = min(BD_angle_vars{i}{n}(:,2));
                    else %empty
                        BD_angle_max(n,i) = 100;
                        BD_angle_min(n,i) = -100;
                    end
                end
            end
        else
            BD_angle_max = 100;
            BD_angle_min = -100;
        end
        if CON_count > 0
            CON_angle_max(1:5,1:CON_count) = nan;
            CON_angle_min(1:5,1:CON_count) = nan;
            for i = 1:CON_count
                for n = 1:5 %5 trials
                    if ~isempty(CON_angle_vars{i}{n})
                        CON_angle_max(n,i) = max(CON_angle_vars{i}{n}(:,2));
                        CON_angle_min(n,i) = min(CON_angle_vars{i}{n}(:,2));
                    else %empty
                        CON_angle_max(n,i) = 100;
                        CON_angle_min(n,i) = -100;
                    end
                end
            end
        else
            CON_angle_max = 100;
            CON_angle_min = -100;
        end

        % select lowest possible common starting angle, highest possible common ending angle
        angle_min = max([max(BD_angle_min) max(CON_angle_min)]);
        angle_max = min([min(BD_angle_max) min(CON_angle_max)]);
        angle_array = (ceil(10*angle_min)/10:0.25:floor(10*angle_max)/10)';

        % reshape and average all isokinetic trials
        if BD_count > 0
            BD_angle_vars_common{5} = [];
            % reshape
            for n = 1:5 %5 trials
                for i = 1:BD_count
                    if ~isempty(BD_angle_vars{i}{n})
                        BD_angle_vars_common{n}(:,i) = spline(BD_angle_vars{i}{n}(:,2), BD_angle_vars{i}{n}(:,1), angle_array);
                        % else - array will remain empty
                    end
                end
            end
            % average
            BD_angle_vars_mean(1:length(angle_array),1:5) = nan;
            BD_angle_vars_SD(1:length(angle_array),1:5) = nan;
            for n=1:5
                if ~isempty(BD_angle_vars_common{n})
                    BD_angle_vars_mean(:,n) = nanmean(BD_angle_vars_common{n}, 2);
                    BD_angle_vars_SD(:,n) = nanstd(BD_angle_vars_common{n},1,2);
                    % else - arrays will remain zero
                end
            end
        end
        if CON_count > 0
            CON_angle_vars_common{5} = [];
            % reshape
            for n = 1:5 %5 trials
                for i = 1:CON_count
                    if ~isempty(CON_angle_vars{i}{n})
                        CON_angle_vars_common{n}(:,i) = spline(CON_angle_vars{i}{n}(:,2), CON_angle_vars{i}{n}(:,1), angle_array);
                        % else - array will remain empty
                    end
                end
            end
            % average
            CON_angle_vars_mean(1:length(angle_array),1:5) = nan;
            CON_angle_vars_SD(1:length(angle_array),1:5) = nan;
            for n=1:5
                if ~isempty(CON_angle_vars_common{n})
                    CON_angle_vars_mean(:,n) = nanmean(CON_angle_vars_common{n}, 2);
                    CON_angle_vars_SD(:,n) = nanstd(CON_angle_vars_common{n},1,2);
                    % else - arrays will remain zero
                end
            end
        end
        % end BD study
        
        
        
    else % intervention
        
        % find common angles across all subjects
        if STR_PRE_count > 0
            STR_PRE_angle_max(1:5,1:STR_PRE_count) = nan;
            STR_PRE_angle_min(1:5,1:STR_PRE_count) = nan;
            for i = 1:STR_PRE_count
                for n = 1:5 %5 trials
                    if ~isempty(STR_PRE_angle_vars{i}{n})
                        STR_PRE_angle_max(n,i) = max(STR_PRE_angle_vars{i}{n}(:,2));
                        STR_PRE_angle_min(n,i) = min(STR_PRE_angle_vars{i}{n}(:,2));
                    else %empty
                        STR_PRE_angle_max(n,i) = 100;
                        STR_PRE_angle_min(n,i) = -100;
                    end
                end
            end
        else
            STR_PRE_angle_max = 100;
            STR_PRE_angle_min = -100;
        end
        if CON_PRE_count > 0
            CON_PRE_angle_max(1:5,1:CON_PRE_count) = nan;
            CON_PRE_angle_min(1:5,1:CON_PRE_count) = nan;
            for i = 1:CON_PRE_count
                for n = 1:5 %5 trials
                    if ~isempty(CON_PRE_angle_vars{i}{n})
                        CON_PRE_angle_max(n,i) = max(CON_PRE_angle_vars{i}{n}(:,2));
                        CON_PRE_angle_min(n,i) = min(CON_PRE_angle_vars{i}{n}(:,2));
                    else %empty
                        CON_PRE_angle_max(n,i) = 100;
                        CON_PRE_angle_min(n,i) = -100;
                    end
                end
            end
        else
            CON_PRE_angle_max = 100;
            CON_PRE_angle_min = -100;
        end
        if STR_POST_count > 0
            STR_POST_angle_max(1:5,1:STR_POST_count) = nan;
            STR_POST_angle_min(1:5,1:STR_POST_count) = nan;
            for i = 1:STR_POST_count
                for n = 1:5 %5 trials
                    if ~isempty(STR_POST_angle_vars{i}{n})
                        STR_POST_angle_max(n,i) = max(STR_POST_angle_vars{i}{n}(:,2));
                        STR_POST_angle_min(n,i) = min(STR_POST_angle_vars{i}{n}(:,2));
                    else %empty
                        STR_POST_angle_max(n,i) = 100;
                        STR_POST_angle_min(n,i) = -100;
                    end
                end
            end
        else
            STR_POST_angle_max = 100;
            STR_POST_angle_min = -100;
        end
        if CON_POST_count > 0
            CON_POST_angle_max(1:5,1:CON_POST_count) = nan;
            CON_POST_angle_min(1:5,1:CON_POST_count) = nan;
            for i = 1:CON_POST_count
                for n = 1:5 %5 trials
                    if ~isempty(CON_POST_angle_vars{i}{n})
                        CON_POST_angle_max(n,i) = max(CON_POST_angle_vars{i}{n}(:,2));
                        CON_POST_angle_min(n,i) = min(CON_POST_angle_vars{i}{n}(:,2));
                    else %empty
                        CON_POST_angle_max(n,i) = 100;
                        CON_POST_angle_min(n,i) = -100;
                    end
                end
            end
        else
            CON_POST_angle_max = 100;
            CON_POST_angle_min = -100;
        end
        
        
        % select lowest possible common starting angle, highest possible common ending angle
        angle_min = max([max(STR_PRE_angle_min) max(CON_PRE_angle_min) max(STR_POST_angle_min) max(CON_POST_angle_min)]);
        angle_max = min([min(STR_PRE_angle_max) min(CON_PRE_angle_max) min(STR_POST_angle_max) min(CON_POST_angle_max)]);
        angle_array = (ceil(10*angle_min)/10:0.25:floor(10*angle_max)/10)';
        
        
        % reshape and average all isokinetic trials
        if STR_PRE_count > 0
            STR_PRE_angle_vars_common{5} = [];
            % reshape
            for n = 1:5 %5 trials
                for i = 1:STR_PRE_count
                    if ~isempty(STR_PRE_angle_vars{i}{n})
                        STR_PRE_angle_vars_common{n}(:,i) = spline(STR_PRE_angle_vars{i}{n}(:,2), STR_PRE_angle_vars{i}{n}(:,1), angle_array);
                        % else - array will remain empty
                    end
                end
            end
            % average
            STR_PRE_angle_vars_mean(1:length(angle_array),1:5) = nan;
            STR_PRE_angle_vars_SD(1:length(angle_array),1:5) = nan;
            for n=1:5
                if ~isempty(STR_PRE_angle_vars_common{n})
                    STR_PRE_angle_vars_mean(:,n) = nanmean(STR_PRE_angle_vars_common{n}, 2);
                    STR_PRE_angle_vars_SD(:,n) = nanstd(STR_PRE_angle_vars_common{n},1,2);
                    % else - arrays will remain zero
                end
            end
        end
        if CON_PRE_count > 0
            CON_PRE_angle_vars_common{5} = [];
            % reshape
            for n = 1:5 %5 trials
                for i = 1:CON_PRE_count
                    if ~isempty(CON_PRE_angle_vars{i}{n})
                        CON_PRE_angle_vars_common{n}(:,i) = spline(CON_PRE_angle_vars{i}{n}(:,2), CON_PRE_angle_vars{i}{n}(:,1), angle_array);
                        % else - array will remain empty
                    end
                end
            end
            % average
            CON_PRE_angle_vars_mean(1:length(angle_array),1:5) = nan;
            CON_PRE_angle_vars_SD(1:length(angle_array),1:5) = nan;
            for n=1:5
                if ~isempty(CON_PRE_angle_vars_common{n})
                    CON_PRE_angle_vars_mean(:,n) = nanmean(CON_PRE_angle_vars_common{n}, 2);
                    CON_PRE_angle_vars_SD(:,n) = nanstd(CON_PRE_angle_vars_common{n},1,2);
                    % else - arrays will remain zero
                end
            end
        end
        if STR_POST_count > 0
            STR_POST_angle_vars_common{5} = [];
            % reshape
            for n = 1:5 %5 trials
                for i = 1:STR_POST_count
                    if ~isempty(STR_POST_angle_vars{i}{n})
                        STR_POST_angle_vars_common{n}(:,i) = spline(STR_POST_angle_vars{i}{n}(:,2), STR_POST_angle_vars{i}{n}(:,1), angle_array);
                        % else - array will remain empty
                    end
                end
            end
            % average
            STR_POST_angle_vars_mean(1:length(angle_array),1:5) = nan;
            STR_POST_angle_vars_SD(1:length(angle_array),1:5) = nan;
            for n=1:5
                if ~isempty(STR_POST_angle_vars_common{n})
                    STR_POST_angle_vars_mean(:,n) = nanmean(STR_POST_angle_vars_common{n}, 2);
                    STR_POST_angle_vars_SD(:,n) = nanstd(STR_POST_angle_vars_common{n},1,2);
                    % else - arrays will remain zero
                end
            end
        end
        if CON_POST_count > 0
            CON_POST_angle_vars_common{5} = [];
            % reshape
            for n = 1:5 %5 trials
                for i = 1:CON_POST_count
                    if ~isempty(CON_POST_angle_vars{i}{n})
                        CON_POST_angle_vars_common{n}(:,i) = spline(CON_POST_angle_vars{i}{n}(:,2), CON_POST_angle_vars{i}{n}(:,1), angle_array);
                        % else - array will remain empty
                    end
                end
            end
            % average
            CON_POST_angle_vars_mean(1:length(angle_array),1:5) = nan;
            CON_POST_angle_vars_SD(1:length(angle_array),1:5) = nan;
            for n=1:5
                if ~isempty(CON_POST_angle_vars_common{n})
                    CON_POST_angle_vars_mean(:,n) = nanmean(CON_POST_angle_vars_common{n}, 2);
                    CON_POST_angle_vars_SD(:,n) = nanstd(CON_POST_angle_vars_common{n},1,2);
                    % else - arrays will remain zero
                end
            end
        end
        % end intervention
        
    end
    

    
    %% OUTPUT individual trial data to file
    % write xls
    if ispc
        filename_output = strcat('data_output/all_strength_output_', datestr(now, 'yyyy-mm-dd HH-MM'));
        
        xlswrite(filename_output, all_strength_output_head, 1, 'A1')
        xlswrite(filename_output, all_strength_isokin_angles_head, 1, 'Y1')
        xlswrite(filename_output, all_strength_output_txt, 1, 'A2')
        xlswrite(filename_output, all_strength_output, 1, 'E2')
        xlswrite(filename_output, all_strength_isokin_angles, 1, 'Y2')
    else
        filename_output = strcat('data_output/all_strength_output_', datestr(now, 'yyyy-mm-dd HH-MM'), '.csv');
        csvwrite(filename_output, all_strength_output)
    end

    
    %% OUTPUT group arrays for STATS, TO FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % variables to export to file
    out_arrays_input_cols = [1 16 31 46 61]; % locations in all_strength_isokin_angles
    if input_project == 1 % BD study
        out_arrays_input_labels = {'DF30' 'PF30' 'PF45' 'PF60' 'PF90' 'PF120'};
    else % intervention
        out_arrays_input_labels = {'DF30' 'PF30' 'PF45' 'PF60' 'PF90'};
    end
    
    % MMM TODO? Some smoothing at some point? avg over 0.5 deg?
    % Flexigym: Torque data were interpolated using a spline function
    % to extract torque at 0.25° intervals. Statistics were applied to torque values at 5° intervals.
    %
    % current data: Torque exists at intervals corresponding to 200 hz. Statistics
    % applied to values at every 2.5° after using spline - no additional averaging/smoothing.
    
    
    
    if input_project == 1 % BD study
        % create table headers (subject numbers)
        out_arrays_headers{1+BD_count+CON_count} = [];
        out_arrays_headers{1} = 'Joint_angle';
        for i=1:BD_count
            out_arrays_headers{i+1} = strcat('BD ', num2str(BD_no(i)));
        end
        for i=1:CON_count
            out_arrays_headers{i+1+BD_count} = strcat('CON ', num2str(CON_no(i)));
        end

        % preallocate output arrays
        cols_abs = length(angles_common2);
        rows = BD_count+CON_count + 1; % adding 1 for column for joint angles
        out_arrays_abs(cols_abs,rows) = zeros;

        % organize and output table for each of the selected variables
        for var = 1:length(out_arrays_input_labels)
            % reset output arrays
            out_arrays_abs(1:cols_abs,1:rows) = zeros;

            % add as first column, joint angles: abs and normalized angles
            out_arrays_abs(:,1) = angles_common2;

            % add BD subjects first
            for subj = 1:BD_count
                if var < 3 % dorsiflexion, PF30
                    out_arrays_abs(:,subj+1) = NaN;
                else
                    correction = 1; % BD subjects' data are one row too far to the right (BD PF60 is in the column of COL PF30)
                    out_arrays_abs(:,subj+1) = all_strength_isokin_angles(subj,out_arrays_input_cols(var-correction):out_arrays_input_cols(var-correction)+14);
                end
            end

            % add CON subjects second
            for subj = 1:CON_count
                if var > 5 % plantarflexion, PF120
                    out_arrays_abs(:,subj+BD_count+1) = NaN;
                else
                    out_arrays_abs(:,subj+BD_count+1) = all_strength_isokin_angles(subj+BD_count,out_arrays_input_cols(var):out_arrays_input_cols(var)+14);
                end
            end

            % create tables and save as file
            out_arrays_abs_table = array2table(out_arrays_abs,'VariableNames',out_arrays_headers);
            filename_output = strcat('data_output/BD_arrays_strength_', out_arrays_input_labels{var}, '_', datestr(now, 'yyyy-mm-dd HH-MM'));
            writetable(out_arrays_abs_table,filename_output,'Delimiter','\t')

            clear out_arrays_abs_table
        end
        % end of BD study
        
    else % intervention
        % create table headers (subject numbers)
        out_arrays_headers{1+STR_PRE_count+CON_PRE_count+STR_POST_count+CON_POST_count} = [];
        out_arrays_headers{1} = 'Joint_angle';
        for i=1:STR_PRE_count
            out_arrays_headers{i+1} = STR_PRE_id{i};
        end
        for i=1:CON_PRE_count
            out_arrays_headers{i+1+STR_PRE_count} = CON_PRE_id{i};
        end
        for i=1:STR_POST_count
            out_arrays_headers{i+1+STR_PRE_count+CON_PRE_count} = STR_POST_id{i};
        end
        for i=1:CON_POST_count
            out_arrays_headers{i+1+STR_PRE_count+CON_PRE_count+STR_POST_count} = CON_POST_id{i};
        end
                
        % organize and output table with TORQUE-ANGLE data per velocity
        % (i.e. one file for DF30, one file for PF45, etc)

        % preallocate output arrays
        cols_abs = length(angles_common2);
        rows = STR_PRE_count+CON_PRE_count+STR_POST_count+CON_POST_count + 1; % adding 1 extra column for holding joint angles
        out_arrays_abs(cols_abs,rows) = zeros;

        for var = 1:length(out_arrays_input_labels)
            
            % reset output arrays
            out_arrays_abs(1:cols_abs,1:rows) = zeros;

            % add as first column, joint angles: abs and normalized angles
            out_arrays_abs(:,1) = angles_common2;

            % add STR PRE subjects first
            for subj = 1:STR_PRE_count
                out_arrays_abs(:,subj+1) = all_strength_isokin_angles(subj,out_arrays_input_cols(var):out_arrays_input_cols(var)+14);
            end
            
            % add CON PRE subjects second
            for subj = 1:CON_PRE_count
                out_arrays_abs(:,subj+STR_PRE_count+1) = all_strength_isokin_angles(subj+STR_PRE_count,out_arrays_input_cols(var):out_arrays_input_cols(var)+14);
            end

            % add STR POST
            for subj = 1:STR_POST_count
                out_arrays_abs(:,subj+STR_PRE_count+CON_PRE_count+1) = all_strength_isokin_angles(subj+STR_PRE_count+CON_PRE_count,out_arrays_input_cols(var):out_arrays_input_cols(var)+14);
            end
            
            % add CON POST
            for subj = 1:CON_POST_count
                out_arrays_abs(:,subj+STR_PRE_count+CON_PRE_count+STR_POST_count+1) = all_strength_isokin_angles(subj+STR_PRE_count+CON_PRE_count+STR_POST_count,out_arrays_input_cols(var):out_arrays_input_cols(var)+14);
            end
            
            % create tables and save as file
            out_arrays_abs_table = array2table(out_arrays_abs,'VariableNames',out_arrays_headers);
            filename_output = strcat('data_output/intervention_arrays_strength_torque-angle', out_arrays_input_labels{var}, '_', datestr(now, 'yyyy-mm-dd HH-MM'));
            writetable(out_arrays_abs_table,filename_output,'Delimiter','\t')

            clear out_arrays_abs_table
        end
               

        % end of intervention
        
    end

    
    
    %% PLOT GROUP FIGURES 
    
    if input_project == 1 % BD study ================================================================
        
        % isometric trials  /////////////////////////////////////////////////////////////////////////
        
        if BD_count > 1 && CON_count > 1 && plot_check
            %variables
            plot_angles = [10, 0, -5 -10, -15];
            plot_isom_torque_min = min([min(BD_data_mean(16:20)) min(CON_data_mean(16:20))]) - max([max(BD_data_SD(16:20)) max(CON_data_SD(16:20))]);
            plot_isom_torque_max = max([max(BD_data_mean(16:20)) max(CON_data_mean(16:20))]) + max([max(BD_data_SD(16:20)) max(CON_data_SD(16:20))]);
            
            plottitle = horzcat('GRP isometric torque-angle, plantar flexion');
            figure('Name',plottitle)
            hold on
            plot(plot_angles(1:3), BD_data_mean(16:18), '-ob', 'MarkerSize',6)
            plot(plot_angles(1:3), CON_data_mean(16:18), '-or', 'MarkerSize',6)
            plot(plot_angles(3:5), BD_data_mean(18:20), '-ob', 'MarkerSize',6)
            plot(plot_angles(3:5), CON_data_mean(18:20), ':or', 'MarkerSize',6)
            errorbar(plot_angles, BD_data_mean(16:20), BD_data_SD(16:20),'-b')
            errorbar(plot_angles, CON_data_mean(16:20), CON_data_SD(16:20),':r')
            %ax = gca;
            %set(ax, 'xdir','reverse')
            axis([-20 15 plot_isom_torque_min plot_isom_torque_max])
            xlabel('Ankle angle (°)')
            ylabel('Isometric peak torque (Nm)')
            title(plottitle,'Interpreter', 'none')
            legend('BD', 'CON', 'Location','Northeast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
        end
        
        
        % isokinetic trials /////////////////////////////////////////////////////////////////////////
        % trial 1 = CON ONLY: DF 30°/s
        % trial 2 = CON: PF 30°/s, BD: PF 45°/s
        % ...
        % trial 5 = CON: PF 90°/s, BD: PF 120°/s
        
        % plot torque-velocity summary - 4 velocities
        if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP ISOKINETIC torque-velocity (all velocities)');
            figure('Name',plottitle)
            hold on
            
            plot_velocity_BD = [45 60 90 120];
            plot_velocity_CON = [30 45 60 90];
            plot(plot_velocity_BD, BD_data_mean(2:5), '-ob', 'MarkerSize',6)
            plot(plot_velocity_CON, CON_data_mean(2:5), '-or', 'MarkerSize',6)
            errorbar(plot_velocity_BD, BD_data_mean(2:5), BD_data_SD(2:5),'-b')
            errorbar(plot_velocity_CON, CON_data_mean(2:5), CON_data_SD(2:5),':r')
            
            axis([20 130 0 150]) %VAR
            xlabel('Velocity (°/s)')
            ylabel('Isokinetic peak torque (Nm)')
            title(plottitle,'Interpreter', 'none')
            legend('BD','CON','location','NorthEast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
        end
        
        % plot torque-velocity summary - common velocities
        if BD_count > 1 && CON_count > 1 && plot_check
            plottitle = horzcat('GRP ISOKINETIC torque-velocity');
            figure('Name',plottitle)
            hold on
            
            plot_velocity_BD = [45 60 90];
            plot_velocity_CON = [45 60 90];
            plot(plot_velocity_BD, BD_data_mean(2:4), '-ob', 'MarkerSize',6)
            plot(plot_velocity_CON, CON_data_mean(3:5), '-or', 'MarkerSize',6)
            errorbar(plot_velocity_BD, BD_data_mean(2:4), BD_data_SD(2:4),'-b')
            errorbar(plot_velocity_CON, CON_data_mean(3:5), CON_data_SD(3:5),':r')
            
            axis([40 100 0 140]) %VAR
            xlabel('Velocity (°/s)')
            ylabel('Isokinetic peak torque (Nm)')
            title(plottitle,'Interpreter', 'none')
            legend('BD','CON','location','NorthEast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
        end
        
        
        
        % isokinetic velocities, 1 figure per plot /////////////////////////////////////////////////////////////////////////
        
        if BD_count > 1 && CON_count > 1 && plot_check
            
            plottitle = horzcat('GRP ISOKINETIC torque-angle, 30°-s, plantar flexion');
            figure('Name',plottitle)
            hold on
            plot(angle_array, CON_angle_vars_mean(:,2),'r','LineWidth',1.5,'LineStyle','--')
            plot(angle_array, CON_angle_vars_mean(:,2)+CON_angle_vars_SD(:,2),'r','LineWidth',0.25,'LineStyle','--')
            plot(angle_array, CON_angle_vars_mean(:,2)-CON_angle_vars_SD(:,2),'r','LineWidth',0.25,'LineStyle','--')
            plot(CON_data_mean(7), CON_data_mean(2), '*r')
            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('CON', 'Location','Northeast')
            title('Isokinetic plantar flexion, 30°/s');
            axis([-10 30 30 150]) %VAR
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
            
            plottitle = horzcat('GRP ISOKINETIC torque-angle, 45°-s, plantar flexion');
            figure('Name',plottitle)
            hold on
            plot(angle_array, BD_angle_vars_mean(:,2),'m','LineWidth',1)
            plot(angle_array, CON_angle_vars_mean(:,3),'m','LineWidth',1.5,'LineStyle','--')
            plot(angle_array, BD_angle_vars_mean(:,2)+BD_angle_vars_SD(:,2),'m','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, CON_angle_vars_mean(:,3)+CON_angle_vars_SD(:,3),'m','LineWidth',0.25,'LineStyle','--')
            plot(angle_array, BD_angle_vars_mean(:,2)-BD_angle_vars_SD(:,2),'m','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, CON_angle_vars_mean(:,3)-CON_angle_vars_SD(:,3),'m','LineWidth',0.25,'LineStyle','--')
            plot(BD_data_mean(7), BD_data_mean(2), '*m')
            plot(CON_data_mean(8), CON_data_mean(3), '*m')
            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('BD', 'CON', 'Location','Northeast')
            title('Isokinetic plantar flexion, 45°/s');
            axis([-10 30 30 150]) %VAR
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
            
            plottitle = horzcat('GRP ISOKINETIC torque-angle, 60°-s, plantar flexion');
            figure('Name',plottitle)
            hold on
            plot(angle_array, BD_angle_vars_mean(:,3),'g','LineWidth',1)
            plot(angle_array, CON_angle_vars_mean(:,4),'g','LineWidth',1.5,'LineStyle','--')
            plot(angle_array, BD_angle_vars_mean(:,3)+BD_angle_vars_SD(:,3),'g','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, CON_angle_vars_mean(:,4)+CON_angle_vars_SD(:,4),'g','LineWidth',0.25,'LineStyle','--')
            plot(angle_array, BD_angle_vars_mean(:,3)-BD_angle_vars_SD(:,3),'g','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, CON_angle_vars_mean(:,4)-CON_angle_vars_SD(:,4),'g','LineWidth',0.25,'LineStyle','--')
            plot(BD_data_mean(8), BD_data_mean(3), '*g')
            plot(CON_data_mean(9), CON_data_mean(4), '*g')
            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('BD', 'CON', 'Location','Northeast')
            title('Isokinetic plantar flexion, 60°/s');
            axis([-10 30 30 150]) %VAR
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
            
            plottitle = horzcat('GRP ISOKINETIC torque-angle, 90°-s, plantar flexion');
            figure('Name',plottitle)
            hold on
            plot(angle_array, BD_angle_vars_mean(:,4),'b','LineWidth',1)
            plot(angle_array, CON_angle_vars_mean(:,5),'b','LineWidth',1.5,'LineStyle','--')
            plot(angle_array, BD_angle_vars_mean(:,4)+BD_angle_vars_SD(:,4),'b','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, CON_angle_vars_mean(:,5)+CON_angle_vars_SD(:,5),'b','LineWidth',0.25,'LineStyle','--')
            plot(angle_array, BD_angle_vars_mean(:,4)-BD_angle_vars_SD(:,4),'b','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, CON_angle_vars_mean(:,5)-CON_angle_vars_SD(:,5),'b','LineWidth',0.25,'LineStyle','--')
            plot(BD_data_mean(9), BD_data_mean(4), '*b')
            plot(CON_data_mean(10), CON_data_mean(5), '*b')
            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('BD', 'CON', 'Location','Northeast')
            title('Isokinetic plantar flexion, 90°/s');
            axis([-10 30 30 150]) %VAR
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
            
            plottitle = horzcat('GRP ISOKINETIC torque-angle, 120°-s, plantar flexion');
            figure('Name',plottitle)
            hold on
            plot(angle_array, BD_angle_vars_mean(:,5),'k','LineWidth',1)
            plot(angle_array, BD_angle_vars_mean(:,5)+BD_angle_vars_SD(:,5),'k','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, BD_angle_vars_mean(:,5)-BD_angle_vars_SD(:,5),'k','LineWidth',0.25,'LineStyle',':')
            plot(BD_data_mean(10), BD_data_mean(5), '*k')
            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('BD', 'Location','Northeast')
            title('Isokinetic plantar flexion, 120°/s');
            axis([-10 30 30 150]) %VAR
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
            
            % subplots
            plottitle = horzcat('GRP ISOKINETIC torque-angle, plantar flexion');
            figure('Name',plottitle)
            
            subplot(2,2,1); % 30 deg/s
            hold on
            plot(angle_array, CON_angle_vars_mean(:,2),'r','LineWidth',1.5,'LineStyle','--')
            plot(angle_array, CON_angle_vars_mean(:,2)+CON_angle_vars_SD(:,2),'r','LineWidth',0.25,'LineStyle','--')
            plot(angle_array, CON_angle_vars_mean(:,2)-CON_angle_vars_SD(:,2),'r','LineWidth',0.25,'LineStyle','--')
            plot(CON_data_mean(7), CON_data_mean(2), '*r')
            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('CON', 'Location','Northeast')
            title('Isokinetic plantar flexion, 30°/s');
            axis([-10 30 30 150]) %VAR
            
            subplot(2,2,2); % 45 deg/s
            hold on
            plot(angle_array, BD_angle_vars_mean(:,2),'m','LineWidth',1)
            plot(angle_array, CON_angle_vars_mean(:,3),'m','LineWidth',1.5,'LineStyle','--')
            plot(angle_array, BD_angle_vars_mean(:,2)+BD_angle_vars_SD(:,2),'m','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, CON_angle_vars_mean(:,3)+CON_angle_vars_SD(:,3),'m','LineWidth',0.25,'LineStyle','--')
            plot(angle_array, BD_angle_vars_mean(:,2)-BD_angle_vars_SD(:,2),'m','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, CON_angle_vars_mean(:,3)-CON_angle_vars_SD(:,3),'m','LineWidth',0.25,'LineStyle','--')
            plot(BD_data_mean(7), BD_data_mean(2), '*m')
            plot(CON_data_mean(8), CON_data_mean(3), '*m')
            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('BD', 'CON', 'Location','Northeast')
            title('Isokinetic plantar flexion, 45°/s');
            axis([-10 30 30 150]) %VAR
            
            subplot(2,2,3); % 60 deg/s
            hold on
            plot(angle_array, BD_angle_vars_mean(:,3),'g','LineWidth',1)
            plot(angle_array, CON_angle_vars_mean(:,4),'g','LineWidth',1.5,'LineStyle','--')
            plot(angle_array, BD_angle_vars_mean(:,3)+BD_angle_vars_SD(:,3),'g','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, CON_angle_vars_mean(:,4)+CON_angle_vars_SD(:,4),'g','LineWidth',0.25,'LineStyle','--')
            plot(angle_array, BD_angle_vars_mean(:,3)-BD_angle_vars_SD(:,3),'g','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, CON_angle_vars_mean(:,4)-CON_angle_vars_SD(:,4),'g','LineWidth',0.25,'LineStyle','--')
            plot(BD_data_mean(8), BD_data_mean(3), '*g')
            plot(CON_data_mean(9), CON_data_mean(4), '*g')
            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('BD', 'CON', 'Location','Northeast')
            title('Isokinetic plantar flexion, 60°/s');
            axis([-10 30 30 150]) %VAR
            
            subplot(2,2,4); % 90 deg/s
            hold on
            plot(angle_array, BD_angle_vars_mean(:,4),'b','LineWidth',1)
            plot(angle_array, CON_angle_vars_mean(:,5),'b','LineWidth',1.5,'LineStyle','--')
            plot(angle_array, BD_angle_vars_mean(:,4)+BD_angle_vars_SD(:,4),'b','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, CON_angle_vars_mean(:,5)+CON_angle_vars_SD(:,5),'b','LineWidth',0.25,'LineStyle','--')
            plot(angle_array, BD_angle_vars_mean(:,4)-BD_angle_vars_SD(:,4),'b','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, CON_angle_vars_mean(:,5)-CON_angle_vars_SD(:,5),'b','LineWidth',0.25,'LineStyle','--')
            plot(BD_data_mean(9), BD_data_mean(4), '*b')
            plot(CON_data_mean(10), CON_data_mean(5), '*b')
            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('BD', 'CON', 'Location','Northeast')
            title('Isokinetic plantar flexion, 90°/s');
            axis([-10 30 30 150]) %VAR
            
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
        end
        
        % isokinetic trials, plot per group /////////////////////////////////////////////////////////////////////////
        
        if BD_count > 1 && plot_check
            plottitle = horzcat('GRP ISOKINETIC torque-angle, plantar flexion, dancers');
            figure('Name',plottitle)
            hold on
            
            % torque x4
            plot(angle_array, BD_angle_vars_mean(:,2),'m','LineWidth',1)
            plot(angle_array, BD_angle_vars_mean(:,3),'g','LineWidth',1)
            plot(angle_array, BD_angle_vars_mean(:,4),'b','LineWidth',1)
            plot(angle_array, BD_angle_vars_mean(:,5),'k','LineWidth',1)
            
            % SD torque + and - x4
            plot(angle_array, BD_angle_vars_mean(:,2)+BD_angle_vars_SD(:,2),'m','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, BD_angle_vars_mean(:,2)-BD_angle_vars_SD(:,2),'m','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, BD_angle_vars_mean(:,3)+BD_angle_vars_SD(:,3),'g','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, BD_angle_vars_mean(:,3)-BD_angle_vars_SD(:,3),'g','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, BD_angle_vars_mean(:,4)+BD_angle_vars_SD(:,4),'b','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, BD_angle_vars_mean(:,4)-BD_angle_vars_SD(:,4),'b','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, BD_angle_vars_mean(:,5)+BD_angle_vars_SD(:,5),'k','LineWidth',0.25,'LineStyle',':')
            plot(angle_array, BD_angle_vars_mean(:,5)-BD_angle_vars_SD(:,5),'k','LineWidth',0.25,'LineStyle',':')
            
            % peak torque, angle of peak torque x4
            plot(BD_data_mean(7), BD_data_mean(2), '*m')
            plot(BD_data_mean(8), BD_data_mean(3), '*g')
            plot(BD_data_mean(9), BD_data_mean(4), '*b')
            plot(BD_data_mean(10), BD_data_mean(5), '*k')
            
            % SD for peak torque, angle of peak torque x4
            errorbar(BD_data_mean(7), BD_data_mean(2), BD_data_SD(2), 'm')
            errorbar(BD_data_mean(8), BD_data_mean(3), BD_data_SD(3), 'g')
            errorbar(BD_data_mean(9), BD_data_mean(4), BD_data_SD(4), 'b')
            errorbar(BD_data_mean(10), BD_data_mean(5), BD_data_SD(5), 'k')
            herrorbar(BD_data_mean(7), BD_data_mean(2), BD_data_SD(7), 'k')
            herrorbar(BD_data_mean(8), BD_data_mean(3), BD_data_SD(8), 'k')
            herrorbar(BD_data_mean(9), BD_data_mean(4), BD_data_SD(9), 'k')
            herrorbar(BD_data_mean(10), BD_data_mean(5), BD_data_SD(10), 'k')
            
            axis([-10 30 30 150]) %VAR
            xlabel('Ankle angle (°)')
            ylabel('Isokinetic torque (Nm)')
            title(plottitle,'Interpreter', 'none')
            labels = {'isokin PF 45' 'isokin PF 60' 'isokin PF 90' 'isokin PF 120'};
            legend(labels, 'Location','Northeast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle),'-dpng')
        end
        
        if CON_count > 1 && plot_check
            plottitle = horzcat('GRP ISOKINETIC torque-angle, plantar flexion, controls');
            figure('Name',plottitle)
            hold on
            plot(angle_array, CON_angle_vars_mean(:,2),'r','LineWidth',1.5)
            plot(angle_array, CON_angle_vars_mean(:,3),'m','LineWidth',1.5)
            plot(angle_array, CON_angle_vars_mean(:,4),'g','LineWidth',1.5)
            plot(angle_array, CON_angle_vars_mean(:,5),'b','LineWidth',1.5)
            
            plot(angle_array, CON_angle_vars_mean(:,2)+CON_angle_vars_SD(:,2),'r','LineWidth',0.25,'LineStyle','--')
            plot(angle_array, CON_angle_vars_mean(:,2)-CON_angle_vars_SD(:,2),'r','LineWidth',0.25,'LineStyle','--')
            plot(angle_array, CON_angle_vars_mean(:,3)+CON_angle_vars_SD(:,3),'m','LineWidth',0.25,'LineStyle','--')
            plot(angle_array, CON_angle_vars_mean(:,3)-CON_angle_vars_SD(:,3),'m','LineWidth',0.25,'LineStyle','--')
            plot(angle_array, CON_angle_vars_mean(:,4)+CON_angle_vars_SD(:,4),'g','LineWidth',0.25,'LineStyle','--')
            plot(angle_array, CON_angle_vars_mean(:,4)-CON_angle_vars_SD(:,4),'g','LineWidth',0.25,'LineStyle','--')
            plot(angle_array, CON_angle_vars_mean(:,5)+CON_angle_vars_SD(:,5),'b','LineWidth',0.25,'LineStyle','--')
            plot(angle_array, CON_angle_vars_mean(:,5)-CON_angle_vars_SD(:,5),'b','LineWidth',0.25,'LineStyle','--')
            
            % peak torque, angle of peak torque x4
            plot(CON_data_mean(7), CON_data_mean(2), '*r')
            plot(CON_data_mean(8), CON_data_mean(3), '*m')
            plot(CON_data_mean(9), CON_data_mean(4), '*g')
            plot(CON_data_mean(10), CON_data_mean(5), '*b')
            
            % SD for peak torque, angle of peak torque x4
            errorbar(CON_data_mean(7), CON_data_mean(2), CON_data_SD(2), 'r')
            errorbar(CON_data_mean(8), CON_data_mean(3), CON_data_SD(3), 'm')
            errorbar(CON_data_mean(9), CON_data_mean(4), CON_data_SD(4), 'g')
            errorbar(CON_data_mean(10), CON_data_mean(5), CON_data_SD(5), 'b')
            herrorbar(CON_data_mean(7), CON_data_mean(2), CON_data_SD(7), 'k')
            herrorbar(CON_data_mean(8), CON_data_mean(3), CON_data_SD(8), 'k')
            herrorbar(CON_data_mean(9), CON_data_mean(4), CON_data_SD(9), 'k')
            herrorbar(CON_data_mean(10), CON_data_mean(5), CON_data_SD(10), 'k')
            
            axis([-10 30 30 150]) %VAR
            xlabel('Ankle angle (°)')
            ylabel('Isokinetic torque (Nm)')
            title(plottitle,'Interpreter', 'none')
            labels = {'isokin PF 30' 'isokin PF 45' 'isokin PF 60' 'isokin PF 90'};
            legend(labels, 'Location','Northeast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
        end
        
        % isokinetic trials, dorsiflex /////////////////////////////////////////////////////////////////////////
        
        if CON_count > 1 && plot_check
            plottitle = horzcat('GRP ISOKINETIC torque-angle, dorsiflexion');
            figure('Name',plottitle)
            hold on
            
            % avg curve and SD
            plot(angle_array, CON_angle_vars_mean(:,1),'b','LineWidth',1,'LineStyle','--')
            plot(angle_array, CON_angle_vars_mean(:,1)+CON_angle_vars_SD(:,1),'b','LineWidth',0.25,'LineStyle','--')
            plot(angle_array, CON_angle_vars_mean(:,1)-CON_angle_vars_SD(:,1),'b','LineWidth',0.25,'LineStyle','--')
            
            % peak torque, angle of peak torque
            plot(CON_data_mean(6), CON_data_mean(1), '*b')
            errorbar(CON_data_mean(6), CON_data_mean(1), CON_data_SD(1), 'b')
            herrorbar(CON_data_mean(6), CON_data_mean(1), CON_data_SD(6), 'k')
            
            ax = gca;
            set(ax, 'xdir','reverse')
            axis([-10 30 0 20]) %VAR
            xlabel('Ankle angle (°)')
            ylabel('Isokinetic torque (Nm)')
            title(plottitle,'Interpreter', 'none')
            legend('isokin DF 30', 'Location','Northeast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
        end
        % end BD study

        
        
        
        
        
        
        
    else % intervention study ================================================================
        
        
        
        % isometric trials /////////////////////////////////////////////////////////////////////////
        
        if STR_PRE_count > 1 && CON_PRE_count > 1 && STR_POST_count > 1 && CON_POST_count > 1 && plot_check
            %variables
            plot_angles = [10, 0, -5 -10, -15];
            plot_isom_torque_min = min([min(STR_PRE_data_mean(16:20)) min(CON_PRE_data_mean(16:20)) min(STR_POST_data_mean(16:20)) min(CON_POST_data_mean(16:20))]) - max([max(STR_PRE_data_SD(16:20)) max(CON_PRE_data_SD(16:20)) max(STR_POST_data_SD(16:20)) max(CON_POST_data_SD(16:20))]);
            plot_isom_torque_max = max([max(STR_PRE_data_mean(16:20)) max(CON_PRE_data_mean(16:20)) max(STR_POST_data_mean(16:20)) max(CON_POST_data_mean(16:20))]) + max([max(STR_PRE_data_SD(16:20)) max(CON_PRE_data_SD(16:20)) max(STR_POST_data_SD(16:20)) max(CON_POST_data_SD(16:20))]);
            
            plottitle = horzcat('GRP isometric torque-angle, plantar flexion');
            figure('Name',plottitle)
            hold on
            
            % straight lines from +10 to -5 degrees (all subjects)
            plot(plot_angles(1:3), STR_PRE_data_mean(16:18), 'Color',col_lightred,'LineStyle','-', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(plot_angles(1:3), STR_POST_data_mean(16:18), 'Color','r','LineStyle','-', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(plot_angles(1:3), CON_PRE_data_mean(16:18), 'Color',col_lightblue,'LineStyle','-', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(plot_angles(1:3), CON_POST_data_mean(16:18), 'Color','b','LineStyle','-', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            % dotted lines from -5 to -15 degrees (some subjects not present)
            plot(plot_angles(3:5), STR_PRE_data_mean(18:20), 'Color',col_lightred,'LineStyle',':', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(plot_angles(3:5), STR_POST_data_mean(18:20), 'Color','r','LineStyle',':', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(plot_angles(3:5), CON_PRE_data_mean(18:20), 'Color',col_lightblue,'LineStyle',':', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(plot_angles(3:5), CON_POST_data_mean(18:20), 'Color','b','LineStyle',':', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            
            errorbar(plot_angles, STR_PRE_data_mean(16:20), STR_PRE_data_SD(16:20), 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
            errorbar(plot_angles, STR_POST_data_mean(16:20), STR_POST_data_SD(16:20), 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
            errorbar(plot_angles, CON_PRE_data_mean(16:20), CON_PRE_data_SD(16:20), 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
            errorbar(plot_angles, CON_POST_data_mean(16:20), CON_POST_data_SD(16:20), 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
            
            %ax = gca;
            %set(ax, 'xdir','reverse')
            axis([-20 15 plot_isom_torque_min plot_isom_torque_max])
            xlabel('Ankle angle (°)')
            ylabel('Isometric peak torque (Nm)')
            title(plottitle,'Interpreter', 'none')
            legend('STR PRE','STR POST','CON PRE','CON POST', 'Location','Northeast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
        end
        
        
        
        % isokinetic trials - INTERVENTION 
        % trial 1 = DF 30°/s
        % trial 2 = PF 30°/s
        % ...
        % trial 5 = PF 90°/s
        
        % plot torque-velocity summary - 4 velocities /////////////////////////////////////////////////////////////////////////
        if STR_PRE_count > 1 && CON_PRE_count > 1 && STR_POST_count > 1 && CON_POST_count > 1 && plot_check
            plottitle = horzcat('GRP ISOKINETIC torque-velocity (all velocities)');
            figure('Name',plottitle)
            hold on
            
            plot_velocity_CON = [30 45 60 90];
            plot(plot_velocity_CON, STR_PRE_data_mean(2:5), 'Color',col_lightred,'LineStyle','-', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(plot_velocity_CON, CON_PRE_data_mean(2:5), 'Color','r','LineStyle','-', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(plot_velocity_CON, STR_POST_data_mean(2:5), 'Color',col_lightblue,'LineStyle','-', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(plot_velocity_CON, CON_POST_data_mean(2:5), 'Color','b','LineStyle','-', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')

            errorbar(plot_velocity_CON, STR_PRE_data_mean(2:5), STR_PRE_data_SD(2:5), 'Color',col_lightred,'LineStyle','-')
            errorbar(plot_velocity_CON, CON_PRE_data_mean(2:5), CON_PRE_data_SD(2:5),':r')
            errorbar(plot_velocity_CON, STR_POST_data_mean(2:5), STR_POST_data_SD(2:5), 'Color',col_lightblue,'LineStyle','-')
            errorbar(plot_velocity_CON, CON_POST_data_mean(2:5), CON_POST_data_SD(2:5),':r')
            
            axis([20 100 0 200]) %VAR
            xlabel('Velocity (°/s)')
            ylabel('Isokinetic peak torque (Nm)')
            title(plottitle,'Interpreter', 'none')
            legend('STR PRE','STR POST','CON PRE','CON POST', 'Location','Northeast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
        end
        
        
        
        % plot torque-angle per velocity /////////////////////////////////////////////////////////////////////////
        if STR_PRE_count > 1 && CON_PRE_count > 1 && STR_POST_count > 1 && CON_POST_count > 1 && plot_check

            plottitle = horzcat('GRP ISOKINETIC torque-angle, 30°-s, dorsiflexion');
            figure('Name',plottitle)
            hold on
            
            % mean lines
            plot(angle_array, STR_PRE_angle_vars_mean(:,1), 'Color',col_lightred,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(angle_array, STR_POST_angle_vars_mean(:,1), 'Color','r','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(angle_array, CON_PRE_angle_vars_mean(:,1), 'Color',col_lightblue,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(angle_array, CON_POST_angle_vars_mean(:,1), 'Color','b','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            
            % +/- SD
            plot(angle_array, STR_PRE_angle_vars_mean(:,1)+STR_PRE_angle_vars_SD(:,1), 'Color',col_lightred,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, STR_PRE_angle_vars_mean(:,1)-STR_PRE_angle_vars_SD(:,1), 'Color',col_lightred,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, STR_POST_angle_vars_mean(:,1)+STR_POST_angle_vars_SD(:,1), 'Color','r','LineStyle',':','LineWidth',0.25)
            plot(angle_array, STR_POST_angle_vars_mean(:,1)-STR_POST_angle_vars_SD(:,1), 'Color','r','LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_PRE_angle_vars_mean(:,1)+CON_PRE_angle_vars_SD(:,1), 'Color',col_lightblue,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_PRE_angle_vars_mean(:,1)-CON_PRE_angle_vars_SD(:,1), 'Color',col_lightblue,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_POST_angle_vars_mean(:,1)+CON_POST_angle_vars_SD(:,1), 'Color','b','LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_POST_angle_vars_mean(:,1)-CON_POST_angle_vars_SD(:,1), 'Color','b','LineStyle',':','LineWidth',0.25)
            
            % peaks
            plot(STR_PRE_data_mean(6), STR_PRE_data_mean(1), 'Color',col_lightred, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(STR_POST_data_mean(6), STR_POST_data_mean(1), 'Color','r', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(CON_PRE_data_mean(6), CON_PRE_data_mean(1), 'Color',col_lightblue, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(CON_POST_data_mean(6), CON_POST_data_mean(1), 'Color','b', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            
            % reverse X axis: start at 30 deg, dorsiflex into -10 deg
            ax = gca;
            set(ax, 'xdir','reverse')

            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('STR PRE','STR POST','CON PRE','CON POST', 'Location','Northeast')
            title('Isokinetic dorsiflexion, 30°/s');
            axis([-10 30 -Inf Inf]) %VAR
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
            
            
            
            plottitle = horzcat('GRP ISOKINETIC torque-angle, 30°-s, plantar flexion');
            figure('Name',plottitle)
            hold on
            
            % mean lines
            plot(angle_array, STR_PRE_angle_vars_mean(:,2), 'Color',col_lightred,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(angle_array, STR_POST_angle_vars_mean(:,2), 'Color','r','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(angle_array, CON_PRE_angle_vars_mean(:,2), 'Color',col_lightblue,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(angle_array, CON_POST_angle_vars_mean(:,2), 'Color','b','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            
            % +/- SD
            plot(angle_array, STR_PRE_angle_vars_mean(:,2)+STR_PRE_angle_vars_SD(:,2), 'Color',col_lightred,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, STR_PRE_angle_vars_mean(:,2)-STR_PRE_angle_vars_SD(:,2), 'Color',col_lightred,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, STR_POST_angle_vars_mean(:,2)+STR_POST_angle_vars_SD(:,2), 'Color','r','LineStyle',':','LineWidth',0.25)
            plot(angle_array, STR_POST_angle_vars_mean(:,2)-STR_POST_angle_vars_SD(:,2), 'Color','r','LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_PRE_angle_vars_mean(:,2)+CON_PRE_angle_vars_SD(:,2), 'Color',col_lightblue,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_PRE_angle_vars_mean(:,2)-CON_PRE_angle_vars_SD(:,2), 'Color',col_lightblue,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_POST_angle_vars_mean(:,2)+CON_POST_angle_vars_SD(:,2), 'Color','b','LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_POST_angle_vars_mean(:,2)-CON_POST_angle_vars_SD(:,2), 'Color','b','LineStyle',':','LineWidth',0.25)
            
            % peaks
            plot(STR_PRE_data_mean(7), STR_PRE_data_mean(2), 'Color',col_lightred, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(STR_POST_data_mean(7), STR_POST_data_mean(2), 'Color','r', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(CON_PRE_data_mean(7), CON_PRE_data_mean(2), 'Color',col_lightblue, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(CON_POST_data_mean(7), CON_POST_data_mean(2), 'Color','b', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            
            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('STR PRE','STR POST','CON PRE','CON POST', 'Location','Northeast')
            title('Isokinetic plantar flexion, 30°/s');
            axis([-10 30 -Inf Inf]) %VAR
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
            
            
            
            plottitle = horzcat('GRP ISOKINETIC torque-angle, 45°-s, plantar flexion');
            figure('Name',plottitle)
            hold on
            
            % mean lines
            plot(angle_array, STR_PRE_angle_vars_mean(:,3), 'Color',col_lightred,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(angle_array, STR_POST_angle_vars_mean(:,3), 'Color','r','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(angle_array, CON_PRE_angle_vars_mean(:,3), 'Color',col_lightblue,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(angle_array, CON_POST_angle_vars_mean(:,3), 'Color','b','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            
            % +/- SD
            plot(angle_array, STR_PRE_angle_vars_mean(:,3)+STR_PRE_angle_vars_SD(:,3), 'Color',col_lightred,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, STR_PRE_angle_vars_mean(:,3)-STR_PRE_angle_vars_SD(:,3), 'Color',col_lightred,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, STR_POST_angle_vars_mean(:,3)+STR_POST_angle_vars_SD(:,3), 'Color','r','LineStyle',':','LineWidth',0.25)
            plot(angle_array, STR_POST_angle_vars_mean(:,3)-STR_POST_angle_vars_SD(:,3), 'Color','r','LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_PRE_angle_vars_mean(:,3)+CON_PRE_angle_vars_SD(:,3), 'Color',col_lightblue,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_PRE_angle_vars_mean(:,3)-CON_PRE_angle_vars_SD(:,3), 'Color',col_lightblue,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_POST_angle_vars_mean(:,3)+CON_POST_angle_vars_SD(:,3), 'Color','b','LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_POST_angle_vars_mean(:,3)-CON_POST_angle_vars_SD(:,3), 'Color','b','LineStyle',':','LineWidth',0.25)
            
            % peaks
            plot(STR_PRE_data_mean(8), STR_PRE_data_mean(3), 'Color',col_lightred, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(STR_POST_data_mean(8), STR_POST_data_mean(3), 'Color','r', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(CON_PRE_data_mean(8), CON_PRE_data_mean(3), 'Color',col_lightblue, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(CON_POST_data_mean(8), CON_POST_data_mean(3), 'Color','b', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            
            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('STR PRE','STR POST','CON PRE','CON POST', 'Location','Northeast')
            title('Isokinetic plantar flexion, 45°/s');
            axis([-10 30 -Inf Inf]) %VAR
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
            
            

            plottitle = horzcat('GRP ISOKINETIC torque-angle, 60°-s, plantar flexion');
            figure('Name',plottitle)
            hold on
            
            % mean lines
            plot(angle_array, STR_PRE_angle_vars_mean(:,4), 'Color',col_lightred,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(angle_array, STR_POST_angle_vars_mean(:,4), 'Color','r','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(angle_array, CON_PRE_angle_vars_mean(:,4), 'Color',col_lightblue,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(angle_array, CON_POST_angle_vars_mean(:,4), 'Color','b','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            
            % +/- SD
            plot(angle_array, STR_PRE_angle_vars_mean(:,4)+STR_PRE_angle_vars_SD(:,4), 'Color',col_lightred,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, STR_PRE_angle_vars_mean(:,4)-STR_PRE_angle_vars_SD(:,4), 'Color',col_lightred,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, STR_POST_angle_vars_mean(:,4)+STR_POST_angle_vars_SD(:,4), 'Color','r','LineStyle',':','LineWidth',0.25)
            plot(angle_array, STR_POST_angle_vars_mean(:,4)-STR_POST_angle_vars_SD(:,4), 'Color','r','LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_PRE_angle_vars_mean(:,4)+CON_PRE_angle_vars_SD(:,4), 'Color',col_lightblue,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_PRE_angle_vars_mean(:,4)-CON_PRE_angle_vars_SD(:,4), 'Color',col_lightblue,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_POST_angle_vars_mean(:,4)+CON_POST_angle_vars_SD(:,4), 'Color','b','LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_POST_angle_vars_mean(:,4)-CON_POST_angle_vars_SD(:,4), 'Color','b','LineStyle',':','LineWidth',0.25)
            
            % peaks
            plot(STR_PRE_data_mean(9), STR_PRE_data_mean(4), 'Color',col_lightred, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(STR_POST_data_mean(9), STR_POST_data_mean(4), 'Color','r', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(CON_PRE_data_mean(9), CON_PRE_data_mean(4), 'Color',col_lightblue, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(CON_POST_data_mean(9), CON_POST_data_mean(4), 'Color','b', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            
            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('STR PRE','STR POST','CON PRE','CON POST', 'Location','Northeast')
            title('Isokinetic plantar flexion, 60°/s');
            axis([-10 30 -Inf Inf]) %VAR
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
            
            

            plottitle = horzcat('GRP ISOKINETIC torque-angle, 90°-s, plantar flexion');
            figure('Name',plottitle)
            hold on
            
            % mean lines
            plot(angle_array, STR_PRE_angle_vars_mean(:,5), 'Color',col_lightred,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(angle_array, STR_POST_angle_vars_mean(:,5), 'Color','r','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(angle_array, CON_PRE_angle_vars_mean(:,5), 'Color',col_lightblue,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(angle_array, CON_POST_angle_vars_mean(:,5), 'Color','b','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            
            % +/- SD
            plot(angle_array, STR_PRE_angle_vars_mean(:,5)+STR_PRE_angle_vars_SD(:,5), 'Color',col_lightred,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, STR_PRE_angle_vars_mean(:,5)-STR_PRE_angle_vars_SD(:,5), 'Color',col_lightred,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, STR_POST_angle_vars_mean(:,5)+STR_POST_angle_vars_SD(:,5), 'Color','r','LineStyle',':','LineWidth',0.25)
            plot(angle_array, STR_POST_angle_vars_mean(:,5)-STR_POST_angle_vars_SD(:,5), 'Color','r','LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_PRE_angle_vars_mean(:,5)+CON_PRE_angle_vars_SD(:,5), 'Color',col_lightblue,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_PRE_angle_vars_mean(:,5)-CON_PRE_angle_vars_SD(:,5), 'Color',col_lightblue,'LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_POST_angle_vars_mean(:,5)+CON_POST_angle_vars_SD(:,5), 'Color','b','LineStyle',':','LineWidth',0.25)
            plot(angle_array, CON_POST_angle_vars_mean(:,5)-CON_POST_angle_vars_SD(:,5), 'Color','b','LineStyle',':','LineWidth',0.25)
            
            % peaks
            plot(STR_PRE_data_mean(10), STR_PRE_data_mean(5), 'Color',col_lightred, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(STR_POST_data_mean(10), STR_POST_data_mean(5), 'Color','r', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(CON_PRE_data_mean(10), CON_PRE_data_mean(5), 'Color',col_lightblue, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(CON_POST_data_mean(10), CON_POST_data_mean(5), 'Color','b', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            
            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('STR PRE','STR POST','CON PRE','CON POST', 'Location','Northeast')
            title('Isokinetic plantar flexion, 90°/s');
            axis([-10 30 -Inf Inf]) %VAR
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
            
        
        
            % same as above, in SUBPLOTS  /////////////////////////////////////////////////////////////////////////
            % not including +/- SD lines here, for legibility
            plottitle = horzcat('GRP ISOKINETIC torque-angle, plantar flexion');
            figure('Name',plottitle)
            
            subplot(2,2,1); % 30 deg/s
            hold on
            plot(angle_array, STR_PRE_angle_vars_mean(:,2), 'Color',col_lightred,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(angle_array, STR_POST_angle_vars_mean(:,2), 'Color','r','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(angle_array, CON_PRE_angle_vars_mean(:,2), 'Color',col_lightblue,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(angle_array, CON_POST_angle_vars_mean(:,2), 'Color','b','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            plot(STR_PRE_data_mean(7), STR_PRE_data_mean(2), 'Color',col_lightred, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(STR_POST_data_mean(7), STR_POST_data_mean(2), 'Color','r', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(CON_PRE_data_mean(7), CON_PRE_data_mean(2), 'Color',col_lightblue, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(CON_POST_data_mean(7), CON_POST_data_mean(2), 'Color','b', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('STR PRE','STR POST','CON PRE','CON POST', 'Location','Northeast')
            title('Isokinetic plantar flexion, 30°/s');
            axis([-10 30 30 150]) %VAR
            
            subplot(2,2,2); % 45 deg/s
            hold on
            plot(angle_array, STR_PRE_angle_vars_mean(:,3), 'Color',col_lightred,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(angle_array, STR_POST_angle_vars_mean(:,3), 'Color','r','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(angle_array, CON_PRE_angle_vars_mean(:,3), 'Color',col_lightblue,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(angle_array, CON_POST_angle_vars_mean(:,3), 'Color','b','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            plot(STR_PRE_data_mean(8), STR_PRE_data_mean(3), 'Color',col_lightred, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(STR_POST_data_mean(8), STR_POST_data_mean(3), 'Color','r', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(CON_PRE_data_mean(8), CON_PRE_data_mean(3), 'Color',col_lightblue, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(CON_POST_data_mean(8), CON_POST_data_mean(3), 'Color','b', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('STR PRE','STR POST','CON PRE','CON POST', 'Location','Northeast')
            title('Isokinetic plantar flexion, 45°/s');
            axis([-10 30 30 150]) %VAR
            
            subplot(2,2,3); % 60 deg/s
            hold on
            plot(angle_array, STR_PRE_angle_vars_mean(:,4), 'Color',col_lightred,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(angle_array, STR_POST_angle_vars_mean(:,4), 'Color','r','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(angle_array, CON_PRE_angle_vars_mean(:,4), 'Color',col_lightblue,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(angle_array, CON_POST_angle_vars_mean(:,4), 'Color','b','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            plot(STR_PRE_data_mean(9), STR_PRE_data_mean(4), 'Color',col_lightred, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(STR_POST_data_mean(9), STR_POST_data_mean(4), 'Color','r', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(CON_PRE_data_mean(9), CON_PRE_data_mean(4), 'Color',col_lightblue, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(CON_POST_data_mean(9), CON_POST_data_mean(4), 'Color','b', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('STR PRE','STR POST','CON PRE','CON POST', 'Location','Northeast')
            title('Isokinetic plantar flexion, 60°/s');
            axis([-10 30 30 150]) %VAR
            
            subplot(2,2,4); % 90 deg/s
            hold on
            plot(angle_array, STR_PRE_angle_vars_mean(:,5), 'Color',col_lightred,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(angle_array, STR_POST_angle_vars_mean(:,5), 'Color','r','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(angle_array, CON_PRE_angle_vars_mean(:,5), 'Color',col_lightblue,'LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(angle_array, CON_POST_angle_vars_mean(:,5), 'Color','b','LineStyle','-') %, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            plot(STR_PRE_data_mean(10), STR_PRE_data_mean(5), 'Color',col_lightred, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightred, 'MarkerEdgeColor',col_lightred)
            plot(STR_POST_data_mean(10), STR_POST_data_mean(5), 'Color','r', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
            plot(CON_PRE_data_mean(10), CON_PRE_data_mean(5), 'Color',col_lightblue, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor',col_lightblue, 'MarkerEdgeColor',col_lightblue)
            plot(CON_POST_data_mean(10), CON_POST_data_mean(5), 'Color','b', 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
            xlabel('Angle (°)')
            ylabel('Torque (Nm)')
            legend('STR PRE','STR POST','CON PRE','CON POST', 'Location','Northeast')
            title('Isokinetic plantar flexion, 90°/s');
            axis([-10 30 30 150]) %VAR
            
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle,'.jpg'),'-dpng')
        end


        
    % isokinetic trials, plot per group /////////////////////////////////////////////////////////////////////////

        if STR_PRE_count > 1 && STR_POST_count >1 && plot_check
            plottitle = horzcat('GRP ISOKINETIC torque-angle, plantar flexion, STRETCHING SIDE, PRE-POST');
            figure('Name',plottitle)
            hold on
            
            % torque x4 PRE
            plot(angle_array, STR_PRE_angle_vars_mean(:,2),'m','LineWidth',1,'LineStyle',':')
            plot(angle_array, STR_PRE_angle_vars_mean(:,3),'g','LineWidth',1,'LineStyle',':')
            plot(angle_array, STR_PRE_angle_vars_mean(:,4),'b','LineWidth',1,'LineStyle',':')
            plot(angle_array, STR_PRE_angle_vars_mean(:,5),'r','LineWidth',1,'LineStyle',':')
            
            % torque x4 POST
            plot(angle_array, STR_POST_angle_vars_mean(:,2),'m','LineWidth',1)
            plot(angle_array, STR_POST_angle_vars_mean(:,3),'g','LineWidth',1)
            plot(angle_array, STR_POST_angle_vars_mean(:,4),'b','LineWidth',1)
            plot(angle_array, STR_POST_angle_vars_mean(:,5),'r','LineWidth',1)

            % peak torque, angle of peak torque x4 PRE
            plot(STR_PRE_data_mean(7), STR_PRE_data_mean(2), '*m')
            plot(STR_PRE_data_mean(8), STR_PRE_data_mean(3), '*g')
            plot(STR_PRE_data_mean(9), STR_PRE_data_mean(4), '*b')
            plot(STR_PRE_data_mean(10), STR_PRE_data_mean(5), '*r')
            
            % peak torque, angle of peak torque x4 POST 
            plot(STR_POST_data_mean(7), STR_POST_data_mean(2), '*m')
            plot(STR_POST_data_mean(8), STR_POST_data_mean(3), '*g')
            plot(STR_POST_data_mean(9), STR_POST_data_mean(4), '*b')
            plot(STR_POST_data_mean(10), STR_POST_data_mean(5), '*r')
            
            % SD for peak torque, angle of peak torque x4 PRE
            errorbar(STR_PRE_data_mean(7), STR_PRE_data_mean(2), STR_PRE_data_SD(2), 'm')
            errorbar(STR_PRE_data_mean(8), STR_PRE_data_mean(3), STR_PRE_data_SD(3), 'g')
            errorbar(STR_PRE_data_mean(9), STR_PRE_data_mean(4), STR_PRE_data_SD(4), 'b')
            errorbar(STR_PRE_data_mean(10), STR_PRE_data_mean(5), STR_PRE_data_SD(5), 'r')
            herrorbar(STR_PRE_data_mean(7), STR_PRE_data_mean(2), STR_PRE_data_SD(7), 'm')
            herrorbar(STR_PRE_data_mean(8), STR_PRE_data_mean(3), STR_PRE_data_SD(8), 'g')
            herrorbar(STR_PRE_data_mean(9), STR_PRE_data_mean(4), STR_PRE_data_SD(9), 'b')
            herrorbar(STR_PRE_data_mean(10), STR_PRE_data_mean(5), STR_PRE_data_SD(10), 'r')
            
            % SD for peak torque, angle of peak torque x4 POST
            errorbar(STR_POST_data_mean(7), STR_POST_data_mean(2), STR_POST_data_SD(2), 'm')
            errorbar(STR_POST_data_mean(8), STR_POST_data_mean(3), STR_POST_data_SD(3), 'g')
            errorbar(STR_POST_data_mean(9), STR_POST_data_mean(4), STR_POST_data_SD(4), 'b')
            errorbar(STR_POST_data_mean(10), STR_POST_data_mean(5), STR_POST_data_SD(5), 'r')
            herrorbar(STR_POST_data_mean(7), STR_POST_data_mean(2), STR_POST_data_SD(7), 'm')
            herrorbar(STR_POST_data_mean(8), STR_POST_data_mean(3), STR_POST_data_SD(8), 'g')
            herrorbar(STR_POST_data_mean(9), STR_POST_data_mean(4), STR_POST_data_SD(9), 'b')
            herrorbar(STR_POST_data_mean(10), STR_POST_data_mean(5), STR_POST_data_SD(10), 'r')
            
            axis([-10 30 -Inf Inf]) %VAR
            xlabel('Ankle angle (°)')
            ylabel('Isokinetic torque (Nm)')
            title(plottitle,'Interpreter', 'none')
            labels = {'PRE isokin PF 30' 'PRE isokin PF 45' 'PRE isokin PF 60' 'PRE isokin PF 90' 'POST isokin PF 30' 'POST isokin PF 45' 'POST isokin PF 60' 'POST isokin PF 90'};
            legend(labels, 'Location','Northeast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle),'-dpng')
        end

        
        
        if CON_PRE_count > 1 && CON_POST_count >1 && plot_check
            plottitle = horzcat('GRP ISOKINETIC torque-angle, plantar flexion, CONTROL SIDE, PRE-POST');
            figure('Name',plottitle)
            hold on
            
            % torque x4 PRE
            plot(angle_array, CON_PRE_angle_vars_mean(:,2),'m','LineWidth',1,'LineStyle',':')
            plot(angle_array, CON_PRE_angle_vars_mean(:,3),'g','LineWidth',1,'LineStyle',':')
            plot(angle_array, CON_PRE_angle_vars_mean(:,4),'b','LineWidth',1,'LineStyle',':')
            plot(angle_array, CON_PRE_angle_vars_mean(:,5),'r','LineWidth',1,'LineStyle',':')
            
            % torque x4 POST
            plot(angle_array, CON_POST_angle_vars_mean(:,2),'m','LineWidth',1)
            plot(angle_array, CON_POST_angle_vars_mean(:,3),'g','LineWidth',1)
            plot(angle_array, CON_POST_angle_vars_mean(:,4),'b','LineWidth',1)
            plot(angle_array, CON_POST_angle_vars_mean(:,5),'r','LineWidth',1)

            % peak torque, angle of peak torque x4 PRE
            plot(CON_PRE_data_mean(7), CON_PRE_data_mean(2), '*m')
            plot(CON_PRE_data_mean(8), CON_PRE_data_mean(3), '*g')
            plot(CON_PRE_data_mean(9), CON_PRE_data_mean(4), '*b')
            plot(CON_PRE_data_mean(10), CON_PRE_data_mean(5), '*r')
            
            % peak torque, angle of peak torque x4 POST 
            plot(CON_POST_data_mean(7), CON_POST_data_mean(2), '*m')
            plot(CON_POST_data_mean(8), CON_POST_data_mean(3), '*g')
            plot(CON_POST_data_mean(9), CON_POST_data_mean(4), '*b')
            plot(CON_POST_data_mean(10), CON_POST_data_mean(5), '*r')
            
            % SD for peak torque, angle of peak torque x4 PRE
            errorbar(CON_PRE_data_mean(7), CON_PRE_data_mean(2), CON_PRE_data_SD(2), 'm')
            errorbar(CON_PRE_data_mean(8), CON_PRE_data_mean(3), CON_PRE_data_SD(3), 'g')
            errorbar(CON_PRE_data_mean(9), CON_PRE_data_mean(4), CON_PRE_data_SD(4), 'b')
            errorbar(CON_PRE_data_mean(10), CON_PRE_data_mean(5), CON_PRE_data_SD(5), 'r')
            herrorbar(CON_PRE_data_mean(7), CON_PRE_data_mean(2), CON_PRE_data_SD(7), 'm')
            herrorbar(CON_PRE_data_mean(8), CON_PRE_data_mean(3), CON_PRE_data_SD(8), 'g')
            herrorbar(CON_PRE_data_mean(9), CON_PRE_data_mean(4), CON_PRE_data_SD(9), 'b')
            herrorbar(CON_PRE_data_mean(10), CON_PRE_data_mean(5), CON_PRE_data_SD(10), 'r')
            
            % SD for peak torque, angle of peak torque x4 POST
            errorbar(CON_POST_data_mean(7), CON_POST_data_mean(2), CON_POST_data_SD(2), 'm')
            errorbar(CON_POST_data_mean(8), CON_POST_data_mean(3), CON_POST_data_SD(3), 'g')
            errorbar(CON_POST_data_mean(9), CON_POST_data_mean(4), CON_POST_data_SD(4), 'b')
            errorbar(CON_POST_data_mean(10), CON_POST_data_mean(5), CON_POST_data_SD(5), 'r')
            herrorbar(CON_POST_data_mean(7), CON_POST_data_mean(2), CON_POST_data_SD(7), 'm')
            herrorbar(CON_POST_data_mean(8), CON_POST_data_mean(3), CON_POST_data_SD(8), 'g')
            herrorbar(CON_POST_data_mean(9), CON_POST_data_mean(4), CON_POST_data_SD(9), 'b')
            herrorbar(CON_POST_data_mean(10), CON_POST_data_mean(5), CON_POST_data_SD(10), 'r')
            
            axis([-10 30 -Inf Inf]) %VAR
            xlabel('Ankle angle (°)')
            ylabel('Isokinetic torque (Nm)')
            title(plottitle,'Interpreter', 'none')
            labels = {'PRE isokin PF 30' 'PRE isokin PF 45' 'PRE isokin PF 60' 'PRE isokin PF 90' 'POST isokin PF 30' 'POST isokin PF 45' 'POST isokin PF 60' 'POST isokin PF 90'};
            legend(labels, 'Location','Northeast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle),'-dpng')
        end

        
        if STR_PRE_count > 1 && STR_POST_count >1 && plot_check
            plottitle = horzcat('GRP ISOKINETIC torque-angle, dorsiflexion, STRETCHING SIDE, PRE-POST');
            figure('Name',plottitle)
            hold on
            
            plot(angle_array, STR_PRE_angle_vars_mean(:,1),'r','LineWidth',1,'LineStyle',':')
            plot(angle_array, STR_POST_angle_vars_mean(:,1),'r','LineWidth',1)
            plot(STR_PRE_data_mean(6), STR_PRE_data_mean(1), '*r')
            plot(STR_POST_data_mean(6), STR_POST_data_mean(1), '*r')
            errorbar(STR_PRE_data_mean(6), STR_PRE_data_mean(1), STR_PRE_data_SD(1), 'r')
            herrorbar(STR_PRE_data_mean(6), STR_PRE_data_mean(1), STR_PRE_data_SD(6), 'r')
            errorbar(STR_POST_data_mean(6), STR_POST_data_mean(1), STR_POST_data_SD(1), 'r')
            herrorbar(STR_POST_data_mean(6), STR_POST_data_mean(1), STR_POST_data_SD(6), 'r')
            
            % reverse X axis: start at 30 deg, dorsiflex into -10 deg
            ax = gca;
            set(ax, 'xdir','reverse')
            
            axis([-10 30 -Inf Inf]) %VAR
            xlabel('Ankle angle (°)')
            ylabel('Isokinetic torque (Nm)')
            title(plottitle,'Interpreter', 'none')
            labels = {'PRE isokin DF 30' 'POST isokin DF 30'};
            legend(labels, 'Location','Northeast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle),'-dpng')
        end
        
        
        
        if CON_PRE_count > 1 && CON_POST_count >1 && plot_check
            plottitle = horzcat('GRP ISOKINETIC torque-angle, dorsiflexion, CONTROL SIDE, PRE-POST');
            figure('Name',plottitle)
            hold on
            
            plot(angle_array, CON_PRE_angle_vars_mean(:,1),'r','LineWidth',1,'LineStyle',':')
            plot(angle_array, CON_POST_angle_vars_mean(:,1),'r','LineWidth',1)
            plot(CON_PRE_data_mean(6), CON_PRE_data_mean(1), '*r')
            plot(CON_POST_data_mean(6), CON_POST_data_mean(1), '*r')
            errorbar(CON_PRE_data_mean(6), CON_PRE_data_mean(1), CON_PRE_data_SD(1), 'r')
            herrorbar(CON_PRE_data_mean(6), CON_PRE_data_mean(1), CON_PRE_data_SD(6), 'r')
            errorbar(CON_POST_data_mean(6), CON_POST_data_mean(1), CON_POST_data_SD(1), 'r')
            herrorbar(CON_POST_data_mean(6), CON_POST_data_mean(1), CON_POST_data_SD(6), 'r')
            
            % reverse X axis: start at 30 deg, dorsiflex into -10 deg
            ax = gca;
            set(ax, 'xdir','reverse')
            
            axis([-10 30 -Inf Inf]) %VAR
            xlabel('Ankle angle (°)')
            ylabel('Isokinetic torque (Nm)')
            title(plottitle,'Interpreter', 'none')
            labels = {'PRE isokin DF 30' 'POST isokin DF 30'};
            legend(labels, 'Location','Northeast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
            print(horzcat('data_plots/',plottitle),'-dpng')
        end
        
        
        
        
    end
    % end of IF BD vs INTERVENTION
    
    
end