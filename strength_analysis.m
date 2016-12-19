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
    plot_conversion = 0;
    plot_norm = 0;



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
%    global dm_CPM_calc_NX dm_CPM_sol_NX
%    global dm_leg_length
    global filepath
    dm_filename = 'data_strength/datamaster_strength.tsv';
    dm_columns = 22; % number of data columns entered per subject % PROJECTSPECIFIC
    linestotal = read_datamaster_strength(dm_filename,dm_columns);
    
    
    
    

    %%% preallocate output arrays
    
    % common arrays for numbers across all subjects:
    all_strength_output = zeros(ceil(linestotal),30); 
    all_strength_output_txt = cell(ceil(linestotal),4);

    % BD-SPECIFIC
    BD_count = 0;
    CON_count = 0;
    
    BD_angle_vars{ceil(linestotal)} = [];
    CON_angle_vars{ceil(linestotal)} = [];
    BD_data = zeros(ceil(linestotal),30);
    CON_data = zeros(ceil(linestotal),30);
    
    % below variables are also deleted + preallocated within the loop
    isokinetic_torque_angle_work(1:5,1:4) = zeros;
    isokinetic_arrays{5} = [];

    


    
    

    
    %%%%%%%%%%%%%%%% LOOP through all lines in datamaster file (except header line)
    for line = 1:linestotal

        
        



        %%% subject/trial identifier
        subjectno = str2double(dm_subjectno{line});
        if subjectno > 100
            filepath = 'data\BD\';
            subject_id = horzcat('dancer ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
            BD_count = BD_count + 1;
            isokinetic_labels = {'isokin DF 30' 'isokin PF 45' 'isokin PF 60' 'isokin PF 90' 'isokin PF 120'};
            %isokinetic_speeds = [30 45 60 90 120];
        else
            filepath = 'data_strength\';
            subject_id = horzcat('control ', dm_subjectno{line}, ' ', dm_side{line}, ' ', dm_timepoint{line}, ' ', dm_trial{line});
            CON_count = CON_count + 1;
            isokinetic_labels = {'isokin DF 30' 'isokin PF 30' 'isokin PF 45' 'isokin PF 60' 'isokin PF 90'};
            %isokinetic_speeds = [30 30 45 60 90];
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
        [convert_norm_angle_a, convert_norm_angle_b, convert_norm_torque_a, convert_norm_torque_b, convert_norm_velocity_a, convert_norm_velocity_b, convert_norm_direction_b] = calculate_norm_constants_complete(subjectno, dm_side{line}, dm_timepoint{line}, line, 'active');
        
        
        
        
        
        
        
        %%% calculations for ISOKINETIC trials
        
        % identify files
        isokinetic_data = {dm_isokinD30{line} dm_isokinP30{line} dm_isokinP45{line} dm_isokinP60{line} dm_isokinP90{line}};
        
        % preallocate
        clear isokinetic_torque_angle_work isokinetic_arrays
        isokinetic_torque_angle_work(1:length(isokinetic_data),1:4) = zeros;
        isokinetic_arrays{length(isokinetic_data)} = [];
        
        % extract data from Noraxon files
        for i = 1:length(isokinetic_data)
            if(strcmpi(isokinetic_data{i}, 'null'))
                % allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case later
                isokinetic_torque_angle_work(i,1) = NaN;
            else 
                [isokinetic_torque_angle_work(i,1), isokinetic_torque_angle_work(i,2), isokinetic_torque_angle_work(i,3), isokinetic_torque_angle_work(i,4), isokinetic_arrays{i}] = extract_isokinetic(isokinetic_data{i}, dm_side{line}, isokinetic_labels{i}, subject_id);
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
            title(plottitle)
            legend(isokinetic_labels)
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
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
             axis([30 130 0 max(isokinetic_torque_angle_work(2:end,1))*1.1]) %VAR
%            ax = gca;
%            set(ax, 'xdir','reverse')
            xlabel('Velocity (°/s)')
            ylabel('Isokinetic peak torque (Nm)')
            title(plottitle)
            legend(isokinetic_labels(1:end),'location','SouthEast')
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
        end
        
        
        
        
        
        
        
        
        
        %%% calculations for ISOMETRIC trials
        
        % identify files
        isometric_data = {dm_isomet_P10_1{line} dm_isomet_P10_2{line} dm_isomet_D00_1{line} dm_isomet_D00_2{line} dm_isomet_D05_1{line} dm_isomet_D05_2{line} dm_isomet_D10_1{line} dm_isomet_D10_2{line} dm_isomet_D15_1{line} dm_isomet_D15_2{line}};
        isometric_labels = {'P10 trial 1' 'P10 trial 2' 'D00 trial 1' 'D00 trial 2' 'D05 trial 1' 'D05 trial 2' 'D10 trial 1' 'D10 trial 2' 'D15 trial 1' 'D15 trial 2'};

        % preallocate
        isometric_torque_angle(1:10,1:2) = zeros;
        isometric_arrays{length(isometric_data)} = [];

        % extract data from Noraxon files
        for i = 1:length(isometric_data)
            if(strcmpi(isometric_data{i}, 'null'))
                % allow for the possibility of discarded trials (null). In that case, just make empty arrays and check for that special case later
                isometric_torque_angle(i,1) = NaN;
                isometric_torque_angle(i,2) = NaN;
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
            ylabel('Isometric torque (Nm)')
            title(plottitle)
            legend(isometric_labels)
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
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
            ax = gca;
            set(ax, 'xdir','reverse')
            xlabel('Ankle angle (°)')
            ylabel('Isometric peak torque (Nm)')
            title(plottitle)
            % legend(isometric_labels)
            saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
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
        
        
        
        
        
        
        

        %%% collect OUTPUT DATA FOR AVERAGE PLOTS
        
        % arrays for averaging isokinetic curves across subjects
        if subjectno > 100
            BD_angle_vars{BD_count} = isokinetic_arrays;
            BD_data(BD_count,1) = isometric_10;
            BD_data(BD_count,2) = isometric_0;
            BD_data(BD_count,3) = isometric_n5;
            BD_data(BD_count,4) = isometric_n10;
            BD_data(BD_count,5) = isometric_n15;
            BD_data(BD_count,6) = isokinetic_torque_angle_work(1,1);
            BD_data(BD_count,7) = isokinetic_torque_angle_work(2,1);
            BD_data(BD_count,8) = isokinetic_torque_angle_work(3,1);
            BD_data(BD_count,9) = isokinetic_torque_angle_work(4,1);
            BD_data(BD_count,10) = isokinetic_torque_angle_work(5,1);
            BD_data(BD_count,11) = isokinetic_torque_angle_work(1,2);
            BD_data(BD_count,12) = isokinetic_torque_angle_work(2,2);
            BD_data(BD_count,13) = isokinetic_torque_angle_work(3,2);
            BD_data(BD_count,14) = isokinetic_torque_angle_work(4,2);
            BD_data(BD_count,15) = isokinetic_torque_angle_work(5,2);
            BD_data(BD_count,16) = isokinetic_torque_angle_work(1,4);
            BD_data(BD_count,17) = isokinetic_torque_angle_work(2,4);
            BD_data(BD_count,18) = isokinetic_torque_angle_work(3,4);
            BD_data(BD_count,19) = isokinetic_torque_angle_work(4,4);
            BD_data(BD_count,20) = isokinetic_torque_angle_work(5,4);
        else
            CON_angle_vars{CON_count} = isokinetic_arrays;
            CON_data(CON_count,1) = isometric_10;
            CON_data(CON_count,2) = isometric_0;
            CON_data(CON_count,3) = isometric_n5;
            CON_data(CON_count,4) = isometric_n10;
            CON_data(CON_count,5) = isometric_n15;
            CON_data(CON_count,6) = isokinetic_torque_angle_work(1,1);
            CON_data(CON_count,7) = isokinetic_torque_angle_work(2,1);
            CON_data(CON_count,8) = isokinetic_torque_angle_work(3,1);
            CON_data(CON_count,9) = isokinetic_torque_angle_work(4,1);
            CON_data(CON_count,10) = isokinetic_torque_angle_work(5,1);
            CON_data(CON_count,11) = isokinetic_torque_angle_work(1,2);
            CON_data(CON_count,12) = isokinetic_torque_angle_work(2,2);
            CON_data(CON_count,13) = isokinetic_torque_angle_work(3,2);
            CON_data(CON_count,14) = isokinetic_torque_angle_work(4,2);
            CON_data(CON_count,15) = isokinetic_torque_angle_work(5,2);
            CON_data(CON_count,16) = isokinetic_torque_angle_work(1,4);
            CON_data(CON_count,17) = isokinetic_torque_angle_work(2,4);
            CON_data(CON_count,18) = isokinetic_torque_angle_work(3,4);
            CON_data(CON_count,19) = isokinetic_torque_angle_work(4,4);
            CON_data(CON_count,20) = isokinetic_torque_angle_work(5,4);
        end


        
        
        
        
        






        %%% OUTPUT final individual data to file
        
        % add data to a common array for all subjects    
        i = 1;
        
        % txt trial ID
        all_strength_output_txt(line,1) = dm_subjectno(line);
        all_strength_output_txt(line,2) = dm_timepoint(line);
        all_strength_output_txt(line,3) = dm_side(line);
        all_strength_output_txt(line,4) = dm_trial(line);
        
        % ----- isokinetic:
        % comes from extract_isokinetic: torque_max, torque_max_angle, torque_max_velocity, work_max, array_output

        % peak torque
        all_strength_output(line,i) = isokinetic_torque_angle_work(1,1);
        i = i+1;
        all_strength_output(line,i) = isokinetic_torque_angle_work(2,1);
        i = i+1;
        all_strength_output(line,i) = isokinetic_torque_angle_work(3,1);
        i = i+1;
        all_strength_output(line,i) = isokinetic_torque_angle_work(4,1);
        i = i+1;
        all_strength_output(line,i) = isokinetic_torque_angle_work(5,1);
        i = i+1;
        % angle of peak torque
        all_strength_output(line,i) = isokinetic_torque_angle_work(1,2);
        i = i+1;
        all_strength_output(line,i) = isokinetic_torque_angle_work(2,2);
        i = i+1;
        all_strength_output(line,i) = isokinetic_torque_angle_work(3,2);
        i = i+1;
        all_strength_output(line,i) = isokinetic_torque_angle_work(4,2);
        i = i+1;
        all_strength_output(line,i) = isokinetic_torque_angle_work(5,2);
        i = i+1;
        % work
        all_strength_output(line,i) = isokinetic_torque_angle_work(1,4);
        i = i+1;
        all_strength_output(line,i) = isokinetic_torque_angle_work(2,4);
        i = i+1;
        all_strength_output(line,i) = isokinetic_torque_angle_work(3,4);
        i = i+1;
        all_strength_output(line,i) = isokinetic_torque_angle_work(4,4);
        i = i+1;
        all_strength_output(line,i) = isokinetic_torque_angle_work(5,4);
        i = i+1;
        
        % ----- isometric:
        
        all_strength_output(line,i) = isometric_10;
        i = i+1;
        all_strength_output(line,i) = isometric_0;
        i = i+1;
        all_strength_output(line,i) = isometric_n5;
        i = i+1;
        all_strength_output(line,i) = isometric_n10;
        i = i+1;
        all_strength_output(line,i) = isometric_n15;

        all_strength_output_head = {'Subject', 'Time', 'Side', 'Trial', ...
            'Peak torque DF 30 (Nm)', 'Peak torque PF 30/45 (Nm)', 'Peak torque PF 45/60 (Nm)', 'Peak torque PF 60/90 (Nm)', 'Peak torque PF 90/120 (Nm)' ...
            'Angle of PT DF 30 (°)', 'Angle of PT PF 30/ (°)', 'Angle of PT PF 45/ (°)', 'Angle of PT PF 60/ (°)', 'Angle of PT PF 90/ (°)' ...
            'Work DF 30 (J)', 'Work PF 30/ (J)', 'Work PF 45/ (J)', 'Work PF 60/ (J)', 'Work PF 90/ (J)' ...
            'Isometric PF torque +10° (Nm)', 'Isometric PF torque 0° (Nm)', 'Isometric PF torque -5° (Nm)', 'Isometric PF torque -10° (Nm)', 'Isometric PF torque -15° (Nm)' ...
            }; % PROJECTSPECIFIC


    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%% LOOP FINISHED
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATIONS ACROSS SUBJECTS - NUMBERS
    
    %%%  mean and stdav of each subject's isometric torque, isokinetic torque/angle/work
    
    % average
    BD_data_mean = mean(BD_data,'omitnan');
    CON_data_mean = mean(CON_data,'omitnan');
    BD_data_SD = std(BD_data,'omitnan');
    CON_data_SD = std(CON_data,'omitnan');
    
    
    if plot_check
        plottitle = horzcat('GRP ISOMETRIC torque-angle, plantar flexion');
        figure('Name',plottitle)
        hold on
        plot_angles = [10, 0, -5 -10, -15];
        plot(plot_angles, BD_data_mean(1:5), 'Marker','o','MarkerSize',6,'MarkerFaceColor','auto','Color','b')
        plot(plot_angles, CON_data_mean(1:5), 'Marker','o','MarkerSize',6,'MarkerFaceColor','auto','Color','r')
        errorbar(plot_angles, BD_data_mean(1:5), BD_data_SD(1:5),'Color','b')
        errorbar(plot_angles, CON_data_mean(1:5), CON_data_SD(1:5),'Color','r')
        ax = gca;
        set(ax, 'xdir','reverse')
        axis([-20 15 min([min(BD_data_mean(1:5)) min(CON_data_mean(1:5))]) max([max(BD_data_mean(1:5)) max(CON_data_mean(1:5))])]) %VAR
        xlabel('Ankle angle (°)')
        ylabel('Isometric torque (Nm)')
        title(plottitle)
        legend('BD', 'CON', 'Location','Northeast')
        saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
     
     
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%% OUTPUT KEY VARIABLES FOR ALL SUBJECTS TO FILE
    % write xls
    if ispc
        filename_output = strcat('data_output/all_strength_output_', datestr(now, 'yyyy-mm-dd HH-MM'));
        
        xlswrite(filename_output, all_strength_output_head, 1, 'A1')
        xlswrite(filename_output, all_strength_output_txt, 1, 'A2')
        xlswrite(filename_output, all_strength_output, 1, 'E2')
    else
        filename_output = strcat('data_output/all_strength_output_', datestr(now, 'yyyy-mm-dd HH-MM'), '.csv');
        csvwrite(filename_output, all_strength_output)
    end

    
    
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATIONS ACROSS SUBJECTS - ARRAYS
    
    %%% average isokinetic arrays
    
    % find common angles across all subjects
    if BD_count > 0
        BD_angle_max(1:5,1:BD_count) = zeros;
        BD_angle_min(1:5,1:BD_count) = zeros;
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
        CON_angle_max(1:5,1:CON_count) = zeros;
        CON_angle_min(1:5,1:CON_count) = zeros;
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
    angle_array = ceil(10*angle_min)/10:0.1:floor(10*angle_max)/10';
    
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
        BD_angle_vars_mean(1:length(angle_array),1:5) = zeros;
        BD_angle_vars_SD(1:length(angle_array),1:5) = zeros;
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
        CON_angle_vars_mean(1:length(angle_array),1:5) = zeros;
        CON_angle_vars_SD(1:length(angle_array),1:5) = zeros;
        for n=1:5
            if ~isempty(CON_angle_vars_common{n})
                CON_angle_vars_mean(:,n) = nanmean(CON_angle_vars_common{n}, 2);
                CON_angle_vars_SD(:,n) = nanstd(CON_angle_vars_common{n},1,2);
                % else - arrays will remain zero
            end
        end
    end
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%% PLOT GROUP FIGURES
    
    
    
    % isokinetic trials

    % MMM TODO GOON: add mea+SD peak torque
    
     % trial 1 = CON ONLY: DF 30°/s
     % trial 2 = CON: PF 30°/s, BD: PF 45°/s
     % trial 5 = CON: PF 90°/s, BD: PF 120°/s

    if BD_count > 1 && CON_count > 1 && plot_check
             plottitle = horzcat('GRP isokinetic torque-angle, plantar flexion');
             figure('Name',plottitle)
             
             subplot(2,2,1);
             hold on
             plot(angle_array, CON_angle_vars_mean(:,2),'r','LineWidth',1.5,'LineStyle','--')
             plot(angle_array, CON_angle_vars_mean(:,2)+CON_angle_vars_SD(:,2),'r','LineWidth',0.25,'LineStyle','--')
             plot(angle_array, CON_angle_vars_mean(:,2)-CON_angle_vars_SD(:,2),'r','LineWidth',0.25,'LineStyle','--')
             xlabel('Angle (°)')
             ylabel('Torque (Nm)')
             legend('CON', 'Location','Northeast')
             title('Isokinetic plantar flexion, 30°/s');
             axis([-10 30 30 150]) %VAR

             subplot(2,2,2);
             hold on
             plot(angle_array, BD_angle_vars_mean(:,2),'m','LineWidth',1)
             plot(angle_array, CON_angle_vars_mean(:,3),'m','LineWidth',1.5,'LineStyle','--')
             plot(angle_array, BD_angle_vars_mean(:,2)+BD_angle_vars_SD(:,2),'m','LineWidth',0.25,'LineStyle',':')
             plot(angle_array, CON_angle_vars_mean(:,3)+CON_angle_vars_SD(:,3),'m','LineWidth',0.25,'LineStyle','--')
             plot(angle_array, BD_angle_vars_mean(:,2)-BD_angle_vars_SD(:,2),'m','LineWidth',0.25,'LineStyle',':')
             plot(angle_array, CON_angle_vars_mean(:,3)-CON_angle_vars_SD(:,3),'m','LineWidth',0.25,'LineStyle','--')
             xlabel('Angle (°)')
             ylabel('Torque (Nm)')
             legend('BD', 'CON', 'Location','Northeast')
             title('Isokinetic plantar flexion, 45°/s');
             axis([-10 30 30 150]) %VAR
             
             subplot(2,2,3);
             hold on
             plot(angle_array, BD_angle_vars_mean(:,3),'g','LineWidth',1)
             plot(angle_array, CON_angle_vars_mean(:,4),'g','LineWidth',1.5,'LineStyle','--')
             plot(angle_array, BD_angle_vars_mean(:,3)+BD_angle_vars_SD(:,3),'g','LineWidth',0.25,'LineStyle',':')
             plot(angle_array, CON_angle_vars_mean(:,4)+CON_angle_vars_SD(:,4),'g','LineWidth',0.25,'LineStyle','--')
             plot(angle_array, BD_angle_vars_mean(:,3)-BD_angle_vars_SD(:,3),'g','LineWidth',0.25,'LineStyle',':')
             plot(angle_array, CON_angle_vars_mean(:,4)-CON_angle_vars_SD(:,4),'g','LineWidth',0.25,'LineStyle','--')
             xlabel('Angle (°)')
             ylabel('Torque (Nm)')
             legend('BD', 'CON', 'Location','Northeast')
             title('Isokinetic plantar flexion, 60°/s');
             axis([-10 30 30 150]) %VAR
             
             subplot(2,2,4);
             hold on
             plot(angle_array, BD_angle_vars_mean(:,4),'b','LineWidth',1)
             plot(angle_array, CON_angle_vars_mean(:,5),'b','LineWidth',1.5,'LineStyle','--')
             plot(angle_array, BD_angle_vars_mean(:,4)+BD_angle_vars_SD(:,4),'b','LineWidth',0.25,'LineStyle',':')
             plot(angle_array, CON_angle_vars_mean(:,5)+CON_angle_vars_SD(:,5),'b','LineWidth',0.25,'LineStyle','--')
             plot(angle_array, BD_angle_vars_mean(:,4)-BD_angle_vars_SD(:,4),'b','LineWidth',0.25,'LineStyle',':')
             plot(angle_array, CON_angle_vars_mean(:,5)-CON_angle_vars_SD(:,5),'b','LineWidth',0.25,'LineStyle','--')
             xlabel('Angle (°)')
             ylabel('Torque (Nm)')
             legend('BD', 'CON', 'Location','Northeast')
             title('Isokinetic plantar flexion, 90°/s');
             axis([-10 30 30 150]) %VAR
             
             saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end

    if BD_count > 1 && plot_check
             plottitle = horzcat('GRP isokinetic torque-angle, plantar flexion, dancers');
             figure('Name',plottitle)
             hold on
             plot(angle_array, BD_angle_vars_mean(:,2),'m','LineWidth',1)
             plot(angle_array, BD_angle_vars_mean(:,3),'g','LineWidth',1)
             plot(angle_array, BD_angle_vars_mean(:,4),'b','LineWidth',1)
             plot(angle_array, BD_angle_vars_mean(:,5),'k','LineWidth',1)

             plot(angle_array, BD_angle_vars_mean(:,2)+BD_angle_vars_SD(:,2),'m','LineWidth',0.25,'LineStyle',':')
             plot(angle_array, BD_angle_vars_mean(:,2)-BD_angle_vars_SD(:,2),'m','LineWidth',0.25,'LineStyle',':')
             plot(angle_array, BD_angle_vars_mean(:,3)+BD_angle_vars_SD(:,3),'g','LineWidth',0.25,'LineStyle',':')
             plot(angle_array, BD_angle_vars_mean(:,3)-BD_angle_vars_SD(:,3),'g','LineWidth',0.25,'LineStyle',':')
             plot(angle_array, BD_angle_vars_mean(:,4)+BD_angle_vars_SD(:,4),'b','LineWidth',0.25,'LineStyle',':')
             plot(angle_array, BD_angle_vars_mean(:,4)-BD_angle_vars_SD(:,4),'b','LineWidth',0.25,'LineStyle',':')
             plot(angle_array, BD_angle_vars_mean(:,5)+BD_angle_vars_SD(:,5),'k','LineWidth',0.25,'LineStyle',':')
             plot(angle_array, BD_angle_vars_mean(:,5)-BD_angle_vars_SD(:,5),'k','LineWidth',0.25,'LineStyle',':')
             
             axis([-10 30 30 150]) %VAR
             xlabel('Ankle angle (°)')
             ylabel('Isokinetic torque (Nm)')
             title(plottitle)
             labels = {'isokin PF 45' 'isokin PF 60' 'isokin PF 90' 'isokin PF 120'};
             legend(labels, 'Location','Northeast')
             saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    if CON_count > 1 && plot_check
             plottitle = horzcat('GRP isokinetic torque-angle, plantar flexion, controls');
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
             
             axis([-10 30 30 150]) %VAR
             xlabel('Ankle angle (°)')
             ylabel('Isokinetic torque (Nm)')
             title(plottitle)
             labels = {'isokin PF 30' 'isokin PF 45' 'isokin PF 60' 'isokin PF 90'};
             legend(labels, 'Location','Northeast')
             saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    
    
    
    if CON_count > 1 && plot_check
             plottitle = horzcat('GRP isokinetic torque-angle, dorsiflexion');
             figure('Name',plottitle)
             hold on
             plot(angle_array, CON_angle_vars_mean(:,1),'b','LineWidth',1,'LineStyle','--')
             plot(angle_array, CON_angle_vars_mean(:,1)+CON_angle_vars_SD(:,1),'b','LineWidth',0.25,'LineStyle','--')
             plot(angle_array, CON_angle_vars_mean(:,1)-CON_angle_vars_SD(:,1),'b','LineWidth',0.25,'LineStyle','--')
             ax = gca;
             set(ax, 'xdir','reverse')
             axis([-10 30 0 20]) %VAR
             xlabel('Ankle angle (°)')
             ylabel('Isokinetic torque (Nm)')
             title(plottitle)
             legend('isokin DF 30', 'Location','Northeast')
             saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
    

    
    
end