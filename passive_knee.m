%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main file for analysis of passive knee extension in Norm
% Marie Moltubakk 16.11.2017
% 
% input argument 1 = plot selection (0 = none, 1 = group plots, 2 = ind plots)
%
% The scripts assume a certain structure for input files from Tracker and
% Noraxon. I.e. number of and order of EMG and other channels. If these
% scripts are to be used for other projects, some of the code must be
% modified. These lines are marked % PROJECTSPECIFIC
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MMM TODO:
%   change group figures to use fits, not raw data
%   cut data at peak torque (reduces rom...)?
%   output start angle?


function [] = passive_knee(input_plot)
    close all
    
    angle_trial_onset = 90; %VAR
    angle_end = 10; %VAR  - max for plots
    
    %% PLOTS - determine which plots to display
    global plot_achilles plot_norm plot_check subject_id plot_individual plot_conversion

    if input_plot >= 1 
        plot_check = 1; % LEVEL 1: group plots
    else
        plot_check = 0;
    end
    if input_plot >= 2
        plot_individual = 1; % LEVEL 1B: individual plots
        plot_norm = 1; % LEVEL 2: main checkpoint plots
    else
        plot_individual = 0;
        plot_norm = 0; % LEVEL 2: main checkpoint plots
    end
    if input_plot >= 3
        plot_conversion = 0; % TMP; % turn on/off plots for data conversion Norm
        plot_achilles = 1; % turn on/off Achilles machine plots
    else
        plot_conversion = 0; % turn on/off plots for data conversion Norm
        plot_achilles = 0; % turn on/off Achilles machine plots
    end


    %% Set constants and globals % PROJECTSPECIFIC

    % sampling frequencies
    global us_zerodispframes noraxonfreq freq_default
    us_zerodispframes = 1; % No of US frames to average as zero displacement
    noraxonfreq = 1500; % sampling frequency of noraxon data
    freq_default = 100; % output frequency for noraxon data without US video (MVC, etc)
%    angle_step_plots = 0.05; % resampled, averaged data extracted every x degrees for PLOTS
    angle_step_stats_abs = 0.5; % every x degrees
%    angle_step_stats_norm = 10;  % every x percent
    
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
    
    
    %% set AXES etc for plots
    col_lightblue = [0.6 0.8 1];
    col_lightred = [1 0.6 0.8];
%    col_orange = [1 0.75 0];
%    col_grey = [0.3 0.3 0.3];
    
    txt_gonio = 'Test angle (°)';
    axis_torque = [angle_end 100 0 Inf];


    %% Read datamaster file, to connect corresponding data files
    % Produces arrays with file names and variables per trial

    global dm_subjectno dm_timepoint dm_side dm_trial 
    global dm_ROM_trial1 dm_ROM_trial2 % used only in sub-scripts: dm_CPM_calc_NX dm_CPM_sol_NX 
    global dm_n_o_trials
    global dm_ROM_ind dm_ROM_common % dm_ROM_trial 
    global filepath filepath2
    dm_filename = 'data/datamaster_knee.tsv';
    linestotal = read_datamaster_knee(dm_filename);
        
    
    %% preallocate output arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % common arrays for all subjects:
    all_passive_output_head = {'Subject', 'Time', 'Side', 'Trial', ...
            'ROM trial (°)', 'ROM subject (L/R/PRE/POST)', 'ROM common (all subjects)', 'ROM @ max torque', 'ROM @ givenangle 50deg', 'starting angle 2trials' ...
            'torque @ trial ROM (N)', 'torque @ subject max ROM', 'torque @ common max ROM', 'torque max', 'torque @ givenangle', ...
            'passive stiffness @ trial ROM (N)', 'stiffness @ subject max ROM', 'stiffness @ common max ROM', 'stiffness @ torque max', 'stiffness @ givenangle', ...
            'stiffness index', ...
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

    STR_PRE_no(ceil(linestotal)) = zeros;
    STR_PRE_angle_vars{ceil(linestotal)} = zeros;
    STR_PRE_angle_vars_mean{ceil(linestotal)} = zeros;
%     STR_PRE_angle_vars_norm{ceil(linestotal)} = zeros;
%     STR_PRE_angle_vars_norm_indlength{ceil(linestotal)} = zeros;
%     STR_PRE_angle_vars_norm_mean{ceil(linestotal)} = zeros;

    STR_POST_no(ceil(linestotal)) = zeros;
    STR_POST_angle_vars{ceil(linestotal)} = zeros;
    STR_POST_angle_vars_mean{ceil(linestotal)} = zeros;
%     STR_POST_angle_vars_norm{ceil(linestotal)} = zeros;
%     STR_POST_angle_vars_norm_indlength{ceil(linestotal)} = zeros;
%     STR_POST_angle_vars_norm_mean{ceil(linestotal)} = zeros;

    CON_PRE_no(ceil(linestotal)) = zeros;
    CON_PRE_angle_vars{ceil(linestotal)} = zeros;
    CON_PRE_angle_vars_mean{ceil(linestotal)} = zeros;
%     CON_PRE_angle_vars_norm{ceil(linestotal)} = zeros;
%     CON_PRE_angle_vars_norm_indlength{ceil(linestotal)} = zeros;
%     CON_PRE_angle_vars_norm_mean{ceil(linestotal)} = zeros;

    CON_POST_no(ceil(linestotal)) = zeros;
    CON_POST_angle_vars{ceil(linestotal)} = zeros;
    CON_POST_angle_vars_mean{ceil(linestotal)} = zeros;
%     CON_POST_angle_vars_norm{ceil(linestotal)} = zeros;
%     CON_POST_angle_vars_norm_indlength{ceil(linestotal)} = zeros;
%     CON_POST_angle_vars_norm_mean{ceil(linestotal)} = zeros;

    
    %% LOOP through all lines in datamaster file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for line = 1:linestotal 
        
        
        %% subject/trial identifier
        trial_timepoint = strcmp(dm_timepoint{line},'POST'); % 0 = PRE, 1 = POST
        trial_leg = strcmp(dm_trial{line},'STR'); % 0 = CON, 1 = STR
        trial_n_o_trials = str2double(dm_n_o_trials{line});
        %out_ROM_trial_max = str2double(dm_ROM_trial{line});
        out_ROM_ind_max = str2double(dm_ROM_ind{line});
        out_ROM_common_max = str2double(dm_ROM_common{line});
        
        filepath = 'data\';
        filepath2 = 'stretcher_knee\';
        subject_id = horzcat('INT_', dm_subjectno{line}, '_', dm_trial{line}, '_', dm_timepoint{line}, '_', dm_side{line});
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
% if CPM (ankle trial) offset data are not appropriate for knee trial
% conversions, consider using default conversions provided by Norm:
%         if dm_side{line} == 'L'
%             convert_norm_angle_a = 1/norm_volt_per_degree;
%         else
%             convert_norm_angle_a = -1/norm_volt_per_degree;
%         end
%         convert_norm_velocity_a = 1/norm_volt_per_velocity;
%         norm_volt_per_nm_a = 0.07490071; % from file "M M M conversion volt-torque 2014-DES NEW DATA.xlsx"
%         % norm_volt_per_nm_b = 0.69283161;
%         % convert_norm_torque_a = 2 * 1/norm_volt_per_nm;
%         convert_norm_torque_a = norm_volt_per_nm_a;
% 
%         convert_norm_velocity_b = 0;
%         convert_norm_angle_b = 0;
%         convert_norm_torque_b = 0;
%         convert_norm_direction_b = 60;

    
        
        %% EXTRACT torque and angle for 2 knee extension trials
        if trial_n_o_trials == 2
            % 2 trials in first file, no second file
            [ROM_angle_1, ROM_torque_1, ROM_angle_2, ROM_torque_2] = extract_force_displ_singletrial_knee(dm_ROM_trial1{line}, dm_side{line}, 2, 'ROMdbl');
        else
            if(strcmpi(dm_ROM_trial1{line}, 'null'))
                % 1 single trial discarded
                ROM_torque_1 = zeros;
                ROM_angle_1 = zeros;
            else
                [ROM_angle_1, ROM_torque_1] = extract_force_displ_singletrial_knee(dm_ROM_trial1{line}, dm_side{line}, 1, 'ROM1');
            end
            if(strcmpi(dm_ROM_trial2{line}, 'null'))
                ROM_torque_2 = zeros;
                ROM_angle_2 = zeros;
            else
                [ROM_angle_2, ROM_torque_2] = extract_force_displ_singletrial_knee(dm_ROM_trial2{line}, dm_side{line}, 1, 'ROM2');
            end
        end
        
        
        %% AVERAGE two passive trials
        if(ROM_torque_1 == 0)
            data_ROM = average_passive_trials_knee(ROM_torque_2, ROM_angle_2);
        elseif(ROM_torque_2 == 0)
            data_ROM = average_passive_trials_knee(ROM_torque_1, ROM_angle_1);
        else % if 2 trials exist
            data_ROM = average_passive_trials_knee(ROM_torque_1, ROM_angle_1, ROM_torque_2, ROM_angle_2);
        end

       
        %% extract TORQUES and ANGLES
        
        % angles are:
        %       out_ROM_trial_max = trial max (different PRE, POST, L, R)
        %       out_ROM_ind_max = subject ind max (lowest PRE, POST, L, R)
        %       out_ROM_common_max = common max (lowest of all subjects)
        
        out_ROM_trial_max = min(data_ROM(:,2));
        %out_ROM_ind_max     - read from file
        %out_ROM_common_max  - read from file
        
        out_ROM_trial_min = max(data_ROM(:,2));
        out_ROM_givenangle = 61; %VAR - 3 subjects less flexible
        
        % corresponding torques (from averaged trials):
        loc_frame = find(data_ROM(:,2) >= out_ROM_trial_max, 1, 'first');
        out_T_trial_max = data_ROM(loc_frame,1);
        loc_frame = find(data_ROM(:,2) >= out_ROM_ind_max, 1, 'first');
        out_T_ind_max = data_ROM(loc_frame,1);
        loc_frame = find(data_ROM(:,2) >= out_ROM_common_max, 1, 'first');
        if isempty(loc_frame) == 0
            out_T_common_max = data_ROM(loc_frame,1);
        else
            out_T_common_max = NaN;
        end
        loc_frame = find(data_ROM(:,2) >= out_ROM_givenangle, 1, 'first');
        if isempty(loc_frame) == 0
            out_T_givenangle = data_ROM(loc_frame,1);
        else
            out_T_givenangle = NaN;
        end
        
        % highest torque
        [out_T_maxtorque,loc_frame] = max(data_ROM(:,1)); % highest T in array
        out_ROM_maxtorque = data_ROM(loc_frame,2);
        
        
        %% PASSIVE STIFFNESS (Nordez 2006) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % passive stiffness = delta torque / delta angle, at various angles

        % data_ROM(:,2) = angle
        % data_ROM(:,1) = torque

        % fit 4th order polynomial to averaged torque-angle curve, using data from STARTING angle 
        fit_ind_max = polyfit(data_ROM(:,2), data_ROM(:,1), 4);

        % extract passive stiffness (derivate of 4th order poly) at:
        out_pstiff_trial_max = (4 * fit_ind_max(1) * out_ROM_trial_max^3) + (3 * fit_ind_max(2) * out_ROM_trial_max^2) + (2 * fit_ind_max(3) * out_ROM_trial_max) + fit_ind_max(4);
        out_pstiff_ind_max = (4 * fit_ind_max(1) * out_ROM_ind_max^3) + (3 * fit_ind_max(2) * out_ROM_ind_max^2) + (2 * fit_ind_max(3) * out_ROM_ind_max) + fit_ind_max(4);
        out_pstiff_common_max = (4 * fit_ind_max(1) * out_ROM_common_max^3) + (3 * fit_ind_max(2) * out_ROM_common_max^2) + (2 * fit_ind_max(3) * out_ROM_common_max) + fit_ind_max(4);
        out_pstiff_maxtorque = (4 * fit_ind_max(1) * out_ROM_maxtorque^3) + (3 * fit_ind_max(2) * out_ROM_maxtorque^2) + (2 * fit_ind_max(3) * out_ROM_maxtorque) + fit_ind_max(4);

        if min(data_ROM(:,2)) > out_ROM_givenangle
            out_pstiff_givenangle = (4 * fit_ind_max(1) * out_ROM_givenangle^3) + (3 * fit_ind_max(2) * out_ROM_givenangle^2) + (2 * fit_ind_max(3) * out_ROM_givenangle) + fit_ind_max(4);
        else
            out_pstiff_givenangle = NaN;
        end
        
        % plotting various methods of curve fit equations:
        if plot_conversion
            fit_ind_max2 = polyfit(data_ROM(:,2), data_ROM(:,1), 3);
            fit_ind_max3 = polyfit(data_ROM(:,2), data_ROM(:,1), 2);
            plot_angle = angle_trial_onset:-0.5:min(data_ROM(:,2));
            plot_torque  = polyval(fit_ind_max,  angle_trial_onset:-0.5:min(data_ROM(:,2)));
            plot_torque2 = polyval(fit_ind_max2, angle_trial_onset:-0.5:min(data_ROM(:,2)));
            plot_torque3 = polyval(fit_ind_max3, angle_trial_onset:-0.5:min(data_ROM(:,2)));
            
            plottitle = horzcat('IND knee torque-angle fit for stiffness, ', subject_id);
            figure('Name',plottitle)
            hold on
            plot(data_ROM(:,2), data_ROM(:,1), 'LineWidth', 1)
            plot(plot_angle, plot_torque, ':', 'LineWidth', 2)
            plot(plot_angle, plot_torque2, '--')
            plot(plot_angle, plot_torque3, '-.')
            %   axis([0 Inf 0 Inf])
            set(gca,'xdir','reverse')
            xlabel(txt_gonio)
            ylabel('Force (N)')
            title(plottitle,'Interpreter', 'none')
            legend('Raw data', 'Fit 4th order','Fit 3rd order','Fit 2nd order','Location','Northwest')
            print(horzcat('data_plots/',plottitle),'-dpng')
        end
        
        %% STIFFNESS INDEX
        
        % fit 2nd order polynomial to averaged torque-angle curve, using data from zero angle to ind max ROM:
        fit_ind_max = polyfit(data_ROM(:,2), data_ROM(:,1), 2);

        % extract stiffness index as 2 * a
        out_pstiff_index = 2 * fit_ind_max(1);
        
        
        %% plot summary of 2 trials: Torque-angle + FIT
        if plot_check
            plot_angle = angle_trial_onset:-0.5:min(data_ROM(:,2));
            plot_torque  = polyval(fit_ind_max,  angle_trial_onset:-0.5:min(data_ROM(:,2)));
            
            plottitle = horzcat('IND knee avg torque-angle, ', subject_id);
            figure('Name',plottitle);
            hold on
            plot(ROM_angle_1,ROM_torque_1,'LineWidth',1) % trial1
            plot(ROM_angle_2,ROM_torque_2,'LineWidth',1) % trial2
            plot(data_ROM(:,2),data_ROM(:,1),'LineWidth',2) % mean data
            plot(plot_angle, plot_torque, ':', 'LineWidth', 2) % fit
            axis(axis_torque)
            set(gca,'Xdir','reverse')
            ylabel('Torque (Nm)')
            xlabel('<-- bent knee -- Angle (deg) -- straight knee -->')
            title(plottitle,'Interpreter', 'none')
            legend('Trial 1','Trial 2', 'Mean', 'Fit', 'Location','Southeast')
            print(horzcat('data_plots/',plottitle),'-dpng')
        end
        
        
        %% PREPARE FOR GROUP ARRAYS (..._angle_vars) - for group plots & stats
        
        % extract angle range common to all trials for current subject
        angle_start = angle_trial_onset + 0.00000001; % tweak for computer storage of numbers, where 8.3000 is stored as 8.2999999999999999 %VAR 
        angle_stop = out_ROM_trial_max;

        % identify locations of start/stop angles in above mentioned arrays
        loc_angle_start = find(data_ROM(:,2) <= angle_start,1,'first');
        loc_angle_stop = find(data_ROM(:,2) <= angle_stop,1,'first');
        if data_ROM(loc_angle_start,2) < angle_start - 1 %VAR
            cprintf('red', 'WARNING: Trial starts after 90 deg:', num2str(data_ROM(loc_angle_start,2)), '.\n')
        end
        
        % contents of below angle_vars arrays / angle_vars contain:
                %   1 angle
                %   2 torque 
                
        if trial_timepoint == 0 && trial_leg == 1 % PRE, STR
            % all data in ONE cell, common angles, RAW data:
            STR_PRE_angle_vars{STR_PRE_count} = [ ...
                data_ROM(loc_angle_start:loc_angle_stop,2) ...  1
                data_ROM(loc_angle_start:loc_angle_stop,1) ...  2
                ];

            %                 % all data in ONE cell, NORMALIZED data:
            %                 STR_PRE_angle_vars_norm_indlength{STR_PRE_count} = STR_PRE_angle_vars{STR_PRE_count};
            %                 STR_PRE_angle_vars_norm_indlength{STR_PRE_count}(:,1) = STR_PRE_angle_vars{1,STR_PRE_count}(:,1)*100/out_ROM_trial_max;                     % 1 angle - normalized to trial max ROM
            %                 STR_PRE_angle_vars_norm_indlength{STR_PRE_count}(:,2) = STR_PRE_angle_vars{1,STR_PRE_count}(:,2)*100/max(STR_PRE_angle_vars{1,STR_PRE_count}(:,2));   % 2 force - to maximal force in trial
            %
            %                 % resample for plots
            %                 % spline normalized data
            %                 STR_PRE_angle_vars_norm{STR_PRE_count} = spline(STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}(:,1)',STR_PRE_angle_vars_norm_indlength{1,STR_PRE_count}',0:angle_step_plots:100)';

        elseif trial_timepoint == 1 && trial_leg == 1 % POST, STR
            STR_POST_angle_vars{STR_POST_count} = [ ...
                data_ROM(loc_angle_start:loc_angle_stop,2) ...  1
                data_ROM(loc_angle_start:loc_angle_stop,1) ...  2
                ];

        elseif trial_timepoint == 0 && trial_leg == 0 % PRE, CON
            CON_PRE_angle_vars{CON_PRE_count} = [ ...
                data_ROM(loc_angle_start:loc_angle_stop,2) ...  1
                data_ROM(loc_angle_start:loc_angle_stop,1) ...  2
                ];

        elseif trial_timepoint == 1 && trial_leg == 0 % POST, CON
            CON_POST_angle_vars{CON_POST_count} = [ ...
                data_ROM(loc_angle_start:loc_angle_stop,2) ...  1
                data_ROM(loc_angle_start:loc_angle_stop,1) ...  2
                ];
        end
        
        
        %% INDIVIDUAL VARIABLES OUTPUT - prepare arrays with single trial data to csv file

        % txt trial ID
        all_passive_output_txt(line,:) = [dm_subjectno(line) dm_timepoint(line) dm_side(line) dm_trial(line)];
        
        % add data to a common array for all subjects    
        all_passive_output(line,:) = [ ...
            out_ROM_trial_max out_ROM_ind_max out_ROM_common_max out_ROM_maxtorque out_ROM_givenangle out_ROM_trial_min ...
        	out_T_trial_max out_T_ind_max out_T_common_max out_T_maxtorque out_T_givenangle ...
        	out_pstiff_trial_max out_pstiff_ind_max out_pstiff_common_max out_pstiff_maxtorque out_pstiff_givenangle ...
            out_pstiff_index ...
            ];

        save all_data_knee_inloop
        close all
    end
    %% LOOP FINISHED  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save all_data_knee
    
    
    %% Truncate angle_vars cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    STR_PRE_angle_vars(STR_PRE_count+1:end) = [];
%    STR_PRE_angle_vars_norm(STR_PRE_count+1:end) = [];
    STR_POST_angle_vars(STR_POST_count+1:end) = [];
%    STR_POST_angle_vars_norm(STR_POST_count+1:end) = [];
    CON_PRE_angle_vars(CON_PRE_count+1:end) = [];
%    CON_PRE_angle_vars_norm(CON_PRE_count+1:end) = [];
    CON_POST_angle_vars(CON_POST_count+1:end) = [];
%    CON_POST_angle_vars_norm(CON_POST_count+1:end) = [];
    
    
    %% OUTPUT individual trial data TO FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % write xls
    if ispc
        filename_output = strcat('data_output/all_passive_knee_output_', datestr(now, 'yyyy-mm-dd HH-MM'), '.xlsx');
        
        xlswrite(filename_output, all_passive_output_head, 1, 'A1')
        xlswrite(filename_output, all_passive_output_txt, 1, 'A2')
        xlswrite(filename_output, all_passive_output, 1, 'E2')
    else
        filename_output = strcat('data_output/all_passive_knee_output_', datestr(now, 'yyyy-mm-dd HH-MM'), '.csv');
        csvwrite(filename_output, all_passive_output)
    end
    
        
    %% OUTPUT group arrays for STATS, TO FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % variables to export to file
        %   1 Angle
        %  2 Torque
    out_arrays_input_cols = [1 2];
    out_arrays_input_labels = {'Angle' 'Torque'};

    if CON_PRE_count > 0 && CON_POST_count > 0 && STR_PRE_count > 0 && STR_POST_count > 0
        
        % preallocate
        STR_PRE_angle_vars_STAT{STR_PRE_count} = zeros;
%        STR_PRE_angle_vars_norm_STAT{STR_PRE_count} = zeros;
        STR_POST_angle_vars_STAT{STR_POST_count} = zeros;
%        STR_POST_angle_vars_norm_STAT{STR_POST_count} = zeros;
        CON_PRE_angle_vars_STAT{CON_PRE_count} = zeros;
%        CON_PRE_angle_vars_norm_STAT{CON_PRE_count} = zeros;
        CON_POST_angle_vars_STAT{CON_POST_count} = zeros;
%        CON_POST_angle_vars_norm_STAT{CON_POST_count} = zeros;
        
        % for absolute arrays: select common angle range = the subject/cell item containing the shortest matrix/ROM
        loc_commonROM = min([cellfun('length',STR_PRE_angle_vars) cellfun('length',STR_POST_angle_vars) cellfun('length',CON_PRE_angle_vars) cellfun('length',CON_POST_angle_vars)]); % location of largest common ROM
        
        
        % resample absolute & normalized arrays - ommit columns with NaN only
        
        for i = 1:STR_PRE_count
            locate_NaN_cols = ~all(isnan(STR_PRE_angle_vars{1,i}(1:5,:)),1);
            resample_axis = (angle_trial_onset:-angle_step_stats_abs:STR_PRE_angle_vars{1,i}(loc_commonROM,1))';
            STR_PRE_angle_vars_STAT{i} = STR_PRE_angle_vars{1,i}(1:numel(resample_axis),:);
            STR_PRE_angle_vars_STAT{i}(:,locate_NaN_cols) = spline(STR_PRE_angle_vars{1,i}(1:loc_commonROM,1)',STR_PRE_angle_vars{1,i}(1:loc_commonROM,locate_NaN_cols).',resample_axis).';
%            resample_axis = 0:angle_step_stats_norm:100;
%            STR_PRE_angle_vars_norm_STAT{i} = STR_PRE_angle_vars_norm{1,i}(1:numel(resample_axis),:);
%            STR_PRE_angle_vars_norm_STAT{i}(:,locate_NaN_cols) = spline(STR_PRE_angle_vars_norm{1,i}(1:loc_commonROM,1)',STR_PRE_angle_vars_norm{1,i}(1:loc_commonROM,locate_NaN_cols).',resample_axis).';
        end
        for i = 1:STR_POST_count
            locate_NaN_cols = ~all(isnan(STR_POST_angle_vars{1,i}(1:5,:)),1);
            resample_axis = (angle_trial_onset:-angle_step_stats_abs:STR_POST_angle_vars{1,i}(loc_commonROM,1))';
            STR_POST_angle_vars_STAT{i} = STR_POST_angle_vars{1,i}(1:numel(resample_axis),:);
            STR_POST_angle_vars_STAT{i}(:,locate_NaN_cols) = spline(STR_POST_angle_vars{1,i}(1:loc_commonROM,1)',STR_POST_angle_vars{1,i}(1:loc_commonROM,locate_NaN_cols).',resample_axis).';
%            resample_axis = 0:angle_step_stats_norm:100;
%            STR_POST_angle_vars_norm_STAT{i} = STR_POST_angle_vars_norm{1,i}(1:numel(resample_axis),:);
%            STR_POST_angle_vars_norm_STAT{i}(:,locate_NaN_cols) = spline(STR_POST_angle_vars_norm{1,i}(1:loc_commonROM,1)',STR_POST_angle_vars_norm{1,i}(1:loc_commonROM,locate_NaN_cols).',resample_axis).';
        end
        for i = 1:CON_PRE_count
            locate_NaN_cols = ~all(isnan(CON_PRE_angle_vars{1,i}(1:5,:)),1);
            resample_axis = (angle_trial_onset:-angle_step_stats_abs:CON_PRE_angle_vars{1,i}(loc_commonROM,1))';
            CON_PRE_angle_vars_STAT{i} = CON_PRE_angle_vars{1,i}(1:numel(resample_axis),:);
            CON_PRE_angle_vars_STAT{i}(:,locate_NaN_cols) = spline(CON_PRE_angle_vars{1,i}(1:loc_commonROM,1)',CON_PRE_angle_vars{1,i}(1:loc_commonROM,locate_NaN_cols).',resample_axis).';
%            resample_axis = 0:angle_step_stats_norm:100;
%            CON_PRE_angle_vars_norm_STAT{i} = CON_PRE_angle_vars_norm{1,i}(1:numel(resample_axis),:);
%            CON_PRE_angle_vars_norm_STAT{i}(:,locate_NaN_cols) = spline(CON_PRE_angle_vars_norm{1,i}(1:loc_commonROM,1)',CON_PRE_angle_vars_norm{1,i}(1:loc_commonROM,locate_NaN_cols).',resample_axis).';
        end
        for i = 1:CON_POST_count
            locate_NaN_cols = ~all(isnan(CON_POST_angle_vars{1,i}(1:5,:)),1);
            resample_axis = (angle_trial_onset:-angle_step_stats_abs:CON_POST_angle_vars{1,i}(loc_commonROM,1))';
            CON_POST_angle_vars_STAT{i} = CON_POST_angle_vars{1,i}(1:numel(resample_axis),:);
            CON_POST_angle_vars_STAT{i}(:,locate_NaN_cols) = spline(CON_POST_angle_vars{1,i}(1:loc_commonROM,1)',CON_POST_angle_vars{1,i}(1:loc_commonROM,locate_NaN_cols).',resample_axis).';
%            resample_axis = 0:angle_step_stats_norm:100;
%            CON_POST_angle_vars_norm_STAT{i} = CON_POST_angle_vars_norm{1,i}(1:numel(resample_axis),:);
%            CON_POST_angle_vars_norm_STAT{i}(:,locate_NaN_cols) = spline(CON_POST_angle_vars_norm{1,i}(1:loc_commonROM,1)',CON_POST_angle_vars_norm{1,i}(1:loc_commonROM,locate_NaN_cols).',resample_axis).';
        end
        
        % create table headers (subject numbers)
        
        % pre and post values
        out_arrays_headers{1+STR_PRE_count+STR_POST_count+CON_PRE_count+CON_POST_count} = [];
        
        out_arrays_headers{1} = 'Joint_angle';
        for i=1:STR_PRE_count
            out_arrays_headers{i+1} = STR_PRE_ID{i};
        end
        for i=1:CON_PRE_count
            out_arrays_headers{i+1+STR_PRE_count} = CON_PRE_ID{i};
        end
        for i=1:STR_POST_count
            out_arrays_headers{i+1+STR_PRE_count+CON_PRE_count} = STR_POST_ID{i};
        end
        for i=1:CON_POST_count
            out_arrays_headers{i+1+STR_PRE_count+CON_PRE_count+STR_POST_count} = CON_POST_ID{i};
        end
        
        % difference (post minus pre) values
        if eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count)
            % output of pre-post DIFFERENCES requires that trials are analysed in systematic order (does not check which subject-numbers to group)
            % - but will not be exported unless equal amount of 4 legs/timepoints - rough coding
            out_arrays_headers_diff{1+STR_PRE_count+CON_PRE_count} = [];
            out_arrays_headers_diff{1} = 'Joint_angle';
            for i=1:STR_PRE_count
                out_arrays_headers_diff{i+1} = strcat('STR_PP_', num2str(STR_PRE_no(i)), '__', num2str(STR_POST_no(i)));
            end
            for i=1:CON_PRE_count
                out_arrays_headers_diff{i+1+STR_PRE_count} = strcat('CON_PP_', num2str(CON_PRE_no(i)), '__',num2str(CON_POST_no(i)));
            end
        end
        
        % preallocate output arrays
        cols_abs = size(STR_PRE_angle_vars_STAT{1},1);
%        cols_norm = size(STR_PRE_angle_vars_norm_STAT{1},1);
        rows = STR_PRE_count+STR_POST_count+CON_PRE_count+CON_POST_count + 1; % adding 1 for column for joint angles
        rows_diff = STR_PRE_count+CON_PRE_count + 1; % adding 1 for column for joint angles
        out_arrays_abs(cols_abs,rows) = zeros;
%        out_arrays_norm(cols_norm,rows) = zeros;
        out_arrays_abs_diff(cols_abs,rows_diff) = zeros;
%        out_arrays_norm_diff(cols_norm,rows_diff) = zeros;

        % organize and output table for each of the selected variables
        for var = 1:length(out_arrays_input_cols)
            % reset output arrays
            out_arrays_abs(1:cols_abs,1:rows) = zeros;
%            out_arrays_norm(cols_norm,rows) = zeros;
            out_arrays_abs_diff(1:cols_abs,1:rows_diff) = zeros;
%            out_arrays_norm_diff(cols_norm,rows_diff) = zeros;
            
            % add as first column, joint angles: abs and normalized angles
            out_arrays_abs(:,1) = STR_PRE_angle_vars_STAT{1}(:,1);
%            out_arrays_norm(:,1) = STR_PRE_angle_vars_norm_STAT{1}(:,1);
            out_arrays_abs_diff(:,1) = STR_PRE_angle_vars_STAT{1}(:,1);
%            out_arrays_norm_diff(:,1) = STR_PRE_angle_vars_norm_STAT{1}(:,1);
            
            % add values: pre and post
            
            % add STR PRE first
            for subj = 1:STR_PRE_count
                % absolute values
                out_arrays_abs(:,subj+1) = STR_PRE_angle_vars_STAT{subj}(:,out_arrays_input_cols(var));
%                % normalized values
%                out_arrays_norm(:,subj+1) = STR_PRE_angle_vars_norm_STAT{subj}(:,out_arrays_input_cols(var));
            end
            
            % add STR POST second
            for subj = 1:STR_POST_count
                % absolute values
                out_arrays_abs(:,subj+STR_PRE_count+1) = STR_POST_angle_vars_STAT{subj}(:,out_arrays_input_cols(var));
%                 % normalized values
%                 out_arrays_norm(:,subj+STR_PRE_count+1) = STR_POST_angle_vars_norm_STAT{subj}(:,out_arrays_input_cols(var));
            end
            
            % add CON PRE
            for subj = 1:CON_PRE_count
                % absolute values
                out_arrays_abs(:,subj+STR_PRE_count+STR_POST_count+1) = CON_PRE_angle_vars_STAT{subj}(:,out_arrays_input_cols(var));
%                % normalized values
%                out_arrays_norm(:,subj+STR_PRE_count+STR_POST_count+1) = CON_PRE_angle_vars_norm_STAT{subj}(:,out_arrays_input_cols(var));
            end
            
            % add CON POST
            for subj = 1:CON_POST_count
                % absolute values
                out_arrays_abs(:,subj+STR_PRE_count+STR_POST_count+CON_PRE_count+1) = CON_POST_angle_vars_STAT{subj}(:,out_arrays_input_cols(var));
%                % normalized values
%                out_arrays_norm(:,subj+STR_PRE_count+STR_POST_count+CON_PRE_count+1) = CON_POST_angle_vars_norm_STAT{subj}(:,out_arrays_input_cols(var));
            end
            
            % add values: difference between PRE-POST
            if eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count)
                % add STR first
                for subj = 1:STR_PRE_count
                    % absolute values
                    out_arrays_abs_diff(:,subj+1) = STR_POST_angle_vars_STAT{subj}(:,out_arrays_input_cols(var)) - STR_PRE_angle_vars_STAT{subj}(:,out_arrays_input_cols(var));
%                    % normalized values
%                    out_arrays_norm_diff(:,subj+1) = STR_POST_angle_vars_norm_STAT{subj}(:,out_arrays_input_cols(var)) - STR_PRE_angle_vars_norm_STAT{subj}(:,out_arrays_input_cols(var));
                end

                % add CON second
                for subj = 1:CON_PRE_count
                    % absolute values
                    out_arrays_abs_diff(:,subj+STR_PRE_count+1) = CON_POST_angle_vars_STAT{subj}(:,out_arrays_input_cols(var)) - CON_PRE_angle_vars_STAT{subj}(:,out_arrays_input_cols(var));
%                    % normalized values
%                    out_arrays_norm_diff(:,subj+STR_PRE_count+1) = CON_POST_angle_vars_norm_STAT{subj}(:,out_arrays_input_cols(var)) - CON_PRE_angle_vars_norm_STAT{subj}(:,out_arrays_input_cols(var));
                end
            end
            
            % create tables and save as file
            % pre and post values
            out_arrays_abs_table = array2table(out_arrays_abs,'VariableNames',out_arrays_headers);
            filename_output = strcat('data_output/intervention_arrays_', out_arrays_input_labels{var} , '_abs_', datestr(now, 'yyyy-mm-dd HH-MM'));
            writetable(out_arrays_abs_table,filename_output,'Delimiter','\t')
            
%             out_arrays_norm_table = array2table(out_arrays_abs,'VariableNames',out_arrays_headers);
%             filename_output = strcat('data_output/intervention_arrays_', out_arrays_input_labels{var} , '_norm_', datestr(now, 'yyyy-mm-dd HH-MM'));
%             writetable(out_arrays_norm_table,filename_output,'Delimiter','\t')
            
            % difference values
            if eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count)
                out_arrays_abs_diff_table = array2table(out_arrays_abs_diff,'VariableNames',out_arrays_headers_diff);
                filename_output = strcat('data_output/intervention_arrays_', out_arrays_input_labels{var} , '_abs_P-P_', datestr(now, 'yyyy-mm-dd HH-MM'));
                writetable(out_arrays_abs_diff_table,filename_output,'Delimiter','\t')

%                 out_arrays_norm_diff_table = array2table(out_arrays_abs_diff,'VariableNames',out_arrays_headers_diff);
%                 filename_output = strcat('data_output/intervention_arrays_', out_arrays_input_labels{var} , '_norm_P-P_', datestr(now, 'yyyy-mm-dd HH-MM'));
%                 writetable(out_arrays_norm_diff_table,filename_output,'Delimiter','\t')
            end
            clear out_arrays_abs_table out_arrays_norm_table out_arrays_abs_diff_table out_arrays_norm_diff_table
        end
    end
        
    
    %% GROUP CALCULATIONS - create variables of MEAN + STDAV for plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  mean and stdav of each subject's INDIVIDUAL MAX ROM, force, elong, EMG, etc

        % STR PRE
        if STR_PRE_count > 0
            n_o_array_elements = length(STR_PRE_angle_vars{1,1}(1,:));
            % preallocate array
            STR_PRE_max(STR_PRE_count,n_o_array_elements) = zeros;
    %        STR_PRE_max_norm(STR_PRE_count,n_o_array_elements) = zeros;
            % collect variables per subject
            for i = 1:STR_PRE_count % per subject
                for j = 1:n_o_array_elements % per element in arrays
                    STR_PRE_max(i,j) = max(STR_PRE_angle_vars{1,i}(end,j));
     %               STR_PRE_max_norm(i,j) = max(STR_PRE_angle_vars_norm{1,i}(end,j));
                end
            end
            % calculate mean and SD of max values across subjects
            STR_PRE_ROM_mean = mean(STR_PRE_max(:,1));
            STR_PRE_ROM_SD = std(STR_PRE_max(:,1));
            STR_PRE_torque_mean = mean(STR_PRE_max(:,2));
            STR_PRE_torque_SD = std(STR_PRE_max(:,2));
            % determine common angle range
            STR_PRE_common_ROM = max(STR_PRE_max(:,1));
        end
        
        % STR POST
        if STR_POST_count > 0
            n_o_array_elements = length(STR_POST_angle_vars{1,1}(1,:));
            % preallocate array
            STR_POST_max(STR_POST_count,n_o_array_elements) = zeros;
     %       STR_POST_max_norm(STR_POST_count,n_o_array_elements) = zeros;
            % collect variables per subject
            for i = 1:STR_POST_count % per subject
                for j = 1:n_o_array_elements % per element in arrays
                    STR_POST_max(i,j) = max(STR_POST_angle_vars{1,i}(end,j));
    %                STR_POST_max_norm(i,j) = max(STR_POST_angle_vars_norm{1,i}(end,j));
                end
            end
            % calculate mean and SD of max values across subjects
            STR_POST_ROM_mean = mean(STR_POST_max(:,1));
            STR_POST_ROM_SD = std(STR_POST_max(:,1));
            STR_POST_torque_mean = mean(STR_POST_max(:,2));
            STR_POST_torque_SD = std(STR_POST_max(:,2));
            % determine common angle range
            STR_POST_common_ROM = max(STR_POST_max(:,1));
        end
        
        % CON PRE
        if CON_PRE_count > 0
            n_o_array_elements = length(CON_PRE_angle_vars{1,1}(1,:));
            % preallocate array
            CON_PRE_max(CON_PRE_count,n_o_array_elements) = zeros;
   %         CON_PRE_max_norm(CON_PRE_count,n_o_array_elements) = zeros;
            % collect variables per subject
            for i = 1:CON_PRE_count % per subject
                for j = 1:n_o_array_elements % per element in arrays
                    CON_PRE_max(i,j) = max(CON_PRE_angle_vars{1,i}(end,j));
    %                CON_PRE_max_norm(i,j) = max(CON_PRE_angle_vars_norm{1,i}(end,j));
                end
            end
            % calculate mean and SD of max values across subjects
            CON_PRE_ROM_mean = mean(CON_PRE_max(:,1));
            CON_PRE_ROM_SD = std(CON_PRE_max(:,1));
            CON_PRE_torque_mean = mean(CON_PRE_max(:,2));
            CON_PRE_torque_SD = std(CON_PRE_max(:,2));
            % determine common angle range
            CON_PRE_common_ROM = max(CON_PRE_max(:,1));
        end

        % CON POST
        if CON_POST_count > 0
            n_o_array_elements = length(CON_POST_angle_vars{1,1}(1,:));
            % preallocate array
            CON_POST_max(CON_POST_count,n_o_array_elements) = zeros;
      %      CON_POST_max_norm(CON_POST_count,n_o_array_elements) = zeros;
            % collect variables per subject
            for i = 1:CON_POST_count % per subject
                for j = 1:n_o_array_elements % per element in arrays
                    CON_POST_max(i,j) = max(CON_POST_angle_vars{1,i}(end,j));
    %                CON_POST_max_norm(i,j) = max(CON_POST_angle_vars_norm{1,i}(end,j));
                end
            end
            % calculate mean and SD of max values across subjects
            CON_POST_ROM_mean = mean(CON_POST_max(:,1));
            CON_POST_ROM_SD = std(CON_POST_max(:,1));
            CON_POST_torque_mean = mean(CON_POST_max(:,2));
            CON_POST_torque_SD = std(CON_POST_max(:,2));
            % determine common angle range
            CON_POST_common_ROM = max(CON_POST_max(:,1));
        end

                         
    %% GROUP CALCULATIONS - create AVERAGE ARRAYS for plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if STR_PRE_count > 0
        %%% average ABSOLUTE arrays, up to all subjects' COMMON MAX ROM, for force, elong, EMG
        
        % preallocate - to shortest length
        len = 10000;
        for i = 1:STR_PRE_count
            if length(STR_PRE_angle_vars{:,i}) < len
                len = length(STR_PRE_angle_vars{:,i});
            end
        end
        STR_PRE_angle_vars_mean_tmp(len,n_o_array_elements,STR_PRE_count) = zeros;
        
        % STR_PRE_angle_vars has same angles (column 1) for all subjects, so max common ROM will be on the same line for all subjects
        loc_end = find(STR_PRE_angle_vars{1,STR_PRE_count}(:,1) <= (STR_PRE_common_ROM + 0.00001), 1, 'first'); % using last subject (STR_PRE_count) - could use any, angle is the same in all
        for i = 1:STR_PRE_count
            STR_PRE_angle_vars_mean_tmp(:,:,i) = STR_PRE_angle_vars{i}(1:loc_end,:);
        end
        STR_PRE_angle_vars_mean = nanmean(STR_PRE_angle_vars_mean_tmp, 3);
        %STR_PRE_angle_vars_SD = nanstd(STR_PRE_angle_vars_mean_tmp,1,3);
        
        %             %%% average NORMALIZED arrays, up to all subjects' INDIVIDUAL ROM, for force, elong, EMG
        %
        %             % preallocate
        %             STR_PRE_angle_vars_norm_mean_tmp(length(STR_PRE_angle_vars_norm{:,1}),n_o_array_elements,STR_PRE_count) = zeros;
        %
        %             % STR_PRE_angle_vars_norm has same angles (column 1) for all subjects, from 0 to 100% for all subjects
        %             for i = 1:STR_PRE_count
        %                 STR_PRE_angle_vars_norm_mean_tmp(:,:,i) = STR_PRE_angle_vars_norm{i}(:,:);
        %             end
        %             STR_PRE_angle_vars_norm_mean = nanmean(STR_PRE_angle_vars_norm_mean_tmp, 3);
        
        %%% clean up
        
        clear STR_PRE_angle_vars_mean_tmp STR_PRE_angle_vars_norm_mean_tmp
    end
    
    
    if STR_POST_count > 0
        %%% average ABSOLUTE arrays, up to all subjects' COMMON MAX ROM, for force, elong, EMG
        
        % preallocate
        len = 10000;
        for i = 1:STR_POST_count
            if length(STR_POST_angle_vars{:,i}) < len
                len = length(STR_POST_angle_vars{:,i});
            end
        end
        STR_POST_angle_vars_mean_tmp(len,n_o_array_elements,STR_POST_count) = zeros;
        
        % STR_POST_angle_vars has same angles (column 1) for all subjects, so max common ROM will be on the same line for all subjects
        loc_end = find(STR_POST_angle_vars{1,STR_POST_count}(:,1) <= (STR_POST_common_ROM + 0.00001), 1, 'first'); % using last subject (STR_POST_count) - could use any, angle is the same in all
        for i = 1:STR_POST_count
            STR_POST_angle_vars_mean_tmp(:,:,i) = STR_POST_angle_vars{i}(1:loc_end,:);
        end
        STR_POST_angle_vars_mean = nanmean(STR_POST_angle_vars_mean_tmp, 3);
        %STR_POST_angle_vars_SD = nanstd(STR_POST_angle_vars_mean_tmp,1,3);
        
        %             %%% average NORMALIZED arrays, up to all subjects' INDIVIDUAL ROM, for force, elong, EMG
        %
        %             % preallocate
        %             STR_POST_angle_vars_norm_mean_tmp(length(STR_POST_angle_vars_norm{:,1}),n_o_array_elements,STR_POST_count) = zeros;
        %
        %             % STR_POST_angle_vars_norm has same angles (column 1) for all subjects, from 0 to 100% for all subjects
        %             for i = 1:STR_POST_count
        %                 STR_POST_angle_vars_norm_mean_tmp(:,:,i) = STR_POST_angle_vars_norm{i}(:,:);
        %             end
        %             STR_POST_angle_vars_norm_mean = nanmean(STR_POST_angle_vars_norm_mean_tmp, 3);
        
        %%% clean up
        
        clear STR_POST_angle_vars_mean_tmp STR_POST_angle_vars_norm_mean_tmp
    end
    
    
    if CON_PRE_count > 0
        %%% average ABSOLUTE arrays, up to all subjects' COMMON MAX ROM, for force, elong, EMG
        
        % preallocate
        len = 10000;
        for i = 1:CON_PRE_count
            if length(CON_PRE_angle_vars{:,i}) < len
                len = length(CON_PRE_angle_vars{:,i});
            end
        end
        CON_PRE_angle_vars_mean_tmp(len,n_o_array_elements,CON_PRE_count) = zeros;
        
        % CON_PRE_angle_vars has same angles (column 1) for all subjects, so max common ROM will be on the same line for all subjects
        loc_end = find(CON_PRE_angle_vars{1,CON_PRE_count}(:,1) <= (CON_PRE_common_ROM + 0.00001), 1, 'first'); % using last subject (CON_PRE_count) - could use any, angle is the same in all
        for i = 1:CON_PRE_count
            CON_PRE_angle_vars_mean_tmp(:,:,i) = CON_PRE_angle_vars{i}(1:loc_end,:);
        end
        CON_PRE_angle_vars_mean = nanmean(CON_PRE_angle_vars_mean_tmp, 3);
        %CON_PRE_angle_vars_SD = nanstd(CON_PRE_angle_vars_mean_tmp,1,3);
        
        %             %%% average NORMALIZED arrays, up to all subjects' INDIVIDUAL ROM, for force, elong, EMG
        %
        %             % preallocate
        %             CON_PRE_angle_vars_norm_mean_tmp(length(CON_PRE_angle_vars_norm{:,1}),n_o_array_elements,CON_PRE_count) = zeros;
        %
        %             % CON_PRE_angle_vars_norm has same angles (column 1) for all subjects, from 0 to 100% for all subjects
        %             for i = 1:CON_PRE_count
        %                 CON_PRE_angle_vars_norm_mean_tmp(:,:,i) = CON_PRE_angle_vars_norm{i}(:,:);
        %             end
        %             CON_PRE_angle_vars_norm_mean = nanmean(CON_PRE_angle_vars_norm_mean_tmp, 3);
        
        %%% clean up
        clear CON_PRE_angle_vars_mean_tmp CON_PRE_angle_vars_norm_mean_tmp
    end
    
    
    if CON_POST_count > 0
        %%% average ABSOLUTE arrays, up to all subjects' COMMON MAX ROM, for force, elong, EMG
        
        % preallocate
        len = 10000;
        for i = 1:CON_POST_count
            if length(CON_POST_angle_vars{:,i}) < len
                len = length(CON_POST_angle_vars{:,i});
            end
        end
        CON_POST_angle_vars_mean_tmp(len,n_o_array_elements,CON_POST_count) = zeros;
        
        % CON_POST_angle_vars has same angles (column 1) for all subjects, so max common ROM will be on the same line for all subjects
        loc_end = find(CON_POST_angle_vars{1,CON_POST_count}(:,1) <= (CON_POST_common_ROM + 0.00001), 1, 'first'); % using last subject (CON_POST_count) - could use any, angle is the same in all
        for i = 1:CON_POST_count
            CON_POST_angle_vars_mean_tmp(:,:,i) = CON_POST_angle_vars{i}(1:loc_end,:);
        end
        CON_POST_angle_vars_mean = nanmean(CON_POST_angle_vars_mean_tmp, 3);
        %CON_POST_angle_vars_SD = nanstd(CON_POST_angle_vars_mean_tmp,1,3);
        
        %             %%% average NORMALIZED arrays, up to all subjects' INDIVIDUAL ROM, for force, elong, EMG
        %
        %             % preallocate
        %             CON_POST_angle_vars_norm_mean_tmp(length(CON_POST_angle_vars_norm{:,1}),n_o_array_elements,CON_POST_count) = zeros;
        %
        %             % CON_POST_angle_vars_norm has same angles (column 1) for all subjects, from 0 to 100% for all subjects
        %             for i = 1:CON_POST_count
        %                 CON_POST_angle_vars_norm_mean_tmp(:,:,i) = CON_POST_angle_vars_norm{i}(:,:);
        %             end
        %             CON_POST_angle_vars_norm_mean = nanmean(CON_POST_angle_vars_norm_mean_tmp, 3);
        
        %%% clean up
        
        clear CON_POST_angle_vars_mean_tmp CON_POST_angle_vars_norm_mean_tmp
        
    end
    
    
    %% PLOT GROUP FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % TORQUE-angle
    
    if plot_check && CON_PRE_count > 0 && CON_POST_count > 0 && STR_PRE_count > 0 && STR_POST_count > 0
        plottitle = horzcat('torque vs angle - 1 GROUP');
        figure('Name',plottitle)
        hold on
        plot(STR_PRE_angle_vars_mean(:,1), STR_PRE_angle_vars_mean(:,2),'Color',col_lightred,'LineStyle','--','LineWidth',1)
        plot(STR_POST_angle_vars_mean(:,1), STR_POST_angle_vars_mean(:,2),'r','LineStyle','-','LineWidth',1)
        plot(CON_PRE_angle_vars_mean(:,1), CON_PRE_angle_vars_mean(:,2),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
        plot(CON_POST_angle_vars_mean(:,1), CON_POST_angle_vars_mean(:,2),'b','LineStyle','-','LineWidth',1)
        
        herrorbar(STR_PRE_ROM_mean, STR_PRE_torque_mean, STR_PRE_ROM_SD, '*m')
        errorbar(STR_PRE_ROM_mean, STR_PRE_torque_mean, STR_PRE_torque_SD, 'Color', col_lightred, 'Marker', '.', 'MarkerFaceColor', col_lightred)
        herrorbar(STR_POST_ROM_mean, STR_POST_torque_mean, STR_POST_ROM_SD, 'r.')
        errorbar(STR_POST_ROM_mean, STR_POST_torque_mean, STR_POST_torque_SD, 'Color', 'r', 'Marker', '.', 'MarkerFaceColor', 'r')
        
        herrorbar(CON_PRE_ROM_mean, CON_PRE_torque_mean, CON_PRE_ROM_SD, '*c')
        errorbar(CON_PRE_ROM_mean, CON_PRE_torque_mean, CON_PRE_torque_SD, 'Color', col_lightblue, 'Marker', '.', 'MarkerFaceColor', col_lightblue)
        herrorbar(CON_POST_ROM_mean, CON_POST_torque_mean, CON_POST_ROM_SD, 'b.')
        errorbar(CON_POST_ROM_mean, CON_POST_torque_mean, CON_POST_torque_SD, 'Color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b')
        
        axis(axis_torque)
        set(gca,'Xdir','reverse')
        xlabel(txt_gonio)
        ylabel('Torque (Nm)')
        title(plottitle,'Interpreter', 'none')
        legend('STR PRE','STR POST','CON PRE','CON POST','S.PRE indmax','Location','Northwest')
        print(horzcat('data_plots/GRP_INT_knee ',plottitle),'-dpng')
    end
    
    if plot_check && STR_PRE_count > 0 && STR_POST_count > 0
        plottitle = horzcat('torque vs angle - 2 STRETCH PRE-POST');
        figure('Name',plottitle)
        hold on
        for i = 1:STR_PRE_count
            plot(STR_PRE_angle_vars{1,i}(:,1),STR_PRE_angle_vars{1,i}(:,2),'LineStyle','--')
        end
        set(gca,'ColorOrderIndex',1)
        for i = 1:STR_POST_count
            plot(STR_POST_angle_vars{1,i}(:,1),STR_POST_angle_vars{1,i}(:,2))
        end
        axis(axis_torque)
        set(gca,'Xdir','reverse')
        xlabel(txt_gonio)
        ylabel('Torque (Nm)')
        title(plottitle,'Interpreter', 'none')
        %legend
        print(horzcat('data_plots/GRP_INT_knee ',plottitle),'-dpng')
    end
    
    if plot_check && CON_PRE_count > 0 && CON_POST_count > 0
        plottitle = horzcat('torque vs angle - 3 CONTROL PRE-POST');
        figure('Name',plottitle)
        hold on
        for i = 1:CON_PRE_count
            plot(CON_PRE_angle_vars{1,i}(:,1),CON_PRE_angle_vars{1,i}(:,2),'LineStyle','--')
        end
        set(gca,'ColorOrderIndex',1)
        for i = 1:CON_POST_count
            plot(CON_POST_angle_vars{1,i}(:,1),CON_POST_angle_vars{1,i}(:,2))
        end
        axis(axis_torque)
        set(gca,'Xdir','reverse')
        xlabel(txt_gonio)
        ylabel('Torque (Nm)')
        title(plottitle,'Interpreter', 'none')
        %legend
        print(horzcat('data_plots/GRP_INT_knee ',plottitle),'-dpng')
    end
    
    % rough coding
    if plot_check && eq(CON_PRE_count, CON_POST_count) && eq(STR_PRE_count, STR_POST_count) && eq(CON_PRE_count,STR_PRE_count) && plot_individual
        for i = 1:CON_PRE_count
            plottitle = horzcat('IND knee torque vs angle - SUBJECT ', CON_PRE_ID{i} ,' PRE-POST');
            figure('Name',plottitle)
            hold on
            
            plot(STR_PRE_angle_vars{1,i}(:,1), STR_PRE_angle_vars{1,i}(:,2),'Color',col_lightred,'LineStyle','--','LineWidth',1)
            plot(STR_POST_angle_vars{1,i}(:,1), STR_POST_angle_vars{1,i}(:,2),'r','LineStyle','-','LineWidth',1)
            plot(CON_PRE_angle_vars{1,i}(:,1), CON_PRE_angle_vars{1,i}(:,2),'Color',col_lightblue,'LineStyle','--','LineWidth',1)
            plot(CON_POST_angle_vars{1,i}(:,1), CON_POST_angle_vars{1,i}(:,2),'b','LineStyle','-','LineWidth',1)
            axis(axis_torque)
            set(gca,'Xdir','reverse')
            xlabel(txt_gonio)
            ylabel('Torque (Nm)')
            title(plottitle,'Interpreter', 'none')
            legend('STR PRE','STR POST','CON PRE','CON POST','Location','Northwest')
            print(horzcat('data_plots/',plottitle),'-dpng')
        end
    end
    
end