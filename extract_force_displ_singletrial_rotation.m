%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract_force_displ_singletrial_rotation
% Marie Moltubakk 2017-12-08
% Read complete, prepared noraon array + prepared us array
% Produce array with correction for ankle rotation upon gonio angle

% used by stiffness analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [time_force_displ_array,maxforce] = extract_force_displ_singletrial_rotation(noraxondata, usdata, usdata_frame, coact_max_torque, coact_max_EMG, at_momentarm, at_rotation_const, side, trial_name)
    global column_achilles column_gonio % noraxonfreq
%    global plot_norm plot_achilles plot_check subject_id % plot_achilles  plot_emg 
    global subject_id
    global filepath
    
    
    %% prepare files
    % Read stiffness trial US data file, determine time stamps, set trigger frame as time = zero
    % Produce US sample frequency + new US array containing time and displacement
    [usdata_prepped,usfreq] = read_us_file(strcat(filepath, usdata, '.txt'), str2double(usdata_frame), trial_name);
    
    % Read stiffness trial Noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
    % Produce a new noraxon data array
    noraxon_prepped = read_noraxon_stiffness(strcat(filepath, noraxondata), usfreq, side, trial_name);
    
    
    %% GONIOMETER - secondary zero offset (primary offset during data collection = only done before 1st of 3 trials)
    gonio_offset = mean(noraxon_prepped(1:3,column_gonio));
    % expexting that gonio should be zero at onset of singletrials
    noraxon_prepped(:,column_gonio) = noraxon_prepped(:,column_gonio) - gonio_offset;
    
    % deal with special cases from raw data
    if strcmp(subject_id, 'INT_4_SOL_PRE_STR_R')
        % reverse goniometer data for all SOL trials
        noraxon_prepped(:,column_gonio) = -noraxon_prepped(:,column_gonio);
    elseif strcmp(subject_id, 'INT_4_GM_PRE_STR_R')
        if strcmp(trial_name, 'MTJ2-rotation')
            % reverse goniometer data for GM trials EXCEPT GM MTJ2
        else
            noraxon_prepped(:,column_gonio) = -noraxon_prepped(:,column_gonio);
        end
    end
    
    
    %% FORCE data - identical to MTJ trial analysis
    % (could be ommitted here)
    
    % Correct TORQUE for co-activation
    noraxon_prepped_coact = correct_coactivation(noraxon_prepped, coact_max_torque, coact_max_EMG, side, trial_name);
    
    % Calculate force by applying internal moment arm
    tendon_force = noraxon_prepped_coact(:,column_achilles) / at_momentarm;
    
    % secondary zero offset for force
    tendon_force_avg = mean(tendon_force(1:5)); % avg of first 5 frames. Previously: 1/20th sec (0,05 s) 
    tendon_force_offset = tendon_force - tendon_force_avg;
    
    
    %%% Force cutting method 2, april 2004
    % prepare for US displacement cutoff for MTJ trials(was: cutoff at x% of derivative)
%    percent_elong = 0.96; % alternative method, keep 96% of total elongation %VAR
%    [max_displ,max_displ_pos] = max(usdata_prepped(:,2));
%    keep_elong = percent_elong * max_displ;
%    elong_cut_index = find(usdata_prepped(:,2)>keep_elong,1,'first');
%    lastframe = elong_cut_index-1;
%    time1 = usdata_prepped(lastframe,1);

    
%    %%% Force cutting method 1, march 2004
%        derivative_cutoff_const = 2.00; % derivative value to cut off at; %VAR
%	derivative_lowest_acceptable = 0; % ignore cutoffs between above value and this value, if haven't reached below %
%	derivative_min_elong = 0.8; % in percent of max elong
%
%        % calculate moving average of displacement data
%        window_mov_avg = round(usfreq / 5);
%        a = 1;
%        b = ones(1,window_mov_avg)/window_mov_avg;
%        usdata_filtfilt = filtfilt(b, a, usdata_prepped(:,2));
%        usdata_derivative(1:length(usdata_filtfilt)) = zeros;
%
%        % derivative of displacement data
%        dt = 1/usfreq;
%        for i = 2:length(usdata_filtfilt)-1
%            usdata_derivative(i) = (usdata_filtfilt(i+1)-usdata_filtfilt(i-1))/dt;
%        end
%
%        % position of X % of derivative
%        [max_deriv_val,max_deriv_pos] = max(usdata_derivative);
%        derivative_cut = derivative_cutoff_const; % * max_deriv_val for percent
%        derivative_cut_index = find(usdata_derivative(max_deriv_pos:length(usdata_derivative))<derivative_cut,1,'first');
%        lastframe = derivative_cut_index+max_deriv_pos;
%        % if not too low low derivative is found very early...
%        while (usdata_filtfilt(lastframe) < max(usdata_filtfilt) * derivative_min_elong && usdata_derivative(lastframe) > derivative_lowest_acceptable)
%            % find NEXT lower derivative
%            derivative_cut_index2 = find(usdata_derivative(lastframe:length(usdata_derivative))<derivative_cut,1,'first');
%            lastframe = lastframe + derivative_cut_index2;
%        end
%        time1 = usdata_prepped(lastframe,1);
%
%        % plot derivative cutoff
%        if plot_check
%            plottitle = horzcat('Displacement cutoff ', subject_id, ' ', trial_name);
%            fig_displcut = figure('Name',plottitle);
%            ylim([-2,inf]);
%            hold on
%            plot(usdata_prepped(:,1),usdata_prepped(:,2),'r')
%            plot(usdata_prepped(:,1),usdata_filtfilt,'g')
%            plot(usdata_prepped(:,1),usdata_derivative,'b')
%            plot([time1 time1], [min(usdata_derivative) max(usdata_derivative)],'k');
%            xlabel('Time (s)'),ylabel('Data (misc)'),title(plottitle);
%            legend('Orig displacement','Filtered displ.','Derivative','Cutoff','Location','Northwest');
%            saveas(fig_displcut, strcat('data_output/displ_cut_', subject_id, '_', trial_name), 'png')
%        end


%    %%% Force cutting method 2, april 2004
%    % cut MTJ usdata at calculated position
%    if strcmp(trial_name(1),'M') == 1  % = MTJ trial - prepare for cutting at max elong
%        usdata_prepped = usdata_prepped(1:lastframe,:); % cutting all columns, in the data without any corrections/filters
%    end
    
    
    %% DISPLACEMENT
    % Correct displacement for ankle rotation
    at_rotation_const = str2double(at_rotation_const);
    
    % length of data
    commonlength = min(length(usdata_prepped),length(noraxon_prepped));
    
    % save time column
    usdata_corrected(:,1) = usdata_prepped(1:commonlength,1);
    
    % for trials to be corrected for rotation, a rotation constant is passed and applied as follows
    usdata_corrected(:,2) = -noraxon_prepped(1:commonlength,column_gonio) * at_rotation_const;


    %% Create outout array: time, force, displacement
    len_array = min(length(tendon_force_offset),length(usdata_prepped));
    time_force_displ_array(:,1) = usdata_prepped(1:len_array,1); % time array
    time_force_displ_array(:,2) = tendon_force_offset(1:len_array); % force array
    time_force_displ_array(:,3) = usdata_corrected(1:len_array,2); % displacement array
    
    % output raw max force from trial
    maxforce = max(tendon_force_offset(1:len_array));

end