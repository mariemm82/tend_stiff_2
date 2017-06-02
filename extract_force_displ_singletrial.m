%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract_force_displ_singletrial
% Marie Moltubakk 11.6.2013
% Read complete, prepared noraon array + prepared us array
% Produce array with time, corrected force, corrected displacement
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [time_force_displ_array,maxforce] = extract_force_displ_singletrial(noraxondata, usdata, usdata_frame, coact_max_torque, coact_max_EMG, at_momentarm, at_rotation_const, side, trial_name)
    global column_achilles column_gonio % noraxonfreq
    global plot_norm plot_us plot_check subject_id % plot_achilles  plot_emg 
    global filepath
    
    
    
    % Read stiffness trial US data file, determine time stamps, set trigger frame as time = zero
    % Produce US sample frequency + new US array containing time and displacement
    [usdata_prepped,usfreq] = read_us_file(strcat(filepath, usdata), str2double(usdata_frame), trial_name);
    
    % Read stiffness trial Noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
    % Produce a new noraxon data array
    noraxon_prepped = read_noraxon_stiffness(strcat(filepath, noraxondata), usfreq, side, trial_name);
        
    
    
    % secondary zero offset for GONIOMETER data (primary offset during data collection, only before 1st of 3 trials)
    gonio_offset = mean(noraxon_prepped(1:3,column_gonio));
    noraxon_prepped(:,column_gonio) = noraxon_prepped(:,column_gonio) - gonio_offset;
    
    
    
    % Correct TORQUE for co-activation
    noraxon_prepped_coact = correct_coactivation(noraxon_prepped, coact_max_torque, coact_max_EMG, side, trial_name);
    
    % Calculate force by applying internal moment arm
    tendon_force = noraxon_prepped_coact(:,column_achilles) / at_momentarm;
    
    % secondary zero offset for force
    tendon_force_avg = mean(tendon_force(1:5)); % avg of first 5 frames. Previously: 1/20th sec (0,05 s) 
    tendon_force_offset = tendon_force - tendon_force_avg;
    
    % output raw max force from trial
    maxforce = max(tendon_force_offset);
    
    
    
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
    
    

    
    % Correct displacement for ankle rotation
    commonlength = min(length(usdata_prepped),length(noraxon_prepped));
    usdata_corrected(:,1) = usdata_prepped(1:commonlength,1); % TIME column
    if at_rotation_const == 0
        % for trials not to be corrected for rotation, 0 is passed (2013, no correction for MTJ - june 2014, no correction for OTJ)
        usdata_corrected(:,2) = usdata_prepped(1:commonlength,2); % use DISPLACEMENT column, without applying correction
    else
        % for trials to be corrected for rotation, a rotation constant is passed and applied as follows
        rotation_correction = noraxon_prepped(1:commonlength,column_gonio) * at_rotation_const;
        usdata_corrected(:,2) = usdata_prepped(1:commonlength,2) + rotation_correction; % use DISPLACEMENT column, WITH applying correction
        % plot US data corrections, both for extmark and rotation
        if plot_us
            plottitle = horzcat('US data correction check for ', subject_id, ' ', trial_name);
            figure('Name',plottitle)
            plot(usdata_prepped(:,1),usdata_prepped(:,3),'r')
            hold on
            plot(usdata_prepped(:,1),usdata_prepped(:,4),'g')
            plot(usdata_corrected(:,1),rotation_correction, 'c')
            plot(usdata_corrected(:,1),usdata_corrected(:,2), 'k','LineWidth',2)
            xlabel('Time (s)'),ylabel('Displacement (mm)'),title(plottitle);
            legend('Orig displacement','Displ. of extmark','Displ. due to rotation','Corrected displ.','Location','Northwest');
        end
    end
    
    
    
    
    
    
    
    
    % Plot synchronization check US vs Norm
    if plot_check && plot_norm
        plottitle = horzcat('SYNC check for ', subject_id, ' ', trial_name);
        figure('Name',plottitle);
        [AX,H1,H2] = plotyy(usdata_corrected(:,1),usdata_corrected(:,2),noraxon_prepped(:,1),tendon_force_offset,'plot');
        hold on
        plot(noraxon_prepped(:,1),noraxon_prepped(:,column_gonio),'m')
%        %%% Force cutting method 2, april 2004
%        if strcmp(trial_name(1),'M') == 1 % MTJ trial - plot cutoff
%            plot([time1 time1], [min(usdata_prepped(:,2)) max(usdata_prepped(:,2))],'k');
%        end
        set(get(AX(1),'Ylabel'),'String','Displacement (mm)')
        set(get(AX(2),'Ylabel'),'String','Force (N)')
        set(AX,{'ycolor'},{'r';'b'})
        set(H1,'color','red')
        set(H2,'color','blue')
        set(AX(1),'YLim',[min(usdata_corrected(:,2)) max(usdata_corrected(:,2))])
        set(AX(1),'YTick',0:1:max(usdata_corrected(:,2)))
        set(AX(1),'box','off')
        set(AX(2),'YLim',[min(tendon_force_offset) max(tendon_force_offset)])
        set(AX(2),'YTick',0:500:max(tendon_force_offset))
        xlabel('Time (s)'),title(plottitle);
        legend('Displ_corr (mm)','Gonio (deg)','Force (N)','Location','Southeast');
    end


    
    
    % output raw force-displacement columns to Excel
    if ispc
        warning('off', 'MATLAB:xlswrite:AddSheet');
        xlswrite(strcat('data_output/stiff_forceoutput_', subject_id, '.xls'),{'Time' 'Displ orig' 'Gonio' 'Displ corr' 'Torque orig' 'Torque coact' 'Force'}, trial_name,'A1')
        xlswrite(strcat('data_output/stiff_forceoutput_', subject_id, '.xls'),noraxon_prepped(:,1), trial_name,'A2')
        xlswrite(strcat('data_output/stiff_forceoutput_', subject_id, '.xls'),usdata_prepped(:,2), trial_name,'B2')
        xlswrite(strcat('data_output/stiff_forceoutput_', subject_id, '.xls'),noraxon_prepped(:,column_gonio), trial_name,'C2')
        xlswrite(strcat('data_output/stiff_forceoutput_', subject_id, '.xls'),usdata_corrected(:,2), trial_name,'D2')
        xlswrite(strcat('data_output/stiff_forceoutput_', subject_id, '.xls'),noraxon_prepped(:,column_achilles), trial_name,'E2')
        xlswrite(strcat('data_output/stiff_forceoutput_', subject_id, '.xls'),noraxon_prepped_coact(:,column_achilles), trial_name,'F2')
        xlswrite(strcat('data_output/stiff_forceoutput_', subject_id, '.xls'),tendon_force_offset, trial_name,'G2')
    else
        %output_csv_head = {'Time' 'Torque orig' 'Torque coact' 'Force' 'Displ orig' 'Gonio' 'Displ corr'};
        output_csv(:,1) = noraxon_prepped(:,1);
        output_csv(:,2) = usdata_prepped(:,2);
        output_csv(:,3) = noraxon_prepped(:,column_gonio);
        output_csv(:,4) = usdata_corrected(:,2);
        output_csv(:,5) = noraxon_prepped(:,column_achilles);
        output_csv(:,6) = noraxon_prepped_coact(:,column_achilles);
        output_csv(:,7) = tendon_force_offset;
        csvwrite(strcat('data_output/stiff_forceoutput_', subject_id, '_', trial_name, '.csv'),output_csv)
    end
        
    % Create outout array: time, force, displacement
    len_array = min(length(tendon_force_offset),length(usdata_prepped));
    time_force_displ_array(:,1) = usdata_prepped(1:len_array,1); % time array
    time_force_displ_array(:,2) = tendon_force_offset(1:len_array); % force array
    time_force_displ_array(:,3) = usdata_corrected(1:len_array,2); % displacement array
end