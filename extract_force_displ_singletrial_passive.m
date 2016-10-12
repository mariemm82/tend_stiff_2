%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract_force_displ_singletrial_passive
% Marie Moltubakk 11.6.2013
% Read complete, prepared noraon array + prepared us array
% Produce array with time, corrected force, corrected displacement
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [force,gonio,angle,displacement,time_us] = extract_force_displ_singletrial_passive(noraxondata, usdata, usdata_frame, max_EMG_TA, max_EMG_GM, max_EMG_GL, max_EMG_SOL, leg_length, side, line, trial_name)
    
    global column_EMG_start column_EMG_end column_l_gm column_r_gm column_l_gl column_r_gl column_l_sol column_r_sol column_l_tibant column_r_tibant column_gonio column_norm_angle column_norm_torque column_norm_velocity column_norm_direction column_achilles
    global plot_norm plot_check subject_id
    global at_momentarm
    global filepath
    
    
    
    % Read US data file, determine time stamps, set trigger frame as time = zero
    % Produce US sample frequency + new US array containing time and displacement
    [usdata_prepped,usfreq] = read_us_file(strcat(filepath, usdata, '.txt'), str2double(usdata_frame), trial_name);
    
    % Read Noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
    % Produce a new noraxon data array
    noraxon_prepped = read_noraxon_passive(strcat(filepath, noraxondata), usfreq, side, trial_name);
    
    
   
    
    
    %%% extract data from big arrays
    % find point where passive torque value is artificially set to zero (= pause at max angle)
    len = find(noraxon_prepped(:,column_norm_torque)==0,1,'first') - 1;
    % if US data does not have the full length of the Norm data:
    if len > length(usdata_prepped)
        cprintf('err',horzcat('WARNING: ', trial_name, 'US data are shorter than Norm data: ', num2str(len), ' Norm frames, ', num2str(length(usdata_prepped)), ' US frames.\n'));
        len = length(usdata_prepped);
    end
    displacement = -usdata_prepped(1:len,2);
    torque = noraxon_prepped(1:len,column_norm_torque);
    % turn dorsiflexion angles from negative to positive for plotting purposes
    angle = -noraxon_prepped(1:len,column_norm_angle);
    gonio = -noraxon_prepped(1:len,column_gonio);
    time = noraxon_prepped(1:len,1);
    time_us = usdata_prepped(1:len,1);
    
    % calculate % EMG activation
    if strcmpi(side,'R') == 1
         emg_GM = noraxon_prepped(1:len,column_r_gm)/max_EMG_GM*100; % percent
         emg_GL = noraxon_prepped(1:len,column_r_gl)/max_EMG_GL*100; % percent
         emg_SOL = noraxon_prepped(1:len,column_r_sol)/max_EMG_SOL*100; % percent
    else % left
         emg_GM = noraxon_prepped(1:len,column_l_gm)/max_EMG_GM*100; % percent
         emg_GL = noraxon_prepped(1:len,column_l_gl)/max_EMG_GL*100; % percent
         emg_SOL = noraxon_prepped(1:len,column_l_sol)/max_EMG_SOL*100; % percent
    end
    cprintf('red', horzcat(trial_name, ' EMG max activation: gm ', num2str(round(max(emg_GM),3)), ' %%, gl ', num2str(round(max(emg_GL),3)), ' %%, sol ', num2str(round(max(emg_SOL),3)), ' %%.\n'));
    
    % convert torque to force
    force = torque / at_momentarm;
    
    % Plot synchronization check US vs Norm
    if plot_check && plot_norm
        plottitle = horzcat('SYNC check for ', subject_id, ' ', trial_name);
        figure('Name',plottitle);
        [AX,H1,H2] = plotyy(time,force/20,time_us,displacement,'plot');
        hold on
        plot(time,gonio,'k')
        plot(time,angle,'g')
        plot(time,emg_GM,'y')
        plot(time,emg_GL,'m')
        plot(time,emg_SOL,'c')
        set(get(AX(2),'Ylabel'),'String','Displacement (mm)')
        set(get(AX(1),'Ylabel'),'String','Force (N)/10 / angle (deg) / EMG')
        set(AX,{'ycolor'},{'b';'r'})
        set(H2,'color','red')
        set(H1,'color','blue')
        set(AX(2),'YLim',[min(displacement) max(displacement)])
        set(AX(2),'YTick',[0:1:1.2*max(displacement)])
        set(AX(2),'box','off')
        set(AX(1),'YLim',[0 max(force)/20])
        set(AX(1),'YTick',[0:round(max(force)/20/10,-1):1.2*max(force)/20])
        set(AX(1),'box','off')
        xlabel('Time (s)'),title(plottitle);
        legend('Force (N) /20','Gonio (deg)','Norm angle (deg)', '%EMG GM', '%EMG GL', '%EMG SOL','Displ (mm)','Location','Northwest');
    end


    
    
% 
% ----------------------- OUTPUT RAW DATA? ------------------
%     
%     % output raw force-displacement columns to Excel
%     if ispc
%         warning('off', 'MATLAB:xlswrite:AddSheet');
%         xlswrite(strcat('data_output/forceoutput_', subject_id, '.xls'),{'Time' 'Displ orig' 'Gonio' 'Displ corr' 'Torque orig' 'Torque coact' 'Force'}, trial_name,'A1')
%         xlswrite(strcat('data_output/forceoutput_', subject_id, '.xls'),noraxon_prepped(:,1), trial_name,'A2')
%         xlswrite(strcat('data_output/forceoutput_', subject_id, '.xls'),usdata_prepped(:,2), trial_name,'B2')
%         xlswrite(strcat('data_output/forceoutput_', subject_id, '.xls'),noraxon_prepped(:,column_gonio), trial_name,'C2')
%         xlswrite(strcat('data_output/forceoutput_', subject_id, '.xls'),usdata_corrected(:,2), trial_name,'D2')
%         xlswrite(strcat('data_output/forceoutput_', subject_id, '.xls'),noraxon_prepped(:,column_achilles), trial_name,'E2')
%         xlswrite(strcat('data_output/forceoutput_', subject_id, '.xls'),noraxon_prepped_coact(:,column_achilles), trial_name,'F2')
%         xlswrite(strcat('data_output/forceoutput_', subject_id, '.xls'),tendon_force_offset, trial_name,'G2')
%     else
%         output_csv_head = {'Time' 'Torque orig' 'Torque coact' 'Force' 'Displ orig' 'Gonio' 'Displ corr'};
%         output_csv(:,1) = noraxon_prepped(:,1);
%         output_csv(:,2) = usdata_prepped(:,2);
%         output_csv(:,3) = noraxon_prepped(:,column_gonio);
%         output_csv(:,4) = usdata_corrected(:,2);
%         output_csv(:,5) = noraxon_prepped(:,column_achilles);
%         output_csv(:,6) = noraxon_prepped_coact(:,column_achilles);
%         output_csv(:,7) = tendon_force_offset;
%         csvwrite(strcat('data_output/forceoutput_', subject_id, '_', trial_name, '.csv'),output_csv)
%     end
end