%%%%%%%%%%%%%%%%%%%%%%%%%%
% correct_coactivation
% Marie Moltubakk 9.6.2013
%%%%%%%%%%%%%%%%%%%%%%%%%%

function noraxon_prepped_coact = correct_coactivation(noraxon_prepped, coact_max_torque, coact_max_EMG, side, trial_name)
    global mute
    global plot_achilles plot_emg plot_check subject_id % plot_norm 
    global column_l_tibant column_r_tibant column_achilles % column_EMG_start column_EMG_end column_l_gm column_r_gm column_l_gl column_r_gl column_l_sol column_r_sol column_norm_angle column_norm_torque column_norm_velocity column_norm_direction 

    % select side
    if strcmpi(side,'L') == 1
        column_tibant = column_l_tibant;
    else
        column_tibant = column_r_tibant;
    end
    
    % calculate TA activation, in percent of MVC, per frame
    percent_activation = noraxon_prepped(:,column_tibant) / coact_max_EMG * 100;
    
    % calculate torque resulting from TA activation, per frame
    TA_torque = percent_activation/100 * coact_max_torque;
    
    % correct external recorded torque (add TA torque)
    corrected_torque = noraxon_prepped(:,column_achilles) + TA_torque;
    
    % add co-activation-corrected achilles torque back into noraxon array
    noraxon_prepped_coact = noraxon_prepped;
    noraxon_prepped_coact(:,column_achilles) = corrected_torque;

    % output coactivation numbers to screen, as text
    if mute == 0
        if strcmp(trial_name,'MTJ1')
            cprintf(horzcat('TA co-activation report: MVC dorsiflexion = ', num2str(coact_max_torque,4), ' Nm. EMG = ', num2str(coact_max_EMG,4),' µV.\n'));
        end
        cprintf(horzcat('    ', trial_name, ' TA max activ. = ', num2str(max(percent_activation),3), '%% -> max torque contrib. = ', num2str(max(TA_torque),3), ' Nm.\n'));
    end
    
    %%% checkpoint plot
    if plot_check && plot_achilles && plot_emg
         plottitle = horzcat('Co-activation correction check for ', subject_id, ' ', trial_name);
         figure('Name',plottitle)
         plot(noraxon_prepped(:,1), noraxon_prepped(:,column_achilles),'b')
         hold on
         plot(noraxon_prepped(:,1), TA_torque,'g')
         plot(noraxon_prepped(:,1), noraxon_prepped_coact(:,column_achilles),'k','LineWidth',2)
         xlabel('Time (s)'),ylabel('Torque (Nm)'),title(plottitle);
         legend('orig ramp torque','coactivation torque','corrected torque','Location','Northwest');
     end

    
end