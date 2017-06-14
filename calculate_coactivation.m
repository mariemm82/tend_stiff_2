%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_coactivation / TIBIALIS ANTERIOR ACTIVATION
% Marie Moltubakk 7.6.2013
% Read complete, prepared noraxon array + number of frames to average (freq * time)
% Produce max torque, max EMG constants
%
% Assumes that the file sent to this method, is the best MVC trial
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [coact_max_torque,coact_max_EMG] = calculate_coactivation(noraxon_prepped, average_frames, side)
    global plot_emg plot_check subject_id
    global column_l_tibant column_r_tibant column_achilles % column_EMG_start column_EMG_end column_l_gm column_r_gm column_l_gl column_r_gl column_l_sol column_r_sol column_norm_angle column_norm_torque column_norm_velocity column_norm_direction 
    
    
    
    % invert torque column, to report dorsiflexion torque as positive
    noraxon_prepped(:,column_achilles) = -noraxon_prepped(:,column_achilles);
    
    
    
    %%% identify highest torque value (average over preset duration)
    averaged_torque = avg_array(average_frames, noraxon_prepped(:,column_achilles));
    % find max EMG and array location of max EMG
    [coact_max_torque, index_torque] = max(averaged_torque(:));
    
    
    
    %%% identify highest EMG RMS value (averate over preset duration)
    
    % select side
    if strcmpi(side,'L') == 1
        column_tibant = column_l_tibant;
    else
        column_tibant = column_r_tibant;
    end
        
    % average array, typically over 500 ms
    averaged_emg = avg_array(average_frames, noraxon_prepped(:,column_tibant));
    % find max EMG and array location of max EMG
    [coact_max_EMG, index_EMG] = max(averaged_emg(:));

    
    
    %%% checkpoint plot
    if plot_check && plot_emg
        plottitle = horzcat('Dorsiflexion MVC averaging check for ', subject_id);
        figure('Name', plottitle)
        plot(noraxon_prepped(:,1), noraxon_prepped(:,column_tibant),'r')
        hold on
        plot(noraxon_prepped(:,1), averaged_emg,'b')
        plot(noraxon_prepped(:,1), noraxon_prepped(:,column_achilles),'g')
        plot(noraxon_prepped(:,1), averaged_torque,'k')
        % vertical lines for placement of max EMG and max torque in array
        plot ([noraxon_prepped(index_EMG,1) noraxon_prepped(index_EMG,1)], [min(noraxon_prepped(:,column_tibant)) max(noraxon_prepped(:,column_tibant))],'b')
        plot ([noraxon_prepped(index_torque,1) noraxon_prepped(index_torque,1)], [min(noraxon_prepped(:,column_tibant)) max(noraxon_prepped(:,column_tibant))],'k')
        xlabel('Time (s)'),ylabel('Data (misc)'),title(plottitle);
        legend('orig EMG','avg EMG','orig torque','avg torque','Location','Northeast');
    end
end




function outputarray = avg_array(window_frames, data_array)
    outputarray(1:length(data_array)) = zeros;
    datalength = length(data_array);
    for i = 1:datalength
        if i <= (window_frames/2) % first X values
            outputarray(i) = mean(data_array(1:i+(window_frames/2))); % first output is value 1 until next 75 values
        elseif i > (datalength-(window_frames/2)) % last X values
            outputarray(i) = mean(data_array(i-(window_frames/2):length(data_array)));
        else % all middle values, where sufficient data is available before and after
            outputarray(i) = mean(data_array(i-(window_frames/2):i+(window_frames/2)));
        end
    end
end