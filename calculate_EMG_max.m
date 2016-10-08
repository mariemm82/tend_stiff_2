%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_EMG_max
% Marie Moltubakk 7.6.2013
% Read complete, prepared noraxon array + number of frames to average (freq * time)
% Produce max torque, max EMG constants
%
% Assumes that the file sent to this method, is the best MVC trial
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [max_torque,max_EMG] = calculate_EMG_max(noraxon_prepped, average_frames, column_emg, invert_DF)
    global plot_emg plot_check subject_id
    global column_achilles
    
    
    
    if invert_DF
        % invert torque column, to report dorsiflexion torque as positive
        noraxon_prepped(:,column_achilles) = -noraxon_prepped(:,column_achilles);
    end
    
    
    
    %%% identify highest torque value (average over preset duration)
    averaged_torque = avg_array(average_frames, noraxon_prepped(:,column_achilles));
    % find max EMG and array location of max EMG
    [max_torque, index_torque] = max(averaged_torque(:));
    
    
    
    %%% identify highest EMG RMS value (averate over preset duration)
    % average array, typically over 500 ms
    averaged_emg = avg_array(average_frames, noraxon_prepped(:,column_emg));
    % find max EMG and array location of max EMG
    [max_EMG, index_EMG] = max(averaged_emg(:));

    
    
    %%% checkpoint plot
    if plot_check && plot_emg
        plottitle = horzcat('MVC averaging check, ', subject_id, ' muscle ', num2str(column_emg));
        figure('Name', plottitle)
        plot(noraxon_prepped(:,1), noraxon_prepped(:,column_emg),'r')
        hold on
        plot(noraxon_prepped(:,1), averaged_emg,'b')
        plot(noraxon_prepped(:,1), noraxon_prepped(:,column_achilles),'g')
        plot(noraxon_prepped(:,1), averaged_torque,'k')
        % vertical lines for placement of max EMG and max torque in array
        plot ([noraxon_prepped(index_EMG,1) noraxon_prepped(index_EMG,1)], [min(noraxon_prepped(:,column_emg)) max(noraxon_prepped(:,column_emg))],'b','linewidth',2) % EMG vertical line 
        plot ([noraxon_prepped(index_torque,1) noraxon_prepped(index_torque,1)], [min(noraxon_prepped(:,column_achilles)) max(noraxon_prepped(:,column_achilles))],'k','linewidth',2) % torque vertical line
        xlabel('Time (s)'),ylabel('Data (misc)'),title(plottitle);
        legend('orig EMG','smooth+max EMG','orig torque','smooth+max torque','Location','Northeast');
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