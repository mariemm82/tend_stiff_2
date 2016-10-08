%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_us_file
% Marie Moltubakk 17.5.2013
% Read US data file, determine time stamps, set trigger frame as time = zero
% Produce US sample frequency, create new US array containing time and displacement
%    usdata_prepped = [time,track_corrected,track_feature,track_extmark];
%%%%%%%%%%%%%%%%%%%%%%%%%%

% NB: Possible rounding off error of time stamp, of  -0,00003 s, as seen when
% comparing timeseries.time and resampled time column in timeseries.data.
% The error does not accumulate/get wors with increasing time.

function [usdata_prepped,usfreq] = read_us_file(usfile, usframe, trial_name)
    global us_zerodispframes % noraxonfreq emg_bandpass emg_rms_ms mvc_window_ms torque_cutoff convert_achilles convert_norm
    % global plot_achilles plot_norm plot_emg plot_check plot_us subject_id

    % import us data
    usdata = importdata(usfile, '\t', 2);



    %%% MANIPULATE TIME/DISPLACEMENT DATA

    % delete frames before trigger click
    % NB: 
    %   US video has 395 frames.
    %   Tracker labels as 0 to 394. Frame 0 = time 0.00, frame 394 = time 21.773
    %   Lichtwark reports 395 data points. First point = time 0.0553, last point = time 21.828
    %   Conclusion: Triggerframe = 0 means no data cutoff, start using line 1 (usframe+1). But Lichtwark timestamps should be lowered by 1 unit.
    %   Norm data start recording at triggerframe, after that, higher sampling frequency and longer duration of recording than the US video
    usdata.data = usdata.data(usframe+1:end,1:6);
    
    % delete last frames, if extmark and feat are not same length
    time_feat = usdata.data(:,1);
    time_ext = usdata.data(:,4);
    time_feat(isnan(time_feat)) = [];
    time_ext(isnan(time_ext)) = [];
    newlength = min(length(time_feat),length(time_ext));
        
    % extract arrays to be used
    time = usdata.data(1:newlength,1);
    track_feature = -usdata.data(1:newlength,2); % PROJECTSPECIFIC
    track_extmark = -usdata.data(1:newlength,5); % PROJECTSPECIFIC

    % set trigger frame as time zero
    time(:) = time(:) - time(1);

    % transform tracks from position to displacement
    track_feature(:) = track_feature(:) - mean(track_feature(1:us_zerodispframes));
    track_extmark(:) = track_extmark(:) - mean(track_extmark(1:us_zerodispframes));

    % correct displacement according to external mark
    track_corrected = track_feature - track_extmark; % NB, no smoothing of US track and extmerk track applied

    % 2014: commented out, US plots will be done for MTJ and OTJ at a later
    % stage in the script. Turn on here for US plots for CPM etc.
%     % plot to verify displacement data
%     if plot_check && plot_us
%         plottitle = horzcat('US data check for ', subject_id, ' ', trial_name);
%         figure('Name',plottitle)
%         plot(time,track_feature,'r')
%         xlabel('Time (s)'),ylabel('Displacement (mm)'),title(plottitle);
%         hold on
%         plot(time,track_extmark,'g')
%         plot(time,track_corrected,'k','LineWidth',2)
%         legend('raw displacement','external marker','corrected displacement','Location','SouthEast');
%     end

    usdata_prepped = [time,track_corrected,track_feature,track_extmark];



    %%% CALCULATE US SAMPLING FREQUENCY
    usfreq = 100/(time(101)-time(1));
end