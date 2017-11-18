%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% read noraxon data file
% Marie Moltubakk 18.5.2013
% Read noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
% Produce a new noraxon data array
% version 21.1.2015, adapted to use individually calculated norm conversion factors

% used by create_angles through extract_force_displ_singletrial_passive.m
% used by passive analyses through extract_force_displ_singletrial_passive_EMG.m
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [noraxon_resampled1, noraxon_resampled2] = read_noraxon_knee_2trials(noraxonfile, finalfreq, side, trial_name)

    %% prepare & import Norm
    global noraxonfreq 
    global angle_cutoff velocity_cutoff torque_cutoff_active % torque_cutoff_bandstop angle_cutoff_active velocity_cutoff_active
    global plot_individual subject_id
    global column_norm_angle column_norm_torque column_norm_velocity column_norm_direction % column_l_tibant column_r_tibant column_achilles
    global convert_norm_angle_a convert_norm_angle_b convert_norm_torque_a convert_norm_torque_b convert_norm_velocity_a convert_norm_velocity_b % convert_norm_direction_b
    
    % import noraxon data
    noraxondata = importdata(noraxonfile, '\t', 1);

    % reset time array to start from 0,0 (trigger click = time zero)
    noraxondata.data(:,1) = noraxondata.data(:,1) - noraxondata.data(1,1);
    time = noraxondata.data(:,1);
    
         
    %% PREPARE AND FILTER ANGLE AND DIRECTION, IDENTIFY ZERO ANGLE

    % convert angle
    if side == 'R'
        angleadjust1 = -180;
    else
        angleadjust1 = 180;
    end
    norm_angle = angleadjust1 + (noraxondata.data(:,column_norm_angle) * convert_norm_angle_a) + convert_norm_angle_b;
    if max(norm_angle) > 340
        norm_angle = norm_angle - 360;
    elseif max(norm_angle) < -180
        norm_angle = norm_angle + 360;
    end
        
    % filter angle
    [B, A] = butter(4, angle_cutoff, 'low');
    norm_angle_filtered = filtfilt(B, A, norm_angle);
    
    
    % prepare direction
    norm_direction = noraxondata.data(:,column_norm_direction);

    
        
    %% PREPARE AND FILTER VELOCITY

    % convert velocity, positive velocity = towards dorsiflexion for both sides
    norm_velocity = -((noraxondata.data(:,column_norm_velocity) * convert_norm_velocity_a) + convert_norm_velocity_b);
    
    % filter velocity
    [B, A] = butter(4, velocity_cutoff, 'low');
    norm_velocity_filtered = filtfilt(B, A, norm_velocity);
    
    % make a rough secondary correction for velocity offset
    % (velocity is only used to determine phases):
    
    % find point where direction changes (=max angle)
    avoidframes = 200; %VAR, # of frames skipped around the point of change in direction, to avoid noise due to the change

    if strcmpi(side,'R') == 1
        %R: the direction data changes only at the END of the stop phase
        % start by finding a "high" direction value (trial onset)
        loc_onset = find(norm_direction>4000,1,'first');
        % find first "low" direction value
        loc_maxROM = loc_onset + find(norm_direction(loc_onset:end)<4000,1,'first');
        % if detected point is just a small spike (value not stable low after 50 frames), try again:
        loc_maxROMcheck = loc_maxROM+50;
        while mean(norm_direction(loc_maxROMcheck:loc_maxROMcheck+20)) > 0
            % find next spike
            loc_maxROMnext = loc_maxROMcheck + find(norm_direction(loc_maxROMcheck:end)<4000,1,'first');
            loc_maxROM = loc_maxROMnext;
            loc_maxROMcheck = loc_maxROM+50;
        end
        % calculate average velocity during one second BEFORE change in direction
        % + avoiding 100 last frames before the change
        zero_velocity = mean(norm_velocity_filtered(loc_maxROM-avoidframes-noraxonfreq/2:loc_maxROM-avoidframes));
    else % left side
        %L: the direction data changes at the BEGINNING of the stop phase
        % find first "high" direction value
        loc_maxROM = find(norm_direction>100,1,'first');
        while mean(norm_direction(loc_maxROM+avoidframes:loc_maxROM+avoidframes+noraxonfreq/2)) < 4500
            % not long enough - find next high point
            loc_maxROM = loc_maxROM+50 + find(norm_direction(loc_maxROM+50:end) > 100, 1,'first');
            % when high value stays long enough for calculating:
        end
        % calculate average velocity during one second AFTER change in direction
        % + avoiding 200 frames at the beginning
        zero_velocity = mean(norm_velocity_filtered(loc_maxROM+avoidframes:loc_maxROM+avoidframes+noraxonfreq/2));
    end
    
    % offset: setting velocity at stop position (maxROM) as zero:
    norm_velocity_filtered = norm_velocity_filtered-zero_velocity;
    % corresponding millivolt correction is:
    zero_velocity_mv = (zero_velocity - convert_norm_velocity_b) * convert_norm_velocity_a;
    % print velocity correction report
    cprintf('magenta', horzcat('Secondary velocity correction: Stop position intially had ', num2str(zero_velocity), ' deg/s --> secondary offset ', num2str(zero_velocity_mv), ' mV.\n'));
    
    
    %% REPLACE ANGLE, DIRECTION, VELOCITY IN DATA

    % swap direction if side is R (should have been done immediately, but
    % some code above uses direction AND corrects if R)
    if side == 'R'
        norm_direction = -noraxondata.data(:,column_norm_direction) + 4000;
    end
    
    % replace angle and direction
    noraxondata.data(:,column_norm_angle) = norm_angle_filtered;
    noraxondata.data(:,column_norm_direction) = norm_direction;
    noraxondata.data(:,column_norm_velocity) = norm_velocity_filtered;

    
    %% DETECT ASCENDING PHASE - DESCENDING PHASE
    
    veloc1 = 1.4; % not 1.7 because some subjects have some spikes dropping LOWER, not higher
    veloc2 = -1.7;
    buffer = finalfreq * 3; % 100 hz * 3 sec = minimum duration of phase
    
    norm_velocity_movingmean = movingmean(norm_velocity_filtered, 400);
    
    % ascending phase ongoing
    changedirection0 = find(norm_velocity_movingmean(1:end) > veloc1, 1, 'first');
    if changedirection0
    else
        changedirection0 = 1;
    end
    
    % find point where ascending phase stops
    changedirection1 = changedirection0 + buffer + find(norm_velocity_movingmean(changedirection0+buffer:end) < 0.5*veloc1, 1, 'first');

    % find point where decending phase starts
    changedirection2 = changedirection1 + find(norm_velocity_movingmean(changedirection1:end) < veloc2, 1, 'first');
    
    % find point where decending phase stops
    change3slower = find(norm_velocity_movingmean(changedirection2+buffer:end) > 0.5*veloc2, 1, 'first');
    change3faster = find(norm_velocity_movingmean(changedirection2+buffer:end) < 2*veloc2, 1, 'first');
    change3find = min([change3slower change3faster]);
    if change3find
        changedirection3 = changedirection2 + buffer + change3find;
    else % velocity steady - use trial until end
        changedirection3 = length(norm_velocity_movingmean);
    end

    % plot to confirm phases of direction changes
    if plot_individual
        plottitle = horzcat('Passive phases check for ', subject_id, ' ', trial_name, '_pt1');
        figure('Name',plottitle)
        plot(time,-norm_direction/100)
        hold on
        plot(time,(norm_angle_filtered),'r') 
        if strcmpi(side,'R') == 1
            plot(time,(noraxondata.data(:,column_norm_torque))/4,'g') %scaled
        else
            plot(time,-(noraxondata.data(:,column_norm_torque))/4,'g') %scaled
        end
        plot(time,(norm_velocity_filtered*10),'c') %scaled
        
        plot ([time(changedirection0) time(changedirection0)], [-100 100],'y')
        plot ([time(changedirection1) time(changedirection1)], [-100 100],'y')
        plot ([time(changedirection2) time(changedirection2)], [-100 100],'y')
        plot ([time(changedirection3) time(changedirection3)], [-100 100],'y')
        
        plot(time,(norm_velocity_movingmean*10),'b') %scaled

        axis([-Inf Inf -60 Inf])
        xlabel('Time (s)')
        ylabel('Misc scaled variables')
        title(plottitle,'Interpreter', 'none')
        legend('Direction/100','Angle','Torque, µV/4','Velocity*10','direction changes','location','northeast')
        print(horzcat('data_plots/IND knee ',plottitle),'-dpng')
    end

    
    
    
    
    
    
    %% UNTIL THIS POINT: TREATING 2 TRIALS TOGETHER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ... but "detecting phases" detects only for first trial
    % first trial ends at changedirection3
    % until here, only manipulating "noraxondata.data"
    
    
    %% AFTER THIS POINT: TREATING ONLY 1ST TRIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    
    
    
    
    
    
    %% PREPARE AND FILTER TORQUE 
    
    % for ascending and descending phase separately:
    
    % convert torque to give positive values for plantar flexion for both sides
    if strcmpi(side,'R') == 1
        norm_torque_ascend = (noraxondata.data(changedirection0:changedirection1,column_norm_torque) * convert_norm_torque_a) + convert_norm_torque_b; 
        norm_torque_descend = (noraxondata.data(changedirection2:changedirection3,column_norm_torque) * convert_norm_torque_a) + convert_norm_torque_b; 
    else % left leg to plantar flexion = positive
        norm_torque_ascend = -((noraxondata.data(changedirection0:changedirection1,column_norm_torque) * convert_norm_torque_a) + convert_norm_torque_b); 
        norm_torque_descend = -((noraxondata.data(changedirection2:changedirection3,column_norm_torque) * convert_norm_torque_a) + convert_norm_torque_b); 
    end

    % filter torque lowpass
    [B, A] = butter(4, torque_cutoff_active, 'low');
    norm_torque_ascend_lowpass = filtfilt(B, A, norm_torque_ascend);
    norm_torque_descend_lowpass = filtfilt(B, A, norm_torque_descend);
    
    % prepare arrays containing times for ascending and descending phase 
    time_torque_ascend = noraxondata.data(changedirection0:changedirection1,1);
    time_torque_descend = noraxondata.data(changedirection2:changedirection3,1);
    
    
    
    % ascending phase:
    
     %%% method 3: spline for torque, from Vidar 2014-05-06
     % spap2 = Least-squares spline approximation 
     % NB, not an array but a struct.....
     w = ones(size(time_torque_ascend)); % w = weights in the error measure
     w([1 end]) = 100; % 100 - forces fit to come very close to the first and last data point    %VAR
     var_a = 1;    %VAR
     var_b = 4; % 2 means least squares straight line, 4 means linear...    %VAR
    
     norm_torque_ascend_spap2 = spap2(var_a, var_b, time_torque_ascend, norm_torque_ascend_lowpass, w);

     %%% method 4: cubic smoothing spline
     p = 0.0005; % 0 = least squares straight line fit, 1 = natural cubic spline interpolant   %VAR
     norm_torque_ascend_csaps = csaps(time_torque_ascend, norm_torque_ascend_lowpass, p, time_torque_ascend); 
     
     %%% method 5: smoothing with lowess (UIO help)
     % lowess = local regression using weighted linear least squares, 1st degree poly model
     % loess = ... 2nd degree poly model
     % rlo(w)ess = robust, lower weight to outliers in the regresion
 %    var_span = 0.9; % (r)lo(w)ess - percent span     %VAR
 %    norm_torque_ascend_smooth = smooth(norm_torque_ascend_lowpass, var_span, 'lowess');    
     
     % chosen filter for futher calculations
     norm_torque_ascend_filtered = norm_torque_ascend_csaps;
    
    if plot_individual
         plottitle = horzcat('Torque filter check ASCEND for ', subject_id, ' ', trial_name, '_pt1');
         figure('Name',plottitle)
         plot(time_torque_ascend,norm_torque_ascend_lowpass,'color','r','linewidth',2)
         hold on
         fnplt(norm_torque_ascend_spap2,'linewidth',1,'color','g')
         plot(time_torque_ascend,norm_torque_ascend_csaps,'b','linewidth',2)
  %       plot(time_torque_ascend,norm_torque_ascend_smooth,'c')
         xlabel('Time (s)')
         ylabel('Torque (Nm)')
         title(plottitle,'Interpreter', 'none')
         legend('initial lowpass','spap2 least-squares spline','csaps cubic smoothing spline','location', 'northwest')
%           print(horzcat('data_plots/IND knee ',plottitle),'-dpng')
    end



    % descending phase:
    
     %%% method 3: spline for torque, from Vidar 2014-05-06
     % spap2 = Least-squares spline approximation 
     % NB, not an array but a struct.....
     w = ones(size(time_torque_descend)); % w = weights in the error measure
     w([1 end]) = 100; % 100 - forces fit to come very close to the first and last data point    %VAR
     var_a = 1;    %VAR
     var_b = 4; % 2 means least squares straight line, 4 means linear...    %VAR
     norm_torque_descend_spap2 = spap2(var_a, var_b, time_torque_descend, norm_torque_descend_lowpass, w);

     %%% method 4: cubic smoothing spline
     p = 0.0005; % 0 = least squares straight line fit, 1 = natural cubic spline interpolant   %VAR
     norm_torque_descend_csaps = csaps(time_torque_descend, norm_torque_descend_lowpass, p, time_torque_descend); 
     
     %%% method 5: smoothing with lowess (UIO help)
     % lowess = local regression using weighted linear least squares, 1st degree poly model
     % loess = ... 2nd degree poly model
     % rlo(w)ess = robust, lower weight to outliers in the regresion
%     var_span = 0.9; % (r)lo(w)ess - percent span     %VAR
 %    norm_torque_descend_smooth = smooth(norm_torque_descend_lowpass, var_span, 'lowess');    
     
     % chosen filter for futher calculations
     norm_torque_descend_filtered = norm_torque_descend_csaps;
    
     
     
     if plot_individual
         plottitle = horzcat('Torque filter check DESCEND for ', subject_id, ' ', trial_name, '_pt1');
         figure('Name',plottitle)
         plot(time_torque_descend,norm_torque_descend_lowpass,'color','r','linewidth',2)
         hold on
         fnplt(norm_torque_descend_spap2,'linewidth',1,'color','g')
         plot(time_torque_descend,norm_torque_descend_csaps,'b','linewidth',2)
   %      plot(time_torque_descend,norm_torque_descend_smooth,'c')
         xlabel('Time (s)')
         ylabel('Torque (Nm)')
         title(plottitle,'Interpreter', 'none')
         legend('initial lowpass','spap2 least-squares spline','csaps cubic smoothing spline','location', 'northeast')
%         print(horzcat('data_plots/IND knee ',plottitle),'-dpng')
     end


    %% INSERT CONVERTED TORQUE into NEW noraxon-variable containing only 1st trial
    noraxon_converted = noraxondata.data(1:changedirection3,:);

    noraxon_converted(:,column_norm_torque) = zeros;
    noraxon_converted(changedirection0:changedirection1,column_norm_torque) = norm_torque_ascend_filtered;
    noraxon_converted(changedirection2:changedirection3,column_norm_torque) = norm_torque_descend_filtered;
    
    
    %% RESAMPLE
    % create time array with old and new timestamps and length
    timearray_orig = make_timearray(noraxonfreq, size(noraxon_converted));
    timearray_new = make_timearray(finalfreq, fix(length(noraxon_converted(:,1))/noraxonfreq*finalfreq));

    % create timeseries type for original data
    noraxon_timeseries = timeseries(noraxon_converted, timearray_orig, 'Name','Noraxon');

    % resample
    noraxon_temp = resample(noraxon_timeseries,timearray_new,'linear');
    noraxon_resampled = noraxon_temp.data;
    
        
    %% CORRECT ANGLE in resampled array (first values "bent" due to filtering)
    % use values #5 to 20 of angle data, to replace values #1-4
    loc_start = 5; %VAR
    loc_stop = 20; %VAR
    angle_correction = noraxon_resampled(1:loc_stop,column_norm_angle);
    
    % calculate linear coeffisients for extrapolation
    
    frames = [ones(length(angle_correction(loc_start:loc_stop)),1) (loc_start:loc_stop)'];
    coeffs = frames\angle_correction(loc_start:loc_stop); % NB backslash --> slope

    % replace first Z values in arrays by using extrapolation
    % y = ax + b ... y = displ, x = angle
    for x = 1:(loc_start-1)
        % place back in original array
        noraxon_resampled(x,column_norm_angle) = ((x*coeffs(2)) + coeffs(1));
    end
    
    
    %% SAVE FIRST TRIAL
    noraxon_resampled1 = noraxon_resampled;
    
        
    %% REMOVE/RESET ARRAYS TO BE REUSED:
    clear noraxon_converted noraxon_resampled noraxon_temp noraxon_timeseries
    clear timearray_orig timearray_resampled
   
    
    
    
    
    


    %% REPEAT ABOVE (lines lines 123... + 190-350) for SECOND TRIAL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % second trial = noraxondata AFTER changedirection3
    % mmm behaving correctly?
    noraxondata.data(1:changedirection3-1,:) = [];
    norm_angle_filtered = noraxondata.data(:,column_norm_angle);
    norm_direction = noraxondata.data(:,column_norm_direction);
    norm_velocity_filtered = noraxondata.data(:,column_norm_velocity);
    time = noraxondata.data(:,1);
    

    %% DETECT ASCENDING PHASE - DESCENDING PHASE
    
    veloc1 = 1.4; % not 1.7 because some subjects have some spikes dropping LOWER, not higher
    veloc2 = -1.7;
    buffer = finalfreq * 3; % 100 hz * 3 sec = minimum duration of phase
    
    norm_velocity_movingmean = movingmean(norm_velocity_filtered, 400);
    
    % ascending phase ongoing
    changedirection0 = find(norm_velocity_movingmean(1:end) > veloc1, 1, 'first');
    if changedirection0
    else
        changedirection0 = 1;
    end
    
    % find point where ascending phase stops
    changedirection1 = changedirection0 + buffer + find(norm_velocity_movingmean(changedirection0+buffer:end) < 0.5*veloc1, 1, 'first');

    % find point where decending phase starts
    changedirection2 = changedirection1 + find(norm_velocity_movingmean(changedirection1:end) < veloc2, 1, 'first');
    
    % find point where decending phase stops
    change3slower = find(norm_velocity_movingmean(changedirection2+buffer:end) > 0.5*veloc2, 1, 'first');
    change3faster = find(norm_velocity_movingmean(changedirection2+buffer:end) < 2*veloc2, 1, 'first');
    change3find = min([change3slower change3faster]);
    if change3find
        changedirection3 = changedirection2 + buffer + change3find;
    else % velocity steady - use trial until end
        changedirection3 = length(norm_velocity_movingmean);
    end

    % plot to confirm phases of direction changes
    if plot_individual
        plottitle = horzcat('Passive phases check for ', subject_id, ' ', trial_name, '_pt2');
        figure('Name',plottitle)
        plot(time,-norm_direction/100)
        hold on
        plot(time,(norm_angle_filtered),'r') 
        if strcmpi(side,'R') == 1
            plot(time,(noraxondata.data(:,column_norm_torque))/4,'g') %scaled
        else
            plot(time,-(noraxondata.data(:,column_norm_torque))/4,'g') %scaled
        end
        plot(time,(norm_velocity_filtered*10),'c') %scaled
        
        plot ([time(changedirection0) time(changedirection0)], [-100 100],'y')
        plot ([time(changedirection1) time(changedirection1)], [-100 100],'y')
        plot ([time(changedirection2) time(changedirection2)], [-100 100],'y')
        plot ([time(changedirection3) time(changedirection3)], [-100 100],'y')
        
        plot(time,(norm_velocity_movingmean*10),'b') %scaled

        axis([-Inf Inf -60 Inf])
        xlabel('Time (s)')
        ylabel('Misc scaled variables')
        title(plottitle,'Interpreter', 'none')
        legend('Direction/100','Angle','Torque, µV/4','Velocity*10','direction changes','location','northeast')
        print(horzcat('data_plots/IND knee ',plottitle),'-dpng')
    end
    
    
    
    %% PREPARE AND FILTER TORQUE 
    
    % for ascending and descending phase separately:
    
    % convert torque to give positive values for plantar flexion for both sides
    if strcmpi(side,'R') == 1
        norm_torque_ascend = (noraxondata.data(changedirection0:changedirection1,column_norm_torque) * convert_norm_torque_a) + convert_norm_torque_b; 
        norm_torque_descend = (noraxondata.data(changedirection2:changedirection3,column_norm_torque) * convert_norm_torque_a) + convert_norm_torque_b; 
    else % left leg to plantar flexion = positive
        norm_torque_ascend = -((noraxondata.data(changedirection0:changedirection1,column_norm_torque) * convert_norm_torque_a) + convert_norm_torque_b); 
        norm_torque_descend = -((noraxondata.data(changedirection2:changedirection3,column_norm_torque) * convert_norm_torque_a) + convert_norm_torque_b); 
    end

    % filter torque lowpass
    [B, A] = butter(4, torque_cutoff_active, 'low');
    norm_torque_ascend_lowpass = filtfilt(B, A, norm_torque_ascend);
    norm_torque_descend_lowpass = filtfilt(B, A, norm_torque_descend);
    
    % prepare arrays containing times for ascending and descending phase 
    time_torque_ascend = noraxondata.data(changedirection0:changedirection1,1);
    time_torque_descend = noraxondata.data(changedirection2:changedirection3,1);
    
    
    
    % ascending phase:
    
     %%% method 3: spline for torque, from Vidar 2014-05-06
     % spap2 = Least-squares spline approximation 
     % NB, not an array but a struct.....
     w = ones(size(time_torque_ascend)); % w = weights in the error measure
     w([1 end]) = 100; % 100 - forces fit to come very close to the first and last data point    %VAR
     var_a = 1;    %VAR
     var_b = 4; % 2 means least squares straight line, 4 means linear...    %VAR
    
     norm_torque_ascend_spap2 = spap2(var_a, var_b, time_torque_ascend, norm_torque_ascend_lowpass, w);

     %%% method 4: cubic smoothing spline
     p = 0.0005; % 0 = least squares straight line fit, 1 = natural cubic spline interpolant   %VAR
     norm_torque_ascend_csaps = csaps(time_torque_ascend, norm_torque_ascend_lowpass, p, time_torque_ascend); 
     
     %%% method 5: smoothing with lowess (UIO help)
     % lowess = local regression using weighted linear least squares, 1st degree poly model
     % loess = ... 2nd degree poly model
     % rlo(w)ess = robust, lower weight to outliers in the regresion
 %    var_span = 0.9; % (r)lo(w)ess - percent span     %VAR
 %    norm_torque_ascend_smooth = smooth(norm_torque_ascend_lowpass, var_span, 'lowess');    
     
     % chosen filter for futher calculations
     norm_torque_ascend_filtered = norm_torque_ascend_csaps;
    
    if plot_individual
         plottitle = horzcat('Torque filter check ASCEND for ', subject_id, ' ', trial_name, '_pt2');
         figure('Name',plottitle)
         plot(time_torque_ascend,norm_torque_ascend_lowpass,'color','r','linewidth',2)
         hold on
         fnplt(norm_torque_ascend_spap2,'linewidth',1,'color','g')
         plot(time_torque_ascend,norm_torque_ascend_csaps,'b','linewidth',2)
  %       plot(time_torque_ascend,norm_torque_ascend_smooth,'c')
         xlabel('Time (s)')
         ylabel('Torque (Nm)')
         title(plottitle,'Interpreter', 'none')
         legend('initial lowpass','spap2 least-squares spline','csaps cubic smoothing spline','location', 'northwest')
         %print(horzcat('data_plots/IND knee ',plottitle),'-dpng')
    end



    % descending phase:
    
     %%% method 3: spline for torque, from Vidar 2014-05-06
     % spap2 = Least-squares spline approximation 
     % NB, not an array but a struct.....
     w = ones(size(time_torque_descend)); % w = weights in the error measure
     w([1 end]) = 100; % 100 - forces fit to come very close to the first and last data point    %VAR
     var_a = 1;    %VAR
     var_b = 4; % 2 means least squares straight line, 4 means linear...    %VAR
     norm_torque_descend_spap2 = spap2(var_a, var_b, time_torque_descend, norm_torque_descend_lowpass, w);

     %%% method 4: cubic smoothing spline
     p = 0.0005; % 0 = least squares straight line fit, 1 = natural cubic spline interpolant   %VAR
     norm_torque_descend_csaps = csaps(time_torque_descend, norm_torque_descend_lowpass, p, time_torque_descend); 
     
     %%% method 5: smoothing with lowess (UIO help)
     % lowess = local regression using weighted linear least squares, 1st degree poly model
     % loess = ... 2nd degree poly model
     % rlo(w)ess = robust, lower weight to outliers in the regresion
%     var_span = 0.9; % (r)lo(w)ess - percent span     %VAR
 %    norm_torque_descend_smooth = smooth(norm_torque_descend_lowpass, var_span, 'lowess');    
     
     % chosen filter for futher calculations
     norm_torque_descend_filtered = norm_torque_descend_csaps;
    
     
     
     if plot_individual
         plottitle = horzcat('Torque filter check DESCEND for ', subject_id, ' ', trial_name, '_pt2');
         figure('Name',plottitle)
         plot(time_torque_descend,norm_torque_descend_lowpass,'color','r','linewidth',2)
         hold on
         fnplt(norm_torque_descend_spap2,'linewidth',1,'color','g')
         plot(time_torque_descend,norm_torque_descend_csaps,'b','linewidth',2)
   %      plot(time_torque_descend,norm_torque_descend_smooth,'c')
         xlabel('Time (s)')
         ylabel('Torque (Nm)')
         title(plottitle,'Interpreter', 'none')
         legend('initial lowpass','spap2 least-squares spline','csaps cubic smoothing spline','location', 'northeast')
         % print(horzcat('data_plots/IND knee ',plottitle),'-dpng')
     end


    %% INSERT CONVERTED TORQUE into NEW noraxon-variable containing only 1st trial
    noraxon_converted = noraxondata.data(1:changedirection3,:);

    noraxon_converted(:,column_norm_torque) = zeros;
    noraxon_converted(changedirection0:changedirection1,column_norm_torque) = norm_torque_ascend_filtered;
    noraxon_converted(changedirection2:changedirection3,column_norm_torque) = norm_torque_descend_filtered;
    
    
    %% RESAMPLE
    % create time array with old and new timestamps and length
    timearray_orig = make_timearray(noraxonfreq, size(noraxon_converted));
    timearray_new = make_timearray(finalfreq, fix(length(noraxon_converted(:,1))/noraxonfreq*finalfreq));

    % create timeseries type for original data
    noraxon_timeseries = timeseries(noraxon_converted, timearray_orig, 'Name','Noraxon');

    % resample
    noraxon_temp = resample(noraxon_timeseries,timearray_new,'linear');
    noraxon_resampled = noraxon_temp.data;
    
        
    %% CORRECT ANGLE in resampled array (first values "bent" due to filtering)
    % use values #5 to 20 of angle data, to replace values #1-4
    loc_start = 5; %VAR
    loc_stop = 20; %VAR
    angle_correction = noraxon_resampled(1:loc_stop,column_norm_angle);
    
    % calculate linear coeffisients for extrapolation
    
    frames = [ones(length(angle_correction(loc_start:loc_stop)),1) (loc_start:loc_stop)'];
    coeffs = frames\angle_correction(loc_start:loc_stop); % NB backslash --> slope

    % replace first Z values in arrays by using extrapolation
    % y = ax + b ... y = displ, x = angle
    for x = 1:(loc_start-1)
        % place back in original array
        noraxon_resampled(x,column_norm_angle) = ((x*coeffs(2)) + coeffs(1));
    end
    
    
    %% SAVE SECOND TRIAL
    noraxon_resampled2 = noraxon_resampled;
    
    
end








%% ----------------------- FUNCTIONS -----------------------------------------


function output = make_timearray(freq,frames)
    output(1:frames,1) = 0.0;
    for i = 1:frames
        output(i) = 1/freq*(i-1);
    end
end




% function outputarray = RMS(emg_rms_frames, EMGarray)
%     outputarray(1:length(EMGarray)) = zeros;
%     datalength = length(EMGarray);
%     for i = 1:datalength
%         if i <= (emg_rms_frames/2) % first X values
%             outputarray(i)=RMSmaths(1, i+(emg_rms_frames/2), EMGarray); % first output is value 1 + next 75 values
%         elseif i > (datalength-(emg_rms_frames/2)) % last X values
%             outputarray(i)=RMSmaths(i-(emg_rms_frames/2), datalength, EMGarray);
%         else % all middle values, where sufficient data is available before and after
%             outputarray(i)=RMSmaths(i-(emg_rms_frames/2), i+(emg_rms_frames/2), EMGarray);
%         end
%     end
% end


% 
% function out = RMSmaths(firstframe,lastframe,EMGarray)
% 
% sum=0;
% for i=firstframe:lastframe
%     sum = sum + EMGarray(i)^2;
% end
% out = sqrt(sum/(lastframe-firstframe+1));
% 
% end