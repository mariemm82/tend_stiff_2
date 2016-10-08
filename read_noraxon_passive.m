%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% read_noraxon_file
% Marie Moltubakk 18.5.2013
% Read noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
% Produce a new noraxon data array

% version 21.1.2015, adapted to use individually calculated norm conversion factors
%%%%%%%%%%%%%%%%%%%%%%%%%%

function noraxon_resampled = read_noraxon_passive(noraxonfile, finalfreq, side, trial_name)
    global us_zerodispframes noraxonfreq emg_bandpass emg_rms_ms mvc_window_ms convert_achilles convert_norm_ind_passive
    global angle_cutoff velocity_cutoff torque_cutoff_bandstop torque_cutoff_active angle_cutoff_active velocity_cutoff_active
    global plot_norm plot_emg plot_check subject_id
    global column_EMG_start column_EMG_end column_l_gm column_r_gm column_l_gl column_r_gl column_l_sol column_r_sol column_l_tibant column_r_tibant column_gonio column_norm_angle column_norm_torque column_norm_velocity column_norm_direction column_achilles
    global convert_norm_angle_a convert_norm_angle_b convert_norm_torque_a convert_norm_torque_b convert_norm_velocity_a convert_norm_velocity_b convert_norm_direction_b
    
    % import noraxon data
    noraxondata = importdata(noraxonfile, '\t', 1);

    % reset time array to start from 0,0 (trigger click = time zero)
    noraxondata.data(:,1) = noraxondata.data(:,1) - noraxondata.data(1,1);
    time = noraxondata.data(:,1);
    
    
    
   
    %%% DATA CONVERSION EMG

    % butterworth filter
    noEMGchannels = column_EMG_end - column_EMG_start + 1;
    [B, A] = butter(4, emg_bandpass); % bandpass because of 2-element vector
    EMGfiltered(length(noraxondata.data),noEMGchannels) = zeros; 
    for i = column_EMG_start:column_EMG_end
       EMGfiltered(:, i-column_EMG_start+1) = filtfilt(B, A, noraxondata.data(:,i));
    end

    % RMS
    % for the very first and last frames (time = emg_rms_ms/2), RMS is
    % calculated based on the material from the first frame until the frame
    % corresponding to emg_rms_ms/2. In F24 project = 0.05 sec.
    emg_rms_frames = (noraxonfreq*(emg_rms_ms/1000));
    if (mod(emg_rms_frames,2)~=0) == 1
        exception = MException('MATLAB:odearguments:InconsistentDataType', 'Change EMG RMS window to a digit dividable by 2 and restart tendstiff.');
        throw(exception);
    end
    EMG_RMS(length(EMGfiltered),noEMGchannels) = zeros;
    for i = 1:length(EMGfiltered(1,:))
       EMG_RMS(:, i) = RMS(emg_rms_frames, EMGfiltered(:,i));
    end

        
    
    
    

     
     
    
    %%% DATA CONVERSION NORM

    

    %%% PREPARE AND FILTER ANGLE AND DIRECTION, IDENTIFY ZERO ANGLE

    % convert angle
    if side == 'R'
        angleadjust1 = -180;
    else
        angleadjust1 = 180;
    end
    norm_angle = angleadjust1 + (noraxondata.data(:,column_norm_angle) * convert_norm_angle_a) + convert_norm_angle_b;
    if max(norm_angle) > 340
        norm_angle = norm_angle - 360;
    elseif max(norm_angle) < -340
        norm_angle = norm_angle + 360;
    end
        
    % filter angle
    [B, A] = butter(4, angle_cutoff, 'low');
    norm_angle_filtered = filtfilt(B, A, norm_angle);
    
    
    
    % prepare direction
    norm_direction = noraxondata.data(:,column_norm_direction);

    % % filter direction - different from the others, because conversion -> +1/0/-1 only
    % [B, A] = butter(4, angle_cutoff, 'low');
    % norm_direction = filtfilt(B, A, noraxondata.data(:,column_norm_direction));
    
    
    
    % define "end of trial"  - data after that point are cut off
    % alt A: use full length of descending phase (past zero degrees)
%    loc_trialend = length(norm_angle_filtered); 
    % alt B: detect point of zero angle
%    loc_2ndhalf = round(length(norm_angle_filtered)/2);
%    loc_trialend = loc_2ndhalf + find(norm_angle_filtered(loc_2ndhalf:end)>0,1,'first');
    % alt C: detect point of starting angle if before zero: symmetric curves ascend and descend
    if norm_angle_filtered(1) > 0 && norm_angle_filtered(1) < norm_angle_filtered(end) % starts and ends with plantar flexion & end more plantar flexed
        loc_2ndhalf = round(length(norm_angle_filtered)/2);
        loc_trialend = loc_2ndhalf + find(norm_angle_filtered(loc_2ndhalf:end) > norm_angle_filtered(1),1,'first'); % use symmetric data
    elseif norm_angle_filtered(1) <= 0 && norm_angle_filtered(end) > 0 % starts in dorsi - use until zero if available
        loc_2ndhalf = round(length(norm_angle_filtered)/2);
        loc_trialend = loc_2ndhalf + find(norm_angle_filtered(loc_2ndhalf:end) > 0,1,'first');
    else
        loc_trialend = length(norm_angle_filtered); % use all
    end
   
    
    
    %%% PREPARE AND FILTER VELOCITY

    % convert velocity, positive velocity = towards dorsiflexion for both sides
    norm_velocity = -((noraxondata.data(:,column_norm_velocity) * convert_norm_velocity_a) + convert_norm_velocity_b);
    
    % filter velocity
    [B, A] = butter(4, velocity_cutoff, 'low');
    norm_velocity_filtered = filtfilt(B, A, norm_velocity);
    
    % make a rough correction for velocity offset:
    
    % find point where direction changes (=max angle)
    avoidframes = 100; %VAR, # of frames skipped around the point of change in direction, to avoid noise due to the change
    if strcmpi(side,'R') == 1
        % R: the direction data changes only at the END of the stop phase
        loc_maxROM = find(norm_direction<4000,1,'first');
        % if detected change in direction is just a small spike (back to high µV value within 100 frames), try again:
        loc_maxROMcheck = loc_maxROM+100;
        while norm_direction(loc_maxROMcheck)>4000
            % find next spike
            loc_maxROMnext = loc_maxROMcheck + find(norm_direction(loc_maxROMcheck:end)<4000,1,'first');
            loc_maxROM = loc_maxROMnext;
            loc_maxROMcheck = loc_maxROM+100;
        end
        % calculate average velocity during one second BEFORE change in direction, avoiding 100 frames around the change
        zero_velocity = mean(norm_velocity_filtered(loc_maxROM-avoidframes-noraxonfreq/2:loc_maxROM-avoidframes));
    else % left side
        %L: the direction data changes at the BEGINNING of the stop phase
        loc_maxROM = find(norm_direction>100,1,'first');
        % calculate average velocity during one second AFTER change in direction, avoiding 100 frames around the change
        zero_velocity = mean(norm_velocity_filtered(loc_maxROM+avoidframes:loc_maxROM+avoidframes+noraxonfreq/2));
    end
    
    
    % mean velocity at maxROM
    norm_velocity_filtered = norm_velocity_filtered-zero_velocity;
    zero_velocity_mv = (zero_velocity - convert_norm_velocity_b) * convert_norm_velocity_a;
    
    % print velocity correction report
    cprintf('magenta', horzcat('Applying conversions: Velocity: ', num2str(zero_velocity), ' deg/s calculated @ stop, corresponding offset ', num2str(zero_velocity_mv), ' µV.\n'));

    
    
    
    
    %%% REPLACE EMG, ANGLE, DIRECTION, VELOCITY IN DATA
    % replace all EMG columns
    for i = column_EMG_start:column_EMG_end
        noraxondata.data(:,i) = EMG_RMS(:,i-column_EMG_start+1);
    end
        
    % replace angle and direction
    noraxondata.data(:,column_norm_angle) = norm_angle_filtered;
    noraxondata.data(:,column_norm_direction) = norm_direction;
    noraxondata.data(:,column_norm_velocity) = norm_velocity_filtered;

    
    
    
    
    %%% DETECT ASCENDING PHASE - DESCENDING PHASE
    
    % detect phases
    vel1 = 2; % velocity in the first/ascending phase   %VAR
    vel2 = -2; % velocity in the second/descending phase   %VAR
    % find point where ascending phase stops
    changedirection = 100 + find(norm_velocity_filtered(100:end)<vel1+(vel2/2),1,'first');
    % find point where decending phase starts
    changedirection2 = changedirection + find(norm_velocity_filtered(changedirection:end)<vel2+(vel1/2),1,'first');
    % find point where decending phase stops - not used - using only until zero angle
    % changedirection3 = changedirection2+500 + find(norm_velocity_filtered(changedirection2+500:end)<vel2+(vel1/2),1,'first');

    % plot to confirm phases of direction changes
    if plot_check && plot_norm
        plottitle = horzcat('Passive phases check for ', subject_id, ' ', trial_name);
        figure('Name',plottitle)
        plot(time,-norm_direction/100)
        hold on
        plot(time,(norm_angle_filtered*5),'r') %scaled
        plot(time,200+(noraxondata.data(:,column_norm_torque))/2,'g') %scaled
        plot(time,(norm_velocity_filtered*20),'c') %scaled
        plot ([time(changedirection) time(changedirection)], [-100 100],'y')
        plot ([time(changedirection2) time(changedirection2)], [-100 100],'y')
        plot ([time(loc_trialend) time(loc_trialend)], [-100 100],'y')
        xlabel('Time (s)'),ylabel('Misc scaled variables'),title(plottitle);
        legend('Direction/100','Angle*5','Torque, µV/2, shift','Velocity*20','direction changes','location','southeast')
    end

    
    
    
    %%% PREPARE AND FILTER TORQUE 
    
    % for ascending and descending phase separately:
    
    % convert torque to give positive values for plantar flexion for both sides
    if strcmpi(side,'R') == 1
        norm_torque_ascend = (noraxondata.data(1:changedirection,column_norm_torque) * convert_norm_torque_a) + convert_norm_torque_b; 
        norm_torque_descend = (noraxondata.data(changedirection2:loc_trialend,column_norm_torque) * convert_norm_torque_a) + convert_norm_torque_b; 
    else % left leg to plantar flexion = positive
        norm_torque_ascend = -((noraxondata.data(1:changedirection,column_norm_torque) * convert_norm_torque_a) + convert_norm_torque_b); 
        norm_torque_descend = -((noraxondata.data(changedirection2:loc_trialend,column_norm_torque) * convert_norm_torque_a) + convert_norm_torque_b); 
    end

    % filter torque lowpass
    [B, A] = butter(4, torque_cutoff_active, 'low');
    norm_torque_ascend_lowpass = filtfilt(B, A, norm_torque_ascend);
    norm_torque_descend_lowpass = filtfilt(B, A, norm_torque_descend);
    
    
    
    % delete data from descending phase (from zero angle and into plantarflexion)

    % deleted data are stored in the following variable:
    % norm_torque_lowpass_timecut = norm_torque_ascend_lowpass(zeroangle:end);
    
    % prepare arrays containing times for ascending and descending phase 
    time_torque_ascend = noraxondata.data(1:changedirection,1);
    time_torque_descend = noraxondata.data(changedirection2:loc_trialend,1);
    
    
    
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
    
    if plot_check && plot_norm
         plottitle = horzcat('Torque filter check ASCEND for ', subject_id, ' ', trial_name);
         figure('Name',plottitle)
         plot(time_torque_ascend,norm_torque_ascend_lowpass,'color','r','linewidth',2)
         hold on
         fnplt(norm_torque_ascend_spap2,'linewidth',1,'color','g')
         plot(time_torque_ascend,norm_torque_ascend_csaps,'b','linewidth',2)
  %       plot(time_torque_ascend,norm_torque_ascend_smooth,'c')
         xlabel('Time (s)'),ylabel('Torque (Nm)'),title(plottitle);
         legend('initial lowpass','spap2 least-squares spline','csaps cubic smoothing spline','location', 'northwest')
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
    
     
     
     if plot_check && plot_norm
         plottitle = horzcat('Torque filter check DESCEND for ', subject_id, ' ', trial_name);
         figure('Name',plottitle)
         plot(time_torque_descend,norm_torque_descend_lowpass,'color','r','linewidth',2)
         hold on
         fnplt(norm_torque_descend_spap2,'linewidth',1,'color','g')
         plot(time_torque_descend,norm_torque_descend_csaps,'b','linewidth',2)
   %      plot(time_torque_descend,norm_torque_descend_smooth,'c')
         xlabel('Time (s)'),ylabel('Torque (Nm)'),title(plottitle);
         legend('initial lowpass','spap2 least-squares spline','csaps cubic smoothing spline','location', 'northwest')
     end





    
    %%% method 1: april 2014, padding and notch filtering    
%    % pad torque signal
%    cc = length(norm_torque_filt1);
%    for i = 1:1000
%        norm_torque_filt1(cc+i) = norm_torque_filt1(cc+i-1) + norm_torque_filt1(cc+i-1) - norm_torque_filt1(cc+i-2);
%    end
%   
%     % filter torque notch ~0,5 hz
%     notchWidth = 0.00015;       %#width of the notch %VAR 0.00015
%     f0 = 1/2.55;                %#notch frequency %VAR 
%     fn = noraxonfreq/2;              %#Nyquist frequency
%     freqRatio = f0/fn;      %#ratio of notch freq. to Nyquist freq.
% 
%     %Compute zeros
%     zerros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];
%     %Compute poles
%     poles = (1-notchWidth) * zerros;
% 
%     b = poly( zerros ); %# Get moving average filter coefficients 
%     a = poly( poles ); %# Get autoregressive filter coefficients
% 
%     %figure;
%     %freqz(b,a,32000,noraxonfreq)
% 
%     % filter signal
%     norm_torque_filtered = filter(b,a,norm_torque_filt1);
% 
%    % cutoff pads
%    norm_torque_filt1 = norm_torque_filt1(1:cc);
%    norm_torque_filtered = norm_torque_filtered(1:cc);
    

    %%% method 2: iirnotch filter
%    wo = f0/(noraxonfreq/2);
%    bw = wo/5;
%    [c,d] = iirnotch(wo,bw);
%    norm_torque_filtered2 = filter(c,d,norm_torque_filt1);

%     % fourier analysis
%     if plot_check
%         y = norm_torque_filt1;
%         T = 1/noraxonfreq;                     % Sample time
%         L = length(y);                     % Length of signal
%         t = (0:L-1)*T;                % Time vector
% 
%         NFFT = 2^nextpow2(L); % Next power of 2 from length of y
%         Y = fft(y,NFFT)/L;
%         f = noraxonfreq/2*linspace(0,1,NFFT/2+1);
%         % Plot single-sided amplitude spectrum.
%         figure,plot(f,2*abs(Y(1:NFFT/2+1))) 
%         title('Single-Sided Amplitude Spectrum of y(t)')
%         xlabel('Frequency (Hz)')
%         ylabel('|Y(f)|')
%     end







    %%% INSERT CONVERTED TORQUE
    noraxon_converted = noraxondata.data(1:loc_trialend,:);

    noraxon_converted(:,column_norm_torque) = zeros;
    noraxon_converted(1:changedirection,column_norm_torque) = norm_torque_ascend_filtered;  % change to chosen filter
    noraxon_converted(changedirection2:loc_trialend,column_norm_torque) = norm_torque_descend_filtered; % change to chosen filter
    

    
    
    
    
    %%% DATA CONVERSION GONIOMETER
    
    % turn left side gonio angles into negative for dorsiflexion
    if side == 'R'
        gonio_raw = noraxondata.data(:,column_gonio); 
    else
        gonio_raw = -noraxondata.data(:,column_gonio);
    end
    
    % smooth gonio
    [B, A] = butter(4, angle_cutoff/2, 'low');
    gonio_filtered = filtfilt(B, A, gonio_raw);
    
    % print warning if recording is started after dorsiflexion is initiated
    % (Norm dorsiflexion = negative)
    if norm_angle_filtered(1) <= -3 %VAR - if more dorsiflexed than -3 
        cprintf('err', horzcat('WARNING: ', trial_name, ' Recording starts after 3!!!! degrees dorsiflexion. DISCARD? Norm starting angle = ', num2str(norm_angle_filtered(1)),'.\n'));
    elseif norm_angle_filtered(1) <= 0  % trial was started too late, does not contain zero
        cprintf('err', horzcat('WARNING: ', trial_name, ' Recording starts after dorsiflexion initiation. Norm starting angle = ', num2str(norm_angle_filtered(1)),'.\n'));
    end % trial was started in time

    %%% gonio and norm angles are aligned at x degrees joint angle
    %%% (= secondary offset)
    angle_at_correction = 0; %VAR - use "+3" for 3 degrees plantar flexion
    
    % find location of correction angle in Norm
    if norm_angle_filtered(1) < angle_at_correction % trial starts later than specified angle - use first frame
        angle_normstart = norm_angle_filtered(1);
        loc_normstart = 1;
    else
        % find location of secondary offset "angle_at_correction" angle in norm data
        angle_normstart = angle_at_correction;
        loc_normstart = find(norm_angle_filtered < angle_at_correction,1,'first') - 1; % select location of last value before going into x degrees and dorsiflexion
    end
        
    % offset gonio at correction angle
    gonio_correction = angle_normstart - gonio_filtered(loc_normstart);
    gonio_corrected = gonio_filtered(:) + gonio_correction;
    
	% replace gonio in array
    noraxon_converted(1:loc_trialend,column_gonio) = gonio_corrected(1:loc_trialend);
    if plot_check && plot_norm
        plottitle = horzcat('Gonio correction check, ', subject_id, ' ', trial_name);
        figure('Name',plottitle)
        plot(time,noraxondata.data(:,column_norm_angle))
        hold on
        plot(time,gonio_filtered)
        plot(time,gonio_corrected)
        xlabel('Time (s)'),ylabel('Angle (deg)'),title(plottitle);
        legend('Norm', 'Gonio orig', 'Gonio corr','location', 'northeast')
    end


    
    
    
    
    %%% RESAMPLE

    % create time array with old and new timestamps and length
    timearray_orig = make_timearray(noraxonfreq, size(noraxon_converted));
    timearray_new = make_timearray(finalfreq, fix(length(noraxon_converted(:,1))/noraxonfreq*finalfreq));

    % create timeseries type for original data
    noraxon_timeseries = timeseries(noraxon_converted, timearray_orig, 'Name','Noraxon');

    % resample
    noraxon_temp = resample(noraxon_timeseries,timearray_new,'linear');
    noraxon_resampled = noraxon_temp.data;
    
    
    
    
    % correct ANGLE in resampled array (first values "bent" due to filtering)
    % use range 5 to 20 of gonio data, to replace values 1-4
    loc_start = 5; %VAR
    loc_stop = 20; %VAR
    angle_correction = noraxon_resampled(1:loc_stop,column_norm_angle);
    
    % calculate linear coeffisients for extrapolation
    
    frames = [ones(length(angle_correction(loc_start:loc_stop)),1) (loc_start:loc_stop)'];
    coeffs = frames\angle_correction(loc_start:loc_stop); % NB backslash --> slope

    % replace first 3 values in arrays by using extrapolation
    % y = ax + b ... y = displ, x = angle
    for x = 1:(loc_start-1)
        % place back in original array
        noraxon_resampled(x,column_norm_angle) = ((x*coeffs(2)) + coeffs(1));
    end
    
    
    
    
    
    
    
    
    %%% CHECKPOINT PLOTS
    
    % EMG triceps surae - for norm strength etc
    if plot_check && plot_emg
        plottitle = horzcat('EMG check, triceps surae, ', subject_id, ' ', trial_name);
        if strcmpi(side,'R') == 1
            plot_GM = noraxon_resampled(:,column_r_gm); 
            plot_GL = noraxon_resampled(:,column_r_gl);
            plot_SOL = noraxon_resampled(:,column_r_sol);
        else %left side
            plot_GM = noraxon_resampled(:,column_l_gm);
            plot_GL = noraxon_resampled(:,column_l_gl);
            plot_SOL = noraxon_resampled(:,column_l_sol);
        end
        figure('Name',plottitle)
        plot(noraxon_resampled(:,1),plot_GM,'y')
        hold on
        plot(noraxon_resampled(:,1),plot_GL,'m')
        plot(noraxon_resampled(:,1),plot_SOL,'c')
        xlabel('Time (s)'),ylabel('EMG (uV)'),title(plottitle);
        legend('GM','GL','SOL');
    end
       
    % plot final data
    if plot_check && plot_norm
        plottitle = horzcat('Norm torque-angle check for ', subject_id, ' ', trial_name);
        figure('Name',plottitle)
        plot(noraxon_resampled(:,column_norm_angle),noraxon_resampled(:,column_norm_torque),'k')
        hold on
        plot(noraxon_resampled(:,column_norm_angle),noraxon_resampled(:,column_norm_velocity),'r')
        set(gca,'xdir','reverse')
        xlabel('Angle (deg)'),ylabel('Torque/velocity'),title(plottitle);
        legend('torque','velocity','Location','Northwest');
    end
end








%%% ----------------------- FUNCTIONS -----------------------------------------




function output = make_timearray(freq,frames)
    output(1:frames,1) = 0.0;
    for i = 1:frames
        output(i) = 1/freq*(i-1);
    end
end




function outputarray = RMS(emg_rms_frames, EMGarray)
    outputarray(1:length(EMGarray)) = zeros;
    datalength = length(EMGarray);
    for i = 1:datalength
        if i <= (emg_rms_frames/2) % first X values
            outputarray(i)=RMSmaths(1, i+(emg_rms_frames/2), EMGarray); % first output is value 1 + next 75 values
        elseif i > (datalength-(emg_rms_frames/2)) % last X values
            outputarray(i)=RMSmaths(i-(emg_rms_frames/2), datalength, EMGarray);
        else % all middle values, where sufficient data is available before and after
            outputarray(i)=RMSmaths(i-(emg_rms_frames/2), i+(emg_rms_frames/2), EMGarray);
        end
    end
end



function out = RMSmaths(firstframe,lastframe,EMGarray)

sum=0;
for i=firstframe:lastframe
    sum = sum + EMGarray(i)^2;
end
out = sqrt(sum/(lastframe-firstframe+1));

end