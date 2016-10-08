%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_noraxon_file
% Marie Moltubakk 18.5.2013
% Read noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
% Produce a new noraxon data array
%%%%%%%%%%%%%%%%%%%%%%%%%%

function noraxon_prepped = read_noraxon_file(noraxonfile, finalfreq, side, trial_name)
    global us_zerodispframes noraxonfreq emg_bandpass emg_rms_ms mvc_window_ms convert_achilles convert_norm_ind_active
    global angle_cutoff velocity_cutoff torque_cutoff_bandstop torque_cutoff_active angle_cutoff_active velocity_cutoff_active
    global plot_achilles plot_norm plot_emg plot_check subject_id
    global column_EMG_start column_EMG_end column_l_gm column_r_gm column_l_gl column_r_gl column_l_sol column_r_sol column_l_tibant column_r_tibant column_gonio column_norm_angle column_norm_torque column_norm_velocity column_norm_direction column_achilles
    
    % import noraxon data
    noraxondata = importdata(noraxonfile, '\t', 1);

    % set trigger frame as time zero
    noraxondata.data(:,1) = noraxondata.data(:,1) - noraxondata.data(1,1);
    

   
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

    


    %%% DATA CONVERSION ACHILLES
    
    % volt to torque for Achilles machine data column
    achilles_torque = noraxondata.data(:,column_achilles) * convert_achilles; 

    % filter achilles torque
    [B, A] = butter(4, torque_cutoff_active, 'low');
    achilles_torque_filtered = filtfilt(B, A, achilles_torque);
    
    
    
    
    %%% DATA CONVERSION GONIOMETER
     if strcmpi(side,'R') == 1
         % right leg towards dorsiflexion = negative
     else % left
         noraxondata.data(:,column_gonio) = -noraxondata.data(:,column_gonio); 
     end
    

    

    
    %%% DATA CONVERSION NORM
    
    % retrieve conversion constants for Norm data
    convert_norm_angle_a = convert_norm_ind_active(1); % PROJECTSPECIFIC
    convert_norm_angle_b = convert_norm_ind_active(2); % PROJECTSPECIFIC
    convert_norm_torque_a = convert_norm_ind_active(3); % PROJECTSPECIFIC
    convert_norm_torque_b = convert_norm_ind_active(4); % PROJECTSPECIFIC
    convert_norm_velocity_a = convert_norm_ind_active(5); % PROJECTSPECIFIC
    convert_norm_velocity_b = convert_norm_ind_active(6); % PROJECTSPECIFIC
    % TODO MMM to be replaced by call to calculate_norm_constants, once conversion testing phase is over. See how it is done in read_noraxon_passive.
%    [convert_norm_angle_a, convert_norm_angle_b, convert_norm_torque_a, convert_norm_torque_b, convert_norm_velocity_a, convert_norm_velocity_b, convert_norm_direction_b] = calculate_norm_constants('active', FP);

    % convert angle
    norm_angle = (noraxondata.data(:,column_norm_angle) * convert_norm_angle_a) + convert_norm_angle_b;
    if min(norm_angle) > 90
        norm_angle = norm_angle - 180;
    else % over 180
        norm_angle = norm_angle + 180;
    end
    
    % filter angle
    [B, A] = butter(4, angle_cutoff_active, 'low');
    norm_angle_filtered = filtfilt(B, A, norm_angle);
    
    % convert torque + positive towards plantar flexion for both sides
    if strcmpi(side,'R') == 1
        norm_torque = (noraxondata.data(:,column_norm_torque) * convert_norm_torque_a) + convert_norm_torque_b; 
    else % left leg to plantar flexion = positive
        norm_torque = -((noraxondata.data(:,column_norm_torque) * convert_norm_torque_a) + convert_norm_torque_b); 
    end
    
    % filter torque
    [B, A] = butter(4, torque_cutoff_active, 'low');
    norm_torque_filtered = filtfilt(B, A, norm_torque);
%    [B, A] = butter(2, torque_cutoff_bandstop, 'stop');
%    norm_torque_filtered = filtfilt(B, A, norm_torque_filtered_tmp);
%     
%     %if plot_check
%         plottitle = horzcat('Torque filter check ', subject_id);
%         figure('Name',plottitle)
%         plot(norm_torque)
%         hold on
%         plot(norm_torque_filtered_tmp)
%         plot(norm_torque_filtered)
%         xlabel('Frame'),ylabel('Torque (Nm)'),title(plottitle);
%         legend('unfiltered', 'lowpass', 'low+bandstop','location','southeast')
%     %end

    
     
	% convert velocity + positive towards plantar flexion for both sides
    if strcmpi(side,'R') == 1
        norm_velocity = (noraxondata.data(:,column_norm_velocity) * convert_norm_velocity_a) + convert_norm_velocity_b;
    else % left leg to plantar flexion = positive
        norm_velocity = -((noraxondata.data(:,column_norm_velocity) * convert_norm_velocity_a) + convert_norm_velocity_b);
    end
    
    % filter velocity
    [B, A] = butter(4, velocity_cutoff_active, 'low');
    norm_velocity_filtered = filtfilt(B, A, norm_velocity);
    
    % prepare direction
    norm_direction_filtered(length(noraxondata.data(:,column_norm_direction))) = zeros;
    
    % filter direction - different from the others, because conversion -> +1/0/-1 only
    [B, A] = butter(4, angle_cutoff_active, 'low'); % TMP confirm correct filter
    norm_direction_filtered = filtfilt(B, A, noraxondata.data(:,column_norm_direction));



    
    
    %%% CREATE NEW NORAXON ARRAY
    noraxon_converted = noraxondata.data;
    
    % replace all EMG columns
    for i = column_EMG_start:column_EMG_end
        noraxon_converted(:,i) = EMG_RMS(:,i-column_EMG_start+1);
    end
        
    % replace Achilles column
    noraxon_converted(:,column_achilles) = achilles_torque_filtered;

    % replace Norm columns
    noraxon_converted(:,column_norm_angle) = norm_angle_filtered;
    noraxon_converted(:,column_norm_torque) = norm_torque_filtered;
    noraxon_converted(:,column_norm_velocity) = norm_velocity_filtered;
    noraxon_converted(:,column_norm_direction) = norm_direction_filtered;
    
    
    
    
    %%% RESAMPLE

    % create time array with old and new timestamps and length
    timearray_orig = make_timearray(noraxonfreq, size(noraxon_converted));
    timearray_new = make_timearray(finalfreq, fix(length(noraxon_converted(:,1))/noraxonfreq*finalfreq));

    % create timeseries type for original data
    noraxon_timeseries = timeseries(noraxon_converted, timearray_orig, 'Name','Noraxon');

    % resample
    noraxon_temp = resample(noraxon_timeseries,timearray_new,'linear');
    noraxon_prepped = noraxon_temp.data;


    

    %%% CHECKPOINT PLOTS

    % EMG tibialis anterior - for tendon stiffness
    if plot_check && plot_achilles && plot_emg 
        plottitle = horzcat('EMG check, TA, ', subject_id, ' ', trial_name);
        if strcmpi(side,'R') == 1
            plot_raw = noraxondata.data(:,column_r_tibant);
            plot_RMS = EMG_RMS(:,column_r_tibant-1);
            plot_resamp = noraxon_prepped(:,column_r_tibant);
        else %left side
            plot_raw = noraxondata.data(:,column_l_tibant);
            plot_RMS = EMG_RMS(:,column_l_tibant-1);
            plot_resamp = noraxon_prepped(:,column_l_tibant);
        end
        figure('Name',plottitle)
        plot(noraxondata.data(:,1),plot_raw,'r')
        hold on
        plot(noraxondata.data(:,1),plot_RMS,'b')
        plot(noraxon_prepped(:,1),plot_resamp,'-ks','MarkerSize',2)
        xlabel('Time (s)'),ylabel('EMG (uV)'),title(plottitle);
        legend('raw EMG','filt+RMS EMG','filt+RMS+resamp EMG');
    end
    
    
    
    % EMG gastroc med - EMG quality check
    if plot_check && plot_norm && plot_emg 
        plottitle = horzcat('EMG check, GM, ', subject_id, ' ', trial_name);
        if strcmpi(side,'R') == 1
            plot_raw = noraxondata.data(:,column_r_gm);
            plot_RMS = EMG_RMS(:,column_r_gm-1);
            plot_resamp = noraxon_prepped(:,column_l_gm);
        else %left side
            plot_raw = noraxondata.data(:,column_l_gm);
            plot_RMS = EMG_RMS(:,column_l_gm-1);
            plot_resamp = noraxon_prepped(:,column_l_gm);
        end
        figure('Name',plottitle)
        plot(noraxondata.data(:,1),plot_raw,'r')
        hold on
        plot(noraxondata.data(:,1),plot_RMS,'b')
        plot(noraxon_prepped(:,1),plot_resamp,'-ks','MarkerSize',2)
        xlabel('Time (s)'),ylabel('EMG (uV)'),title(plottitle);
        legend('raw EMG','filt+RMS EMG','filt+RMS+resamp EMG');
    end
    
    % EMG triceps surae - for norm strength etc
    if plot_check && plot_norm && plot_emg 
        plottitle = horzcat('EMG check, triceps surae, ', subject_id, ' ', trial_name);
        if strcmpi(side,'R') == 1
            plot_GM = noraxon_prepped(:,column_r_gm); 
            plot_GL = noraxon_prepped(:,column_r_gl);
            plot_SOL = noraxon_prepped(:,column_r_sol);
        else %left side
            plot_GM = noraxon_prepped(:,column_l_gm);
            plot_GL = noraxon_prepped(:,column_l_gl);
            plot_SOL = noraxon_prepped(:,column_l_sol);
        end
        figure('Name',plottitle)
        plot(noraxon_prepped(:,1),plot_GM,'y')
        hold on
        plot(noraxon_prepped(:,1),plot_GL,'m')
        plot(noraxon_prepped(:,1),plot_SOL,'c')
        xlabel('Time (s)'),ylabel('EMG (uV)'),title(plottitle);
        legend('GM','GL','SOL');
    end

%     % Achilles torque
%      if plot_check && plot_achilles
%          plottitle = horzcat('Achilles torque check for ', subject_id, ' ', trial_name);
%          figure('Name',plottitle)
%          plot(noraxondata.data(:,1),achilles_torque,'r')
%          hold on
%          plot(noraxondata.data(:,1),achilles_torque_filtered,'b')
%          plot(noraxon_prepped(:,1),noraxon_prepped(:,column_achilles),'-ks','MarkerSize',2)
%          xlabel('Time (s)'),ylabel('Achilles torque (Nm)'),title(plottitle);
%          legend('raw torque','filt torque','filt+resamp torque','Location','Southeast');
%      end

    % Norm torque, raw and resampled
    if plot_check && plot_norm
        plottitle = horzcat('Norm torque check for ', subject_id, ' ', trial_name);
        figure('Name',plottitle)
        plot(noraxondata.data(:,1),norm_torque,'r')
        hold on
        plot(noraxondata.data(:,1),norm_torque_filtered,'b')
        plot(noraxon_prepped(:,1),noraxon_prepped(:,column_norm_torque),'-ks','MarkerSize',2)
        xlabel('Time (s)'),ylabel('Norm torque (Nm)'),title(plottitle);
        legend('raw torque','filt torque','filt+resamp torque','Location','Southeast');
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