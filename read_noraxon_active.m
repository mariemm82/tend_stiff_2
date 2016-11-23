%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% read_noraxon_active
% Marie Moltubakk 2.11.2016
% Read noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
% Produce a new noraxon data array
%%%%%%%%%%%%%%%%%%%%%%%%%%

function noraxon_resampled = read_noraxon_active(noraxonfile, finalfreq, side, trial_name)
    global us_zerodispframes noraxonfreq emg_bandpass emg_rms_ms mvc_window_ms convert_achilles convert_norm_ind_passive
    global angle_cutoff velocity_cutoff torque_cutoff_bandstop torque_cutoff_active angle_cutoff_active velocity_cutoff_active
    global plot_check plot_individual subject_id
    global column_EMG_start column_EMG_end column_l_gm column_r_gm column_l_gl column_r_gl column_l_sol column_r_sol column_l_tibant column_r_tibant column_gonio column_norm_angle column_norm_torque column_norm_velocity column_norm_direction column_achilles
    global convert_norm_angle_a convert_norm_angle_b convert_norm_torque_a convert_norm_torque_b convert_norm_velocity_a convert_norm_velocity_b convert_norm_direction_b
    
    % import noraxon data
    noraxondata = importdata(noraxonfile, '\t', 1);

    % reset time array to start from 0,0 (trigger click = time zero)
    noraxondata.data(:,1) = noraxondata.data(:,1) - noraxondata.data(1,1);
    
    
    


%     %%% DATA CONVERSION EMG
% 
%     % butterworth filter
%     noEMGchannels = column_EMG_end - column_EMG_start + 1;
%     [B, A] = butter(4, emg_bandpass); % bandpass because of 2-element vector
%     EMGfiltered(length(noraxondata.data),noEMGchannels) = zeros; 
%     for i = column_EMG_start:column_EMG_end
%        EMGfiltered(:, i-column_EMG_start+1) = filtfilt(B, A, noraxondata.data(:,i));
%     end
% 
%     % RMS
%     % for the very first and last frames (time = emg_rms_ms/2), RMS is
%     % calculated based on the material from the first frame until the frame
%     % corresponding to emg_rms_ms/2. In F24 project = 0.05 sec.
%     emg_rms_frames = (noraxonfreq*(emg_rms_ms/1000));
%     if (mod(emg_rms_frames,2)~=0) == 1
%         exception = MException('MATLAB:odearguments:InconsistentDataType', 'Change EMG RMS window to a digit dividable by 2 and restart tendstiff.');
%         throw(exception);
%     end
%     EMG_RMS(length(EMGfiltered),noEMGchannels) = zeros;
%     for i = 1:length(EMGfiltered(1,:))
%        EMG_RMS(:, i) = RMS(emg_rms_frames, EMGfiltered(:,i));
%     end

        
    
    
    

     
     
    
    

    %%% PREPARE AND FILTER ANGLE AND DIRECTION, IDENTIFY ZERO ANGLE

    % convert angle
    if side == 'R'
        angleadjust1 = -180;
    else
        angleadjust1 = 180;
    end
    norm_angle = angleadjust1 + (noraxondata.data(:,column_norm_angle) * convert_norm_angle_a) + convert_norm_angle_b;
    if max(norm_angle) > 325
        norm_angle = norm_angle - 360;
    elseif max(norm_angle) < -325
        norm_angle = norm_angle + 360;
    end
        
    % filter angle
    [B, A] = butter(4, angle_cutoff_active, 'low');  % TODO MMM change all cutoff to active
    norm_angle_filtered = filtfilt(B, A, norm_angle);
    
    % prepare direction
    norm_direction = noraxondata.data(:,column_norm_direction);

    % % filter direction - different from the others, because conversion -> +1/0/-1 only
    % [B, A] = butter(4, angle_cutoff, 'low');
    % norm_direction = filtfilt(B, A, noraxondata.data(:,column_norm_direction));
    
    
    
    
    
    
    %%% PREPARE AND FILTER VELOCITY

    % convert velocity, positive velocity = towards dorsiflexion for both sides
    norm_velocity = -((noraxondata.data(:,column_norm_velocity) * convert_norm_velocity_a) + convert_norm_velocity_b);
    
    % filter velocity
    [B, A] = butter(4, velocity_cutoff_active, 'low');
    norm_velocity_filtered = filtfilt(B, A, norm_velocity);
    
    % make a rough correction for velocity offset: - MMM TODO ???
    

    
    
    
    
    
    %%% REPLACE EMG, ANGLE, DIRECTION, VELOCITY IN DATA
%     % replace all EMG columns
%     for i = column_EMG_start:column_EMG_end
%         noraxondata.data(:,i) = EMG_RMS(:,i-column_EMG_start+1);
%     end
        
    % replace angle and direction
    noraxondata.data(:,column_norm_angle) = norm_angle_filtered;
    noraxondata.data(:,column_norm_direction) = norm_direction;
    noraxondata.data(:,column_norm_velocity) = norm_velocity_filtered;

    
    
    
    
    
    %%% PREPARE AND FILTER TORQUE 
    
    % convert torque to give positive values for plantar flexion for both sides
    if strcmpi(side,'R') == 1
        norm_torque = (noraxondata.data(:,column_norm_torque) * convert_norm_torque_a) + convert_norm_torque_b; 
    else % left leg to plantar flexion = positive
        norm_torque = -((noraxondata.data(:,column_norm_torque) * convert_norm_torque_a) + convert_norm_torque_b); 
    end

    % filter torque lowpass
    [B, A] = butter(4, torque_cutoff_active, 'low');
    norm_torque_lowpass = filtfilt(B, A, norm_torque);



    %%% INSERT CONVERTED TORQUE
    noraxondata.data(:,column_norm_torque) = norm_torque_lowpass;
    
    
    
    
    
    
    
    
    %%% RESAMPLE

    noraxon_converted = noraxondata.data;

    % create time array with old and new timestamps and length
    timearray_orig = make_timearray(noraxonfreq, size(noraxon_converted));
    timearray_new = make_timearray(finalfreq, fix(length(noraxon_converted(:,1))/noraxonfreq*finalfreq));

    % create timeseries type for original data
    noraxon_timeseries = timeseries(noraxon_converted, timearray_orig, 'Name','Noraxon');

    % resample
    noraxon_temp = resample(noraxon_timeseries,timearray_new,'linear');
    noraxon_resampled = noraxon_temp.data;
    
    
    
    
    
    
    
    %%% CHECKPOINT PLOTS

       
    % plot final data
    if plot_individual
        plottitle = horzcat('Norm torque check for ', subject_id, ' ', trial_name);
        figure('Name',plottitle)
        plot(noraxon_resampled(:,1),noraxon_resampled(:,column_norm_torque),'k','LineWidth',2)
        hold on
        plot(noraxon_resampled(:,1),noraxon_resampled(:,column_norm_angle),'b')
        plot(noraxon_resampled(:,1),noraxon_resampled(:,column_norm_velocity),'r')
        xlabel('Time (s)'),ylabel('Torque/angle/velocity'),title(plottitle);
        legend('Torque','Angle','Velocity','Location','Northwest');
    end
end








%%% ----------------------- FUNCTIONS -----------------------------------------




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
% 
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