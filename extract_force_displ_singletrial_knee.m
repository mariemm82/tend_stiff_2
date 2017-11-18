%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract_force_displ_singletrial_knee
% Marie Moltubakk 16.11.2017
% Read complete, prepared noraon array (no US for knee)
% Produce array with time, corrected force, corrected displacement
% 
% used by passive analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%



function varargout = extract_force_displ_singletrial_knee(noraxondata, side, n_o_trials, trialname)
    
    %% prepare
%    nOutputs = nargout;
%    varargout = cell(1,nOutputs);
    global column_norm_angle column_norm_torque
    global filepath filepath2 freq_default
        
    
    %% Read Noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
    % Produce a new noraxon data array
    if n_o_trials == 1
        noraxon_prepped = read_noraxon_knee(strcat(filepath, filepath2, noraxondata), freq_default, side, trialname);
    else % 2 trials
        [noraxon_prepped1, noraxon_prepped2] = read_noraxon_knee_2trials(strcat(filepath, filepath2, noraxondata), freq_default, side, trialname);
    end
   
       
    %% extract data from big arrays
    
    if n_o_trials == 1
        % locate endpoint of ascending phase
        % Norm torque = zero during stop at maximal ROM - use this to identify ascending phase (0 to max ROM).
        loc_ascending_ongoing = find(noraxon_prepped(:,column_norm_torque) > 0, 1, 'first');
        loc_ascending_end = loc_ascending_ongoing + find(noraxon_prepped(loc_ascending_ongoing:end,column_norm_torque)==0,1,'first') - 2;

        angle = noraxon_prepped(loc_ascending_ongoing:loc_ascending_end,column_norm_angle);
        torque = noraxon_prepped(loc_ascending_ongoing:loc_ascending_end,column_norm_torque);
        varargout = {angle torque};
    else
        % locate endpoint of ascending phase
        % Norm torque = zero during stop at maximal ROM - use this to identify ascending phase (0 to max ROM).
        loc_ascending_ongoing = find(noraxon_prepped1(:,column_norm_torque) > 0, 1, 'first');
        loc_ascending_end = loc_ascending_ongoing + find(noraxon_prepped1(loc_ascending_ongoing:end,column_norm_torque)==0,1,'first') - 2;

        angle1 = noraxon_prepped1(loc_ascending_ongoing:loc_ascending_end,column_norm_angle);
        torque1 = noraxon_prepped1(loc_ascending_ongoing:loc_ascending_end,column_norm_torque);

        loc_ascending_ongoing = find(noraxon_prepped2(:,column_norm_torque) > 0, 1, 'first');
        loc_ascending_end = loc_ascending_ongoing + find(noraxon_prepped2(loc_ascending_ongoing:end,column_norm_torque)==0,1,'first') - 2;

        angle2 = noraxon_prepped2(loc_ascending_ongoing:loc_ascending_end,column_norm_angle);
        torque2 = noraxon_prepped2(loc_ascending_ongoing:loc_ascending_end,column_norm_torque);
        
        varargout = {angle1 torque1 angle2 torque2};
    end
end