%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_momentarm
% Marie Moltubakk 11.6.2013
% Read complete, prepared noraon array + prepared us array
% Produce AT moment arm constant in meter, by tendon excursion method
%%%%%%%%%%%%%%%%%%%%%%%%%%
    


function at_momentarm = calculate_momentarm(noraxon_prepped, usdata_prepped, leg_length)
%    global dm_leg_length line

    
    
    %%% METHOD 1 (ORIGINALLY PREFERRED) - from ankle rotation in Norm
    % Tendon excursion method
    % Olivier reference: Maganaris 2004: Imaging-based estimates of moment arm length in intact human muscle-tendons
    % Original reference: An KN, Takahashi K, Harrigan TP, Chao EY (1984) Determinations of muscle orientations and moment arms. J Biomech Eng 106:280–282
    % Fath: Direct comparison of in vivo Achilles tendon moment arms obtained from ultrasound and MR scans. Describing in detail the methods for calculating AT moment arm from US scans of the MTJ
    
    % global plot_achilles plot_norm plot_emg plot_check subject_id
    % global column_EMG_start column_EMG_end column_l_gm column_r_gm column_l_gl column_r_gl column_l_sol column_r_sol column_l_tibant column_r_tibant column_norm_angle column_norm_torque column_norm_velocity column_norm_direction column_achilles    
    
    % this description MIGHT be from Fath, but not sure!
    
    % Identify movement phases (see angle constant calculation)
    
    % Apply mathematics
    %    A: At angle zero: First derivative of tendon displacement with respect to ankle angle. Use angular interval of +-5 degrees, spline fitted to find values at exact angles.
    %    B: Second and third order polyn to tendon displ data. Analytical differentiation at ankle angle of 0 deg.
    
    % Separate into dorsiflexed-to-plantarflexed and opposite
    
    % moment arm equals the derivative of tendon travel with respect to joint angulation


    

    
    %%% METHOD 2 - anthopometry method: Spoor 1990: Estimation of instantaneous moment arms of lower-leg muscles.
    percent_leg_length = 10.978444;
    at_momentarm = percent_leg_length/100 * str2double(leg_length)/100; % leg length is in cm
 



    


    % Output moment arm report
    cprintf('blue', horzcat('AT moment arm = ', num2str(at_momentarm), ' m.\n'));

end
