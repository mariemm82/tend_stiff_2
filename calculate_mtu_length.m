%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_mtu_length
% Marie Moltubakk 12.2.2015
% Read length of lower leg + ankle angle array
% Produce array with leg length across angle array
%%%%%%%%%%%%%%%%%%%%%%%%%%



function [MTU_length_array, MTU_elong_array, MTU_strain_array, resting_vars] = calculate_mtu_length(angle_displ_SOL, angle_displ_GMMTJ, angle_displ_GMFAS, angle_GM_Fukunaga, in_resting_at_SOL_length, in_resting_at_GM_length, standing_calf_length, resting_GM_pennation, resting_GM_faslen, angle_common, resting_ankle_angle)
    global at_momentarm  % subject_id

        
    %% create array of angles common to all 6 trials
    % angle_common is read from file from create_angles_passive, to avoid potential tiny discrepancies between various calculations of max
    angle_array = (0:0.05:(angle_common+0.01))'; %VAR  - adding 0.01 to ensure that the actual common angle is included - data stored as 10.9999999 
    
       
    %% LENGTH of whole MTU (ELONGATION further down is based on prone/resting length)

    %%% Grieve (1978)
    % calculate percent change in MTU length, across angle array (degrees)
    % Grieve = gastroc MTU length, from tendon to tendon
    % angle array = full ROM test (angles in degrees)
    % my data: ankle in neutral position as 0 degrees, dorsiflexed as positive
    % Grieve: ankle in neutral position = 90 degrees, dorsiflexed as more than 90 degrees  MMM 
    
    % Grieve constants for ankle
    a0 = -22.18468;
    a1 = 0.30141;
    a2 = -0.00061;
    
    % calf/MTU elongation - in % change from standing length
    MTU_length_Grieve_rel = (a0 + (a1*(90+angle_array)) + (a2*(90+angle_array).^2))/100;

    % calf/MTU length - in mm, across angle array
    MTU_length_Grieve_abs = standing_calf_length * (1 + MTU_length_Grieve_rel);
        
    % resting MTU length (prone position)
    resting_calf_length_rel = (a0 + (a1*(90+resting_ankle_angle)) + (a2*(90+resting_ankle_angle).^2))/100;
    resting_calf_length_abs = standing_calf_length * (1 + resting_calf_length_rel);
    
    
    
    %%% Hawkins & Hull (1990)
    % H&H calculates "normalized muscle tendon length" = normalized to shank = lateral malleolus to lateral epicondyle
    % base equation: calf_length_array_HH = C0 + (C1*aH) + (C2*aK) + (C3*(aK^2)) + (C4*aA);
    aH = 90; % hip angle in degrees (value is not used since C1=0)
    aK = 0; % knee angle, straight knee = 0
    aA = (angle_array+90); % ankle angle: full plantar flex = 0, range = 60 to 120deg flexions
    C1 = 0;
    C3 = 0;
    % constants for GM:
    C0 = 0.9;
    C2 = -0.00062;
    C4 = 0.00214;
    
    % calculating elongation (% of H&H 0deg measured length), across angle array (from 0deg to max ROM)
    MTU_length_HH_GM_rel = C0 + (C1*aH) + (C2*aK) + (C3*(aK^2)) + (C4*aA);
    
%     % constants for GL:
%     C0 = 0.894;
%     C2 = -0.0005;
%     C4 = 0.00214;
%     MTU_length_HH_GL_rel = C0 + (C1*aH) + (C2*aK) + (C3*(aK^2)) + (C4*aA);

    % constants for SOL:
    C0 = 0.563;
    C2 = 0;
    C4 = 0.00193;
    
    MTU_length_HH_SOL_rel = C0 + (C1*aH) + (C2*aK) + (C3*(aK^2)) + (C4*aA);
    % SOL MTU length at rest per H&H method:
    aA = resting_ankle_angle+90;
    MTU_length_HH_SOL_rel_rest = C0 + (C1*aH) + (C2*aK) + (C3*(aK^2)) + (C4*aA);
    
    % correction for H&H: 
    % we MEASURED 0deg MTU length (calc insert to lat epi), 
    % while H&H calculates 0deg MTU length to ~109% of their measured length (lat mall to lat epi)
    
    % using our 0deg length (calc-epi) + H&H equation to calculate corresponding H&H 0deg length (mall-epi)
    standing_leg_length_HH = standing_calf_length/MTU_length_HH_GM_rel(1);
    
    % applying H&H percent wise elongation to above H&H 0deg length
    MTU_length_HH_GM_abs = standing_leg_length_HH * MTU_length_HH_GM_rel;
    
    MTU_length_HH_SOL_abs = standing_leg_length_HH * MTU_length_HH_SOL_rel;
    MTU_length_HH_SOL_abs_rest = standing_leg_length_HH * MTU_length_HH_SOL_rel_rest;
    
    
    
    %%% elongation from ankle rotation and moment arm (geometry)
    % convert momentarm from m to mm (* 1000)
    MTU_length_geometry = standing_calf_length + (at_momentarm * 1000 * sind(angle_array));

    
    
%     %%% plot method comparison
%     plottitle = horzcat('Comparison of MTU length models, ', subject_id);
%     figure('Name',plottitle)
%     hold on
%     plot(angle_array,MTU_length_Grieve_abs)
%     plot(angle_array,MTU_length_HH_GM_abs)
%     plot(angle_array,MTU_length_geometry)
%     legend('Grieve','Hawkin&Hull','Geometry', 'Location', 'NorthWest')
%     xlabel('Gonio angle (°)')
%     ylabel('Length (mm)')
%     title(plottitle)
    


    %%% choice of method for GM MTU length:
    MTU_GM_length = MTU_length_Grieve_abs; % alternatives: MTU_length_geometry; MTU_length_HH_GM_abs
    
        
    %% Length of free Achilles tendon

    % Initial length = measured seated in Achilles machine, ankle 0 degrees, US still images of "railway" across lower leg. Length in mm.
    resting_AT_SOL_length = str2double(in_resting_at_SOL_length);
    % reshape displacement:
    AT_SOL_displ = spline(angle_displ_SOL(:,1), angle_displ_SOL(:,2), angle_array);
    % convert to length + reshape:
    % SOL MTJ displacement = elongation of proximal portion.
    % AT elongation = elongation of MTJ - elongation proximal portion.
    % AT elongation = MTU elongation - SOL MTJ displacement.
    % AT length = MTU elongation - SOL MTJ displacement + AT length at zero.
    tend_AT_length_array = (MTU_GM_length-MTU_GM_length(1)) + resting_AT_SOL_length - AT_SOL_displ;
    
    
    %% Length of Achilles tendon up to GM insertion ("GM tendon")
    resting_AT_GM_length = str2double(in_resting_at_GM_length);
    % reshape displacement:
    AT_GM_displ = spline(angle_displ_GMMTJ(:,1), angle_displ_GMMTJ(:,2), angle_array);
    % convert to length + reshape
    tend_GM_length_array = (MTU_GM_length-MTU_GM_length(1)) + resting_AT_GM_length - AT_GM_displ;
    
        
    %% GMFAS displacement: spline to fit same angle array
    GMFAS_displ_array = spline(angle_displ_GMFAS(:,1), angle_displ_GMFAS(:,2), angle_array);
   
       
    %% GM msc and apo length array
    msc_GM_length_array = MTU_GM_length - tend_GM_length_array;
    apo_GM_array = tend_GM_length_array - tend_AT_length_array;
    
    
    %% SOL msc length array
    msc_SOL_length_array = MTU_length_HH_SOL_abs - tend_AT_length_array;
    
    
    %% GM muscle and SEE elongation from anthropometry (Lichtwark/Fukunaga)
    % angle_GM_Fukunaga contains:
        %   averaged angle (currently calculated from gonio)
        %   averaged fasicle length
        %   averaged pennation angle
        %   averaged fascicle elongation
        %   averaged fascicle strain
        
    resting_GM_msc_len = resting_GM_faslen * cosd(resting_GM_pennation);

    if angle_GM_Fukunaga == 0 % fascicle data not existing
        msc_GM_elong_Fukunaga(1:length(angle_array),1) = zeros;
    else
%         % version 1: lower leg = tend + msc
%         loc_frame = find(angle_GM_Fukunaga(:,1)>=0,1,'first');
%         loc_frame2 = find(angle_GM_Fukunaga(:,1)>=angle_common,1,'first');
%         % change to correct angle array
%         msc_GM_elong_Fukunaga = angle_GM_Fukunaga(loc_frame:loc_frame2,2);
%         % calculate other elongations, lengths, strains
%         msc_GM_length_Fukunaga = (calf_length - AT_GM_length) + msc_GM_elong_Fukunaga;
%         tend_GM_length_Fukunaga = MTU_GM_length - msc_GM_length_Fukunaga;
%         tend_GM_elong_Fukunaga = tend_GM_length_Fukunaga - tend_GM_length_Fukunaga(1);

        % version 2: msc = fascicle horiz displ, SEE = lower leg minus msc
        loc_frame = find(angle_GM_Fukunaga(:,1)>=0,1,'first');
        loc_frame2 = find(angle_GM_Fukunaga(:,1)>=angle_common,1,'first');
        
        % MUSCLE: length (longitudinal axis) from Fukunaga fasicles/pennation:
        msc_GM_length_Fukunaga = angle_GM_Fukunaga(loc_frame:loc_frame2,2) .* cosd(angle_GM_Fukunaga(loc_frame:loc_frame2,3));
        % MUSCLE: strain and elongation - according to GM msc length at prone/resting ankle angle
        msc_GM_elong_Fukunaga = msc_GM_length_Fukunaga - resting_GM_msc_len;
        msc_GM_strain_Fukunaga = msc_GM_elong_Fukunaga / resting_GM_msc_len * 100;
        
        % SEE: remainder of MTU length (e.g. Grieve) = SEE length:
        SEE_length_Fukunaga = MTU_GM_length - msc_GM_length_Fukunaga;
        % SEE: strain and elongation - according to SEE length at prone/resting ankle angle
        resting_SEE_len = resting_calf_length_abs - resting_GM_msc_len;
        SEE_elong_Fukunaga = SEE_length_Fukunaga - resting_SEE_len; % elong = change from resting SEE length (resting calf minus resting msc)
        SEE_strain_Fukunaga = SEE_elong_Fukunaga / resting_SEE_len * 100; % strain = in percent of resting SEE length
    end
    
    
    %% Final data
    
    MTU_length_array = [ ...
        angle_array ...
        tend_AT_length_array ...    % = free AT
        tend_GM_length_array ...    % = GM tend, from calc to GM insert
        MTU_GM_length ...           % = from calc to knee
        zeros(length(angle_array),1) ...  % GMFAS displacement has no length
        apo_GM_array ...            % = from end of free AT to GM ins
        msc_GM_length_array ...     % = from GM ins to knee
        msc_SOL_length_array ...    % = H&H SOL length minus free AT
        msc_GM_length_Fukunaga ...  % = GM muscle length based on GM faslen + penn.ang. (Lichtwark/Fukunaga)
        SEE_length_Fukunaga ...     % = SEE length based on GM faslen + penn.ang. (Lichtwark/Fukunaga)
        ]; 

    % calculate additional resting lengths:
    resting_apo_length = resting_AT_GM_length - resting_AT_SOL_length;
    resting_msc_GM_length = resting_calf_length_abs - resting_AT_GM_length;
    resting_msc_SOL_length = MTU_length_HH_SOL_abs_rest - resting_AT_SOL_length;
    
    MTU_elong_array = [...
        angle_array ...
        tend_AT_length_array-resting_AT_SOL_length ...% based on resting free AT length           OLD: tend_AT_length_array-tend_AT_length_array(1)
        tend_GM_length_array-resting_AT_GM_length ... % based on resting whole AT length          OLD: tend_GM_length_array-tend_GM_length_array(1)
        MTU_GM_length-resting_calf_length_abs ...        % based on resting leg length
        GMFAS_displ_array ...                            % DISPLACEMENT from GMFAS tracking
        apo_GM_array-resting_apo_length ...              % based on resting length                 OLD: apo_GM_array-apo_GM_array(1)
        msc_GM_length_array-resting_msc_GM_length...     % based on resting length                  OLD: msc_GM_length_array-msc_GM_length_array(1)
        msc_SOL_length_array-resting_msc_SOL_length...   % based on resting_msc_SOL_length               OLD: (like above)
        msc_GM_elong_Fukunaga ...                        % based on RESTING GM msc length
        SEE_elong_Fukunaga ...                           % based on RESTING SEE length
        ]; 
    
    MTU_strain_array = [ ...
        angle_array ...
        (tend_AT_length_array-resting_AT_SOL_length) / resting_AT_SOL_length*100 ...
        (tend_GM_length_array-resting_AT_GM_length) / resting_AT_GM_length*100 ...     
        (MTU_GM_length-resting_calf_length_abs) / resting_calf_length_abs*100 ...
        zeros(length(angle_array),1) ...                                                  % GMFAS has no length --> no strain
        (apo_GM_array-resting_apo_length) / resting_apo_length*100 ...
        (msc_GM_length_array-resting_msc_GM_length) /resting_msc_GM_length*100 ...
        (msc_SOL_length_array-resting_msc_SOL_length) / resting_msc_SOL_length*100 ...
        msc_GM_strain_Fukunaga ...
    	SEE_strain_Fukunaga ...
        ]; 
    
    resting_vars = [resting_ankle_angle resting_AT_SOL_length resting_AT_GM_length resting_calf_length_abs resting_msc_GM_length 0 resting_apo_length resting_msc_SOL_length resting_GM_msc_len resting_SEE_len resting_GM_pennation resting_GM_faslen];
    % matching the following col references:
%         col_angle_elong = 1;
%         col_AT = 2;
%         col_GMtend = 3;
%         col_leg = 4;
%         col_GMFAS = 5;
%         col_GMapo = 6;
%         col_GMmsc = 7;
%         col_SOLmsc = 8;
%         col_GMmsc_Fukunaga = 9;
%         col_SEE_Fukunaga = 10;
%         pennation angle = 11;
%         faslen = 12;
    
end