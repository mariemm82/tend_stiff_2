%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_mtu_length
% Marie Moltubakk 12.2.2015
% Read length of lower leg + ankle angle array
% Produce array with leg length across angle array
%%%%%%%%%%%%%%%%%%%%%%%%%%



function [MTU_length_array, MTU_elong_array, MTU_strain_array] = calculate_mtu_length(angle_displ_SOL, angle_displ_GMMTJ, angle_displ_GMFAS, initial_at_SOL_length, initial_at_GM_length, initial_calf_length, angle_common)



    %%%%%%%% create array of angles common to all 6 trials
    % angle_common is read from file from create_angles_passive, to avoid potential tiny discrepancies between various calculations of max
    angle_array = (0:0.05:(angle_common+0.01))'; %VAR  - adding 0.01 to ensure that the actual common angle is included - data stored as 10.9999999 



    %%%%%%%% Length of whole MTU

    %%% extract initial calf length in mm
    % my initial calf length = centimeters from lateral epicondyle to calc insertion
    calf_length = str2double(initial_calf_length) * 10;
    
    %%% Grieve 1978:
    % calculate percent change in MTU length, across angle array (degrees)
    % Grieve = gastroc MTU length, from tendon to tendon
    % angle array = full ROM test (angles in degrees)
    % my data report ankle in neutral position as 0 degrees, Grieve: ankle in neutral position = 90 degrees
    a0 = -22.18468;
    a1 = 0.30141;
    a2 = -0.00061;
    % calf/MTU length in % change from initial length:
    MTU_length_Grieve_rel = (a0 + (a1*(90+angle_array)) + (a2*(90+angle_array).^2))/100;

    % absolute calf/MTU length (mm), across angle array
    MTU_length_Grieve_abs = calf_length * (1 + MTU_length_Grieve_rel);

    %%% Hawkins & Hull (1990)
    % "normalized muscle tendon length" = normalized to shank = lateral malleolus to lateral epicondyle
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
    MTU_length_HH_GM_rel = C0 + (C1*aH) + (C2*aK) + (C3*(aK^2)) + (C4*aA);
%     % constants for GL:
%     C0 = 0.894;
%     C2 = -0.0005;
%     C4 = 0.00214;
%     calf_length_array_HH_GL = C0 + (C1*aH) + (C2*aK) + (C3*(aK^2)) + (C4*aA);
    % constants for SOL:
    C0 = 0.563;
    C2 = 0;
    C4 = 0.00193;
    MTU_length_HH_SOL_rel = C0 + (C1*aH) + (C2*aK) + (C3*(aK^2)) + (C4*aA);
    
    % correction: 
    % we KNOW the exact initial MTU length = calc insert to lat epi = 100% length (instead of H&H which says MTU length = ~108% of lat mall to lat epi)
    % first data point = 0 degrees ankle = 100% length
    % extracting length change from this point (% of calf/MTU length)
    MTU_length_HH_GM_abs = calf_length * (1+(MTU_length_HH_GM_rel - MTU_length_HH_GM_rel(1)));
    MTU_length_HH_SOL_abs = calf_length * MTU_length_HH_SOL_rel;
    
    %%% choice of method for GM MTU length:
    MTU_GM_length = MTU_length_HH_GM_abs; %alternative MTU_length_Grieve_abs
    
    
    
    
    %%%%%%% Length of free Achilles tendon

    % Initial length = measured seated in Achilles machine, ankle 0 degrees, US still images of "railway" across lower leg. Length in mm.
    AT_SOL_length = str2double(initial_at_SOL_length);
    % reshape displacement:
    AT_SOL_displ = spline(angle_displ_SOL(:,1), angle_displ_SOL(:,2), angle_array);
    % convert to length + reshape
    % OLD, error:    AT_SOL_length_array = spline(angle_displ_SOL(:,1), (AT_SOL_length + angle_displ_SOL(:,2)), angle_array);
    AT_SOL_length_array = (MTU_GM_length-MTU_GM_length(1)) + AT_SOL_length - AT_SOL_displ;
    
    
    
    %%%%%%% Length of Achilles tendon up to GM insertion ("GM tendon")
    AT_GM_length = str2double(initial_at_GM_length);
    % reshape displacement:
    AT_GM_displ = spline(angle_displ_GMMTJ(:,1), angle_displ_GMMTJ(:,2), angle_array);
    % convert to length + reshape
    % OLD, error:     AT_GM_length_array = spline(angle_displ_GMMTJ(:,1), (AT_GM_length + angle_displ_GMMTJ(:,2)), angle_array);
    AT_GM_length_array = (MTU_GM_length-MTU_GM_length(1)) + AT_GM_length - AT_GM_displ;
    
    
    
    %%%%% GMFAS displacement, spline to fit same angle array
    GMFAS_displ_array = spline(angle_displ_GMFAS(:,1), angle_displ_GMFAS(:,2), angle_array);
   

    
    
    %%%% GM msc length array
    msc_GM_length_array = MTU_GM_length - AT_GM_length_array;
    apo_GM_array = AT_GM_length_array - AT_SOL_length_array;
    
    

    %%%% SOL msc length array
    msc_SOL_length_array = MTU_length_HH_SOL_abs - AT_SOL_length_array;
    
    
    
    
    %%%%%% Final data
    MTU_length_array = [ ...
        angle_array ...
        AT_SOL_length_array ...% = free AT
        AT_GM_length_array ... % = GM length from calc to GM insert
        MTU_GM_length ...      % = from calc to knee
        zeros(length(angle_array),1) ...                          % GMFAS track has no length
        apo_GM_array ...       % = from SOL to GM ins
        msc_GM_length_array ... % = from GM ins to knee
        msc_SOL_length_array];  % = H&H SOL length minus free AT
    
    MTU_elong_array = [...
        angle_array ...
        AT_SOL_length_array-AT_SOL_length_array(1) ...
        AT_GM_length_array-AT_GM_length_array(1) ...    % ((at_GM_length_array-at_SOL_length_array) - (at_GM_length_array(1)-at_SOL_length_array(1))    ) ... % = GM apo, from end of free AT to GM insert
        MTU_GM_length-MTU_GM_length(1) ...
        GMFAS_displ_array ...                           % GMFAS track has displacement 
        apo_GM_array-apo_GM_array(1) ...
        msc_GM_length_array-msc_GM_length_array(1)...
        msc_SOL_length_array-msc_SOL_length_array(1)];
    
    MTU_strain_array = [ ...
        angle_array ...
        (AT_SOL_length_array-AT_SOL_length_array(1)) / AT_SOL_length_array(1)*100 ...
        (AT_GM_length_array-AT_GM_length_array(1)) / AT_GM_length_array(1)*100 ...       % ((at_GM_length_array-at_SOL_length_array) - (at_GM_length_array(1)-at_SOL_length_array(1))) / (at_GM_length_array(1)-at_SOL_length_array(1))*100     ...
        (MTU_GM_length-MTU_GM_length(1)) / MTU_GM_length(1)*100 ...
        zeros(length(angle_array),1) ...                          % GMFAS track has no length
        (apo_GM_array-apo_GM_array(1)) / apo_GM_array(1)*100 ...
        (msc_GM_length_array-msc_GM_length_array(1)) / msc_GM_length_array(1)*100 ...
        (msc_SOL_length_array-msc_SOL_length_array(1)) / msc_SOL_length_array(1)*100 ];
    
    
end