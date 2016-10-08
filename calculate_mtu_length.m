%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_mtu_length
% Marie Moltubakk 12.2.2015
% Read length of lower leg + ankle angle array
% Produce array with leg length across angle array
%%%%%%%%%%%%%%%%%%%%%%%%%%



function [MTU_length_array] = calculate_mtu_length(angle_displ_SOL, angle_displ_GMMTJ, angle_displ_GMFAS, initial_at_SOL_length, initial_at_GM_length, initial_calf_length, angle_common)



    %%%%%%%% calculate angle range common to all 6 trials
    % OLD:
    % angle_common = min([max(angle_displ_SOL(:,1)) max(angle_displ_GMMTJ(:,1)) max(angle_displ_GMFAS(:,1))]);
    % NEW:
    % angle_common is read from file from create_angles_passive, to avoid potential tiny discrepancies between various calculations of max
    angle_array = (0:0.05:(angle_common+0.01))'; %VAR  - adding 0.01 to ensure that the actual common angle is included - data stored as 10.9999999 



    %%%%%%%% Length of whole MTU - Grieve 1978

    %%% extract calf length in mm
    % initial calf length = centimeters from lateral epicondyle to calc insertion
    calf_length = str2double(initial_calf_length) * 10;

    %%% calculate percent change in MTU length, across angle array (degrees)
    % angle array = full ROM test (degrees)
    % my data report neutral position = 0 degrees, Grieve data report neutral position = 90 degrees
    a0 = -22.18468;
    a1 = 0.30141;
    a2 = -0.00061;
    percent_change_mtu = (a0 + (a1*(90+angle_array)) + (a2*(90+angle_array).^2))/100;

    %%% calf/MTU length (mm), across angle array
    calf_length_array = calf_length * (1 + percent_change_mtu);



    %%%%%%% Length of free Achilles tendon

    % Initial length = measured seated in Achilles machine, ankle 0 degrees, US still images of "railway" across lower leg. Length in mm.
    at_SOL_length = str2double(initial_at_SOL_length);
    
    % convert to length + reshape
    at_SOL_length_array = spline(angle_displ_SOL(:,1), (at_SOL_length + angle_displ_SOL(:,2)), angle_array);



    %%%%%%% Length of Achilles tendon up to GM insertion ("GM tendon")
    at_GM_length = str2double(initial_at_GM_length);
    
    % convert to length + reshape
    at_GM_length_array = spline(angle_displ_GMMTJ(:,1), (at_GM_length + angle_displ_GMMTJ(:,2)), angle_array);
    
    
    
    %%%%% GMFAS displacement to same angle array
    GMFAS_displ_array = spline(angle_displ_GMFAS(:,1), angle_displ_GMFAS(:,2), angle_array);
   
   
    %%%%%% Final data
    MTU_length_array = [angle_array at_SOL_length_array at_GM_length_array calf_length_array GMFAS_displ_array];
    
    
    
end