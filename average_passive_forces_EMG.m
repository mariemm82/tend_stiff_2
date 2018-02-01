%%%%%%%%%%%%%%%%%%%%%%%%%%
% average_passive_forces + EMG
% Marie Moltubakk + REB + M 25.5.2015
% Read 3 trials of passive force and gonio angle
% Produce array with average force/angle
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [output_array] = average_passive_forces_EMG(trials_SOL, trials_GMMTJ, trials_GMFAS, subject_nr)
    global subject_id
    
    angle_step = 0.05; % VAR - reshaped, averaged data extracted every x degrees
    smoother = 10; %VAR 
    
    % no need to extrapolate data to zero degrees in this method, because the input variables (data_SOL etc are already extrapolated)
    
    placement_force = 1;
    placement_gonio = 2;
    % placement_displ = 3;
    placement_torque = 4;
    placement_emg_gm = 5;
    placement_emg_gl = 6;
    placement_emg_sol = 7;
    
    %% common angles, spline, average 
    
    % select the smallest range of the gonio data as basis for angle array
    common_angle_start = max([min(trials_SOL(:,placement_gonio)) min(trials_GMMTJ(:,placement_gonio)) min(trials_GMFAS(:,placement_gonio))]);
    common_angle_stop = min([max(trials_SOL(:,placement_gonio)) max(trials_GMMTJ(:,placement_gonio)) max(trials_GMFAS(:,placement_gonio))]);
    average_angle_array = (common_angle_start:angle_step:common_angle_stop)';

    % reshape and average FORCE across common angle array
    common_force_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_force), average_angle_array);
    common_force_gonio_GMMTJ = spline(trials_GMMTJ(:,placement_gonio), trials_GMMTJ(:,placement_force), average_angle_array);
    common_force_gonio_GMFAS = spline(trials_GMFAS(:,placement_gonio), trials_GMFAS(:,placement_force), average_angle_array);
    average_force_gonio = (smooth(common_force_gonio_SOL,smoother) + smooth(common_force_gonio_GMMTJ,smoother) + smooth(common_force_gonio_GMFAS,smoother)) / 3;
    
    % reshape and average TORQUE across common angle array
    common_torque_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_torque), average_angle_array);
    common_torque_gonio_GMMTJ = spline(trials_GMMTJ(:,placement_gonio), trials_GMMTJ(:,placement_torque), average_angle_array);
    common_torque_gonio_GMFAS = spline(trials_GMFAS(:,placement_gonio), trials_GMFAS(:,placement_torque), average_angle_array);
    average_torque_gonio = (smooth(common_torque_gonio_SOL,smoother) + smooth(common_torque_gonio_GMMTJ,smoother) + smooth(common_torque_gonio_GMFAS,smoother)) / 3;
    
    
    %% handle subjects with erroraneous / missing EMG data
    
    % subject without EMG
    
    if strcmp(subject_id,'INT_15_CON_PRE_L') || strcmp(subject_id,'INT_31_CON_PRE_L') || strcmp(subject_id,'INT_31_STR_PRE_R')
        average_emg_sol_gonio(1:length(average_angle_array),1) = NaN;
        average_emg_gm_gonio(1:length(average_angle_array),1) = NaN;
        average_emg_gl_gonio(1:length(average_angle_array),1) = NaN;
        
%       elseif subject_nr == 110 % BD
%         average_emg_sol_gonio(1:length(average_angle_array),1) = NaN;
%         average_emg_gm_gonio(1:length(average_angle_array),1) = NaN;
%         average_emg_gl_gonio(1:length(average_angle_array),1) = NaN;
          
    % subjects with partial EMG 
        
    elseif strcmp(subject_id,'INT_29_CON_PRE_L') 
        % 29 - use only SOL trials
        % reshape and average EMG SOL across common angle array
        common_emg_sol_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_emg_sol), average_angle_array);
        average_emg_sol_gonio = smooth(common_emg_sol_gonio_SOL,smoother);
        % reshape and average EMG GM across common angle array
        common_emg_gm_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_emg_gm), average_angle_array);
        average_emg_gm_gonio = smooth(common_emg_gm_gonio_SOL,smoother);
        % reshape and average EMG GL across common angle array
        common_emg_gl_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_emg_gl), average_angle_array);
        average_emg_gl_gonio = smooth(common_emg_gl_gonio_SOL,smoother);
        
    elseif strcmp(subject_id,'INT_22_CON_PRE_R') 
        % 22	PRE	R	CON - use only GMMTJ EMG
        % reshape and average EMG SOL across common angle array
        common_emg_sol_gonio_GMMTJ = spline(trials_GMMTJ(:,placement_gonio), trials_GMMTJ(:,placement_emg_sol), average_angle_array);
        average_emg_sol_gonio = smooth(common_emg_sol_gonio_GMMTJ,smoother);
        % reshape and average EMG GM across common angle array
        common_emg_gm_gonio_GMMTJ = spline(trials_GMMTJ(:,placement_gonio), trials_GMMTJ(:,placement_emg_gm), average_angle_array);
        average_emg_gm_gonio = smooth(common_emg_gm_gonio_GMMTJ,smoother);
        % reshape and average EMG GL across common angle array
        common_emg_gl_gonio_GMMTJ = spline(trials_GMMTJ(:,placement_gonio), trials_GMMTJ(:,placement_emg_gl), average_angle_array);
        average_emg_gl_gonio = smooth(common_emg_gl_gonio_GMMTJ,smoother);
        
    elseif strcmp(subject_id,'INT_22_STR_PRE_L') 
        % 22	PRE	L - don't use EMG from GMFAS
        % reshape and average EMG SOL across common angle array
        common_emg_sol_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_emg_sol), average_angle_array);
        common_emg_sol_gonio_GMMTJ = spline(trials_GMMTJ(:,placement_gonio), trials_GMMTJ(:,placement_emg_sol), average_angle_array);
        average_emg_sol_gonio = (smooth(common_emg_sol_gonio_SOL,smoother) + smooth(common_emg_sol_gonio_GMMTJ,smoother)) / 2;

        % reshape and average EMG GM across common angle array
        common_emg_gm_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_emg_gm), average_angle_array);
        common_emg_gm_gonio_GMMTJ = spline(trials_GMMTJ(:,placement_gonio), trials_GMMTJ(:,placement_emg_gm), average_angle_array);
        average_emg_gm_gonio = (smooth(common_emg_gm_gonio_SOL,smoother) + smooth(common_emg_gm_gonio_GMMTJ,smoother)) / 2;

        % reshape and average EMG GL across common angle array
        common_emg_gl_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_emg_gl), average_angle_array);
        common_emg_gl_gonio_GMMTJ = spline(trials_GMMTJ(:,placement_gonio), trials_GMMTJ(:,placement_emg_gl), average_angle_array);
        average_emg_gl_gonio = (smooth(common_emg_gl_gonio_SOL,smoother) + smooth(common_emg_gl_gonio_GMMTJ,smoother)) / 2;
    
    elseif strcmp(subject_id,'INT_21_STR_PRE_R') 
        % INT_21_STR_PRE_R - don't use GMMTJ
        % reshape and average EMG SOL across common angle array
        common_emg_sol_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_emg_sol), average_angle_array);
        common_emg_sol_gonio_GMFAS = spline(trials_GMFAS(:,placement_gonio), trials_GMFAS(:,placement_emg_sol), average_angle_array);
        average_emg_sol_gonio = (smooth(common_emg_sol_gonio_SOL,smoother) + smooth(common_emg_sol_gonio_GMFAS,smoother)) / 2;

        % reshape and average EMG GM across common angle array
        common_emg_gm_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_emg_gm), average_angle_array);
        common_emg_gm_gonio_GMFAS = spline(trials_GMFAS(:,placement_gonio), trials_GMFAS(:,placement_emg_gm), average_angle_array);
        average_emg_gm_gonio = (smooth(common_emg_gm_gonio_SOL,smoother) + smooth(common_emg_gm_gonio_GMFAS,smoother)) / 2;

        % reshape and average EMG GL across common angle array
        common_emg_gl_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_emg_gl), average_angle_array);
        common_emg_gl_gonio_GMFAS = spline(trials_GMFAS(:,placement_gonio), trials_GMFAS(:,placement_emg_gl), average_angle_array);
        average_emg_gl_gonio = (smooth(common_emg_gl_gonio_SOL,smoother) + smooth(common_emg_gl_gonio_GMFAS,smoother)) / 2;
    
        
%     elseif subject_nr == 106 % BD
%         % 106: use only GMFAS
%         % reshape and average EMG SOL across common angle array
%         common_emg_sol_gonio_GMFAS = spline(trials_GMFAS(:,placement_gonio), trials_GMFAS(:,placement_emg_sol), average_angle_array);
%         average_emg_sol_gonio = smooth(common_emg_sol_gonio_GMFAS,smoother);
%         % reshape and average EMG GM across common angle array
%         common_emg_gm_gonio_GMFAS = spline(trials_GMFAS(:,placement_gonio), trials_GMFAS(:,placement_emg_gm), average_angle_array);
%         average_emg_gm_gonio = smooth(common_emg_gm_gonio_GMFAS,smoother);
%         % reshape and average EMG GL across common angle array
%         common_emg_gl_gonio_GMFAS = spline(trials_GMFAS(:,placement_gonio), trials_GMFAS(:,placement_emg_gl), average_angle_array);
%         average_emg_gl_gonio = smooth(common_emg_gl_gonio_GMFAS,smoother);

        
      else
        
    % normal subjects
    
        % reshape and average EMG SOL across common angle array
        common_emg_sol_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_emg_sol), average_angle_array);
        common_emg_sol_gonio_GMMTJ = spline(trials_GMMTJ(:,placement_gonio), trials_GMMTJ(:,placement_emg_sol), average_angle_array);
        common_emg_sol_gonio_GMFAS = spline(trials_GMFAS(:,placement_gonio), trials_GMFAS(:,placement_emg_sol), average_angle_array);
        average_emg_sol_gonio = (smooth(common_emg_sol_gonio_SOL,smoother) + smooth(common_emg_sol_gonio_GMMTJ,smoother) + smooth(common_emg_sol_gonio_GMFAS,smoother)) / 3;

        % reshape and average EMG GM across common angle array
        common_emg_gm_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_emg_gm), average_angle_array);
        common_emg_gm_gonio_GMMTJ = spline(trials_GMMTJ(:,placement_gonio), trials_GMMTJ(:,placement_emg_gm), average_angle_array);
        common_emg_gm_gonio_GMFAS = spline(trials_GMFAS(:,placement_gonio), trials_GMFAS(:,placement_emg_gm), average_angle_array);
        average_emg_gm_gonio = (smooth(common_emg_gm_gonio_SOL,smoother) + smooth(common_emg_gm_gonio_GMMTJ,smoother) + smooth(common_emg_gm_gonio_GMFAS,smoother)) / 3;

        % reshape and average EMG GL across common angle array
        common_emg_gl_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_emg_gl), average_angle_array);
        common_emg_gl_gonio_GMMTJ = spline(trials_GMMTJ(:,placement_gonio), trials_GMMTJ(:,placement_emg_gl), average_angle_array);
        common_emg_gl_gonio_GMFAS = spline(trials_GMFAS(:,placement_gonio), trials_GMFAS(:,placement_emg_gl), average_angle_array);
        average_emg_gl_gonio = (smooth(common_emg_gl_gonio_SOL,smoother) + smooth(common_emg_gl_gonio_GMMTJ,smoother) + smooth(common_emg_gl_gonio_GMFAS,smoother)) / 3;
    
    end
    
    %% output variable
    output_array = [average_force_gonio average_angle_array average_torque_gonio average_emg_gm_gonio average_emg_gl_gonio average_emg_sol_gonio];
end
