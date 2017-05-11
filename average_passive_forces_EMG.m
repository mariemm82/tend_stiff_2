%%%%%%%%%%%%%%%%%%%%%%%%%%
% average_passive_forces + EMG
% Marie Moltubakk + REB + M 25.5.2015
% Read 3 trials of passive force and gonio angle
% Produce array with average force/angle
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [output_array] = average_passive_forces_EMG(trials_SOL, trials_GMMTJ, trials_GMFAS, subject_nr)

    %global plot_achilles plot_norm plot_emg plot_check plot_us subject_id
    angle_step = 0.05; % VAR - reshaped, averaged data extracted every x degrees
    smoother = 10; %VAR 
    
    % no need to extrapolate data to zero degrees in this method, because the input variables (data_SOL etc are already extrapolated)
    
    placement_force = 1;
    placement_gonio = 2;
    % placement_displ = 3;
    placement_emg_gm = 4;
    placement_emg_gl = 5;
    placement_emg_sol = 6;

    % select the smallest range of the gonio data as basis for angle array
    common_angle_start = max([min(trials_SOL(:,placement_gonio)) min(trials_GMMTJ(:,placement_gonio)) min(trials_GMFAS(:,placement_gonio))]);
    common_angle_stop = min([max(trials_SOL(:,placement_gonio)) max(trials_GMMTJ(:,placement_gonio)) max(trials_GMFAS(:,placement_gonio))]);
    average_angle_array = (common_angle_start:angle_step:common_angle_stop)';

    % reshape and average FORCE across common angle array
    common_force_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_force), average_angle_array);
    common_force_gonio_GMMTJ = spline(trials_GMMTJ(:,placement_gonio), trials_GMMTJ(:,placement_force), average_angle_array);
    common_force_gonio_GMFAS = spline(trials_GMFAS(:,placement_gonio), trials_GMFAS(:,placement_force), average_angle_array);
    average_force_gonio = (smooth(common_force_gonio_SOL,smoother) + smooth(common_force_gonio_GMMTJ,smoother) + smooth(common_force_gonio_GMFAS,smoother)) / 3;
    
    % handle subjects with erroraneous / missing EMG data
    
    % if subject without EMG (31, 110)
    if subject_nr == 31 || subject_nr == 110
        average_emg_sol_gonio(1:length(average_angle_array),1) = NaN;
        average_emg_gm_gonio(1:length(average_angle_array),1) = NaN;
        average_emg_gl_gonio(1:length(average_angle_array),1) = NaN;
        
    % if subject with partial EMG 
    elseif subject_nr == 106 % (106 - use only GMFAS)
        % reshape and average EMG SOL across common angle array
        common_emg_sol_gonio_GMFAS = spline(trials_GMFAS(:,placement_gonio), trials_GMFAS(:,placement_emg_sol), average_angle_array);
        average_emg_sol_gonio = smooth(common_emg_sol_gonio_GMFAS,smoother);
        % reshape and average EMG GM across common angle array
        common_emg_gm_gonio_GMFAS = spline(trials_GMFAS(:,placement_gonio), trials_GMFAS(:,placement_emg_gm), average_angle_array);
        average_emg_gm_gonio = smooth(common_emg_gm_gonio_GMFAS,smoother);
        % reshape and average EMG GL across common angle array
        common_emg_gl_gonio_GMFAS = spline(trials_GMFAS(:,placement_gonio), trials_GMFAS(:,placement_emg_gl), average_angle_array);
        average_emg_gl_gonio = smooth(common_emg_gl_gonio_GMFAS,smoother);
    
    elseif subject_nr == 29 % (29 - use only SOL trials)
        % reshape and average EMG SOL across common angle array
        common_emg_sol_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_emg_sol), average_angle_array);
        average_emg_sol_gonio = smooth(common_emg_sol_gonio_SOL,smoother);
        % reshape and average EMG GM across common angle array
        common_emg_gm_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_emg_gm), average_angle_array);
        average_emg_gm_gonio = smooth(common_emg_gm_gonio_SOL,smoother);
        % reshape and average EMG GL across common angle array
        common_emg_gl_gonio_SOL = spline(trials_SOL(:,placement_gonio), trials_SOL(:,placement_emg_gl), average_angle_array);
        average_emg_gl_gonio = smooth(common_emg_gl_gonio_SOL,smoother);
    
    else
    % if normal subject
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
    
    output_array = [average_force_gonio average_angle_array average_emg_gm_gonio average_emg_gl_gonio average_emg_sol_gonio];
end
