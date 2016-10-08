%%%%%%%%%%%%%%%%%%%%%%%%%%
% average_passive_forces + EMG
% Marie Moltubakk + REB + M 25.5.2015
% Read 3 trials of passive force and gonio angle
% Produce array with average force/angle
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [output_array] = average_passive_forces_EMG(data_SOL, data_GMMTJ, data_GMFAS)

    %global plot_achilles plot_norm plot_emg plot_check plot_us subject_id
    
    % no need to extrapolate data to zero degrees in this method, because the input variables (data_SOL etc are already extrapolated)
    
    placement_force = 1;
    placement_gonio = 2;
    % placement_displ = 3;
    placement_emg_gm = 4;
    placement_emg_gl = 5;
    placement_emg_sol = 6;

    % select the smallest range of the gonio data as basis for angle array
    common_angle_start = max([min(data_SOL(:,placement_gonio)) min(data_GMMTJ(:,placement_gonio)) min(data_GMFAS(:,placement_gonio))]);
    common_angle_stop = min([max(data_SOL(:,placement_gonio)) max(data_GMMTJ(:,placement_gonio)) max(data_GMFAS(:,placement_gonio))]);
    average_angle_array = common_angle_start:0.05:common_angle_stop;

    % reshape and average FORCE across common angle array
    common_force_gonio_SOL = spline(data_SOL(:,placement_gonio), data_SOL(:,placement_force), average_angle_array);
    common_force_gonio_GMMTJ = spline(data_GMMTJ(:,placement_gonio), data_GMMTJ(:,placement_force), average_angle_array);
    common_force_gonio_GMFAS = spline(data_GMFAS(:,placement_gonio), data_GMFAS(:,placement_force), average_angle_array);
    average_force_gonio = (common_force_gonio_SOL + common_force_gonio_GMMTJ + common_force_gonio_GMFAS) / 3;
    
    % reshape and average EMG sol across common angle array
    common_emg_sol_gonio_SOL = spline(data_SOL(:,placement_gonio), data_SOL(:,placement_emg_sol), average_angle_array);
    common_emg_sol_gonio_GMMTJ = spline(data_GMMTJ(:,placement_gonio), data_GMMTJ(:,placement_emg_sol), average_angle_array);
    common_emg_sol_gonio_GMFAS = spline(data_GMFAS(:,placement_gonio), data_GMFAS(:,placement_emg_sol), average_angle_array);
    average_emg_sol_gonio = (common_emg_sol_gonio_SOL + common_emg_sol_gonio_GMMTJ + common_emg_sol_gonio_GMFAS) / 3;
    
    % reshape and average EMG sol across common angle array
    common_emg_gm_gonio_SOL = spline(data_SOL(:,placement_gonio), data_SOL(:,placement_emg_gm), average_angle_array);
    common_emg_gm_gonio_GMMTJ = spline(data_GMMTJ(:,placement_gonio), data_GMMTJ(:,placement_emg_gm), average_angle_array);
    common_emg_gm_gonio_GMFAS = spline(data_GMFAS(:,placement_gonio), data_GMFAS(:,placement_emg_gm), average_angle_array);
    average_emg_gm_gonio = (common_emg_gm_gonio_SOL + common_emg_gm_gonio_GMMTJ + common_emg_gm_gonio_GMFAS) / 3;
    
    % reshape and average EMG sol across common angle array
    common_emg_gl_gonio_SOL = spline(data_SOL(:,placement_gonio), data_SOL(:,placement_emg_gl), average_angle_array);
    common_emg_gl_gonio_GMMTJ = spline(data_GMMTJ(:,placement_gonio), data_GMMTJ(:,placement_emg_gl), average_angle_array);
    common_emg_gl_gonio_GMFAS = spline(data_GMFAS(:,placement_gonio), data_GMFAS(:,placement_emg_gl), average_angle_array);
    average_emg_gl_gonio = (common_emg_gl_gonio_SOL + common_emg_gl_gonio_GMMTJ + common_emg_gl_gonio_GMFAS) / 3;
    
    output_array = rot90([average_emg_sol_gonio; average_emg_gl_gonio; average_emg_gm_gonio; average_angle_array; average_force_gonio],3);
    % rotate 270 degrees --> (force, gonio, (no displ), gm, gl, sol)
end
