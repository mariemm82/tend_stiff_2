%%%%%%%%%%%%%%%%%%%%%%%%%%
% average_passive_forces
% Marie Moltubakk + REB + M 25.5.2015
% Read 3 trials of passive force and gonio angle
% Produce array with average force/angle
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [output_array] = average_passive_forces(data_sol, data_gmmtj, data_gmfas)

    % no need to extrapolate data to zero degrees in this method, because the input variables (data_SOL etc are already extrapolated)

    placement_force = 1;
    placement_gonio = 2;

    % select the smallest range of the gonio data as basis for angle array
    common_angle_start = max([min(data_sol(:,placement_gonio)) min(data_gmmtj(:,placement_gonio)) min(data_gmfas(:,placement_gonio))]);
    common_angle_stop = min([max(data_sol(:,placement_gonio)) max(data_gmmtj(:,placement_gonio)) max(data_gmfas(:,placement_gonio))]);
    average_angle_array = common_angle_start:0.05:common_angle_stop;

    % reshape and average force across common angle array
    common_force_gonio_sol = spline(data_sol(:,placement_gonio), data_sol(:,placement_force), average_angle_array);
    common_force_gonio_gmmtj = spline(data_gmmtj(:,placement_gonio), data_gmmtj(:,placement_force), average_angle_array);
    common_force_gonio_gmfas = spline(data_gmfas(:,placement_gonio), data_gmfas(:,placement_force), average_angle_array);
    average_force_gonio = (common_force_gonio_sol + common_force_gonio_gmmtj + common_force_gonio_gmfas) / 3;

    output_array = [average_force_gonio; average_angle_array]';

end
