%%%%%%%%%%%%%%%%%%%%%%%%%%
% average_passive_trials KNEE
% Marie Moltubakk 6.2.2015 + 16-11-2017
% Read two trials of torque, norm angle
% Produce array with average data
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [output_array] = average_passive_trials_knee(varargin) % (torque1, gonio1, torque2, gonio2)

if nargin == 4 % two trials submitted, to be averaged
    torque1 = varargin{1};
    gonio1 = varargin{2};
    torque2 = varargin{3};
    gonio2 = varargin{4};

    
    %% if gonio angle starts after zero degrees, extrapolate data
%     
%     % trial 1
%     if min(gonio1) > 0.0
%         % select first 2 degrees of existing gonio data
%         angle_end = min(gonio1) + 2; %VAR extrapolate using first 2 degrees of the gonio movement
%         loc_angle_end = find(gonio1 >= angle_end,1,'first');
%         % calculate linear coeffisients
%         gonio_modified = [ones(length(gonio1(1:loc_angle_end)),1) gonio1(1:loc_angle_end)];
%         coeffs_torque = gonio_modified\torque1(1:loc_angle_end); % NB backslash --> slope
%         coeffs_displ = gonio_modified\displ1(1:loc_angle_end); % NB backslash --> slope
%         % enlarge arrays by adding values at zero angle to the front
%         angle_half = min(gonio1)/2;
%         gonio1 = [0; angle_half; gonio1];
%         torque1 = [coeffs_torque(1); (angle_half*coeffs_torque(2)) + coeffs_torque(1); torque1];
%         displ1 = [coeffs_displ(1); (angle_half*coeffs_displ(2)) + coeffs_displ(1); displ1];
%     end
%     
%     % repeat for trial 2
%     if min(gonio2) > 0.0
%         % select first 2 degrees of existing gonio data
%         angle_end = min(gonio2) + 2; %VAR extrapolate using first 2 degrees of the gonio movement
%         loc_angle_end = find(gonio2 >= angle_end,1,'first');
%         % calculate linear coeffisients
%         gonio_modified = [ones(length(gonio2(1:loc_angle_end)),1) gonio2(1:loc_angle_end)];
%         coeffs_torque = gonio_modified\torque2(1:loc_angle_end); % NB backslash --> slope
%         coeffs_displ = gonio_modified\displ2(1:loc_angle_end); % NB backslash --> slope
%         % enlarge arrays by adding values at zero angle to the front
%         angle_half = min(gonio2)/2;
%         gonio2 = [0; angle_half; gonio2];
%         torque2 = [coeffs_torque(1); (angle_half*coeffs_torque(2)) + coeffs_torque(1); torque2];
%         displ2 = [coeffs_displ(1); (angle_half*coeffs_displ(2)) + coeffs_displ(1); displ2];
% 
%     end
    
        
    %% create array with common angles (.05 intervals) for resampling
    
    % select the smallest range of the gonio data as basis for angle array
    common_angle_stop = max([min(gonio1) min(gonio2)]);
    common_angle_start = min([max(gonio1) max(gonio2)]);
    % round to .05 
    common_angle_start_rounded = round(common_angle_start/0.05)*0.05;
    % create array of angles
    average_angle_array = common_angle_start_rounded:-0.05:round(common_angle_stop/0.05)*0.05;

    % chop off ends if rounded away from max values
    if average_angle_array(1) > common_angle_start
        average_angle_array(1) = [];
    end
    if average_angle_array(end) < common_angle_stop
        average_angle_array(end) = [];
    end
    
        
    %% reshape, resample, average
    
    % delete potential duplicate values from gonio arrays
    % http://www.mathworks.com/matlabcentral/fileexchange/8354-consolidator
    [gonio1unique,torque1unique,ind] = consolidator(gonio1,torque1,'mean',0);
    [gonio2unique,torque2unique,ind] = consolidator(gonio2,torque2,'mean',0);

    % reshape and average TORQUE across common angle array
    common_torque1_gonio = spline(gonio1unique, torque1unique, average_angle_array);
    common_torque2_gonio = spline(gonio2unique, torque2unique, average_angle_array);
    average_torque_gonio = (common_torque1_gonio + common_torque2_gonio) / 2;

    
    
    
else %nargin 2 - one trial
    torque1 = varargin{1};
    gonio1 = varargin{2};
   

    %% if gonio angle starts after zero degrees, extrapolate data
%     % trial 1
%     if min(gonio1) > 0.0
%         % select first 2 degrees of existing gonio data
%         angle_end = min(gonio1) + 2; %VAR extrapolate using first 2 degrees of the gonio movement
%         loc_angle_end = find(gonio1 >= angle_end,1,'first');
%         % calculate linear coeffisients
%         gonio_modified = [ones(length(gonio1(1:loc_angle_end)),1) gonio1(1:loc_angle_end)];
%         coeffs_torque = gonio_modified\torque1(1:loc_angle_end); % NB backslash --> slope
%         coeffs_displ = gonio_modified\displ1(1:loc_angle_end); % NB backslash --> slope
%         % enlarge arrays by adding values at zero angle to the front
%         gonio1 = [0; gonio1];
%         torque1 = [coeffs_torque(1); torque1];
%         displ1 = [coeffs_displ(1); displ1];
%     end

    
    %% create array with common angles (.05 intervals) for resampling
    
    % select the smallest range of the gonio data as basis for angle array
    common_angle_start = max(gonio1); % high angle = very bent knee
    common_angle_stop = min(gonio1); % lower angle = extended knee
    % round to .05 
    common_angle_start_rounded = round(common_angle_start/0.05)*0.05;
    % create array of angles
    average_angle_array = common_angle_start_rounded:-0.05:round(common_angle_stop/0.05)*0.05;

    % chop off ends if rounded away from max values
    if average_angle_array(1) > common_angle_start
        average_angle_array(1) = [];
    end
    if average_angle_array(length(average_angle_array)) < common_angle_stop
        average_angle_array(length(average_angle_array)) = [];
    end

    
    %% reshape, resample, average
    
    % delete potential duplicate values from gonio arrays
    % http://www.mathworks.com/matlabcentral/fileexchange/8354-consolidator
    [gonio1unique,torque1unique,~] = consolidator(gonio1,torque1,'mean',0);

    % reshape and average TORQUE across common angle array
    average_torque_gonio = spline(gonio1unique, torque1unique, average_angle_array);

end


output_array = [average_torque_gonio; average_angle_array]';

end
