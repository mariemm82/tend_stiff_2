%%%%%%%%%%%%%%%%%%%%%%%%%%
% average_passive_trials_licht
% Marie Moltubakk 1.10.2015
% Read two trials of lichtwark data + Norm data from torque, gonio angle, norm angle, displacement
% Produce array with average data
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [output_array] = average_passive_trials_licht(varargin) % (gonio1, angle1, faslen1, pennation1, gonio2, angle2, faslen2, pennation2, description)

global plot_check subject_id plot_licht

if nargin == 9 % two trials submitted, to be averaged
    gonio1 = varargin{1};
    %angle1 = varargin{2};
    faslen1 = varargin{3};
    pennation1 = varargin{4};
    gonio2 = varargin{5};
    %angle2 = varargin{6};
    faslen2 = varargin{7};
    pennation2 = varargin{8};
    description = varargin{9};
    
    
    
    %%% if gonio angle starts after zero degrees, extrapolate data
    
    % trial 1
    if min(gonio1) > 0.0
        % select first 2 degrees of existing gonio data
        angle_end = min(gonio1) + 2; %VAR extrapolate using first 2 degrees of the gonio movement
        loc_angle_end = find(gonio1 >= angle_end,1,'first');
        % calculate linear coeffisients
        gonio_modified = [ones(length(gonio1(1:loc_angle_end)),1) gonio1(1:loc_angle_end)];
        coeffs_faslen = gonio_modified\faslen1(1:loc_angle_end); % NB backslash --> slope
        coeffs_pennation = gonio_modified\pennation1(1:loc_angle_end); % NB backslash --> slope
        % enlarge arrays by adding values at zero angle to the front
        angle_half = min(gonio1)/2;
        gonio1 = [0; angle_half; gonio1];
        faslen1 = [coeffs_faslen(1); (angle_half*coeffs_faslen(2)) + coeffs_faslen(1); faslen1];
        pennation1 = [coeffs_pennation(1); (angle_half*coeffs_pennation(2)) + coeffs_pennation(1); pennation1];
    end
    
    % repeat for trial 2
    if min(gonio2) > 0.0
        % select first 2 degrees of existing gonio data
        angle_end = min(gonio2) + 2; %VAR extrapolate using first 2 degrees of the gonio movement
        loc_angle_end = find(gonio2 >= angle_end,1,'first');
        % calculate linear coeffisients
        gonio_modified = [ones(length(gonio2(1:loc_angle_end)),1) gonio2(1:loc_angle_end)];
        coeffs_faslen = gonio_modified\faslen2(1:loc_angle_end); % NB backslash --> slope
        coeffs_pennation = gonio_modified\pennation2(1:loc_angle_end); % NB backslash --> slope
        % enlarge arrays by adding values at zero angle to the front
        angle_half = min(gonio2)/2;
        gonio2 = [0; angle_half; gonio2];
        faslen2 = [coeffs_faslen(1); (angle_half*coeffs_faslen(2)) + coeffs_faslen(1); faslen2];
        pennation2 = [coeffs_pennation(1); (angle_half*coeffs_pennation(2)) + coeffs_pennation(1); pennation2];
    end
    
    
    
    %%% create array with common angles (.05 intervals) for resampling
    
    % select the smallest range of the gonio data as basis for angle array
    common_angle_start = max([min(gonio1) min(gonio2)]);
    common_angle_stop = min([max(gonio1) max(gonio2)]);
    % round to .05 
    common_angle_start_rounded = round(common_angle_start/0.05)*0.05;
    % create array of angles
    average_angle_array = common_angle_start_rounded:0.05:round(common_angle_stop/0.05)*0.05;

    % chop off ends if rounded away from max values
    if average_angle_array(1) < common_angle_start
        average_angle_array(1) = [];
    end
    if average_angle_array(length(average_angle_array)) > common_angle_stop
        average_angle_array(length(average_angle_array)) = [];
    end

    
    
    %%% reshape, resample, average
    
    % delete potential duplicate values from gonio arrays
    % http://www.mathworks.com/matlabcentral/fileexchange/8354-consolidator
    [~,faslen1unique,~] = consolidator(gonio1,faslen1,'mean',0);
    [~,faslen2unique,~] = consolidator(gonio2,faslen2,'mean',0);
    [gonio1unique,pennation1unique,~] = consolidator(gonio1,pennation1,'mean',0);
    [gonio2unique,pennation2unique,~] = consolidator(gonio2,pennation2,'mean',0);

    % reshape FASLEN across common angle array
    % if necessary, will extrapolate to include zero degrees gonio angle
    common_faslen1_gonio = spline(gonio1unique, faslen1unique, average_angle_array);
    common_faslen2_gonio = spline(gonio2unique, faslen2unique, average_angle_array);
    
    % average FASLEN across common angle array
    average_faslen_gonio = (common_faslen1_gonio + common_faslen2_gonio) / 2;

    % reshape PENNATION across common angle array
    % if necessary, will extrapolate to include zero degrees gonio angle
    common_pennation1_gonio = spline(gonio1unique, pennation1unique, average_angle_array);
    common_pennation2_gonio = spline(gonio2unique, pennation2unique, average_angle_array);
    
    % average PENNATION across common angle array
    average_pennation_gonio = (common_pennation1_gonio + common_pennation2_gonio) / 2;
    
    
    
    %%% plot original and averaged data
    if plot_check && plot_licht
        plottitle = horzcat('Lichtwark averaging, ', subject_id, ' - ', description);
        figure('Name',plottitle);
        max_faslen = max([max(faslen1) max(faslen2)]);
        min_faslen = min([min(faslen1) min(faslen2)]);
        max_pennation = max([max(pennation1) max(pennation2)]);
        min_pennation = min([min(pennation1) min(pennation2)]);
        % top panel = fas len
        AXa = subplot(2,1,1);
        plot(gonio1,faslen1)
        hold on
        plot(gonio2,faslen2)
        plot(average_angle_array,average_faslen_gonio,'LineWidth',2)
        set(get(AXa,'Ylabel'),'String','Fascicle length (mm)')
        set(AXa,'YLim',[min_faslen*.9 max_faslen*1.1])
        title(plottitle)
        % bottom panel = pennation angle
        AXc = subplot(2,1,2);
        plot(gonio1,pennation1)
        hold on
        plot(gonio2,pennation2)
        plot(average_angle_array,average_pennation_gonio,'LineWidth',2)
        set(get(AXc,'Ylabel'),'String','Pennation angle (deg)')
        set(AXc,'YLim',[min_pennation*.9 max_pennation*1.1])
        xlabel('Gonio angle (deg)');
        legend('Trial 1','Trial 2','Spline','Location','Southeast');
    end
    
    
    
    
    
    
else % nargin == 5, e.g. only one trial submitted
    
    gonio1 = varargin{1};
    %angle1 = varargin{2};
    faslen1 = varargin{3};
    pennation1 = varargin{4};
    description = varargin{5};
    
    
    
    %%% if gonio angle starts after zero degrees, extrapolate data
    
    % trial 1
    if min(gonio1) > 0.0
        % select first 2 degrees of existing gonio data
        angle_end = min(gonio1) + 2; %VAR extrapolate using first 2 degrees of the gonio movement
        loc_angle_end = find(gonio1 >= angle_end,1,'first');
        % calculate linear coeffisients
        gonio_modified = [ones(length(gonio1(1:loc_angle_end)),1) gonio1(1:loc_angle_end)];
        coeffs_faslen = gonio_modified\faslen1(1:loc_angle_end); % NB backslash --> slope
        coeffs_pennation = gonio_modified\pennation1(1:loc_angle_end); % NB backslash --> slope
        % enlarge arrays by adding values at zero angle to the front
        angle_half = min(gonio1)/2;
        gonio1 = [0; angle_half; gonio1];
        faslen1 = [coeffs_faslen(1); (angle_half*coeffs_faslen(2)) + coeffs_faslen(1); faslen1];
        pennation1 = [coeffs_pennation(1); (angle_half*coeffs_pennation(2)) + coeffs_pennation(1); pennation1];
    end

    
    
    %%% create array with common angles (.05 intervals) for resampling
    
    % select the smallest range of the gonio data as basis for angle array
    common_angle_start = min(gonio1);
    common_angle_stop = max(gonio1);
    % round to .05 
    common_angle_start_rounded = round(common_angle_start/0.05)*0.05;
    % create array of angles
    average_angle_array = common_angle_start_rounded:0.05:round(common_angle_stop/0.05)*0.05;

    % chop off ends if rounded away from max values
    if average_angle_array(1) < common_angle_start
        average_angle_array(1) = [];
    end
    if average_angle_array(length(average_angle_array)) > common_angle_stop
        average_angle_array(length(average_angle_array)) = [];
    end

    
    
    %%% reshape, resample, average
    
    % delete potential duplicate values from gonio arrays
    % http://www.mathworks.com/matlabcentral/fileexchange/8354-consolidator
    [~,faslen1unique,~] = consolidator(gonio1,faslen1,'mean',0);
    [gonio1unique,pennation1unique,~] = consolidator(gonio1,pennation1,'mean',0);

    % reshape FASLEN across common angle array
    average_faslen_gonio = spline(gonio1unique, faslen1unique, average_angle_array);

    % reshape PENNATION across common angle array
    average_pennation_gonio = spline(gonio1unique, pennation1unique, average_angle_array);
    
    % plot original and averaged data
    if plot_check && plot_licht
        plottitle = horzcat('Lichtwark averaging, ', subject_id, ' - ', description);
        figure('Name',plottitle);
        max_faslen = max(faslen1);
        min_faslen = min(faslen1);
        max_pennation = max(pennation1);
        min_pennation = min(pennation1);
        % top panel = fas len
        AXa = subplot(2,1,1);
        plot(gonio1,faslen1)
        hold on
        plot(average_angle_array,average_faslen_gonio,'LineWidth',2)
        set(get(AXa,'Ylabel'),'String','Fascicle length (mm)')
        set(AXa,'YLim',[min_faslen*.9 max_faslen*1.1])
        title(plottitle)
        AXc = subplot(2,1,2);
        plot(gonio1,pennation1)
        hold on
        plot(average_angle_array,average_pennation_gonio,'LineWidth',2)
        set(get(AXc,'Ylabel'),'String','Pennation angle (deg)')
        set(AXc,'YLim',[min_pennation*.9 max_pennation*1.1])
        xlabel('Gonio angle (deg)');
        legend('Trial','Spline','Location','Southeast');
    end
    
end





output_array = rot90([average_pennation_gonio; average_faslen_gonio; average_angle_array],3);
% rotate 270 degrees --> angle, faslen, pennation

end






% OLD

%         %%% create array with common angles (.05 intervals) for resampling
% 
%         % select the smallest range of the gonio data as basis for angle array
%         common_angle_start = max([min(gonio1) min(gonio2)]);
%         common_angle_stop = min([max(gonio1) max(gonio2)]);
%         % round to .05 & create array of angles
%         if common_angle_start > 0
%             % if zero angle does not exist, spline must create it, so array must start from zero
%             common_angle_start_rounded = 0;
%         else
%             common_angle_start_rounded = round(common_angle_start/0.05)*0.05;
%         end
%         average_angle_array = common_angle_start_rounded:0.05:round(common_angle_stop/0.05)*0.05;
%         % old version starting after zero if data don't exist from zero:
%     %     average_angle_array = round(common_angle_start/0.05)*0.05:0.05:round(common_angle_stop/0.05)*0.05;
% 
%         % chop off ends if rounded away from max values
%     %    if average_angle_array(1) < common_angle_start
%     %        average_angle_array(1) = [];
%     %    end
%         if average_angle_array(length(average_angle_array)) > common_angle_stop
%             average_angle_array(length(average_angle_array)) = [];
%         end