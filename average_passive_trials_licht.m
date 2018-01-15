%%%%%%%%%%%%%%%%%%%%%%%%%%
% average_passive_trials_licht
% Marie Moltubakk 1.10.2015
% Read two trials of lichtwark data + Norm data from torque, gonio angle, norm angle, displacement
% Produce array with average data
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [output_array] = average_passive_trials_licht(varargin) % (gonio1, angle1, faslen1, pennation1, time1, gonio2, angle2, faslen2, pennation2, time2, angle_common,, description)

global plot_check subject_id plot_licht
angle_step = 0.05; % VAR - reshaped, averaged data extracted every x degrees
smoother = 10; %VAR 
errormsg_extrapol = 0.5; 

if nargin == 12 % two trials submitted, to be averaged
    gonio1 = varargin{1};
    %angle1 = varargin{2};
    faslen1 = varargin{3};
    pennation1 = varargin{4};
    time1 = varargin{5};
    gonio2 = varargin{6};
    %angle2 = varargin{7};
    faslen2 = varargin{8};
    pennation2 = varargin{9};
    time2 = varargin{10};
    angle_common = varargin{11};
    description = varargin{12};
    
        
    %% curve fitting GONIO data to TIME array
    
    % method 1: fit 4th order polynomial to averaged gonio-angle curve
    fit_gonio1 = polyfit(time1, gonio1, 4);
    fit_gonio2 = polyfit(time2, gonio2, 4);
    
    % create array across times
    gonio_new1(1:length(time1),1) = zeros;
    for ang = 1:length(time1)
        gonio_new1(ang) = (fit_gonio1(1) * time1(ang)^4) + (fit_gonio1(2) * time1(ang)^3) + (fit_gonio1(3) * time1(ang)^2) + (fit_gonio1(4) * time1(ang)) + fit_gonio1(5);
    end
    gonio_new2(1:length(time2),1) = zeros;
    for ang = 1:length(time2)
        gonio_new2(ang) = (fit_gonio2(1) * time2(ang)^4) + (fit_gonio2(2) * time2(ang)^3) + (fit_gonio2(3) * time2(ang)^2) + (fit_gonio2(4) * time2(ang)) + fit_gonio2(5);
    end
    
    %    % method 2: cubic smoothing spline
    %    p = 0.0005; % 0 = least squares straight line fit, 1 = natural cubic spline interpolant   %VAR
    %    gonio_new = csaps(time1, gonio1, p, time_torque_ascend); 
    
        
    %% if new gonio angle is smaller than original max, extrapolate all data
    
    common_angle_orig = angle_common + 0.0001; % round(max([max(gonio1) max(gonio2)]),2);
    
    % trial 1
    if max(gonio_new1) < common_angle_orig
        if common_angle_orig - max(gonio_new1) > errormsg_extrapol
            cprintf('red', horzcat('WARNING: Extrapolating max ROM1 by ', num2str(common_angle_orig - max(gonio_new1)), ' degrees.\n'))
        end
        % extrapolate using last X degrees of existing data
        angle_end = max(gonio_new1) - 2; %VAR
        loc_angle_end = find(gonio_new1 >= angle_end,1,'first');
        % calculate linear coeffisients
        gonio_modifiedEND = [ones(length(gonio_new1(loc_angle_end:end)),1) gonio_new1(loc_angle_end:end)];
        coeffs_faslenEND = gonio_modifiedEND\faslen1(loc_angle_end:end); % NB backslash --> slope
        coeffs_pennationEND = gonio_modifiedEND\pennation1(loc_angle_end:end);
        % enlarge arrays by adding values at the end
        gonio_new1 = [gonio_new1; common_angle_orig];
        faslen1 = [faslen1; (common_angle_orig*coeffs_faslenEND(2)) + coeffs_faslenEND(1)];
        pennation1 = [pennation1; (common_angle_orig*coeffs_pennationEND(2)) + coeffs_pennationEND(1)];
    end
    
    % repeat for trial 2
    if max(gonio_new2) < common_angle_orig
        if common_angle_orig - max(gonio_new1) > errormsg_extrapol
            cprintf('red', horzcat('WARNING: Extrapolating max ROM2 by ', num2str(common_angle_orig - max(gonio_new1)), ' degrees.\n'))
        end
        % extrapolate using last X degrees of existing data
        angle_end = max(gonio_new2) - 2; %VAR
        loc_angle_end = find(gonio_new2 >= angle_end,1,'first');
        % calculate linear coeffisients
        gonio_modifiedEND = [ones(length(gonio_new2(loc_angle_end:end)),1) gonio_new2(loc_angle_end:end)];
        coeffs_faslenEND = gonio_modifiedEND\faslen2(loc_angle_end:end); % NB backslash --> slope
        coeffs_pennationEND = gonio_modifiedEND\pennation2(loc_angle_end:end);
        % enlarge arrays by adding values at the end
        gonio_new2 = [gonio_new2; common_angle_orig];
        faslen2 = [faslen2; (common_angle_orig*coeffs_faslenEND(2)) + coeffs_faslenEND(1)];
        pennation2 = [pennation2; (common_angle_orig*coeffs_pennationEND(2)) + coeffs_pennationEND(1)];
    end
        
    
    %% if gonio angle starts after zero degrees, extrapolate data
    
    % trial 1
    if min(gonio_new1) > 0.0
        % extrapolate using first X degrees of existing data
        angle_end = min(gonio_new1) + 5; %VAR
        loc_angle_end = find(gonio_new1 >= angle_end,1,'first');
        % calculate linear coeffisients
        gonio_modified = [ones(length(gonio_new1(1:loc_angle_end)),1) gonio_new1(1:loc_angle_end)];
        coeffs_faslen = gonio_modified\faslen1(1:loc_angle_end); % NB backslash --> slope
        coeffs_pennation = gonio_modified\pennation1(1:loc_angle_end); % NB backslash --> slope
        % enlarge arrays by adding values at zero angle to the front
        angle_half = min(gonio_new1)/2;
        gonio_new1 = [0; angle_half; gonio_new1];
        faslen1 = [coeffs_faslen(1); (angle_half*coeffs_faslen(2)) + coeffs_faslen(1); faslen1];
        pennation1 = [coeffs_pennation(1); (angle_half*coeffs_pennation(2)) + coeffs_pennation(1); pennation1];
    end
    
    % repeat for trial 2
    if min(gonio_new2) > 0.0
        % extrapolate using first X degrees of existing data
        angle_end = min(gonio_new2) + 5; %VAR
        loc_angle_end = find(gonio_new2 >= angle_end,1,'first');
        % calculate linear coeffisients
        gonio_modified = [ones(length(gonio_new2(1:loc_angle_end)),1) gonio_new2(1:loc_angle_end)];
        coeffs_faslen = gonio_modified\faslen2(1:loc_angle_end); % NB backslash --> slope
        coeffs_pennation = gonio_modified\pennation2(1:loc_angle_end); % NB backslash --> slope
        % enlarge arrays by adding values at zero angle to the front
        angle_half = min(gonio_new2)/2;
        gonio_new2 = [0; angle_half; gonio_new2];
        faslen2 = [coeffs_faslen(1); (angle_half*coeffs_faslen(2)) + coeffs_faslen(1); faslen2];
        pennation2 = [coeffs_pennation(1); (angle_half*coeffs_pennation(2)) + coeffs_pennation(1); pennation2];
    end
    
        
    %% create array with common angles (.05 intervals) for resampling
    
    % select the smallest range of the gonio data as basis for angle array
    common_angle_start = max([min(gonio_new1) min(gonio_new2)]);
    common_angle_stop = min([max(gonio_new1) max(gonio_new2)]);
    % create array of angles
    average_angle_array = (ceil(common_angle_start/angle_step)*angle_step:angle_step:floor(common_angle_stop/angle_step)*angle_step)';

        
    %% reshape, resample, average
    
    % reshape FASLEN across common angle array
    % if necessary, will extrapolate to include zero degrees gonio angle
    common_faslen1_gonio = spline(gonio_new1, faslen1, average_angle_array);
    common_faslen2_gonio = spline(gonio_new2, faslen2, average_angle_array);
    
    % average FASLEN across common angle array
    average_faslen_gonio = (smooth(common_faslen1_gonio,smoother) + smooth(common_faslen2_gonio,smoother)) / 2;

    % reshape PENNATION across common angle array
    % if necessary, will extrapolate to include zero degrees gonio angle
    common_pennation1_gonio = spline(gonio_new1, pennation1, average_angle_array);
    common_pennation2_gonio = spline(gonio_new2, pennation2, average_angle_array);
    
    % average PENNATION across common angle array
    average_pennation_gonio = (smooth(common_pennation1_gonio,smoother) + smooth(common_pennation2_gonio,smoother)) / 2;
        
    
    %% plot original and averaged data
    if plot_check && plot_licht
        plottitle = horzcat('IND Lichtwark averaging, ', subject_id, ' - ', description);
        figure('Name',plottitle);
        max_faslen = max([max(faslen1) max(faslen2)]);
        min_faslen = min([min(faslen1) min(faslen2)]);
        max_pennation = max([max(pennation1) max(pennation2)]);
        min_pennation = min([min(pennation1) min(pennation2)]);
        % top panel = fas len
        AXa = subplot(2,1,1);
        plot(gonio_new1,faslen1)
        hold on
        plot(gonio_new2,faslen2)
        plot(average_angle_array,average_faslen_gonio,'LineWidth',2)
        set(get(AXa,'Ylabel'),'String','Fascicle length (mm)')
        set(AXa,'YLim',[min_faslen*.9 max_faslen*1.1])
        title(plottitle)
        % bottom panel = pennation angle
        AXc = subplot(2,1,2);
        plot(gonio_new1,pennation1)
        hold on
        plot(gonio_new2,pennation2)
        plot(average_angle_array,average_pennation_gonio,'LineWidth',2)
        set(get(AXc,'Ylabel'),'String','Pennation angle (deg)')
        set(AXc,'YLim',[min_pennation*.9 max_pennation*1.1])
        xlabel('Gonio angle (deg)');
        legend('Trial 1','Trial 2','Spline','Location','Northeast');
        saveas(gcf, horzcat('data_plots/',plottitle,'.jpg'))
    end
    
        
else % nargin == 7, e.g. only one trial submitted
    
    gonio1 = varargin{1};
    %angle1 = varargin{2};
    faslen1 = varargin{3};
    pennation1 = varargin{4};
    time1 = varargin{5};
    angle_common = varargin{6};
    description = varargin{7};
        
    
    %% curve fitting gonio data

    % method 1: fit 4th order polynomial to averaged gonio-angle curve
    fit_gonio1 = polyfit(time1, gonio1, 4);

    % create array across angles
    gonio_new1(1:length(time1),1) = zeros;
    for ang = 1:length(time1)
        gonio_new1(ang) = (fit_gonio1(1) * time1(ang)^4) + (fit_gonio1(2) * time1(ang)^3) + (fit_gonio1(3) * time1(ang)^2) + (fit_gonio1(4) * time1(ang)) + fit_gonio1(5);
    end
    
    
    %% if new gonio angle is smaller than original max, extrapolate all data
    
    common_angle_orig = angle_common + 0.0001; % round(max(gonio1),2);
    
    % trial 1
    if max(gonio_new1) < common_angle_orig
        if common_angle_orig - max(gonio_new1) > errormsg_extrapol
            cprintf('red', horzcat('WARNING: Extrapolating max ROM1 by ', num2str(common_angle_orig - max(gonio_new1)), ' degrees.\n'))
        end
        % extrapolate using last X degrees of existing data
        angle_end = max(gonio_new1) - 2; %VAR
        loc_angle_end = find(gonio_new1 >= angle_end,1,'first');
        % calculate linear coeffisients
        gonio_modifiedEND = [ones(length(gonio_new1(loc_angle_end:end)),1) gonio_new1(loc_angle_end:end)];
        coeffs_faslenEND = gonio_modifiedEND\faslen1(loc_angle_end:end); % NB backslash --> slope
        coeffs_pennationEND = gonio_modifiedEND\pennation1(loc_angle_end:end);
        % enlarge arrays by adding values at the end
        gonio_new1 = [gonio_new1; common_angle_orig];
        faslen1 = [faslen1; (common_angle_orig*coeffs_faslenEND(2)) + coeffs_faslenEND(1)];
        pennation1 = [pennation1; (common_angle_orig*coeffs_pennationEND(2)) + coeffs_pennationEND(1)];
    end
    
    
    %% if gonio angle starts after zero degrees, extrapolate data
    
    % trial 1
    if min(gonio_new1) > 0.0
        % select first 2 degrees of existing gonio data
        angle_end = min(gonio_new1) + 2; %VAR extrapolate using first 2 degrees of the gonio movement
        loc_angle_end = find(gonio_new1 >= angle_end,1,'first');
        % calculate linear coeffisients
        gonio_modified = [ones(length(gonio_new1(1:loc_angle_end)),1) gonio_new1(1:loc_angle_end)];
        coeffs_faslen = gonio_modified\faslen1(1:loc_angle_end); % NB backslash --> slope
        coeffs_pennation = gonio_modified\pennation1(1:loc_angle_end); % NB backslash --> slope
        % enlarge arrays by adding values at zero angle to the front
        angle_half = min(gonio_new1)/2;
        gonio_new1 = [0; angle_half; gonio_new1];
        faslen1 = [coeffs_faslen(1); (angle_half*coeffs_faslen(2)) + coeffs_faslen(1); faslen1];
        pennation1 = [coeffs_pennation(1); (angle_half*coeffs_pennation(2)) + coeffs_pennation(1); pennation1];
    end

    
    %% create array with common angles (.05 intervals) for resampling
    
    % select the smallest range of the gonio data as basis for angle array
    common_angle_start = min(gonio_new1);
    common_angle_stop = max(gonio_new1);
    % create array of angles
    average_angle_array = (ceil(common_angle_start/angle_step)*angle_step:angle_step:floor(common_angle_stop/angle_step)*angle_step)';
    
    
    %% reshape, resample, average
    
    % reshape FASLEN across common angle array
    average_faslen_gonio = smooth(spline(gonio_new1, faslen1, average_angle_array),smoother);

    % reshape PENNATION across common angle array
    average_pennation_gonio = smooth(spline(gonio_new1, pennation1, average_angle_array),smoother);
    
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
        plot(gonio_new1,faslen1)
        hold on
        plot(average_angle_array,average_faslen_gonio,'LineWidth',2)
        set(get(AXa,'Ylabel'),'String','Fascicle length (mm)')
        set(AXa,'YLim',[min_faslen*.9 max_faslen*1.1])
        title(plottitle)
        AXc = subplot(2,1,2);
        plot(gonio_new1,pennation1)
        hold on
        plot(average_angle_array,average_pennation_gonio,'LineWidth',2)
        set(get(AXc,'Ylabel'),'String','Pennation angle (deg)')
        set(AXc,'YLim',[min_pennation*.9 max_pennation*1.1])
        xlabel('Gonio angle (deg)');
        legend('Trial','Spline','Location','Northeast');
    end
    
end

output_array = [average_angle_array average_faslen_gonio average_pennation_gonio];

end