%%%%%%%%%%%%%%%%%%%%%%%%%%
% final_stiffness
% Marie Moltubakk 16.6.2013
% Read 3+3 trials of MTJ and OTJ scans
% Produce stiffness coeffisients and output arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%



function [fitresult, gof, force_elong_array, final_cutoff_force] = final_stiffness(time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3, time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3, forceintervals, force_cutoff_manual, force_max_trials) %new2014-04-14

global plot_achilles subject_id %plot_norm plot_emg plot_check

    %% calculate tendon elongation
    
    % print maximal force from each trial to the screen
    cprintf('black', horzcat('Ramp trials, ind max force = '))
    cprintf('*blue', horzcat(num2str(round(min(force_max_trials))), ' N'))
    cprintf('black', horzcat(', trials = ', num2str(round(force_max_trials(1))), ' - ', num2str(round(force_max_trials(2))), ' - ', num2str(round(force_max_trials(3))), ' - ', num2str(round(force_max_trials(4))), ' - ', num2str(round(force_max_trials(5))), ' - ', num2str(round(force_max_trials(6))), ' N.\n'))

    % find lowest common force between all 6 trials
    commonforce = min(force_max_trials);

    % create array of force values to use for averaging
    force_array = (0:forceintervals:commonforce)';

    % cut arrays at common force level (force across 3x MTJ + 3x OTJ)
    mtjindex = 1;
    otjindex = 1;
    
    if(isempty(time_force_displ_mtj1))
        % Check for empty arrays: do nothing if empty (=discarded trial)
        time_force_displ_mtj1(1,1:3) = zeros();
        plot_displ_mtj(:,1) = zeros();
    else
        % detect position of max common force in array
        del_mtj1 = find(time_force_displ_mtj1(:,2)>force_array(end), 1, 'first');

        % delete array contents after commonforce
        % keeping 10 data points after commonforce, in order to do proper interpolation around the size of commonforce
        keep_datapoints = 10; %VAR
        time_force_displ_mtj1(del_mtj1+keep_datapoints:end,:) = [];

        % use interpolation to extract displacement values at common force levels
        displ_mtj(:,mtjindex) = interp1(time_force_displ_mtj1(:,2),time_force_displ_mtj1(:,3),force_array,'linear'); % TMP 2014 MMM change to spline?
        plot_displ_mtj(:,1) = displ_mtj(:,mtjindex);
        mtjindex = mtjindex+1;
    end

    if(isempty(time_force_displ_mtj2))
        time_force_displ_mtj2(1,1:3) = zeros();
        plot_displ_mtj(:,2) = zeros();
    else
        del_mtj2 = find(time_force_displ_mtj2(:,2)>force_array(end), 1, 'first');
        time_force_displ_mtj2(del_mtj2+keep_datapoints:end,:) = [];
        displ_mtj(:,mtjindex) = interp1(time_force_displ_mtj2(:,2),time_force_displ_mtj2(:,3),force_array,'linear');
        plot_displ_mtj(:,2) = displ_mtj(:,mtjindex);
        mtjindex = mtjindex+1;
    end

    if(isempty(time_force_displ_mtj3))
        time_force_displ_mtj3(1,1:3) = zeros();
        plot_displ_mtj(:,3) = zeros();
    else
        del_mtj3 = find(time_force_displ_mtj3(:,2)>force_array(end), 1, 'first');
        time_force_displ_mtj3(del_mtj3+keep_datapoints:end,:) = [];
        displ_mtj(:,mtjindex) = interp1(time_force_displ_mtj3(:,2),time_force_displ_mtj3(:,3),force_array,'linear');
        plot_displ_mtj(:,3) = displ_mtj(:,mtjindex);
    end

    if(isempty(time_force_displ_otj1))
        time_force_displ_otj1(1,1:3) = zeros();
        plot_displ_otj(:,1) = zeros();
    else
        del_otj1 = find(time_force_displ_otj1(:,2)>force_array(end), 1, 'first');
        time_force_displ_otj1(del_otj1+keep_datapoints:end,:) = [];
        displ_otj(:,otjindex) = interp1(time_force_displ_otj1(:,2),time_force_displ_otj1(:,3),force_array,'linear');
        plot_displ_otj(:,1) = displ_otj(:,otjindex);
        otjindex = otjindex+1;
    end

    if(isempty(time_force_displ_otj2))
        time_force_displ_otj2(1,1:3) = zeros();
        plot_displ_otj(:,2) = zeros();
    else
        del_otj2 = find(time_force_displ_otj2(:,2)>force_array(end), 1, 'first');
        time_force_displ_otj2(del_otj2+keep_datapoints:end,:) = [];
        displ_otj(:,otjindex) = interp1(time_force_displ_otj2(:,2),time_force_displ_otj2(:,3),force_array,'linear');
        plot_displ_otj(:,2) = displ_otj(:,otjindex);
        otjindex = otjindex+1;
    end

    if(isempty(time_force_displ_otj3))
        time_force_displ_otj3(1,1:3) = zeros();
        plot_displ_otj(:,3) = zeros();
    else
        del_otj3 = find(time_force_displ_otj3(:,2)>force_array(end), 1, 'first');
        time_force_displ_otj3(del_otj3+keep_datapoints:end,:) = [];
        displ_otj(:,otjindex) = interp1(time_force_displ_otj3(:,2),time_force_displ_otj3(:,3),force_array,'linear');
        plot_displ_otj(:,3) = displ_otj(:,otjindex);
    end

    % average displacement trials, calculate elongation
    displ_mtj_mean = mean(displ_mtj,2);
    displ_otj_mean = mean(displ_otj,2);
    tend_elong = displ_mtj_mean - displ_otj_mean;
    
    if plot_achilles
        % MTJ force-disp x 3
        plottitle = horzcat('MTJ stiffness check for ', subject_id);
        figure('Name',plottitle)
        plot(time_force_displ_mtj1(:,3),time_force_displ_mtj1(:,2),'b')
        hold on
        plot(time_force_displ_mtj2(:,3),time_force_displ_mtj2(:,2),'r')
        plot(time_force_displ_mtj3(:,3),time_force_displ_mtj3(:,2),'g')
        plot(mean(displ_mtj,2), force_array, 'k','LineWidth',2)
        plot(plot_displ_mtj(:,1), force_array,'b.')
        plot(plot_displ_mtj(:,2), force_array,'r.')
        plot(plot_displ_mtj(:,3), force_array,'g.')
        xlabel('Displacement (mm)'),ylabel('Tendon force (N)'),title(plottitle);
        legend('trial 1','trial 2','trial 3','mean','Location','Southeast');
        saveas(gcf, strcat('data_plots_stiff/IND_stiff_MTJ_avg_', subject_id, 'png'))
    end

    if plot_achilles
        % OTJ force-disp x 3
        plottitle = horzcat('OTJ stiffness check for ', subject_id);
        figure('Name',plottitle)
        plot(time_force_displ_otj1(:,3),time_force_displ_otj1(:,2),'b')
        hold on
        plot(time_force_displ_otj2(:,3),time_force_displ_otj2(:,2),'r')
        plot(time_force_displ_otj3(:,3),time_force_displ_otj3(:,2),'g')
        plot(mean(displ_otj,2), force_array, 'k','LineWidth',2)
        plot(plot_displ_otj(:,1), force_array,'b.')
        plot(plot_displ_otj(:,2), force_array,'r.')
        plot(plot_displ_otj(:,3), force_array,'g.')
        xlabel('Displacement (mm)'),ylabel('Tendon force (N)'),title(plottitle);
        legend('trial 1','trial 2','trial 3','mean','Location','Southeast');
        saveas(gcf, strcat('data_plots_stiff/IND_stiff_OTJ_avg_', subject_id, 'png'))
    end
    
    
    %% cut force (method 3, June 2014)
    
    % prepare for cutoff at 90% of max elongation
    percent_elong = 0.90; %VAR
    [max_elong,loc_max_elong] = max(tend_elong);
    keep_elong = percent_elong * max_elong;
    loc_elong_cut = find(tend_elong>keep_elong,1,'first') - 1; % finds the first value larger than 90% limit - we will use the largest BEFORE 90% limit

    % write cutoff to screen
    cprintf('black', horzcat('Max elong = ', num2str(round(max_elong,2)), ' mm @ ', num2str(force_array(loc_max_elong)), ' N, ', num2str(percent_elong*100), '%% elong -> cutoff at '))
    cprintf('*blue', horzcat(num2str(force_array(loc_elong_cut)), ' N.\n'))

    % cut off even earlier if datamaster is set to cut at lower value
    if force_array(loc_elong_cut) > str2double(force_cutoff_manual)
        cprintf('red',horzcat('Force cutoff is lowered manually (datamaster), to ', force_cutoff_manual, ' N.\n'));
        loc_elong_cut = find(force_array>str2double(force_cutoff_manual),1,'first') - 1;
    end

    % save the applied max force level, for output to file
    final_cutoff_force = force_array(loc_elong_cut);
    
    
    %% curve fitting for stiffness
    [fitresult, gof] = fit_stiffness(tend_elong, force_array, loc_elong_cut, displ_mtj_mean, displ_otj_mean);

    % write stiffness coefficients to screen
    cprintf('*blue', horzcat('Stiffness coeffs = ', regexprep(num2str(coeffvalues(fitresult)),' +',' '), '. R2 = ', num2str(gof.rsquare), '\n'))

    % create stiffness data array for later use
    force_elong_array = [tend_elong(1:loc_elong_cut) force_array(1:loc_elong_cut)];

end