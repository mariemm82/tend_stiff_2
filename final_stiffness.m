%%%%%%%%%%%%%%%%%%%%%%%%%%
% final_stiffness
% Marie Moltubakk 16.6.2013
% Read ...
% Produce ...
%%%%%%%%%%%%%%%%%%%%%%%%%%



function [fitresult,gof,stiff_frames,stiff_usedforce,commonforce] = final_stiffness(time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3, time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3, forceintervals, cutforce, rampmax) %new2014-04-14
    global plot_achilles plot_norm plot_emg plot_check subject_id 
    
    % Check for discarded trials (NULL / empty arrays) and set max so high it will not be chosen.
    if(isempty(time_force_displ_mtj1))
        max_time_force_displ_mtj1 = 100000;
    else
        max_time_force_displ_mtj1 = max(time_force_displ_mtj1(:,2));
    end
    
    if(isempty(time_force_displ_mtj2))
        max_time_force_displ_mtj2 = 100000;
    else
        max_time_force_displ_mtj2 = max(time_force_displ_mtj2(:,2));
    end
    
    if(isempty(time_force_displ_mtj3))
        max_time_force_displ_mtj3 = 100000;
    else
        max_time_force_displ_mtj3 = max(time_force_displ_mtj3(:,2));
    end
    
    if(isempty(time_force_displ_otj1))
        max_time_force_displ_otj1 = 100000;
    else
        max_time_force_displ_otj1 = max(time_force_displ_otj1(:,2));
    end
    
    if(isempty(time_force_displ_otj2))
        max_time_force_displ_otj2 = 100000;
    else
        max_time_force_displ_otj2 = max(time_force_displ_otj2(:,2));
    end
    
    if(isempty(time_force_displ_otj3))
        max_time_force_displ_otj3 = 100000;
    else
        max_time_force_displ_otj3 = max(time_force_displ_otj3(:,2));
    end
    
    % print maximal force from each trial to the screen
    report = sprintf(horzcat('Ramp trials uncut force:\n Common max = ', num2str(round(min(rampmax))), ' N, trials = ', num2str(round(rampmax(1))), ' - ', num2str(round(rampmax(2))), ' - ', num2str(round(rampmax(3))), ' - ', num2str(round(rampmax(4))), ' - ', num2str(round(rampmax(5))), ' - ', num2str(round(rampmax(6))), ' N.'));
    disp(report)
%    report = sprintf(horzcat('Ramp trials cut force:\n Max = ', num2str(round(max_time_force_displ_mtj1)), ' - ', num2str(round(max_time_force_displ_mtj2)), ' - ', num2str(round(max_time_force_displ_mtj3)), ' - ', num2str(round(max_time_force_displ_otj1)), ' - ', num2str(round(max_time_force_displ_otj2)), ' - ', num2str(round(max_time_force_displ_otj3)), ' N.'));
%    disp(report)
    
    % find lowest common force between all 6 trials
    commonforce = min(rampmax);

%%% moved elsewhere    
%     % use excel-determined common force instead of lowest force from the trials, print warning if excel-determined force is too high
%     if(commonforce < str2double(cutforce)) %new2014-04-12
%         % datamaster force level is set too high, print warning
%         report = sprintf(horzcat(' Common max force = ', num2str(round(commonforce)), ' N. Ignoring preset threshold of ', cutforce, ' N.'));
%         useforce = commonforce;
%     else
%         report = sprintf(horzcat(' Common max force = ', num2str(round(commonforce)), ' N. Cutoff at preset threshold = ', cutforce, ' N.'));
%         useforce = str2double(cutforce);
%     end
%     disp(report)
    
    % determine array of force values to use for averaging
    nmz_force = (0:forceintervals:commonforce)';
    
    % cut arrays at common force level - 3x MTJ, 3x OTJ
    mtjindex = 1;
    otjindex = 1;
    if(isempty(time_force_displ_mtj1))
        % Check for empty arrays: do nothing if empty (=discarded trial)
        time_force_displ_mtj1(1,1:3) = zeros();
        plot_nmz_disp_mtj(:,1) = zeros();
    else
        % detect position of max common force in array
        del_mtj1 = find(time_force_displ_mtj1(:,2)>max(nmz_force), 1, 'first');
        
        % delete array contents after commonforce
        % keeping 10 data points after commonforce, in order to do proper interpolation around the size of commonforce 
        keep_datapoints = 10; %VAR
        time_force_displ_mtj1(del_mtj1+keep_datapoints:length(time_force_displ_mtj1),:) = [];
        
        % use interpolation to extract displacement values at common force levels
        nmz_disp_mtj(:,mtjindex) = interp1(time_force_displ_mtj1(:,2),time_force_displ_mtj1(:,3),nmz_force,'linear'); % TMP 2014 MMM change to spline?
        plot_nmz_disp_mtj(:,1) = nmz_disp_mtj(:,mtjindex);
        mtjindex = mtjindex+1;
    end
    
    if(isempty(time_force_displ_mtj2))
        time_force_displ_mtj2(1,1:3) = zeros();
        plot_nmz_disp_mtj(:,2) = zeros();
    else
        del_mtj2 = find(time_force_displ_mtj2(:,2)>max(nmz_force), 1, 'first');
        time_force_displ_mtj2(del_mtj2+keep_datapoints:length(time_force_displ_mtj2),:) = [];
        nmz_disp_mtj(:,mtjindex) = interp1(time_force_displ_mtj2(:,2),time_force_displ_mtj2(:,3),nmz_force,'linear');
        plot_nmz_disp_mtj(:,2) = nmz_disp_mtj(:,mtjindex);
        mtjindex = mtjindex+1;
    end
    
    if(isempty(time_force_displ_mtj3))
        time_force_displ_mtj3(1,1:3) = zeros();
        plot_nmz_disp_mtj(:,3) = zeros();
    else
        del_mtj3 = find(time_force_displ_mtj3(:,2)>max(nmz_force), 1, 'first');
        time_force_displ_mtj3(del_mtj3+keep_datapoints:length(time_force_displ_mtj3),:) = [];
        nmz_disp_mtj(:,mtjindex) = interp1(time_force_displ_mtj3(:,2),time_force_displ_mtj3(:,3),nmz_force,'linear');
        plot_nmz_disp_mtj(:,3) = nmz_disp_mtj(:,mtjindex);
    end
    
    if(isempty(time_force_displ_otj1))
        time_force_displ_otj1(1,1:3) = zeros();
        plot_nmz_disp_otj(:,1) = zeros();
    else
        del_otj1 = find(time_force_displ_otj1(:,2)>max(nmz_force), 1, 'first');
        time_force_displ_otj1(del_otj1+keep_datapoints:length(time_force_displ_otj1),:) = [];
        nmz_disp_otj(:,otjindex) = interp1(time_force_displ_otj1(:,2),time_force_displ_otj1(:,3),nmz_force,'linear');
        plot_nmz_disp_otj(:,1) = nmz_disp_otj(:,otjindex);
        otjindex = otjindex+1;
    end
    
    if(isempty(time_force_displ_otj2))
        time_force_displ_otj2(1,1:3) = zeros();
        plot_nmz_disp_otj(:,2) = zeros();
    else
        del_otj2 = find(time_force_displ_otj2(:,2)>max(nmz_force), 1, 'first');
        time_force_displ_otj2(del_otj2+keep_datapoints:length(time_force_displ_otj2),:) = [];
        nmz_disp_otj(:,otjindex) = interp1(time_force_displ_otj2(:,2),time_force_displ_otj2(:,3),nmz_force,'linear');
        plot_nmz_disp_otj(:,2) = nmz_disp_otj(:,otjindex);
        otjindex = otjindex+1;
    end
    
    if(isempty(time_force_displ_otj3))
        time_force_displ_otj3(1,1:3) = zeros();
        plot_nmz_disp_otj(:,3) = zeros();
    else
        del_otj3 = find(time_force_displ_otj3(:,2)>max(nmz_force), 1, 'first');
        time_force_displ_otj3(del_otj3+keep_datapoints:length(time_force_displ_otj3),:) = [];
        nmz_disp_otj(:,otjindex) = interp1(time_force_displ_otj3(:,2),time_force_displ_otj3(:,3),nmz_force,'linear');
        plot_nmz_disp_otj(:,3) = nmz_disp_otj(:,otjindex);
    end
    
    % calculate and plot averaged displacement
    nmz_disp_mtj_mean = mean(nmz_disp_mtj');
    nmz_disp_otj_mean = mean(nmz_disp_otj');
	nmz_disp = nmz_disp_mtj_mean - nmz_disp_otj_mean;
    
    % rotate arrays
    nmz_disp_mtj_mean = nmz_disp_mtj_mean';
    nmz_disp_otj_mean = nmz_disp_otj_mean';
    nmz_disp = nmz_disp';
    
    
    
    if plot_achilles
        % MTJ force-disp x 3
        plottitle = horzcat('MTJ stiffness check for ', subject_id);
        figure('Name',plottitle)
        plot(time_force_displ_mtj1(:,3),time_force_displ_mtj1(:,2),'b')
        hold on
        plot(time_force_displ_mtj2(:,3),time_force_displ_mtj2(:,2),'r')
        plot(time_force_displ_mtj3(:,3),time_force_displ_mtj3(:,2),'g')
        plot(mean(nmz_disp_mtj'), nmz_force, 'k','LineWidth',2)
        plot(plot_nmz_disp_mtj(:,1), nmz_force,'b.')
        plot(plot_nmz_disp_mtj(:,2), nmz_force,'r.')
        plot(plot_nmz_disp_mtj(:,3), nmz_force,'g.')
        xlabel('Displacement (mm)'),ylabel('Tendon force (N)'),title(plottitle);
        legend('trial 1','trial 2','trial 3','mean','Location','Southeast');
    end

    if plot_achilles
        % OTJ force-disp x 3
        plottitle = horzcat('OTJ stiffness check for ', subject_id);
        figure('Name',plottitle)
        plot(time_force_displ_otj1(:,3),time_force_displ_otj1(:,2),'b')
        hold on
        plot(time_force_displ_otj2(:,3),time_force_displ_otj2(:,2),'r')
        plot(time_force_displ_otj3(:,3),time_force_displ_otj3(:,2),'g')
        plot(mean(nmz_disp_otj'), nmz_force, 'k','LineWidth',2)
        plot(plot_nmz_disp_otj(:,1), nmz_force,'b.')
        plot(plot_nmz_disp_otj(:,2), nmz_force,'r.')
        plot(plot_nmz_disp_otj(:,3), nmz_force,'g.')
        xlabel('Displacement (mm)'),ylabel('Tendon force (N)'),title(plottitle);
        legend('trial 1','trial 2','trial 3','mean','Location','Southeast');
    end
    
    
    
    %%% Force cutting method 3, june 2014
    % prepare for cutoff at 90% of max elongation
    percent_elong = 0.90; %VAR
    [max_elong,max_elong_pos] = max(nmz_disp);
    keep_elong = percent_elong * max_elong;
    elong_cut_index = find(nmz_disp>keep_elong,1,'first');
% commenting out: will include the first force level that is HIGHER than 90% of displacement
%    elong_cut_index = elong_cut_index-1;

    % write cutoff to screen
    report = sprintf(horzcat('Max elong ', num2str(max_elong), ' mm @ ', num2str(round(nmz_force(max_elong_pos))), ' N, ', num2str(percent_elong*100), '%% elong -> ', num2str(round(nmz_force(elong_cut_index))), ' N.'));
    disp(report)
    
    % cut off even earlier if datamaster is set to cut at lower value
    force_cut_index = find(nmz_force>str2double(cutforce),1,'first');
    force_cut_index = force_cut_index-1;
    if force_cut_index < elong_cut_index
        cut_index = force_cut_index;
        report = sprintf(horzcat('Datamaster cut force is set lower = ', num2str(nmz_force(cut_index)), ' N.'));
        disp(report)
    else
        cut_index = elong_cut_index;
    end
    
    % saving the applied max force level, for output to file
    stiff_usedforce = nmz_force(cut_index);

    % curve fitting for stiffness
    [fitresult, gof] = fit_stiffness(nmz_disp, nmz_force, cut_index, nmz_disp_mtj_mean, nmz_disp_otj_mean);
    
    % extract string of stiffness equation coefficients
    coeffstring = '';
    coeffvector = coeffvalues(fitresult);
    for i = 1:numcoeffs(fitresult)
        coeffstring = horzcat(coeffstring, ' ', num2str(coeffvector(i)));
    end
    
    % write stiffness coefficients to screen
    report = sprintf(horzcat('Stiffness report:\n Stiffness coeff =', coeffstring, '. R2 = ', num2str(gof.rsquare)));
    disp(report)    

    % create stiffness data array for later use
    stiff_frames(:,1) = nmz_disp(1:elong_cut_index);
    stiff_frames(:,2) = nmz_force(1:elong_cut_index);

end