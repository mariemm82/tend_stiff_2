%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract_force_displ_singletrial_passive
% Marie Moltubakk 11.6.2013, EMG added 1.6.2015
% Read complete, prepared noraon array + prepared us array
% Produce array with time, corrected force, corrected displacement

% used by passive analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [force,gonio,angle,displ_final,emg_GM,emg_GL,emg_SOL,time_us] = extract_force_displ_singletrial_passive_EMG(noraxondata, usdata, usdata_frame, max_EMG_TA, max_EMG_GM, max_EMG_GL, max_EMG_SOL, leg_length, side, line, trial_name)
    
    global column_l_gm column_r_gm column_l_gl column_r_gl column_l_sol column_r_sol column_gonio column_norm_angle column_norm_torque % column_l_tibant column_r_tibant  column_norm_velocity column_norm_direction column_achilles column_EMG_start column_EMG_end 
    global plot_check subject_id plot_norm % plot_emg plot_us 
    global at_momentarm
    global filepath
    
    
    
    % Read US data file, determine time stamps, set trigger frame as time = zero
    % Produce US sample frequency + new US array containing time and displacement
    % aa = strcat(filepath, usdata, '.txt') % 
    [usdata_prepped,usfreq] = read_us_file(strcat(filepath, usdata, '.txt'), str2double(usdata_frame), trial_name);
    
    
    
    % Read Noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
    % Produce a new noraxon data array
    noraxon_prepped = read_noraxon_passive(strcat(filepath, noraxondata), usfreq, side, trial_name);
    
    
   
    
    
    
    %%%%%%%%%%%%%%%%%%%%% extract data from big arrays
    
    %%% locate endpoint of ascending phase
    % Norm torque = zero during stop at maximal ROM - use this to identify ascending phase (0 to max ROM).
    loc_ascending_end = find(noraxon_prepped(:,column_norm_torque)==0,1,'first') - 1;
    % if US data does not have the full length of the Norm data:
    if loc_ascending_end > length(usdata_prepped)
        gonio_usmax = -noraxon_prepped(length(usdata_prepped),column_gonio);
        gonio_normmax = -noraxon_prepped(loc_ascending_end,column_gonio);
        if gonio_normmax-gonio_usmax > 2 %VAR
            cprintf('*err',horzcat('WARNING - ', trial_name, ': ', num2str(loc_ascending_end), ' Norm frames, ', num2str(length(usdata_prepped)), ' US frames. Not recorded/tracked until max ROM? Last ', num2str(round(gonio_normmax-gonio_usmax,3)), '° cut off. DISCARD?\n'));
        else
            cprintf('err',horzcat('WARNING - ', trial_name, ': ', num2str(loc_ascending_end), ' Norm frames, ', num2str(length(usdata_prepped)), ' US frames. Not recorded/tracked until max ROM? Last ', num2str(round(gonio_normmax-gonio_usmax,3)), '° cut off.\n'));
        end
        loc_ascending_end = length(usdata_prepped);
    end
    
    
    
    %%% time and angles
    
    time = noraxon_prepped(1:loc_ascending_end,1);
    time_us = usdata_prepped(1:loc_ascending_end,1);
    
    % turn dorsiflexion angles from negative to positive for plotting purposes
    angle = -noraxon_prepped(1:loc_ascending_end,column_norm_angle);
    gonio = -noraxon_prepped(1:loc_ascending_end,column_gonio);
    
    % permanent extrapolation of gonio and other data, in cases where gonio
    % starts after zero degrees, happens in "average_passive_trials" (after
    % all data are synchronized)
    
    
    
    
    %%% displacement
    % when gonio angle is zero, set displacement to zero
    
    % offset displacement: displacement = zero when gonio angle = zero
    displ_raw = -usdata_prepped(1:loc_ascending_end,2);
    
    % find last gonio angle BEFORE crossing zero angle
    loc_gonio_zero = find(gonio>=0,1,'first')-1; % if first value fulfills requirement, zero angle does not exist. Subtracting 1 --> loc_ = 0
    
    if loc_gonio_zero > 0 % zero angle exists in data
        % interpolation to find displacement offset at exactly zero:
        ang1 = gonio(loc_gonio_zero);
        ang2 = gonio(loc_gonio_zero+1);
        ang_ratio = -ang1/(ang2-ang1);
        displ1 = usdata_prepped(loc_gonio_zero,2);
        displ2 = usdata_prepped(loc_gonio_zero+1,2);
        displ_offset = displ1 + ang_ratio*(displ2-displ1);
        
    else % gonio does not start from zero - extrapolation is needed
        % extrapolation to find displacement offset if gonio angle were zero

        % select first X degrees of the existing gonio data
        range = 4; %VAR
        ang_interpol_end = gonio(1) + range;
        loc_ang_interpol_end = find(gonio >= ang_interpol_end,1,'first');
        
        % calculate linear coeffisients for extrapolation
        gonio_modified = [ones(length(gonio(1:loc_ang_interpol_end)),1) gonio(1:loc_ang_interpol_end)];
        coeffs_displ = gonio_modified\displ_raw(1:loc_ang_interpol_end); % NB backslash --> slope
        
        % enlarge arrays by adding values at zero angle to the front
        % y = ax + b ... y = displ, x = angle
        x = 0;
        % displacement if angle were zero (= value to be offset):
        displ_offset = - ((x*coeffs_displ(2)) + coeffs_displ(1));
    end
    displ_final = displ_raw + displ_offset;
    
    
    
    %%% EMG: calculate % EMG activation
    if strcmpi(side,'R') == 1
         emg_GM = noraxon_prepped(1:loc_ascending_end,column_r_gm)/max_EMG_GM*100; % percent
         emg_GL = noraxon_prepped(1:loc_ascending_end,column_r_gl)/max_EMG_GL*100; % percent
         emg_SOL = noraxon_prepped(1:loc_ascending_end,column_r_sol)/max_EMG_SOL*100; % percent
    else % left
         emg_GM = noraxon_prepped(1:loc_ascending_end,column_l_gm)/max_EMG_GM*100; % percent
         emg_GL = noraxon_prepped(1:loc_ascending_end,column_l_gl)/max_EMG_GL*100; % percent
         emg_SOL = noraxon_prepped(1:loc_ascending_end,column_l_sol)/max_EMG_SOL*100; % percent
    end
    cprintf('green', horzcat(trial_name, ' EMG max activation: gm ', num2str(round(max(emg_GM),3)), ' %%, gl ', num2str(round(max(emg_GL),3)), ' %%, sol ', num2str(round(max(emg_SOL),3)), ' %%.\n'));
    
    
    
    %%% force / torque 
    torque = noraxon_prepped(1:loc_ascending_end,column_norm_torque);
    force = torque / at_momentarm;
    
    
    
    %%% Plot synchronization check US vs Norm
    if plot_check && plot_norm
        plottitle = horzcat('SYNC check for ', subject_id, ' ', trial_name);
        figure('Name',plottitle);
        [AX,H1,H2] = plotyy(time,force/20,time_us,displ_final,'plot');
        hold on
        plot(time,gonio,'k')
        plot(time,angle,'g')
        plot(time,emg_GM,'y')
        plot(time,emg_GL,'m')
        plot(time,emg_SOL,'c')
        set(get(AX(2),'Ylabel'),'String','Displacement (mm)')
        set(get(AX(1),'Ylabel'),'String','Force (N)/10 / angle (deg) / EMG')
        set(AX,{'ycolor'},{'b';'r'})
        set(H2,'color','red')
        set(H1,'color','blue')
        set(AX(2),'YLim',[min(displ_final) max(displ_final)])
        set(AX(2),'YTick',(0:1:1.2*max(displ_final)))
        set(AX(2),'box','off')
        set(AX(1),'YLim',[0 max(force)/20])
        set(AX(1),'YTick',(0:round(max(force)/20/10,-1):1.2*max(force)/20))
        set(AX(1),'box','off')
        xlabel('Time (s)'),title(plottitle);
        legend('Force (N) /20','Gonio (deg)','Norm angle (deg)', '%EMG GM', '%EMG GL', '%EMG SOL','Displ (mm)','Location','Northwest');
    end


end