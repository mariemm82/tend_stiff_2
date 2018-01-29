%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract_force_displ_singletrial_passive
% Marie Moltubakk 7.12.2017
% Read complete, prepared noraon array + prepared us array
% Produce gonio-elong data for ankle rotation fits
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [coeffvals, rsquare, gonio, angle, displ] = extract_rot_corr(noraxondata, usdata, usdata_frame, side, trial_name)
    angle_end = -3; %VAR fit data from start of trial up to gonio angle 3

    global mute
    global column_gonio column_norm_angle %column_l_gm column_r_gm column_l_gl column_r_gl column_l_sol column_r_sol  column_norm_torque % column_l_tibant column_r_tibant  column_norm_velocity column_norm_direction column_achilles column_EMG_start column_EMG_end 
    %global plot_check subject_id plot_norm % plot_emg plot_us 
    %global at_momentarm
    global filepath
    
    %% gather files
    
    % Read US data file, determine time stamps, set trigger frame as time = zero
    % Produce US sample frequency + new US array containing time and displacement
    % aa = strcat(filepath, usdata, '.txt') % 
    [usdata_prepped,usfreq] = read_us_file(strcat(filepath, usdata, '.txt'), str2double(usdata_frame), trial_name);
    
        
    % Read Noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
    % Produce a new noraxon data array
    noraxon_prepped = read_noraxon_passive(strcat(filepath, noraxondata), usfreq, side, trial_name);
    
    
    %% extract data
    gonio_raw = noraxon_prepped(:,column_gonio);
    
    % find gonio angle of 5 degrees dorsiflex 
    % (first point where angle gets more negative than 5 is INCLUDED)
    loc_gonio_end = find(gonio_raw <= angle_end,1,'first');
    if loc_gonio_end == 1
        coeffvals = [NaN NaN];
        rsquare = NaN;
        gonio = NaN;
        angle = NaN;
        displ = NaN;
    else
        gonio = gonio_raw(1:loc_gonio_end);

        % turn dorsiflexion Norm angles from negative to positive for plotting purposes
        angle = -noraxon_prepped(1:loc_gonio_end,column_norm_angle);

        displ = -usdata_prepped(1:loc_gonio_end,2);


        %% calculate correction factors
        % create linear equation
        [fitresult, gof] = fit_ankle_rotation(gonio, displ, horzcat('Stretch ', trial_name));
        coeffvals = coeffvalues(fitresult);
        at_rotation_const = coeffvals(1);
        rsquare = gof.rsquare;

        if at_rotation_const > 0 %VAR
            cprintf('red',horzcat('Ankle rotation ', trial_name, ': ', num2str(at_rotation_const), ' mm/deg (start gonio angle ', num2str(-gonio(1)), ', norm angle ', num2str(-angle(1)),')', '.\n'));
        elseif at_rotation_const < -0.25 %VAR
            cprintf('red',horzcat('Ankle rotation ', trial_name, ': ', num2str(at_rotation_const), ' mm/deg (start gonio angle ', num2str(-gonio(1)), ', norm angle ', num2str(-angle(1)),')', '.\n'));
        else
            if mute == 0
                cprintf('blue',horzcat('Ankle rotation ', trial_name, ': ', num2str(at_rotation_const), ' mm/deg (start gonio angle ', num2str(-gonio(1)), ', norm angle ', num2str(-angle(1)),')', '.\n'));
            end
        end
    end
    
end