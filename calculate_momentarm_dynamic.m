%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_momentarm_dynamic
% Marie Moltubakk 15.01.2018
% Produce AT moment arm array respecting goniometer angle (Fath paper)
%%%%%%%%%%%%%%%%%%%%%%%%%%
    


function at_momentarm_dynamic = calculate_momentarm_dynamic(angle_array)
    global at_momentarm
    
    % Fath, COR rest (mm --> m)
    fath90 = 51.7 / 1000;
    fath105 = 55.4 / 1000;
    fath120 = 56.7 / 1000;
    
    % Fath, TE rest
    % 30.8
    % 34.3
    % 37.9
    % 41.1

    % Fath, COR MVC
    % 46.5
    % 53.8
    % 58.5
    % 61.8

    
    % percent increase in moment arm, per degree:
    fath_percent_1 = ((fath105-fath90)/15)/fath90;
    fath_percent_2 = ((fath120-fath105)/15)/fath105;
    
    % angle where to change from first to second % increase:
    fath_threshold = 15;
    
    % do we need to use both percent increase variables?
    loc_threshold = find(angle_array > fath_threshold,1,'first');
    if isempty(loc_threshold)
        %only first 
        at_momentarm_dynamic = at_momentarm * (1 + angle_array * fath_percent_1);
    else
        at_momentarm_dynamic(1:loc_threshold-1,1) = at_momentarm * (1 + angle_array(1:loc_threshold-1) * fath_percent_1);
        at_momentarm_dynamic(loc_threshold:length(angle_array),1) = at_momentarm_dynamic(loc_threshold-1) * (1 + (angle_array(loc_threshold:end)-fath_threshold) * fath_percent_2);
    end

    % Output moment arm report
    cprintf('magenta', horzcat('AT moment arm dynamic = ', num2str(1000*min(at_momentarm_dynamic)), ' to ', num2str(1000*max(at_momentarm_dynamic)), ' mm.\n'));

end
