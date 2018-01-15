%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_momentarm_dynamic
% Marie Moltubakk 15.01.2018
% Produce AT moment arm array respecting goniometer angle (Fath paper)
%%%%%%%%%%%%%%%%%%%%%%%%%%
    


function at_momentarm_dynamic = calculate_momentarm_dynamic(angle_array)
    global at_momentarm
    
    fath90 = 51.7;
    fath105 = 55.4;
    fath120 = 56.7;
    
    % percent increase in moment arm, per degree:
    fath_percent_1 = ((fath105-fath90)/15)/fath90;
    fath_percent_2 = ((fath120-fath105)/15)/fath105;
    
    % angle where to change from first to second % increase:
    fath_threshold = 15;
    
    % do we need to use both percent increase variables?
    loc_angle_zero = find(angle_array >= angle_zero,1,'first');
    angle_array;
    at_momentarm_dynamic;


    % Output moment arm report
    cprintf('magenta', horzcat('AT moment arm dynamic = ', num2str(min(at_momentarm_dynamic)), ' to ', num2str(max(at_momentarm_dynamic)), ' m.\n'));

end
