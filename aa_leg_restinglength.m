function resting_calf_length_abs = aa_leg_restinglength(resting_ankle_angle_pos, standing_calf_length_cm)
    a0 = -22.18468;
    a1 = 0.30141;
    a2 = -0.00061;

    resting_ankle_angle = -resting_ankle_angle_pos;
    standing_calf_length = standing_calf_length_cm * 10;
    
    % calf/MTU elongation - in % change from standing length
    resting_calf_length_rel = (a0 + (a1*(90+resting_ankle_angle)) + (a2*(90+resting_ankle_angle).^2))/100;
    % calf/MTU length - in mm, across angle array
    resting_calf_length_abs = standing_calf_length * (1 + resting_calf_length_rel);

end