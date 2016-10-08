%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_rotation_correction
% Marie Moltubakk 17.6.2013
% Read complete, prepared noraon array + prepared us array
% Produce ???? AT moment arm constant, by tendon excursion method
%%%%%%%%%%%%%%%%%%%%%%%%%%

function at_rotation_const = calculate_rotation_correction(noraxon_rot, usdata_rot)

   % angle to which gonio data should be collected (find NEXT frame where gonio goes above x degrees, in plantarflexion series)
   stopangle = 6; %VAR



    global plot_conversion plot_check subject_id
    global column_gonio column_norm_angle 
    
    % extract Norm machine angle series
    norm_angle_filtered = noraxon_rot(:,column_norm_angle);
    
    % Identify movement phases and stop phases, from NORM machine angles
    angle_goingdown = find(norm_angle_filtered<-2.0,1,'first'); %find more than 2.5 degrees into dorsi
    angle_goingup = find(norm_angle_filtered(angle_goingdown:end)>-1.5,1,'first'); %find closer to zero than 1,5 degrees dorsi
    angle_goingup = angle_goingdown + angle_goingup - 1;
    
    ten1_start = find(norm_angle_filtered>9.9,1,'first');

    ten1_stop = find(norm_angle_filtered(ten1_start:end)<=9.9,1,'first');
    ten1_stop = ten1_start + ten1_stop - 1;
    
    five_start = find(norm_angle_filtered(ten1_stop:end)<=-4.9,1,'first');
    five_start = five_start + ten1_stop - 1;
    
    five_stop = find(norm_angle_filtered(five_start:end)>=-4.9,1,'first');
    five_stop = five_stop + five_start - 1;
    
    ten2_start = find(norm_angle_filtered(five_stop:end)>=9.9,1,'first');
    ten2_start = ten2_start + five_stop - 1;

    ten2_stop = find(norm_angle_filtered(ten2_start:end)<=9.9,1,'first');
    ten2_stop = ten2_stop + ten2_start - 1;

    % extract data from FIRST movement phase plantarflexion
    displ0 = usdata_rot(angle_goingup:ten1_start,2);
    angle0 = noraxon_rot(angle_goingup:ten1_start,column_gonio);

    % extract data from movement phase DORSIFLEXION (used only for plot)
    displ1 = usdata_rot(ten1_stop:five_start,2);
    angle1 = noraxon_rot(ten1_stop:five_start,column_gonio);
    
    % extract data from movement phase PLANTARFLEXION
    displ2 = usdata_rot(five_stop:ten2_start,2);
    angle2 = noraxon_rot(five_stop:ten2_start,column_gonio);
    
   if plot_check && plot_conversion
        plottitle = horzcat('Ankle rotation correction check 2 for ', subject_id);
        figure('Name',plottitle);
        plot(angle0, displ0, 'g.'); % plantar direction first phase
        hold on
        plot(angle1, displ1, 'r.'); % dorsi direction
        plot(angle2, displ2, 'b.'); % plantar direction second (default) phase
        xlabel('Goniometer ankle angle (deg)'),ylabel('Calcaneus displacement (mm)'),title(plottitle);
        legend('Plantarflex startup', 'Dorsiflex', 'Plantarflex second phase', 'Location','Northeast');
   end
    
   
   % calculation based on second plantar flexion phase (default)
   % find first frame where gonio passes zero degrees, in plantarflexion series
   zeroangle = find(angle2>=0,1,'first');
   % we want to start from the frame right BEFORE zero, not after zero
   zeroangle = zeroangle-1;
   
   % find NEXT frame where gonio goes above 8 degrees, in plantarflexion series
   fiveangle = find(angle2(zeroangle:end)>stopangle,1,'first');
   if isempty(fiveangle) % if goniometer doesn't go to 5 degrees
       fiveangle = length(angle2); % then use maximal gonio angle
   else
       fiveangle = fiveangle + zeroangle-1;
   end

   % create linear equation for plantarflex direction, starting from zero angle
   [fitresult, gof] = fit_ankle_rotation(angle2(zeroangle:fiveangle), displ2(zeroangle:fiveangle),'2nd');
   coeffvals = coeffvalues(fitresult);
   

   first_plantar_usable = 1;
   if angle_goingup > ten1_start % 1st plantarflex phase does not start from dorsiflex angle (no neg angle)
       first_plantar_usable = 0;
   end
   
   if first_plantar_usable
       % calculation based on FIRST plantar flexion phase (crisis)
       % find first frame where gonio passes zero degrees, in plantarflexion series
       zeroangle = find(angle0>=0,1,'first');
       % we want to start from the frame right BEFORE zero, not after zero
       if zeroangle > 1 %don't go to frame 0 if the gonio does not start lower than zero
           zeroangle = zeroangle-1;
       end

       % find NEXT frame where gonio goes above 8 degrees, in plantarflexion series
       fiveangle = find(angle0(zeroangle:end)>stopangle,1,'first');
       if isempty(fiveangle) % if goniometer doesn't go to X degrees
           fiveangle = length(angle0); % then use maximal gonio angle
       else
           fiveangle = fiveangle + zeroangle-1;
       end

       % create linear equation for plantarflex direction, starting from zero angle
       [fitresult, gof] = fit_ankle_rotation(angle0(zeroangle:fiveangle), displ0(zeroangle:fiveangle),'1st');
       coeffvals0 = coeffvalues(fitresult);
   else % will not calculate any realistic coeffvals:
       coeffvals0(1) = 10000;
   end
   

   
   % select the best rotation constant
   at_rotation_const = min(coeffvals0(1),coeffvals(1));
   
   
    
    
   
   % plot norm angle, gonio angle, displacement, zone lines
   if plot_check && plot_conversion
       plottitle = horzcat('SYNC check ankle rotation for ', subject_id);
       fig_anklerot = figure('Name', plottitle);
       plot(noraxon_rot(1:ten2_stop,1),norm_angle_filtered(1:ten2_stop),'b'); % norm angle
       hold on
       plot(noraxon_rot(1:ten2_stop,1),noraxon_rot(1:ten2_stop,column_gonio),'m'); % goniometer 
       plot(usdata_rot(1:ten2_stop,1),(8*usdata_rot(1:ten2_stop,2)),'r'); % displacement MAGNIFIED
       plot([noraxon_rot(ten1_start,1) noraxon_rot(ten1_start,1)], [min(norm_angle_filtered(1:ten2_stop)) max(norm_angle_filtered(1:ten2_stop))],'g');
       plot([noraxon_rot(ten1_stop,1) noraxon_rot(ten1_stop,1)], [min(norm_angle_filtered(1:ten2_stop)) max(norm_angle_filtered(1:ten2_stop))],'g');
       plot([noraxon_rot(five_start,1) noraxon_rot(five_start,1)], [min(norm_angle_filtered(1:ten2_stop)) max(norm_angle_filtered(1:ten2_stop))],'g');
       plot([noraxon_rot(five_stop,1) noraxon_rot(five_stop,1)], [min(norm_angle_filtered(1:ten2_stop)) max(norm_angle_filtered(1:ten2_stop))],'g');
       plot([noraxon_rot(ten2_start,1) noraxon_rot(ten2_start,1)], [min(norm_angle_filtered(1:ten2_stop)) max(norm_angle_filtered(1:ten2_stop))],'g');
       xlabel('Time (s)'),ylabel('Data (misc)'),title(plottitle);
       legend('Norm angle', 'Goniometer', 'Calc displ MAGNIFIED', 'Phases','Location','Northeast');
       text(0.2, 10.2, horzcat('Displ/deg = ', num2str(at_rotation_const)), 'Color', 'k');
       saveas(fig_anklerot, strcat('data_output/ankle_rot_', subject_id), 'png')
   end
   
   % output angle conversion numbers to screen, as text
   report = sprintf(horzcat('Ankle rotation: Displacement per degree of rotation = ', num2str(at_rotation_const), '.'));
   disp(report)
    
end