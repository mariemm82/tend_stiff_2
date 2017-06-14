%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_rotation_correction
% Marie Moltubakk 17.6.2013
% Read complete, prepared noraon array + prepared us array
% Produce ???? AT moment arm constant, by tendon excursion method
%%%%%%%%%%%%%%%%%%%%%%%%%%

function at_rotation_const = calculate_rotation_correction(noraxon_rot, usdata_rot)

   % angle to which gonio data should be collected (find NEXT frame where gonio goes above x degrees, in plantarflexion series)
   stopangle = 6; %VAR



    global plot_norm plot_check plot_achilles subject_id
    global column_gonio column_norm_angle 
    
    
    
    
    % secondary zero offset for GONIOMETER data (primary offset during data collection, only before 1st of 3 trials)
    %     this is normally done right after read_noraxon_stiffness, in
    %     extract_force_displ_singletrial - but for noraxon_CPM trials,
    %     extr_f_displ is not used --> running here instead 
    gonio_offset = mean(noraxon_rot(1:3,column_gonio)) - mean(noraxon_rot(1:3,column_norm_angle));
    noraxon_rot(:,column_gonio) = noraxon_rot(:,column_gonio) - gonio_offset;
    
    
    
    
    % extract Norm machine angle series
    norm_angle_filtered = noraxon_rot(:,column_norm_angle);
    
    % Identify movement phases and stop phases, from NORM machine angles
    angle_goingdown = find(norm_angle_filtered<-2.0,1,'first'); %find moment of large dorsi - more than 2 degrees
    angle_goingup = find(norm_angle_filtered(angle_goingdown:end)>-1.5,1,'first'); %find subsequent movement towards plantar
    angle_goingup = angle_goingdown + angle_goingup - 1;
    
    % below limits are only used to identify peaks: 
    %   fit between angle and displacement is done from 0 to 6 degrees of
    %   plantar flexion (6 = stopangle)
    peak_df = 9; %VAR
    peak_pf = -4; %VAR
    
    ten_plateau_1_start = find(norm_angle_filtered(angle_goingup:end)>peak_df,1,'first');
    ten_plateau_1_start = angle_goingup + ten_plateau_1_start - 1;

    ten_plateau_1_stop = find(norm_angle_filtered(ten_plateau_1_start:end)<=peak_df,1,'first');
    ten_plateau_1_stop = ten_plateau_1_start + ten_plateau_1_stop - 1;
    
    five_plateau_start = find(norm_angle_filtered(ten_plateau_1_stop:end)<=peak_pf,1,'first');
    five_plateau_start = five_plateau_start + ten_plateau_1_stop - 1;
    
    five_plateau_stop = find(norm_angle_filtered(five_plateau_start:end)>=peak_pf,1,'first');
    five_plateau_stop = five_plateau_stop + five_plateau_start - 1;
    
    ten_plateau_2_start = find(norm_angle_filtered(five_plateau_stop:end)>=peak_df,1,'first');
    ten_plateau_2_start = ten_plateau_2_start + five_plateau_stop - 1;

    ten_plateau_2_stop = find(norm_angle_filtered(ten_plateau_2_start:end)<=peak_df,1,'first');
    ten_plateau_2_stop = ten_plateau_2_stop + ten_plateau_2_start - 1;

    % extract data from FIRST movement phase plantarflexion
    displ0 = usdata_rot(angle_goingup:ten_plateau_1_start,2);
    angle0 = noraxon_rot(angle_goingup:ten_plateau_1_start,column_gonio);

    % extract data from movement phase DORSIFLEXION (used only for plot)
    displ1 = usdata_rot(ten_plateau_1_stop:five_plateau_start,2);
    angle1 = noraxon_rot(ten_plateau_1_stop:five_plateau_start,column_gonio);
    
    % extract data from movement phase PLANTARFLEXION
    displ2 = usdata_rot(five_plateau_stop:ten_plateau_2_start,2);
    angle2 = noraxon_rot(five_plateau_stop:ten_plateau_2_start,column_gonio);
    
   if plot_check && plot_norm
        plottitle = horzcat('Ankle rotation correction check 2 for ', subject_id);
        figure('Name',plottitle);
        plot(angle0, displ0, 'g.'); % plantar direction first phase
        hold on
        plot(angle1, displ1, 'r.'); % dorsi direction
        plot(angle2, displ2, 'b.'); % plantar direction second phase (default) 
        xlabel('<--dorsiflex --- Goniometer ankle angle (deg) --- plantarflex-->'),ylabel('Calcaneus displacement (mm)'),title(plottitle);
        legend('Plantarflex onset (-2 to 10)', 'Dorsiflex (10 to -5)', 'Plantarflex (-5 to 10)', 'Location','Northeast');
   end
    
   
   % calculation based on second plantar flexion phase (default)
   % find first frame where gonio passes zero degrees, in plantarflexion series
   zeroangle = find(angle2>=0,1,'first');
   % we want to start from the frame right BEFORE zero, not after zero
   zeroangle = zeroangle-1;
   
   % find NEXT frame where gonio goes above X degrees, in plantarflexion series
   fiveangle = find(angle2(zeroangle:end)>stopangle,1,'first');
   if isempty(fiveangle) % if goniometer doesn't go to 5 degrees
       fiveangle = length(angle2); % then use maximal gonio angle
   else
       fiveangle = fiveangle + zeroangle-1;
   end

   % create linear equation for plantarflex direction, starting from zero angle
   [fitresult, ~] = fit_ankle_rotation(angle2(zeroangle:fiveangle), displ2(zeroangle:fiveangle),'2nd');
   coeffvals_2nd = coeffvalues(fitresult);
   
    % check if the first plantar flexion phase (startup) may also be used
   first_plantar_usable = 1;
   if angle_goingup > ten_plateau_1_start % 1st plantarflex phase does not start from dorsiflex angle (no neg angle)
       first_plantar_usable = 0;
   end
   
   if first_plantar_usable
       % calculation based on FIRST plantar flexion phase (to be used in case of crisis)
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
       [fitresult, ~] = fit_ankle_rotation(angle0(zeroangle:fiveangle), displ0(zeroangle:fiveangle),'1st');
       coeffvals_1st = coeffvalues(fitresult);
   else % will not calculate any realistic coeffvals:
       coeffvals_1st(1) = 10000;
   end
   

   
   % select the best rotation constant: the lowest one (note: should
   % preferably be negative - some movements result in positive constants)
   at_rotation_const = min(coeffvals_1st(1),coeffvals_2nd(1));
   
    
    
   
   % plot norm angle, gonio angle, displacement, zone lines
   if plot_check && plot_achilles
       % check for matlab version, do not use yyaxis on home computer (2015)
       mat_version = version('-release');
       
       plottitle = horzcat('SYNC check ankle rotation for ', subject_id);
       fig_anklerot = figure('Name', plottitle);
       plot(noraxon_rot(1:ten_plateau_2_stop,1),norm_angle_filtered(1:ten_plateau_2_stop),'b'); % norm angle
       hold on
       % left axis
       if strcmp(mat_version,'2015b') == 0
           yyaxis left
       end
       ylabel('Angle (deg)')
       plot(noraxon_rot(1:ten_plateau_2_stop,1),noraxon_rot(1:ten_plateau_2_stop,column_gonio),'m'); % goniometer 
       % right axis if possible
       if strcmp(mat_version,'2015b') == 0
           yyaxis right
           ylabel('Calc insertion displacement (mm)')
           set(gca,'YDir','Reverse')
       end
       plot(usdata_rot(1:ten_plateau_2_stop,1),usdata_rot(1:ten_plateau_2_stop,2)); % displacement
       % left axis
       if strcmp(mat_version,'2015b') == 0
          yyaxis left
       end
       plot([noraxon_rot(ten_plateau_1_start,1) noraxon_rot(ten_plateau_1_start,1)], [min(norm_angle_filtered(1:ten_plateau_2_stop)) max(norm_angle_filtered(1:ten_plateau_2_stop))],'g');
       plot([noraxon_rot(ten_plateau_1_stop,1) noraxon_rot(ten_plateau_1_stop,1)], [min(norm_angle_filtered(1:ten_plateau_2_stop)) max(norm_angle_filtered(1:ten_plateau_2_stop))],'g');
       plot([noraxon_rot(five_plateau_start,1) noraxon_rot(five_plateau_start,1)], [min(norm_angle_filtered(1:ten_plateau_2_stop)) max(norm_angle_filtered(1:ten_plateau_2_stop))],'g');
       plot([noraxon_rot(five_plateau_stop,1) noraxon_rot(five_plateau_stop,1)], [min(norm_angle_filtered(1:ten_plateau_2_stop)) max(norm_angle_filtered(1:ten_plateau_2_stop))],'g');
       plot([noraxon_rot(ten_plateau_2_start,1) noraxon_rot(ten_plateau_2_start,1)], [min(norm_angle_filtered(1:ten_plateau_2_stop)) max(norm_angle_filtered(1:ten_plateau_2_stop))],'g');
       
       xlabel('Time (s)')
       ylabel('Angle (deg)')
       title(plottitle);
       legend('Norm angle', 'Goniometer', 'Phases','Location','Northeast');
       text(0.2, 8.2, horzcat('Displ/deg = ', num2str(at_rotation_const)), 'Color', 'k');
       saveas(fig_anklerot, strcat('data_plots_stiff/IND_ankle_rot_', subject_id), 'png')
   end
   
   % output angle conversion numbers to screen, as text
   cprintf('blue',horzcat('Ankle rotation: Displacement per degree of rotation = ', num2str(at_rotation_const), ' (elim ', num2str(max(coeffvals_1st(1),coeffvals_2nd(1))), ')', '.\n'));
   
    
end