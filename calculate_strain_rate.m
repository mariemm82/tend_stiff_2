%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_strain_rate
% Marie Moltubakk 16.6.2013
% Read 6 trials containing time-force-displacement data
% Produce one variable (strain_rate) with calculated strain rate for each of the 6 trials
%%%%%%%%%%%%%%%%%%%%%%%%%%



function strain_rate = calculate_strain_rate(time_force_displ_mtj1, time_force_displ_mtj2, time_force_displ_mtj3, time_force_displ_otj1, time_force_displ_otj2, time_force_displ_otj3)
    global plot_achilles plot_norm plot_emg plot_check subject_id

    % TODOLATER (IF STRAIN RATE NEEDED):
    % change to averaged MTJ displacement minus averaged OTJ displacement - what to use for time?
    % subract start time when displacement changes?

    
    % largest deformation 
    commonforce = min([max(time_force_displ_mtj1(:,2)) max(time_force_displ_mtj2(:,2)) max(time_force_displ_mtj3(:,2)) max(time_force_displ_otj1(:,2)) max(time_force_displ_otj2(:,2)) max(time_force_displ_otj3(:,2))]);
    
    % detect array contents after commonforce
    del_mtj1 = find(time_force_displ_mtj1(:,2)>=commonforce, 1, 'first');
    del_mtj2 = find(time_force_displ_mtj2(:,2)>=commonforce, 1, 'first');
    del_mtj3 = find(time_force_displ_mtj3(:,2)>=commonforce, 1, 'first');
    del_otj1 = find(time_force_displ_otj1(:,2)>=commonforce, 1, 'first');
    del_otj2 = find(time_force_displ_otj2(:,2)>=commonforce, 1, 'first');
    del_otj3 = find(time_force_displ_otj3(:,2)>=commonforce, 1, 'first');
    
    % delete array contents after commonforce
    time_force_displ_mtj1(del_mtj1:length(time_force_displ_mtj1),:) = [];
    time_force_displ_mtj2(del_mtj2:length(time_force_displ_mtj2),:) = [];
    time_force_displ_mtj3(del_mtj3:length(time_force_displ_mtj3),:) = [];
    time_force_displ_otj1(del_otj1:length(time_force_displ_otj1),:) = [];
    time_force_displ_otj2(del_otj2:length(time_force_displ_otj2),:) = [];
    time_force_displ_otj3(del_otj3:length(time_force_displ_otj3),:) = [];
    
    % find maximal deformation per trial, within commonforce
    max_displ1 = max(time_force_displ_mtj1(:,3));
    index_max_mtj1 = find(time_force_displ_mtj1(:,3)==max_displ1, 1, 'first');
    max_displ2 = max(time_force_displ_mtj2(:,3));
    index_max_mtj2 = find(time_force_displ_mtj2(:,3)==max_displ2, 1, 'first');
    max_displ3 = max(time_force_displ_mtj3(:,3));
    index_max_mtj3 = find(time_force_displ_mtj3(:,3)==max_displ3, 1, 'first');
    max_displ4 = max(time_force_displ_otj1(:,3));
    index_max_otj1 = find(time_force_displ_otj1(:,3)==max_displ4, 1, 'first');
    max_displ5 = max(time_force_displ_otj2(:,3));
    index_max_otj2 = find(time_force_displ_otj2(:,3)==max_displ5, 1, 'first');
    max_displ6 = max(time_force_displ_otj3(:,3));
    index_max_otj3 = find(time_force_displ_otj3(:,3)==max_displ6, 1, 'first');
    
    % strain rate = max deformation / time
    strain_rate(1) = time_force_displ_mtj1(index_max_mtj1,3) / time_force_displ_mtj1(index_max_mtj1,2);
    strain_rate(2) = time_force_displ_mtj2(index_max_mtj2,3) / time_force_displ_mtj2(index_max_mtj2,2);
    strain_rate(3) = time_force_displ_mtj3(index_max_mtj3,3) / time_force_displ_mtj3(index_max_mtj3,2);
    strain_rate(4) = time_force_displ_otj1(index_max_otj1,3) / time_force_displ_otj1(index_max_otj1,2);
    strain_rate(5) = time_force_displ_otj2(index_max_otj2,3) / time_force_displ_otj2(index_max_otj2,2);
    strain_rate(6) = time_force_displ_otj3(index_max_otj3,3) / time_force_displ_otj3(index_max_otj3,2);

end