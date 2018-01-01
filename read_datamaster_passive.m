%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_datamaster_passive - for combining passive dorsiflexion with US data
% Marie Moltubakk 4.2.2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = read_datamaster_passive(file)
    global dm_subjectno dm_timepoint dm_side dm_trial 
    global dm_ROM_sol1_NX dm_ROM_sol1_US dm_ROM_sol1_US_frame dm_ROM_sol2_NX dm_ROM_sol2_US dm_ROM_sol2_US_frame
    global dm_ROM_gmmtj1_NX dm_ROM_gmmtj1_US dm_ROM_gmmtj1_US_frame dm_ROM_gmmtj2_NX dm_ROM_gmmtj2_US dm_ROM_gmmtj2_US_frame
    global dm_ROM_gmfas1_NX dm_ROM_gmfas1_US dm_ROM_gmfas1_US_frame dm_ROM_gmfas2_NX dm_ROM_gmfas2_US dm_ROM_gmfas2_US_frame dm_ROM_gmfas1_licht dm_ROM_gmfas2_licht
    global dm_MVC_PF dm_MVC_DF dm_CPM_calc_NX dm_CPM_calc_US dm_CPM_calc_US_frame dm_CPM_sol_NX dm_CPM_sol_US dm_CPM_sol_US_frame 
    global dm_leg_length dm_at_SOL_length dm_at_GM_length dm_GMmsc_penn dm_GMmsc_faslen dm_ankle_angle_rest
    global dm_at_SOL_zero dm_at_GM_zero

    % import datamaster file
    data = fileread(file);
    datamaster = textscan(data, '%q','Delimiter','\t');
    
    dm_columns = 40; % number of data columns entered per subject % PROJECTSPECIFIC
    
    % restructure imported data into multiple columns
    % n o lines = 1 header + 1 per subject entry
    nolines = length(datamaster{1,1})/dm_columns;

    for i = 2:nolines
        delta = (i-1)*dm_columns;
        dm_subjectno{i-1,1}=datamaster{1,1}{delta+1,1};
        dm_timepoint{i-1,1}=datamaster{1,1}{delta+2,1};
        dm_side{i-1,1}=datamaster{1,1}{delta+3,1};
        dm_trial{i-1,1}=datamaster{1,1}{delta+4,1};
        dm_ROM_sol1_NX{i-1,1}=datamaster{1,1}{delta+5,1};
        dm_ROM_sol1_US{i-1,1}=datamaster{1,1}{delta+6,1};
        dm_ROM_sol1_US_frame{i-1,1}=datamaster{1,1}{delta+7,1};
        dm_ROM_sol2_NX{i-1,1}=datamaster{1,1}{delta+8,1};
        dm_ROM_sol2_US{i-1,1}=datamaster{1,1}{delta+9,1};
        dm_ROM_sol2_US_frame{i-1,1}=datamaster{1,1}{delta+10,1};
        dm_ROM_gmmtj1_NX{i-1,1}=datamaster{1,1}{delta+11,1};
        dm_ROM_gmmtj1_US{i-1,1}=datamaster{1,1}{delta+12,1};
        dm_ROM_gmmtj1_US_frame{i-1,1}=datamaster{1,1}{delta+13,1};
        dm_ROM_gmmtj2_NX{i-1,1}=datamaster{1,1}{delta+14,1};
        dm_ROM_gmmtj2_US{i-1,1}=datamaster{1,1}{delta+15,1};
        dm_ROM_gmmtj2_US_frame{i-1,1}=datamaster{1,1}{delta+16,1};
        dm_ROM_gmfas1_NX{i-1,1}=datamaster{1,1}{delta+17,1};
        dm_ROM_gmfas1_US{i-1,1}=datamaster{1,1}{delta+18,1};
        dm_ROM_gmfas1_US_frame{i-1,1}=datamaster{1,1}{delta+19,1};
        dm_ROM_gmfas2_NX{i-1,1}=datamaster{1,1}{delta+20,1};
        dm_ROM_gmfas2_US{i-1,1}=datamaster{1,1}{delta+21,1};
        dm_ROM_gmfas2_US_frame{i-1,1}=datamaster{1,1}{delta+22,1};
        % added for Licht
        dm_ROM_gmfas1_licht{i-1,1}=datamaster{1,1}{delta+23,1};
        dm_ROM_gmfas2_licht{i-1,1}=datamaster{1,1}{delta+24,1};
        % needed after applying Licht
        dm_MVC_PF{i-1,1}=datamaster{1,1}{delta+25,1};
        dm_MVC_DF{i-1,1}=datamaster{1,1}{delta+26,1};
        dm_CPM_calc_NX{i-1,1}=datamaster{1,1}{delta+27,1};
        dm_CPM_calc_US{i-1,1}=datamaster{1,1}{delta+28,1};
        dm_CPM_calc_US_frame{i-1,1}=datamaster{1,1}{delta+29,1};
        dm_CPM_sol_NX{i-1,1}=datamaster{1,1}{delta+30,1};
        dm_CPM_sol_US{i-1,1}=datamaster{1,1}{delta+31,1};
        dm_CPM_sol_US_frame{i-1,1}=datamaster{1,1}{delta+32,1};
        dm_leg_length{i-1,1}=datamaster{1,1}{delta+33,1};
        dm_at_SOL_length{i-1,1}=datamaster{1,1}{delta+34,1};
        dm_at_GM_length{i-1,1}=datamaster{1,1}{delta+35,1};
        dm_GMmsc_penn{i-1,1}=datamaster{1,1}{delta+36,1};
        dm_GMmsc_faslen{i-1,1}=datamaster{1,1}{delta+37,1};
        dm_ankle_angle_rest{i-1,1}=datamaster{1,1}{delta+38,1};
        dm_at_SOL_zero{i-1,1}=datamaster{1,1}{delta+39,1};
        dm_at_GM_zero{i-1,1}=datamaster{1,1}{delta+40,1};
    end

    out = nolines-1; % lines in datamaster to be analysed, minus header
end