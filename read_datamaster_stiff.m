%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_datamaster file for TENDON STIFFNESS
% Marie Moltubakk 17.5.2013
% 
% NB 2017-10-15: BD study has only 29? columns. If going back, must add dummy columns:
%      dm_group = leg's group (stretching or control)
%      dm_tendonlength = tendon length at 0 deg ankle angle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = read_datamaster_stiff(file)
    global dm_subjectno dm_timepoint dm_side dm_trial dm_group
    global dm_stiff1_NX dm_stiff1_US dm_stiff1_US_frame dm_stiff2_NX dm_stiff2_US dm_stiff2_US_frame dm_stiff3_NX dm_stiff3_US dm_stiff3_US_frame 
    global dm_heel1_NX dm_heel1_US dm_heel1_US_frame dm_heel2_NX dm_heel2_US dm_heel2_US_frame dm_heel3_NX dm_heel3_US dm_heel3_US_frame
    global dm_MVC_PF dm_MVC_DF dm_CPM_calc_NX dm_CPM_calc_US dm_CPM_calc_US_frame dm_CPM_sol_NX dm_leg_length
    global dm_tendonlength dm_cutforce %new2014-04-14 + 2017-11-15
    global dm_rot_const

    
    % import datamaster file
    data = fileread(file);
    datamaster = textscan(data, '%q','Delimiter','\t');
    
    datamaster_columns = 33; % number of data columns entered per subject % PROJECTSPECIFIC
    
    % restructure imported data into multiple columns
    % n o lines = 1 header + 1 per subject entry
    nolines = length(datamaster{1,1})/datamaster_columns;

    for i = 2:nolines
        delta = (i-1)*datamaster_columns;
        dm_subjectno{i-1,1}=datamaster{1,1}{delta+1,1};
        dm_timepoint{i-1,1}=datamaster{1,1}{delta+2,1};
        dm_side{i-1,1}=datamaster{1,1}{delta+3,1};
        dm_trial{i-1,1}=datamaster{1,1}{delta+4,1};
        dm_group{i-1,1}=datamaster{1,1}{delta+5,1}; % new 2017-10-15
        
        dm_stiff1_NX{i-1,1}=datamaster{1,1}{delta+6,1};
        dm_stiff1_US{i-1,1}=datamaster{1,1}{delta+7,1};
        dm_stiff1_US_frame{i-1,1}=datamaster{1,1}{delta+8,1};
        dm_stiff2_NX{i-1,1}=datamaster{1,1}{delta+9,1};
        dm_stiff2_US{i-1,1}=datamaster{1,1}{delta+10,1};
        dm_stiff2_US_frame{i-1,1}=datamaster{1,1}{delta+11,1};
        dm_stiff3_NX{i-1,1}=datamaster{1,1}{delta+12,1};
        dm_stiff3_US{i-1,1}=datamaster{1,1}{delta+13,1};
        dm_stiff3_US_frame{i-1,1}=datamaster{1,1}{delta+14,1};
        
        dm_heel1_NX{i-1,1}=datamaster{1,1}{delta+15,1};
        dm_heel1_US{i-1,1}=datamaster{1,1}{delta+16,1};
        dm_heel1_US_frame{i-1,1}=datamaster{1,1}{delta+17,1};
        dm_heel2_NX{i-1,1}=datamaster{1,1}{delta+18,1};
        dm_heel2_US{i-1,1}=datamaster{1,1}{delta+19,1};
        dm_heel2_US_frame{i-1,1}=datamaster{1,1}{delta+20,1};
        dm_heel3_NX{i-1,1}=datamaster{1,1}{delta+21,1};
        dm_heel3_US{i-1,1}=datamaster{1,1}{delta+22,1};
        dm_heel3_US_frame{i-1,1}=datamaster{1,1}{delta+23,1};
        
        dm_MVC_PF{i-1,1}=datamaster{1,1}{delta+24,1};
        dm_MVC_DF{i-1,1}=datamaster{1,1}{delta+25,1};
        dm_CPM_calc_NX{i-1,1}=datamaster{1,1}{delta+26,1};
        dm_CPM_calc_US{i-1,1}=datamaster{1,1}{delta+27,1};
        dm_CPM_calc_US_frame{i-1,1}=datamaster{1,1}{delta+28,1};

        dm_CPM_sol_NX{i-1,1}=datamaster{1,1}{delta+29,1};
        
        dm_leg_length{i-1,1}=datamaster{1,1}{delta+30,1};
        dm_tendonlength{i-1,1}=datamaster{1,1}{delta+31,1}; % new 2017-11-15
        dm_cutforce{i-1,1}=datamaster{1,1}{delta+32,1}; %new2014-04-14
        dm_rot_const{i-1,1}=datamaster{1,1}{delta+33,1}; %new2017-12-08
    end

    out = nolines-1; % lines in datamaster to be analysed, minus header
end