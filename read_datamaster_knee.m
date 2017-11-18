%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_datamaster_knee - for passive knee extension
% Marie Moltubakk 16.11.2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = read_datamaster_knee(file)
    global dm_subjectno dm_timepoint dm_side dm_trial 
    global dm_ROM_trial1 dm_ROM_trial2 dm_CPM_calc_NX dm_CPM_sol_NX 
    global dm_n_o_trials
    global dm_ROM_trial dm_ROM_ind dm_ROM_common

    % import datamaster file
    data = fileread(file);
    datamaster = textscan(data, '%q','Delimiter','\t');
    
    dm_columns = 12; % number of data columns entered per subject % PROJECTSPECIFIC
    
    % restructure imported data into multiple columns
    % n o lines = 1 header + 1 per subject entry
    nolines = length(datamaster{1,1})/dm_columns;

    for i = 2:nolines
        delta = (i-1)*dm_columns;
        dm_subjectno{i-1,1}=datamaster{1,1}{delta+1,1};
        dm_timepoint{i-1,1}=datamaster{1,1}{delta+2,1};
        dm_side{i-1,1}=datamaster{1,1}{delta+3,1};
        dm_trial{i-1,1}=datamaster{1,1}{delta+4,1};
        dm_ROM_trial1{i-1,1}=datamaster{1,1}{delta+5,1};
        dm_ROM_trial2{i-1,1}=datamaster{1,1}{delta+6,1};
        dm_CPM_calc_NX{i-1,1}=datamaster{1,1}{delta+7,1};
        dm_CPM_sol_NX{i-1,1}=datamaster{1,1}{delta+8,1};
        dm_n_o_trials{i-1,1}=datamaster{1,1}{delta+9,1};
        dm_ROM_trial{i-1,1}=datamaster{1,1}{delta+10,1};
        dm_ROM_ind{i-1,1}=datamaster{1,1}{delta+11,1};
        dm_ROM_common{i-1,1}=datamaster{1,1}{delta+12,1};
    end

    out = nolines-1; % lines in datamaster to be analysed, minus header
end