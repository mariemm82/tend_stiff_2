%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_datamaster_strength - for multiple isometric and isokinetic trials
% Marie Moltubakk 26.10.2016
%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = read_datamaster_strength(file,datamaster_columns)
    global dm_subjectno dm_timepoint dm_side dm_trial 
    global dm_isomet_P10_1 dm_isomet_P10_2 dm_isomet_D00_1 dm_isomet_D00_2 dm_isomet_D05_1 dm_isomet_D05_2 dm_isomet_D10_1 dm_isomet_D10_2 dm_isomet_D15_1 dm_isomet_D15_2
    global dm_isokinD30 dm_isokinP30 dm_isokinP45 dm_isokinP60 dm_isokinP90
    global dm_CPM_calc_NX dm_CPM_sol_NX
    global dm_leg_length

    % import datamaster file
    data = fileread(file);
    datamaster = textscan(data, '%q','Delimiter','\t');

    % restructure imported data into multiple columns
    % n o lines = 1 header + 1 per subject entry
    nolines = length(datamaster{1,1})/datamaster_columns;

    for i = 2:nolines
        delta = (i-1)*datamaster_columns;
        dm_subjectno{i-1,1}=datamaster{1,1}{delta+1,1};
        dm_timepoint{i-1,1}=datamaster{1,1}{delta+2,1};
        dm_side{i-1,1}=datamaster{1,1}{delta+3,1};
        dm_trial{i-1,1}=datamaster{1,1}{delta+4,1};
        dm_isomet_P10_1{i-1,1}=datamaster{1,1}{delta+5,1};
        dm_isomet_P10_2{i-1,1}=datamaster{1,1}{delta+6,1};
        dm_isomet_D00_1{i-1,1}=datamaster{1,1}{delta+7,1};
        dm_isomet_D00_2{i-1,1}=datamaster{1,1}{delta+8,1};
        dm_isomet_D05_1{i-1,1}=datamaster{1,1}{delta+9,1};
        dm_isomet_D05_2{i-1,1}=datamaster{1,1}{delta+10,1};
        dm_isomet_D10_1{i-1,1}=datamaster{1,1}{delta+11,1};
        dm_isomet_D10_2{i-1,1}=datamaster{1,1}{delta+12,1};
        dm_isomet_D15_1{i-1,1}=datamaster{1,1}{delta+13,1};
        dm_isomet_D15_2{i-1,1}=datamaster{1,1}{delta+14,1};
        dm_isokinD30{i-1,1}=datamaster{1,1}{delta+15,1};
        dm_isokinP30{i-1,1}=datamaster{1,1}{delta+16,1};
        dm_isokinP45{i-1,1}=datamaster{1,1}{delta+17,1};
        dm_isokinP60{i-1,1}=datamaster{1,1}{delta+18,1};
        dm_isokinP90{i-1,1}=datamaster{1,1}{delta+19,1};
        dm_CPM_calc_NX{i-1,1}=datamaster{1,1}{delta+20,1};
        dm_CPM_sol_NX{i-1,1}=datamaster{1,1}{delta+21,1};
        dm_leg_length{i-1,1}=datamaster{1,1}{delta+22,1};
    end

    out = nolines-1; % lines in datamaster to be analysed, minus header
end