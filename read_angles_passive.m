%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_angles_passive
% Marie Moltubakk + REB + M 25.5.2015
% 
% used by passiveUS to import pre-generated angle data across subjects
% (generated in create_angles_passive)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_angles_passive(file)
%    global ang_subjectno 
%    global ang_pre_r ang_pre_l ang_post_r ang_post_l ang_ind_max ang_common_max
    global input_gon_pre_r input_gon_pre_l input_gon_post_r input_gon_post_l input_gon_ind_max input_gon_common_max input_gon_ind_rmax input_gon_ind_lmax
    
    % import datamaster file
    data = fileread(file);
    datamaster = textscan(data, '%q');

    % restructure imported data into multiple columns
    datamaster_columns = 15; % number of columns in file %VAR
    nolines = length(datamaster{1,1})/datamaster_columns;

    for i = 1:nolines
        delta = (i-1)*datamaster_columns;
%         ang_subjectno{i,1}=datamaster{1,1}{delta+1,1};
%         ang_pre_r{i,1}=datamaster{1,1}{delta+2,1};
%         ang_pre_l{i,1}=datamaster{1,1}{delta+3,1};
%         ang_post_r{i,1}=datamaster{1,1}{delta+4,1};
%         ang_post_l{i,1}=datamaster{1,1}{delta+5,1};
%         ang_ind_max{i,1}=datamaster{1,1}{delta+6,1};
%         ang_common_max{i,1}=datamaster{1,1}{delta+7,1};
        input_gon_pre_r{i,1}=datamaster{1,1}{delta+8,1};
        input_gon_pre_l{i,1}=datamaster{1,1}{delta+9,1};
        input_gon_post_r{i,1}=datamaster{1,1}{delta+10,1};
        input_gon_post_l{i,1}=datamaster{1,1}{delta+11,1};
        input_gon_ind_max{i,1}=datamaster{1,1}{delta+12,1};
        input_gon_common_max{i,1}=datamaster{1,1}{delta+13,1};
        input_gon_ind_rmax{i,1}=datamaster{1,1}{delta+14,1};
        input_gon_ind_lmax{i,1}=datamaster{1,1}{delta+15,1};
    end
end