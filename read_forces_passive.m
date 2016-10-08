%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_forces_passive
% Marie Moltubakk 5.8.2016
%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_forces_passive(file)
    %global input_for_subjectno
    global input_for_pre_r input_for_pre_l input_for_post_r input_for_post_l input_for_ind_max input_for_common_max input_for_ind_rmax input_for_ind_lmax
    
    % import datamaster file
    data = fileread(file);
    datamaster = textscan(data, '%q');

    % restructure imported data into multiple columns
    datamaster_columns = 9; % number of columns in file %VAR
    nolines = length(datamaster{1,1})/datamaster_columns;

    for i = 1:nolines
        delta = (i-1)*datamaster_columns;
        %input_for_subjectno{i,1}=datamaster{1,1}{delta+1,1};
        input_for_pre_r{i,1}=datamaster{1,1}{delta+2,1};
        input_for_pre_l{i,1}=datamaster{1,1}{delta+3,1};
        input_for_post_r{i,1}=datamaster{1,1}{delta+4,1};
        input_for_post_l{i,1}=datamaster{1,1}{delta+5,1};
        input_for_ind_max{i,1}=datamaster{1,1}{delta+6,1};
        input_for_common_max{i,1}=datamaster{1,1}{delta+7,1};
        input_for_ind_rmax{i,1}=datamaster{1,1}{delta+8,1};
        input_for_ind_lmax{i,1}=datamaster{1,1}{delta+9,1};
    end
end