%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main file for passive trial analysis
% Marie Moltubakk 29.4.2014
% 
% Note 1:
% The scripts assume a certain structure for input files from Tracker and
% Noraxon. I.e. number of and order of EMG and other channels. If these
% scripts are to be used for other projects, some of the code must be
% modified. These lines are marked % PROJECTSPECIFIC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all

global plot_achilles plot_norm plot_emg plot_check plot_us subject_id

plot_check = 1; % turn on/off main checkpoint plots
plot_achilles = 0; % turn on/off troubleshoot plots
plot_norm = 0; 
plot_emg = 0;  % RMS 3 EMG channels per trial
plot_us = 0;



%%% Set constants % PROJECTSPECIFIC
global us_zerodispframes noraxonfreq freq_default emg_bandpass emg_rms_ms mvc_window_ms torque_cutoff angle_cutoff convert_achilles convert_norm

us_zerodispframes = 1; % No of US frames to average as zero displacement

noraxonfreq = 1500; % sampling frequency of noraxon data
freq_default = 100; % output frequency for noraxon data without US video (MVC, etc)

emg_bandpass = [10/(noraxonfreq/2) 500/(noraxonfreq/2)]; % cutoff frequencies for EMG butterworth bandpass filter
emg_rms_ms = 100; % milliseconds RMS window, must be divisible by 4 - ref Basmajian 1985 = 50 ms, ref Aagaard paper Passive tensile stress = 200 ms
mvc_window_ms = 500; % milliseconds window for determining MVC torque and EMG

angle_cutoff = 10/(noraxonfreq/2); % cutoff freq, Norm filtering angle
torque_cutoff = [0.36/(noraxonfreq/2) 0.42/(noraxonfreq/2)]; % cutoff freq, Norm filtering - ref Winter 1990 = 15 hz. Kongsgaard = 8hz

% convertion factors Norm, µV -> angle/torque/velocity: 
%convert_norm(1) = 0.09918;      %convert_norm_angle_a
%convert_norm(2) = -168.714827;  %convert_norm_angle_b
%convert_norm(3) = 3518.257824066; %convert_norm_angle_c
convert_norm(5) = 0.076965;     %convert_norm_torque_a
convert_norm(6) = 7.862109;     %convert_norm_torque_b
convert_norm(7) = 0.20535;      %convert_norm_velocity_a
convert_norm(8) = 10.92228;     %convert_norm_velocity_b

% column placement in Noraxon data
global column_EMG_start column_EMG_end column_l_gm column_r_gm column_l_gl column_r_gl column_l_sol column_r_sol column_l_tibant column_r_tibant column_gonio column_norm_angle column_norm_torque column_norm_velocity column_norm_direction column_achilles
column_EMG_start = 2; 
column_EMG_end = 9;
column_l_gm = 2;
column_r_gm = 3;
column_l_gl = 4;
column_r_gl = 5;
column_l_sol = 6;
column_r_sol = 7;
column_l_tibant = 8;
column_r_tibant = 9;
column_gonio = 10;
column_norm_angle = 11;
column_norm_torque = 12;
column_norm_velocity = 13;
column_norm_direction = 14;
column_achilles = 15;





%%% Loop through all files in folder

allROMFiles = dir('data_passive\');
allROMNames = {allROMFiles.name};
NOFiles = length(allROMNames);



% prepare variables to check for matching up trials
working_subject = '0';
working_timepoint = '0';
working_side = '0';
new_subject = 1;
new_timepoint = 1;
new_side = 1;
set_no = 0;

for line = 1:NOFiles-3; % filename nr 1 and 2 are '.' and '..', last filename is CPM folder - ignore these three
    
    %%% input data for file to be analysed
    filename = allROMNames{line+2}; % filename nr 1 and 2 are '.' and '..'
    inputfile = horzcat('data_passive/', filename);
    subject_nr = filename(1:2); % return 2 first characters
    bb = length(filename);
    cc = strfind(filename, 'ROM');
    timepoint = filename(cc+4); % return TIMEPOINT part of filename
    side = filename(cc+5); % return SIDE part of filename
    trialname = filename(cc+7:bb-4); % return trial part of filename
    
    
    
    %%% check if data belong to same or other subject / timepoint / side
    if strcmp(working_subject,subject_nr)
        new_subject = 0;
    else
        new_subject = 1;
        working_subject = subject_nr;
    end
    if strcmp(working_timepoint,timepoint)
        new_timepoint = 0;
    else
        new_timepoint = 1;
        working_timepoint = timepoint;
    end
    if strcmp(working_side,side)
        new_side = 0;
    else
        new_side = 1;
        working_side = side;
    end
    if new_subject || new_timepoint || new_side
        set_no = set_no+1;
    end
    
    
    
    %%% if new trial, print header
    if (new_subject || new_timepoint || new_side)
        subject_id = horzcat('subject ', num2str(subject_nr), ' ', timepoint, side, ' ', trialname);
        disp(sprintf(horzcat('----------------', subject_id, '------------------')))
    end
    
    
    
    %%% if new trial:
    %%% Calculate individual Norm conversion factors (a,b) for ANGLE
    if (new_subject || new_timepoint || new_side)
    
        % identify correct CPM file
        if timepoint == '2'
            timetxt = 'PRE';
        elseif timepoint == '5'
            timetxt = 'POST';
        else
            timetxt = 'ERROR';
        end
        CPMfilename = horzcat('data_passive\CPM\', subject_nr, '_', timetxt, '_CPM_', timepoint, side, '.tsv');

        % Read raw torque data from Norm, filter
        % Produce individual conversion factors for angle
        [convert_ind_angle_a, convert_ind_angle_b] = calculate_angle_constants(angle_cutoff, CPMfilename);

        % place individually calculated constant into array with common constants 
        convert_norm(1) = convert_ind_angle_a;  
        convert_norm(2) = convert_ind_angle_b;  
        convert_norm(3) = 0; 
    end    
    
    
    
    %%% for all trials:
    %%% Read stiffness trial Noraxon data file, set first frame as time = zero, EMG+torque data treatment, resample
    % Produce a new noraxon data array
    noraxon_prepped = read_noraxon_passive(inputfile, freq_default, side, trialname); %output freq, side, trialname
    
    
    
    
    %%% for all trials: 
    %%% extract and print variables individually
    
    % max range of motion
    max_ROM = min(noraxon_prepped(:,column_norm_angle));
    % torque at max ROM
    torque_max_ROM = max(noraxon_prepped(:,column_norm_torque));
    report = sprintf(horzcat('Passive torque at ', num2str(round(max_ROM)), ' = ', num2str(torque_max_ROM), ' Nm (max ROM).'));
    disp(report)
    
    % extract torque at X and Y degrees
    % (negative = dorsiflex)
    if new_subject %set new angles
        array_maxangles = [1 -32.9; 2 -24.3; 3 -21.9; 4 -40.7; 5 -18.6; 6 -15.1; 7 -24.1; 8 -11.7; 9 -31.9; 10 -15.3; 11 -21.2; 13 -13.9; 15 -14.3; 16 -11.2; 18 -18.2; 19 -11.9; 20 -16.1; 21 -10.4; 22 -10.7; 24 -11.2; 25 -17.5; 26 -18.6; 28 -14.1; 29 -6.0; 30 -9.6; 31 -6.0]; %VAR
        row = find(array_maxangles(:,1)==str2num(subject_nr));
        if row
            angle1 = array_maxangles(row,2); %VAR 
            angle2 = -6; %VAR
        else
            report = sprintf(horzcat('Error, subject number ', subject_nr, ' not found in array - using default angles'));
            disp(report)
            angle1 = -15; %VAR 
            angle2 = -6; %VAR
        end
    end
    
    if max_ROM <= angle1
        loc_angle1 = find(noraxon_prepped(:,column_norm_angle)<angle1,1,'first');
        torque_angle1 = noraxon_prepped(loc_angle1,column_norm_torque);
        report = sprintf(horzcat('Passive torque at ', num2str(angle1), ' = ', num2str(torque_angle1), ' Nm.'));
    else
        torque_angle1 = -1000;
        report = sprintf(horzcat('Subject does not reach ', num2str(angle1), ' degrees.'));
    end
    disp(report)
    if max_ROM <= angle2
        loc_angle2 = find(noraxon_prepped(:,column_norm_angle)<angle2,1,'first');
        torque_angle2 = noraxon_prepped(loc_angle2,column_norm_torque);
        report = sprintf(horzcat('Passive torque at  ', num2str(angle2), ' = ', num2str(torque_angle2), ' Nm.'));
    else
        torque_angle2 = -1000;
        report = sprintf(horzcat('Subject does not reach ', num2str(angle2), ' degrees.'));
    end
    disp(report)
    
    
    
    %%% collect data belonging to one subject ... NEW JUNE 2014
    
    
    
    
    
    %%% if at the end of a subject / trial / side: 
    %%% Output final data to array
    
    % add data to a common array for all subjects    
    all_output_txt(line,1) = {subject_nr};
    all_output_txt(line,2) = {timepoint};
    all_output_txt(line,3) = {side};
    all_output_txt(line,4) = {trialname};
    all_output_txt(line,5) = {max_ROM};
    all_output_txt(line,6) = {torque_max_ROM};
    all_output_txt(line,7) = {angle1};
    all_output_txt(line,8) = {torque_angle1};
    all_output_txt(line,9) = {angle2};
    all_output_txt(line,10) = {torque_angle2};
end



%%% Output key variables for all subjects to file

% write xls
if ispc
    filename_output = strcat('data_output/passive_', datestr(now, 'yyyy-mm-dd HH-MM'));
    all_output_head = {'Subject', 'Timepoint', 'Side', 'Trial', 'Max ROM', 'Max torque', 'Angle1', 'Torque a1', 'Angle2', 'Torque a2'}; % PROJECTSPECIFIC
    xlswrite(filename_output, all_output_head, 1, 'A1')
    xlswrite(filename_output, all_output_txt, 1, 'A2')
else
    filename_output = strcat('data_output/passive_', datestr(now, 'yyyy-mm-dd HH-MM'), '.csv');
    csvwrite(filename_output, all_output_txt)
end






%%% TODO
% more appropriate velocity conversion favtors