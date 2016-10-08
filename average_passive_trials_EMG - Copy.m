%%%%%%%%%%%%%%%%%%%%%%%%%%
% average_passive_trials WITH EMG
% Marie Moltubakk 6.2.2015
% Read two trials of torque, gonio angle, norm angle, displacement + EMG
% Produce array with average data
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [output_array] = average_passive_trials_EMG(varargin) % (torque1, gonio1, angle1, displ1, emg_gm1, emg_gl1, emg_sol1, torque2, gonio2, angle2, displ2, emg_gm2, emg_gl2, emg_sol2)

method = 2; % MMM TMP - two different ways of consolidating gonio data




if nargin == 14 % two trials submitted, to be averaged %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    torque1 = varargin{1};
    gonio1 = varargin{2};
    angle1 = varargin{3};
    displ1 = varargin{4};
    emg_gm1 = varargin{5};
    emg_gl1 = varargin{6};
    emg_sol1 = varargin{7};
    torque2 = varargin{8};
    gonio2 = varargin{9};
    angle2 = varargin{10};
    displ2 = varargin{11};
    emg_gm2 = varargin{12};
    emg_gl2 = varargin{13};
    emg_sol2 = varargin{14};
    
   
    
    if method == 1
        %%% if gonio angle starts after zero degrees, extrapolate all data

        % trial 1
        if min(gonio1) > 0.0
            % select first 2 degrees of existing gonio data
            angle_end = min(gonio1) + 5; %VAR extrapolate using first X degrees of the gonio movement
            loc_angle_end = find(gonio1 >= angle_end,1,'first');
            % calculate linear coeffisients
            gonio_modified = [ones(length(gonio1(1:loc_angle_end)),1) gonio1(1:loc_angle_end)];
            coeffs_torque = gonio_modified\torque1(1:loc_angle_end); % NB backslash --> slope
            coeffs_displ = gonio_modified\displ1(1:loc_angle_end);
            coeffs_emg_gm = gonio_modified\emg_gm1(1:loc_angle_end);
            coeffs_emg_gl = gonio_modified\emg_gl1(1:loc_angle_end);
            coeffs_emg_sol = gonio_modified\emg_sol1(1:loc_angle_end);
            % enlarge arrays by adding values at zero angle to the front
            angle_half = min(gonio1)/2;
            gonio1 = [0; angle_half; gonio1];
            torque1 = [coeffs_torque(1); (angle_half*coeffs_torque(2)) + coeffs_torque(1); torque1];
            displ1 = [coeffs_displ(1); (angle_half*coeffs_displ(2)) + coeffs_displ(1); displ1];
            emg_gm1 = [coeffs_emg_gm(1); (angle_half*coeffs_emg_gm(2)) + coeffs_emg_gm(1); emg_gm1];
            emg_gl1 = [coeffs_emg_gl(1); (angle_half*coeffs_emg_gl(2)) + coeffs_emg_gl(1); emg_gl1];
            emg_sol1 = [coeffs_emg_sol(1); (angle_half*coeffs_emg_sol(2)) + coeffs_emg_sol(1); emg_sol1];
        end

        % repeat for trial 2
        if min(gonio2) > 0.0
            % select first 2 degrees of existing gonio data
            angle_end = min(gonio2) + 5; %VAR extrapolate using first X degrees of the gonio movement
            loc_angle_end = find(gonio2 >= angle_end,1,'first');
            % calculate linear coeffisients
            gonio_modified = [ones(length(gonio2(1:loc_angle_end)),1) gonio2(1:loc_angle_end)];
            coeffs_torque = gonio_modified\torque2(1:loc_angle_end); % NB backslash --> slope
            coeffs_displ = gonio_modified\displ2(1:loc_angle_end); % NB backslash --> slope
            coeffs_emg_gm = gonio_modified\emg_gm2(1:loc_angle_end);
            coeffs_emg_gl = gonio_modified\emg_gl2(1:loc_angle_end);
            coeffs_emg_sol = gonio_modified\emg_sol2(1:loc_angle_end);
            % enlarge arrays by adding values at zero angle to the front
            angle_half = min(gonio2)/2;
            gonio2 = [0; angle_half; gonio2];
            torque2 = [coeffs_torque(1); (angle_half*coeffs_torque(2)) + coeffs_torque(1); torque2];
            displ2 = [coeffs_displ(1); (angle_half*coeffs_displ(2)) + coeffs_displ(1); displ2];
            emg_gm2 = [coeffs_emg_gm(1); (angle_half*coeffs_emg_gm(2)) + coeffs_emg_gm(1); emg_gm2];
            emg_gl2 = [coeffs_emg_gl(1); (angle_half*coeffs_emg_gl(2)) + coeffs_emg_gl(1); emg_gl2];
            emg_sol2 = [coeffs_emg_sol(1); (angle_half*coeffs_emg_sol(2)) + coeffs_emg_sol(1); emg_sol2];
        end
        


        %%% create array with common angles (.05 intervals) for resampling

        % select the smallest range of the gonio data as basis for angle array
        common_angle_start = max([min(gonio_1) min(gonio_2)]);
        common_angle_stop = min([max(gonio_1) max(gonio_2)]);
        % round to .05 
        common_angle_start_rounded = ceil(common_angle_start/0.05)*0.05;
        % create array of angles
        average_angle_array = common_angle_start_rounded:0.05:floor(common_angle_stop/0.05)*0.05;


        
        %%% reshape, resample, average

        % if there are duplicate in gonio arrays, consolidate data and delete duplicate
        % http://www.mathworks.com/matlabcentral/fileexchange/8354-consolidator
        [~,torque1unique,~] = consolidator(gonio1,torque1,'mean',0);
        [~,torque2unique,~] = consolidator(gonio2,torque2,'mean',0);
        [~,displ1unique,~] = consolidator(gonio1,displ1,'mean',0);
        [~,displ2unique,~] = consolidator(gonio2,displ2,'mean',0);
        [~,emg_sol1unique,~] = consolidator(gonio1,emg_sol1,'mean',0);
        [~,emg_sol2unique,~] = consolidator(gonio2,emg_sol2,'mean',0);
        [~,emg_gm1unique,~] = consolidator(gonio1,emg_gm1,'mean',0);
        [~,emg_gm2unique,~] = consolidator(gonio2,emg_gm2,'mean',0);
        [gonio1unique,emg_gl1unique,~] = consolidator(gonio1,emg_gl1,'mean',0);
        [gonio2unique,emg_gl2unique,~] = consolidator(gonio2,emg_gl2,'mean',0);

        % reshape and average TORQUE across common angle array
        common_torque1_gonio = spline(gonio1unique, torque1unique, average_angle_array);
        common_torque2_gonio = spline(gonio2unique, torque2unique, average_angle_array);
        average_torque_gonio = (common_torque1_gonio + common_torque2_gonio) / 2;

        % reshape and average DISPLACEMENT across common angle array
        common_displ1_gonio = spline(gonio1unique, displ1unique, average_angle_array);
        common_displ2_gonio = spline(gonio2unique, displ2unique, average_angle_array);
        average_displ_gonio = (common_displ1_gonio + common_displ2_gonio) / 2;

        % reshape and average EMG across common angle array
        common_emg_gm1_gonio = spline(gonio1unique, emg_gm1unique, average_angle_array);  % do EMG_gm1 contents need to be UNIQUE, as is done for torque and displacement?
        common_emg_gm2_gonio = spline(gonio2unique, emg_gm2unique, average_angle_array);
        average_emg_gm_gonio = (common_emg_gm1_gonio + common_emg_gm2_gonio) / 2;
        common_emg_gl1_gonio = spline(gonio1unique, emg_gl1unique, average_angle_array);  % do EMG_gl1 contents need to be UNIQUE, as is done for torque and displacement?
        common_emg_gl2_gonio = spline(gonio2unique, emg_gl2unique, average_angle_array);
        average_emg_gl_gonio = (common_emg_gl1_gonio + common_emg_gl2_gonio) / 2;
        common_emg_sol1_gonio = spline(gonio1unique, emg_sol1unique, average_angle_array);  % do EMG_sol1 contents need to be UNIQUE, as is done for torque and displacement?
        common_emg_sol2_gonio = spline(gonio2unique, emg_sol2unique, average_angle_array);
        average_emg_sol_gonio = (common_emg_sol1_gonio + common_emg_sol2_gonio) / 2;

        
        
        
    elseif method == 2

        %%% curve fitting gonio data - GONIO 1
        
        % method 1: fit 4th order polynomial to averaged gonio-angle curve
        fit_gonio1 = polyfit(angle1, gonio1, 4);
        
        % create array across angles
        gonio_new1(1:length(angle1),1) = zeros;
        for ang = 1:length(angle1)
            gonio_new1(ang) = (fit_gonio1(1) * angle1(ang)^4) + (fit_gonio1(2) * angle1(ang)^3) + (fit_gonio1(3) * angle1(ang)^2) + (fit_gonio1(4) * angle1(ang)) + fit_gonio1(5);
        end
        
        %    % method 2: cubic smoothing spline
        %    p = 0.0005; % 0 = least squares straight line fit, 1 = natural cubic spline interpolant   %VAR
        %    gonio_new = csaps(angle1, gonio1, p, time_torque_ascend); 
        
        %%% curve fitting gonio data - GONIO 2
        
        % method 1: fit 4th order polynomial to averaged gonio-angle curve
        fit_gonio2 = polyfit(angle2, gonio2, 4);
                
        % create array across angles
        gonio_new2(1:length(angle2),1) = zeros;
        for ang = 1:length(angle2)
            gonio_new2(ang) = (fit_gonio2(1) * angle2(ang)^4) + (fit_gonio2(2) * angle2(ang)^3) + (fit_gonio2(3) * angle2(ang)^2) + (fit_gonio2(4) * angle2(ang)) + fit_gonio2(5);
        end
        
        % plot check gonio curve fit
        figure,plot(angle1,gonio1) % TMP MMM
        hold on
        plot(angle1,gonio_new1)
        plot(angle2,gonio2)
        plot(angle2,gonio_new2)




        %%% create array with common angles (.05 intervals) for resampling

        % select the smallest range of the gonio data as basis for angle array
        common_angle_start = max([min(gonio_new1) min(gonio_new2)]);
        common_angle_stop = min([max(gonio_new1) max(gonio_new2)]);
        % round to .05 
        common_angle_start_rounded = ceil(common_angle_start/0.05)*0.05;
        % create array of angles
        average_angle_array = common_angle_start_rounded:0.05:floor(common_angle_stop/0.05)*0.05;




        %%% reshape, resample, average

        % reshape and average TORQUE across common angle array
        common_torque1_gonio = spline(gonio_new1, torque1, average_angle_array);
        common_torque2_gonio = spline(gonio_new2, torque2, average_angle_array);
        average_torque_gonio = (common_torque1_gonio + common_torque2_gonio) / 2;

        % reshape and average DISPLACEMENT across common angle array
        common_displ1_gonio = spline(gonio_new1, displ1, average_angle_array);
        common_displ2_gonio = spline(gonio_new2, displ2, average_angle_array);
        average_displ_gonio = (common_displ1_gonio + common_displ2_gonio) / 2;

        % reshape and average EMG across common angle array
        common_emg_gm1_gonio = spline(gonio_new1, emg_gm1, average_angle_array);
        common_emg_gm2_gonio = spline(gonio_new2, emg_gm2, average_angle_array);
        average_emg_gm_gonio = (common_emg_gm1_gonio + common_emg_gm2_gonio) / 2;
        common_emg_gl1_gonio = spline(gonio_new1, emg_gl1, average_angle_array);
        common_emg_gl2_gonio = spline(gonio_new2, emg_gl2, average_angle_array);
        average_emg_gl_gonio = (common_emg_gl1_gonio + common_emg_gl2_gonio) / 2;
        common_emg_sol1_gonio = spline(gonio_new1, emg_sol1, average_angle_array);
        common_emg_sol2_gonio = spline(gonio_new2, emg_sol2, average_angle_array);
        average_emg_sol_gonio = (common_emg_sol1_gonio + common_emg_sol2_gonio) / 2;



    end

    
    
    
    
    
else % nargin == 7, one trial submitted %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    torque1 = varargin{1};
    gonio1 = varargin{2};
    angle1 = varargin{3};
    displ1 = varargin{4};
    emg_gm1 = varargin{5};
    emg_gl1 = varargin{6};
    emg_sol1 = varargin{7};


    
    if method == 1

        %%% if gonio angle starts after zero degrees, extrapolate data

        % trial 1
        if min(gonio1) > 0.0
            % select first 2 degrees of existing gonio data
            angle_end = min(gonio1) + 5; %VAR extrapolate using first 5 degrees of the gonio movement
            loc_angle_end = find(gonio1 >= angle_end,1,'first');
            % calculate linear coeffisients
            gonio_modified = [ones(length(gonio1(1:loc_angle_end)),1) gonio1(1:loc_angle_end)];
            coeffs_torque = gonio_modified\torque1(1:loc_angle_end); % NB backslash --> slope
            coeffs_displ = gonio_modified\displ1(1:loc_angle_end); % NB backslash --> slope
            coeffs_emg_gm = gonio_modified\emg_gm1(1:loc_angle_end);
            coeffs_emg_gl = gonio_modified\emg_gl1(1:loc_angle_end);
            coeffs_emg_sol = gonio_modified\emg_sol1(1:loc_angle_end);
            % enlarge arrays by adding values at zero angle to the front
            angle_half = min(gonio1)/2;
            gonio1 = [0; angle_half; gonio1];
            torque1 = [coeffs_torque(1); (angle_half*coeffs_torque(2)) + coeffs_torque(1); torque1];
            displ1 = [coeffs_displ(1); (angle_half*coeffs_displ(2)) + coeffs_displ(1); displ1];
            emg_gm1 = [coeffs_emg_gm(1); (angle_half*coeffs_emg_gm(2)) + coeffs_emg_gm(1); emg_gm1];
            emg_gl1 = [coeffs_emg_gl(1); (angle_half*coeffs_emg_gl(2)) + coeffs_emg_gl(1); emg_gl1];
            emg_sol1 = [coeffs_emg_sol(1); (angle_half*coeffs_emg_sol(2)) + coeffs_emg_sol(1); emg_sol1];
        end
        
        
        
        %%% create array with common angles (.05 intervals) for resampling

        % select the smallest range of the gonio data as basis for angle array
        common_angle_start = min(gonio1);
        common_angle_stop = max(gonio1);
        % round to .05 
        common_angle_start_rounded = round(common_angle_start/0.05)*0.05;
        % create array of angles
        average_angle_array = common_angle_start_rounded:0.05:round(common_angle_stop/0.05)*0.05;

        % chop off ends if rounded away from max values
        if average_angle_array(1) < common_angle_start
            average_angle_array(1) = [];
        end
        if average_angle_array(length(average_angle_array)) > common_angle_stop
            average_angle_array(length(average_angle_array)) = [];
        end
        


        %%% reshape, resample, average

        % delete potential duplicate values from gonio arrays
        % http://www.mathworks.com/matlabcentral/fileexchange/8354-consolidator
        [~,torque1unique,~] = consolidator(gonio1,torque1,'mean',0);
        [~,displ1unique,~] = consolidator(gonio1,displ1,'mean',0);
        [~,emg_sol1unique,~] = consolidator(gonio1,emg_sol1,'mean',0);
        [~,emg_gm1unique,~] = consolidator(gonio1,emg_gm1,'mean',0);
        [gonio1unique,emg_gl1unique,~] = consolidator(gonio1,emg_gl1,'mean',0);

        % reshape and average TORQUE across common angle array
        average_torque_gonio = spline(gonio1unique, torque1unique, average_angle_array);

        % reshape and average DISPLACEMENT across common angle array
        average_displ_gonio = spline(gonio1unique, displ1unique, average_angle_array);

        % reshape and average EMG across common angle array
        average_emg_gm_gonio = spline(gonio1unique, emg_gm1unique, average_angle_array);  % do EMG_gm1 contents need to be UNIQUE, as is done for torque and displacement?
        average_emg_gl_gonio = spline(gonio1unique, emg_gl1unique, average_angle_array);  % do EMG_gl1 contents need to be UNIQUE, as is done for torque and displacement?
        average_emg_sol_gonio = spline(gonio1unique, emg_sol1unique, average_angle_array);  % do EMG_sol1 contents need to be UNIQUE, as is done for torque and displacement?


        
    else
        %%% curve fitting gonio data - GONIO 1
        
        % method 1: fit 4th order polynomial to averaged gonio-angle curve
        fit_gonio1 = polyfit(angle1, gonio1, 4);
        
        % create array with common angles (.05 intervals)
        average_angle_array = (0:0.05:max(angle1))';
        
        % create array across angles
        gonio_new1(1:length(average_angle_array),1) = zeros;
        for ang = 1:length(average_angle_array)
            gonio_new1(ang) = (fit_gonio1(1) * average_angle_array(ang)^4) + (fit_gonio1(2) * average_angle_array(ang)^3) + (fit_gonio1(3) * average_angle_array(ang)^2) + (fit_gonio1(4) * average_angle_array(ang)) + fit_gonio1(5);
        end
          
    
        
        
        
        %%% MMM GOON
        
        
        
        
        
    end

end



output_array = rot90([average_emg_sol_gonio; average_emg_gl_gonio; average_emg_gm_gonio; average_displ_gonio; average_angle_array; average_torque_gonio],3);
% rotate 270 degrees --> torque, gonio, displ, emgx3

end
