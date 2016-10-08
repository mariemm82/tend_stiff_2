%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_us_licht
% Marie Moltubakk & R & M 25.8.2015
% Read US data file from LICHTWARK matlab script
% Produce array containing fascicle length and pennation angle from Lichtwark analyses:
%   usdata_licht = 
%   1 = time
%   2 = GM fascicle length
%   3 = GM pennation angle
%   4 = time (dupicate)     - IF existing
%   5 = SOL fascicle length - IF existing
%   6 = SOL pennation angle - IF existing
%%%%%%%%%%%%%%%%%%%%%%%%%%


function [usdata_licht] = read_us_licht(usfile, usframe, trial_name)
    
    % import us data
    usdata = importdata(usfile, '\t', 1);
    the_size = size(usdata.data);
    number_of_columns = the_size(2);
    
    
    % cut off data until triggerframe
    % NB: 
    %   US video has 395 frames.
    %   Tracker labels as 0 to 394. Frame 0 = time 0.00, frame 394 = time 21.773
    %   Lichtwark reports 395 data points. First point = time 0.0553, last point = time 21.828
    %   Conclusion: Triggerframe = 0 means no data cutoff, start using line 1. But Lichtwark timestamps should be lowered by 1 unit.
    %   Norm data start recording at triggerframe, after that, higher sampling frequency and longer duration of recording than the US video
    usdata_offset = usdata.data(usframe+1:end,:);

    % check number and type of fascicles tracked, average duplicates
    if number_of_columns == 3 % contains 1 GM trial
        usdata_offset(:,1) = usdata_offset(:,1) - usdata_offset(1,1); % correct lichtwark time error (set first frame = 0.00)
        % keep all 3 columns as as
        usdata_licht = usdata_offset;
        
    elseif number_of_columns == 5  % contains 2 GM trials
        usdata_offset(:,1) = usdata_offset(:,1) - usdata_offset(1,1); % correct lichtwark time error (set first frame = 0.00)
        % average 2 GM
        usdata_licht(:,1) = usdata_offset(:,1);
        usdata_licht(:,2) = mean([usdata_offset(:,2) usdata_offset(:,4)], 2);
        usdata_licht(:,3) = mean([usdata_offset(:,3) usdata_offset(:,5)], 2);
        
    elseif number_of_columns == 6 % contains 1 GM trial + 1 SOL trial
        usdata_offset(:,1) = usdata_offset(:,1) - usdata_offset(1,1); % correct lichtwark time error (set first frame = 0.00)
        usdata_offset(:,4) = usdata_offset(:,4) - usdata_offset(1,4); % correct lichtwark time error (set first frame = 0.00)
        % keep all 6 columns as is: 1 GM + 1 SOL
        usdata_licht = usdata_offset;
        
    elseif number_of_columns == 8  % contains 2 of one muscle, 1 of the other
        % check data headers to decide which columns to average
        str = usdata.textdata{1};
        mode = str2num(str(end-1));
        if mode == 1 % Average 2 GM facicles + 1 sol
            usdata_offset(:,1) = usdata_offset(:,1) - usdata_offset(1,1); % correct lichtwark time error (set first frame = 0.00)
            usdata_offset(:,6) = usdata_offset(:,6) - usdata_offset(1,6); % correct lichtwark time error (set first frame = 0.00)
            usdata_licht(:,1) = usdata_offset(:,1);
            usdata_licht(:,2) = mean([usdata_offset(:,2) usdata_offset(:,4)], 2);
            usdata_licht(:,3) = mean([usdata_offset(:,3) usdata_offset(:,5)], 2);
            usdata_licht(:,4) = usdata_offset(:,6);
            usdata_licht(:,5) = usdata_offset(:,7);
            usdata_licht(:,6) = usdata_offset(:,8);
        else % 1 GM + average 2 sol
            usdata_offset(:,1) = usdata_offset(:,1) - usdata_offset(1,1); % correct lichtwark time error (set first frame = 0.00)
            usdata_offset(:,4) = usdata_offset(:,4) - usdata_offset(1,4); % correct lichtwark time error (set first frame = 0.00)
            usdata_licht(:,1) = usdata_offset(:,1);
            usdata_licht(:,2) = usdata_offset(:,2);
            usdata_licht(:,3) = usdata_offset(:,3);
            usdata_licht(:,4) = usdata_offset(:,4);
            usdata_licht(:,5) = mean([usdata_offset(:,5) usdata_offset(:,7)], 2);
            usdata_licht(:,6) = mean([usdata_offset(:,6) usdata_offset(:,8)], 2);
        end
        
    elseif number_of_columns == 10  % contains 2 GM trials + 2 SOL trials
        usdata_offset(:,1) = usdata_offset(:,1) - usdata_offset(1,1); % correct lichtwark time error (set first frame = 0.00)
        usdata_offset(:,7) = usdata_offset(:,7) - usdata_offset(1,7); % correct lichtwark time error (set first frame = 0.00)
        % average 2 GM & 2 SOL columns
        usdata_licht(:,1) = usdata_offset(:,1);
        usdata_licht(:,2) = mean([usdata_offset(:,2) usdata_offset(:,4)], 2);
        usdata_licht(:,3) = mean([usdata_offset(:,3) usdata_offset(:,5)], 2);
        usdata_licht(:,4) = usdata_offset(:,6);
        usdata_licht(:,5) = mean([usdata_offset(:,7) usdata_offset(:,9)], 2);
        usdata_licht(:,6) = mean([usdata_offset(:,8) usdata_offset(:,10)], 2);
    end
    
    
    % convert rad to deg
    usdata_licht(:,3) = usdata_licht(:,3) * 180 / pi;
    if(length(usdata_licht(1,:))) == 6
        usdata_licht(:,6) = usdata_licht(:,6) * 180 / pi;
    end
    
    % print report for Lichtwark data
    cprintf('blue', horzcat(trial_name, ': GM fascicle length ', num2str(round(min(usdata_licht(:,2)),1)), '-', num2str(round(max(usdata_licht(:,2)),1)), ' mm, pennation = ', num2str(round(max(usdata_licht(:,3)),1)), '-', num2str(round(min(usdata_licht(:,3)),1)), ' deg. (NB: not exact zero angle)\n'));

    
    

end



