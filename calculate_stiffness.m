%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_stiffness
% Marie Moltubakk 26.6.2013
% Read stiffness equation and max force
% Produce stiffness for a range, i.e. 80-100%, of ind max force
%%%%%%%%%%%%%%%%%%%%%%%%%%

function stiffness = calculate_stiffness(stiff_eq, force100, percent_start, percent_stop)
    
    force_start = force100 * percent_start;
    force_stop = force100 * percent_stop;
    
    % displacement at start % of individual max range
    objective = @(x) stiff_eq(x) - force_start;
    displ_start = fzero(objective, 6);
    
    % displacement at stop % of individual max range
    objective = @(x) stiff_eq(x) - force_stop;
    displ_stop = fzero(objective, 6);
    
    stiffness = (force_stop-force_start)/(displ_stop-displ_start);
    
    cprintf('blue', horzcat(' Stiffness ', num2str(percent_start*100), '-', num2str(percent_stop*100), ' ind = ', num2str(stiffness,4), ' N/mm.\n'))
end