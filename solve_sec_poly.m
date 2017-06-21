%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve_sec_poly
% Marie Moltubakk 21.6.2017
% Input second order polynomial equation and desired Y values
% Output array of corresponding X values
%%%%%%%%%%%%%%%%%%%%%%%%%%

function x_array = solve_sec_poly(a, b, c, y_start, y_end, y_intervals, x_guess)
% a, b, c = coeffficients of 2nd order equation
% y_start, end, intervals = e.g. from 0 N to 2500 N, with 100 N intervals
% x_guess = a rough estimate of the range of x values (e.g. 4 mm), to select the appropriate root of the 2nd order equation

    f = fittype('a*x^2+b*x+c');
    equation = cfit(f, a, b, c);

    y_array = (y_start:y_intervals:y_end)';

    x_array(length(y_array),1) = NaN;

    for i = 1:length(y_array)
        equation.c = c - y_array(i);
        x_array(i) = fzero(equation,x_guess);
    end
    
%     % create plot if desired:
%     equation.c = c;
%     figure,plot(x_array,y_array,'.')
%     hold on
%     plot(equation)
%     legend('output X/Y', 'equation','location','Northwest')
end