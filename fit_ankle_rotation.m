function [fitresult, gof] = fit_ankle_rotation(angle2, displ2, phasename)
global plot_norm subject_id

%CREATEFIT(DISPL2,ANGLE2)
%  Create a fit.
%
%  Data for 'untitled fit 4' fit:
%      X Input : displ2
%      Y Output: angle2
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 17-Jun-2013 13:00:05


%% Fit: 'untitled fit 4'.
[xData, yData] = prepareCurveData( angle2, displ2 );

% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( ft );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

if plot_norm
    % Plot fit with data.
    plottitle = horzcat('FIT PLOT for ankle rotation correction for ', subject_id);
    fignavn = figure( 'Name', plottitle);
    h = plot( fitresult, xData, yData );
    legend( h, horzcat('angle vs. displ ', phasename), horzcat('Linear fit: displ/deg = ', num2str(coeffvalues(fitresult))), 'Location', 'NorthEast' );
    % Label axes
    xlabel( 'Goniometer ankle angle (deg)' );
    ylabel( 'Calcaneus displacement (mm)' );
    title(plottitle);
    grid on
    saveas(fignavn, strcat('data_plots_stiff/IND_anklerot_fit_', subject_id, phasename), 'png')
end