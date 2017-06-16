function [fitresult, gof] = fit_stiffness(elongation, force, loc_elong_cut, displ_MTJ, displ_OTJ)
global subject_id

%CREATEFIT(NMZ_DISP,NMZ_FORCE)
%      X Input : nmz_disp
%      Y Output: nmz_force
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.


%% Fit: 'Stiffness fit'.
[xData, yData] = prepareCurveData( elongation(1:loc_elong_cut), force(1:loc_elong_cut) );

% Set up fittype and options.
ft = fittype( 'poly2' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf -Inf]; % last variable 0 <- require fit through origo
opts.Upper = [Inf Inf Inf]; % last variable 0

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
coeffvals = coeffvalues(fitresult);

% Plot fit with data.
plottitle = horzcat('FIT PLOT, tendon stiffness ', subject_id);
fignavn = figure('Name', plottitle);
hold on
h = plot(fitresult, xData, yData);
plot(displ_MTJ,force,'g.','DisplayName','MTJ displacement')
plot(displ_OTJ,force,'g.','DisplayName','OTJ displacement')
plot(elongation(loc_elong_cut+1:end),force(loc_elong_cut+1:end),'y.','DisplayName','Elongation in cut off range')
legend( h, 'Force-elongation', 'Stiffness fit', 'Location', 'SouthEast' );
% Label axes
xlabel( 'Tendon elongation (mm)' );
ylabel( 'Tendon force (N)' );
title(plottitle);
text(3.2, 600, horzcat('Y = ', num2str(coeffvals(1)), 'x^2 + ', num2str(coeffvals(2)), 'x + ', num2str(coeffvals(3))), 'Color', 'r')

saveas(fignavn, strcat('data_plots_stiff/IND_stiff_fit_', subject_id), 'png')
