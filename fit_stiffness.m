function [fitresult, gof] = fit_stiffness(nmz_disp, nmz_force, elong_cut_index, nmz_disp_mtj_mean, nmz_disp_otj_mean)
global subject_id

%CREATEFIT(NMZ_DISP,NMZ_FORCE)
%  Create a fit.
%
%  Data for 'Stiffness fit' fit:
%      X Input : nmz_disp
%      Y Output: nmz_force
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 25-Jun-2013 09:18:40


%% Fit: 'Stiffness fit'.
[xData, yData] = prepareCurveData( nmz_disp(1:elong_cut_index), nmz_force(1:elong_cut_index) );

% Set up fittype and options.
ft = fittype( 'poly2' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf 0];
opts.Upper = [Inf Inf 0];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
coeffvals = coeffvalues(fitresult);

% Plot fit with data.
plottitle = horzcat('FIT PLOT, tendon stiffness ', subject_id);
fignavn = figure( 'Name', plottitle);
h = plot( fitresult, xData, yData );
hold on
plot(nmz_disp_mtj_mean,nmz_force,'g.','DisplayName','MTJ displacement')
plot(nmz_disp_otj_mean,nmz_force,'g.','DisplayName','OTJ displacement')
plot(nmz_disp(elong_cut_index+1:end),nmz_force(elong_cut_index+1:end),'y.','DisplayName','Elongation in cut off range')
legend( h, 'Force-elongation', 'Stiffness fit', 'Location', 'SouthEast' );
% Label axes
xlabel( 'Tendon elongation (mm)' );
ylabel( 'Tendon force (N)' );
title(plottitle);
text(3.2, 400, horzcat('Y = ', num2str(coeffvals(1)), 'x^2 + ', num2str(coeffvals(2)), 'x + 0'), 'Color', 'r')

saveas(fignavn, strcat('data_output/stiff_fit_', subject_id), 'png')
