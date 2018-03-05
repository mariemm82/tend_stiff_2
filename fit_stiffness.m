function [fitresult, gof] = fit_stiffness(elongation, force, loc_cut, displ_MTJ, displ_OTJ)
global subject_id plot_check plot_achilles
% 0 = use fit through zero
% 1 = use fit with free beginning
choice_of_fit = 0; %VAR



%CREATEFIT(NMZ_DISP,NMZ_FORCE)
%      X Input : nmz_disp
%      Y Output: nmz_force
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.


%% Fit: 'Stiffness fit'.
[xData, yData] = prepareCurveData( elongation(1:loc_cut), force(1:loc_cut) );

% Set up fittype and options.
ft = fittype( 'poly2' );
opts = fitoptions( ft );



%% version 1, through zero
opts.Lower = [-Inf -Inf 0]; % last variable 0 <- require fit through origo
opts.Upper = [Inf Inf 0]; % last variable 0

% Fit model to data.
[fitresult0, gof0] = fit( xData, yData, ft, opts );
coeffvals = coeffvalues(fitresult0);

if plot_check %TMP plot_achilles
% Plot fit with data.
    plottitle = horzcat('FIT PLOT, tendon stiffness through ZERO, ', subject_id);
    fignavn = figure('Name', plottitle);
    hold on
    h = plot(fitresult0, xData, yData);
    plot(displ_MTJ,force,'g.','DisplayName','MTJ displacement')
    plot(displ_OTJ,force,'g.','DisplayName','OTJ displacement')
    plot(elongation(loc_cut+1:end),force(loc_cut+1:end),'y.','DisplayName','Elongation in cut off range')
    legend( h, 'Force-elongation', 'Stiffness fit', 'Location', 'SouthEast' );
    % Label axes
    axis([-1 Inf 0 Inf])
    xlabel( 'Tendon elongation (mm)' );
    ylabel( 'Tendon force (N)' );
    title(plottitle,'Interpreter', 'none');
    loc_text_x = max(displ_MTJ)*0.6;
    loc_text_y = force(end)*0.25;
    text(loc_text_x, loc_text_y, horzcat('Y = ', num2str(coeffvals(1)), 'x^2 + ', num2str(coeffvals(2)), 'x + ', num2str(coeffvals(3))), 'Color', 'r')
    
    saveas(fignavn, strcat('data_plots_stiff/IND_stiff_FIT_', subject_id, '_zero'), 'png')
end

%% version 2, NOT through zero
opts.Lower = [-Inf -Inf -Inf]; % last variable 0 <- require fit through origo
opts.Upper = [Inf Inf Inf]; % last variable 0

% Fit model to data.
[fitresult1, gof1] = fit( xData, yData, ft, opts );
coeffvals = coeffvalues(fitresult1);

if plot_achilles
% Plot fit with data.
    plottitle = horzcat('FIT PLOT, tendon stiffness FREE onset, ', subject_id);
    fignavn = figure('Name', plottitle);
    hold on
    h = plot(fitresult1, xData, yData);
    plot(displ_MTJ,force,'g.','DisplayName','MTJ displacement')
    plot(displ_OTJ,force,'g.','DisplayName','OTJ displacement')
    plot(elongation(loc_cut+1:end),force(loc_cut+1:end),'y.','DisplayName','Elongation in cut off range')
    legend( h, 'Force-elongation', 'Stiffness fit', 'Location', 'SouthEast' );
    % Label axes
    axis([-1 Inf 0 Inf])
    xlabel( 'Tendon elongation (mm)' );
    ylabel( 'Tendon force (N)' );
    title(plottitle,'Interpreter', 'none');
    loc_text_x = max(displ_MTJ)*0.6;
    loc_text_y = force(end)*0.25;
    text(loc_text_x, loc_text_y, horzcat('Y = ', num2str(coeffvals(1)), 'x^2 + ', num2str(coeffvals(2)), 'x + ', num2str(coeffvals(3))), 'Color', 'r')

    saveas(fignavn, strcat('data_plots_stiff/IND_stiff_FIT_', subject_id, '_free'), 'png')
    close
end


%% final selection of data
if choice_of_fit == 1
    fitresult = fitresult1;
    gof = gof1;
else
    fitresult = fitresult0;
    gof = gof0;
end