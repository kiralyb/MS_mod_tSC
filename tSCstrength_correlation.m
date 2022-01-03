function [rho,p] = tSCstrength_correlation(x,strength,xlims)
%TSCSTRENGTH_CORRELATION Correlation of a property with tSC strength
%   [RHO,P] = TSCSTRENGTH_CORRELATION(X,STRENGTH) Testing the correlation of
%   a theta cycle property (X) with the strength of different tSCs. The 
%   propetry is plotted against tSC strength and the spearman correlation 
%   is examined. The spearman correlation coefficient (RHO) and the P-value
%   for testing the hypothesis of no correlation against the alternative
%   hypothesis of a nonzero correlation is returned.
%
%   Required input arguments:
%       X: Theta cycle property
%       STRENGTH: tSC strengths
%       XLIMS: X-axis limits of the scatter plots

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   31-Oct-2021

scatter(x,strength,'.','k','MarkerEdgeAlpha',0.05)
[rho,p] = corr(double(x),strength,'Type','Spearman');
coef = polyfit(double(x), strength, 1);
refline(coef(1), coef(2));
setmyplot_balazs;
set(gca,'ylim',[-0.05,0.10])
ylabel('tSC strength')
set(gca,'xlim',xlims)
ylim = get(gca,'ylim');
text(xlims(1),ylim(2) * 0.95,sprintf('r=%f',rho))
