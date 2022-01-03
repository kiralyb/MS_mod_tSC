function tsc_barplot(y)
%TSC_BARPLOT Bar plots for tSC
%   TSC_BARPLOT(Y) Bar plots with means and error bars to compare theta
%   cycles with different tSCs.
%
%   Required input arguments:
%       Y: Data for theta cycles with different tSCs.

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   31-Oct-2021
global Colors
mean_Y = nanmean(y);
tSC_num = length(mean_Y);
barColors = Colors;
if tSC_num == 2
    barColors = [[1,1,1];[0,0,0]];
end

for c = 1:tSC_num
    B(c) = bar(c,mean_Y(c));
    hold on
    set(B(c),'FaceColor',barColors(c,:));
end
errorbar(1:tSC_num,nanmean(y(:,1:tSC_num)),zeros(1,tSC_num),nanstd(y(:,1:tSC_num)),'Color','k','LineStyle','none');
set(gca,'xtick',1:tSC_num)
if tSC_num == 5
    set(gca,'xticklabel',({'tSC1','tSC2','tSC3','tSC4','tSC5'}))
elseif tSC_num == 2
    set(gca,'xticklabel',({'tSC','no tSC'}))
end
set(gca,'XTickLabelRotation',45)
setmyplot_balazs;
xlim([0,tSC_num + 1])