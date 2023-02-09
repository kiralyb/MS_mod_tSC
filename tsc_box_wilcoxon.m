function p = tsc_box_wilcoxon(prop,non_tsc_coupled,tsc_coupled,propdim)

boxplot([prop(non_tsc_coupled,propdim);prop(tsc_coupled,propdim)],[zeros(1,length(non_tsc_coupled)),ones(1,length(tsc_coupled))]')
h=findobj(gca,'tag','Outliers');
delete(h);
ylim auto
set(gca,'xticklabel',{'not tSC coupled','tSC coupled'});
p = ranksum(prop(non_tsc_coupled,propdim),prop(tsc_coupled,propdim));
sigstar([1,2],p)
setmyplot_balazs
set(gca,'XTickLabelRotation',45)

end