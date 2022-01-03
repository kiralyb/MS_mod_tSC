function [p] = tsc_box_anova(y,y_norm,filt)
%TSC_BOX_ANOVA Box plot and ANOVA test for tSC 
%   [P] = TSC_BOX_ANOVA(Y,Y_NORM,FILT) Performs repeated measures ANOVA
%   test followed by Tukey's post hoc test for multiple comparisons between
%   theta cycles with different tSC on measure Y and creates boxplots. The
%   p-values of Tukey's post hoc test for each comparison is returned.
%
%   Required input arguments:
%       Y: Data for theta cycles with different tSC or without tSCs.
%       Y_NORM: normalized value of Y.
%       FILT: filters on which values of Y should be used.

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   31-Oct-2021
global Colors

% Perform repeated meassures ANOVA test followed by Tukey's post hoc test
n = size(y_norm,2);
[~,~,stats] = rmanova(y(filt,1:n));
[results,~] = multcompare(stats,'Display','off');
p = results(:,[1,2,6]);

% plotting
boxplot(y_norm(filt,:),'Colors',Colors,'symbol','','Widths',0.5) 
h=findobj(gca,'tag','Outliers');
delete(h);
ylim auto
if ~isequal(y,y_norm)
    hold on
    line([0,n + 1],[0,0],'Color','red','LineWidth',1,'LineStyle','--')
end
set(gca,'xtick',1:n)
try
    set(gca,'xticklabel',({'tSC1','tSC2','tSC3','tSC4','tSC5','no'}))
catch    
end
set(gca,'XTickLabelRotation',45)
setmyplot_balazs;


function [F_exp,p,stats] = rmanova(Y)

% delete subjects with missing data
Y(sum(isnan(Y),2)>0,:)=[];

% Calculate means
m = mean2(Y); 
m_experiment = mean(Y,1);
m_subject = mean(Y,2); 

% calculate the mean squared differences
Subject_ssqrd = 0;
Experiment_ssqrd = 0;
ssqrd_ER = 0;
[n_subjects,n_exp] = size(Y);
for subject = 1:n_subjects
    for experiment = 1:n_exp
        Subject_ssqrd = Subject_ssqrd + (m_subject(subject) - m)^2;
        Experiment_ssqrd = Experiment_ssqrd + (m_experiment(experiment) - m)^2;
        ssqrd_ER  = ssqrd_ER + (Y(subject,experiment) - m_experiment(experiment) - m_subject(subject) + m)^2;
    end
end
MS_subject = Subject_ssqrd / (n_subjects - 1);
MS_experiment = Experiment_ssqrd  / (n_exp - 1);
MS_ER = ssqrd_ER / ((n_subjects - 1) * (n_exp - 1));

% calculate the F-statistics
F_exp = MS_experiment / MS_ER;
F_subject = MS_subject / MS_ER;

% calculate p-values
p_experiment  = 1 - fcdf(F_exp, n_exp - 1, (n_subjects - 1) * (n_exp - 1));
p_subject = 1 - fcdf(F_subject, n_subjects - 1, (n_subjects - 1) * (n_exp - 1));
p = [p_experiment, p_subject];

% fill stats output 
stats.gnames = num2str((1:n_exp)');
stats.n = n_subjects * ones(1,n_exp);
stats.source = 'anova1';
stats.means = m_experiment';
stats.df = n_subjects * n_exp - n_subjects - n_exp + 1;
stats.s = sqrt(MS_ER);