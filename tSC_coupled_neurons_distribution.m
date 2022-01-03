function tSC_coupled_neurons_distribution(significance_matrix,tsc_num,pRayleigh_theta,rythmgroupnum)
%TSC_COUPLED_NEURONS_DISTRIBUTION tSC copuling distibution plots 
%   TSC_COUPLED_NEURONS_DISTRIBUTION(SIGNIFICANCE_MATRIX,TSC_NUM,PRAYLEIGH_THETA,RYTHMGROUPNUM)
%   plots 4 tSC coupling distribution plots. First, the number of
%   neurons coupled to each tSC. Second, the histogram of the
%   number of tSCs the neurons are coupled to. Third, the number and ratio
%   of neurons coupled to tSCs in different rytmicity groups. Fourth, the 
%   number and ratio of tSC coupled neurons among theta coupled and not
%   coupled neurons. 
%
%   Required input arguments:
%       SIGNIFICANCE_MATRIX: logic matrix showing whether each neuron is
%       coupled to different tSCs.
%       TSC_NUM: number of tSCs plotted
%       PRAYLEIGH_THETA: results of Rayleigh's test for theta coupling
%       RYTHMGROUPNUM: array showing the rythmicity group of the neurons
%
%   See also  TSC_NEURON_FIRINGPROP and TSC_COUPLING_TEST

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   31-Oct-2021
global Colors

% Plot 'number of tSC coupled to' histogram
figure
subplot(1,tsc_num,1)
tsc_coup = sum(significance_matrix(1:tsc_num,:))';
tsc_coup(tsc_coup == 0 & significance_matrix(end,:)' == 0) = NaN;
histogram(tsc_coup,-0.5:1:tsc_num + 0.5,'FaceColor','w','EdgeColor','k','EdgeAlpha',1);
set(gca,'xtick',0:tsc_num);
setmyplot_balazs;
xlim([-1,tsc_num + 1])
xlabel('Number of tSCs coupled to')

% Plot 'number of tSC coupled neurons' barplot
subplot(1,tsc_num,2)
for tSC = 1:tsc_num + 1
    if tSC == tsc_num + 1
        inx = find(sum(significance_matrix(1:tsc_num,:)) == 0);
        x_pos = tSC + 1;
    else
        inx = find(significance_matrix(tSC,:) == 1);
        x_pos = tSC;
    end
    hold on
    B = bar(x_pos,length(inx));
    set(B,'FaceColor',Colors(x_pos,:));
    set(gca,'xtick',[1:tsc_num,tsc_num + 2]);
    set(gca,'xticklabel',{'tSC1','tSC2','tSC3','tSC4','no'});
    xlim([0,tsc_num + 3])
    ylabel('Number of neurons')
    setmyplot_balazs;
    set(gca,'XTickLabelRotation',45)
end

% Plot the number and ratio of tSC coupled neurons for different rythmicity
% groups.
subplot(1,tsc_num,3)
ind = tsc_coup > 0;
rythmgroupnum_filt = NaN(1,length(rythmgroupnum));
rythmgroupnum_filt(rythmgroupnum == 1 & ~isnan(tsc_coup)) = 1;
rythmgroupnum_filt(rythmgroupnum == 13 & ~isnan(tsc_coup)) = 2;
rythmgroupnum_filt(rythmgroupnum == 2 & ~isnan(tsc_coup)) = 3;
rythmgroupnum_filt(rythmgroupnum == 17 & ~isnan(tsc_coup)) = 4;
histogram(rythmgroupnum_filt,'FaceColor','k','FaceAlpha',1);
hold on
histogram(rythmgroupnum_filt(ind),'FaceColor','w','FaceAlpha',1,'EdgeColor','k','EdgeAlpha',1);
xlim([0,5]);
set(gca,'xtick',1:4);
set(gca,'xticklabel',{'pacemaker','theta rythmic','tonic','not rythmic'});
setmyplot_balazs;
set(gca,'XTickLabelRotation',45)

% Plot the number and ratio of tSC coupled neurons in the function of theta
% coupling
subplot(1,tsc_num,4)
bar(1,sum(pRayleigh_theta(:,end) < 0.01 & ~isnan(tsc_coup)),'FaceColor','k');
hold on
bar(1,sum(pRayleigh_theta(:,end) < 0.01 & tsc_coup ~= 0 & ~isnan(tsc_coup)),'FaceColor','w')
bar(2,sum(pRayleigh_theta(~isnan(tsc_coup),end) > 0.01) + sum(isnan(pRayleigh_theta(~isnan(tsc_coup),end))),'FaceColor','k')
bar(2,sum(pRayleigh_theta(:,end) > 0.01 & tsc_coup ~= 0 & ~isnan(tsc_coup)),'FaceColor','w')
xlim([0,3])
set(gca,'xtick',[1,2]);
set(gca,'xticklabel',{'theta coupled','not theta coupled'});
set(gca,'XTickLabelRotation',45)
setmyplot_balazs;
