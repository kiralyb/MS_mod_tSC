function tSC_coupling_epprops_comp
%Compar differnt electrophysiological properties of tSC coupled and not
%coupled neurons from preprocessed porperty data.
% Properties: firing rate, theta coupling, intra-burst spike
% number, intra-burst frequency, burst-skip ratio, Inter spike interval histgorams,
% speed correlation, 

tsc_coupled = find(sum(significance_matrix(1:4,:))>0);
non_tsc_coupled = find(significance_matrix(5,:)==1);

figure
subplot(1,5,1)
tsc_box_wilcoxon(frate,non_tsc_coupled,tsc_coupled,7)
ylabel('firing rate')
subplot(1,5,2)
tsc_box_wilcoxon(hmvl_theta,non_tsc_coupled,tsc_coupled,7)
ylabel('theta phase coupling strength')
subplot(1,5,3)
tsc_box_wilcoxon(burstlength,non_tsc_coupled,tsc_coupled,7)
ylabel('intra-burst spike number')
subplot(1,5,4)
tsc_box_wilcoxon(burstrate,non_tsc_coupled,tsc_coupled,7)
ylabel('intra-burst frequency')
subplot(1,5,5)
tsc_box_wilcoxon(skipratio,non_tsc_coupled,tsc_coupled,7)
ylabel('burst-skip ratio')

% Inter spike intervall histgorams
figure
isihist =  vertcat(Matrix(tsc_coupled).isihist);
errorshade(logspace(0,4,100),mean(isihist),std(isihist)/sqrt(n),'LineColor',[1,0,0],'ShadeColor',[1,0,0]);
hold on
isihist =  vertcat(Matrix(non_tsc_coupled).isihist);
errorshade(logspace(0,4,100),mean(isihist),std(isihist)/sqrt(n),'LineColor',Colors(end,:),'ShadeColor',Colors(end,:));
setmyplot_balazs
set(gca,'xscale','log')
ylabel('Average probability')
xlabel('Inter spike interval (ms)')

% psuedo continous spike train power
figure
pow =  vertcat(Matrix(tsc_coupled).smooth_pow);
errorshade(1:size(pow,2),squeeze(nanmean(pow(:,:))),squeeze(nanstd(pow(:,:),0,1)),'LineColor',[1,0,0],'ShadeColor',[1,0,0]);
hold on
pow =  vertcat(Matrix(non_tsc_coupled).smooth_pow);
errorshade(1:size(pow,2),squeeze(nanmean(pow(:,:))),squeeze(nanstd(pow(:,:),0,1)),'LineColor',Colors(end,:),'ShadeColor',Colors(end,:));
set(gca,'xtick',[4,12,22,35,54,80,169])
setmyplot_balazs
ylabel('Average normalized power')
xlabel('Frequency (Hz)')

% speed correlation
p_szig = 0.01;
frspeed_rho = vertcat(Matrix(tsc_coupled).frspeed_rho);
frspeed_p = vertcat(Matrix(tsc_coupled).frspeed_p);
figure
subplot(1,2,1)
pie([sum(frspeed_p<p_szig & frspeed_rho>0),sum(frspeed_p<p_szig & frspeed_rho<0),sum(frspeed_p>p_szig)])
ap = [sum(frspeed_p<p_szig & frspeed_rho>0),sum(frspeed_p<p_szig & frspeed_rho<0),sum(frspeed_p>p_szig)];

subplot(1,2,2)
hist(frspeed_rho(frspeed_p<p_szig))
frspeed_rho = vertcat(Matrix(non_tsc_coupled).frspeed_rho);
frspeed_p = vertcat(Matrix(non_tsc_coupled).frspeed_p);
figure
subplot(1,2,1)
pie([sum(frspeed_p<p_szig & frspeed_rho>0),sum(frspeed_p<p_szig & frspeed_rho<0),sum(frspeed_p>p_szig)])
ad = [sum(frspeed_p<p_szig & frspeed_rho>0),sum(frspeed_p<p_szig & frspeed_rho<0),sum(frspeed_p>p_szig)];

[h,p]=prop_test(ap,ad)

subplot(1,2,2)
hist(frspeed_rho(frspeed_p<p_szig))
setmyplot_balazs

% theta coupling
figure
subplot(1,2,1)
pie([sum(pRayleigh_theta(tsc_coupled,7) < 0.01),sum(pRayleigh_theta(tsc_coupled,7) > 0.01)])
ap = [sum(pRayleigh_theta(tsc_coupled,7) < 0.01),sum(pRayleigh_theta(tsc_coupled,7) > 0.01)];
subplot(1,2,2)
pie([sum(pRayleigh_theta(non_tsc_coupled,7) < 0.01),sum(pRayleigh_theta(non_tsc_coupled,7) > 0.01)])
ad = [sum(pRayleigh_theta(non_tsc_coupled,7) < 0.01),sum(pRayleigh_theta(non_tsc_coupled,7) > 0.01)];
[h,p]=prop_test(ap,ad)

figure;
subplot(1,2,1)
rose(hang_theta(sum(significance_matrix(1:4,:))>0 & pRayleigh_theta(:,7)' < 0.01,7),12)
[mu,kappa,Value,p,Rsquare]=b_watson(hang_theta(sum(significance_matrix(1:4,:))>0 & pRayleigh_theta(:,7)' < 0.01,7));
subplot(1,2,2); rose(hang_theta(significance_matrix(5,:)>0 & pRayleigh_theta(:,7)' < 0.01,7),12)
[mu,kappa2,Value,p,Rsquare]=b_watson(hang_theta(significance_matrix(5,:)>0 & pRayleigh_theta(:,7)' < 0.01,7));
[u2 p] = b_watsontwo(hang_theta(sum(significance_matrix(1:4,:))>0 & pRayleigh_theta(:,7)' < 0.01,7),hang_theta(significance_matrix(5,:)>0 & pRayleigh_theta(:,7)' < 0.01,7))

[param2,err] = watsonkfit(hang_theta(sum(significance_matrix(1:4,:))>0 & pRayleigh_theta(:,7)' < 0.01,7),2,1)
% figure
% pie([sum(SWR_fr(tsc_coupled)>All_fr(tsc_coupled))/length(tsc_coupled),sum(SWR_fr(tsc_coupled)<All_fr(tsc_coupled))/length(tsc_coupled)])
% pie([sum(SWR_fr(non_tsc_coupled)>All_fr(non_tsc_coupled))/length(non_tsc_coupled),sum(SWR_fr(non_tsc_coupled)<All_fr(non_tsc_coupled))/length(non_tsc_coupled)])

% speed correlation strength
tsc_coupled = find(sum(significance_matrix(1:4,[1:20,46:end]))>0);
non_tsc_coupled = find(significance_matrix(5,[1:20,46:end])==1);
figure
tsc_box_wilcoxon(abs(frspeed_rho),non_tsc_coupled,tsc_coupled,1)
ylabel('speed correlation coefficient')
