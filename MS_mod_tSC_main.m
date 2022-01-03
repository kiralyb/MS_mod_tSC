function MS_mod_tSC_main
%MS_MOD_TSC_MAIN Main wrapper for the figures presented in the Király et
%al. manuscript. 
%   This code is for further analysis after running
%   TSC_RUN_SESSION_ANALYSES on the sessions of simultaneous septal unit
%   and hippocampal LFP recordings.

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   23-Dec-2021

%% Load awake mouse data
mainpath='L:\Balint\tsc\awake_mouse\';
load([mainpath, 'ses_Matrix.mat'],'ses_Matrix');
n_ses = length(ses_Matrix);
load([mainpath, 'Matrix.mat'],'Matrix');
n = length(Matrix);

% define global variables
global Colors Samplingrate
Colors = [[0.2,0.5,1];[0.4,0.8,0.4];[0.5,0,1];[1,0.5,0];[0.5,0.5,0.5];[0,0,0]]; % define tSC colors, and the 'no tSC' color
Samplingrate = 1000; % sampling rate in Hz
%% Figure 1 and S1 
% Supra-theta spectral components in the MS and the hippocampal CA1
% radiatum layer of freely moving mice.
%--------------------------------------------------------------------------
% Panel A
tSC_raster(mainpath,'20161695','1824','27')
xlim([502153,503153])

% Panel B
% example ISI
logISI(mainpath,'20161988', '101102', 1, 16)
% average ISI histogram
figure
isihist =  vertcat(Matrix.isihist);
errorshade(logspace(0,4,100),mean(isihist),std(isihist)/sqrt(n),'LineColor',Colors(end,:),'ShadeColor',Colors(end,:));
setmyplot_balazs
set(gca,'xscale','log')
ylabel('Average probability')
xlabel('Inter spike interval (ms)')

% Panel C
figure
pow =  vertcat(Matrix.smooth_pow);
% example spectrum 
cellID=189;
plot(1:size(pow,2),squeeze(pow(cellID,:))) 
set(gca,'xtick',[4,12,22,35,54,80,169])
setmyplot_balazs
ylabel('Normalized power')
xlabel('Frequency (Hz)')
% average spectrum
figure
errorshade(1:size(pow,2),squeeze(nanmean(pow(:,:))),squeeze(nanstd(pow(:,:),0,1))/sqrt(n),'LineColor',Colors(end,:),'ShadeColor',Colors(end,:));
set(gca,'xtick',[4,12,22,35,54,80,169])
setmyplot_balazs
ylabel('Average normalized power')
xlabel('Frequency (Hz)')

% Panel D
neuron_spectra(mainpath,'20161989','139140','3',1,7);

% Panel E-F 
% tSCs were extracted in Phyton with the 'Analysis framework for the
% extraction of theta-nested spectral components' package by
% Lopes-dos-Santos et al., 2018. Figure were created with the
% 'plotMeantSCcycle' function.

% Panel G
figure
exID = 6;
tSC_num = 5;
tSC_spect = horzcat(ses_Matrix.tSCs);
tSC_spect = reshape(tSC_spect,[size(tSC_spect,1),tSC_num,n_ses]);
tSC_spectrum_plotter(tSC_spect,exID)

% Panel H
figure
avg_cycle_freq = 1 ./ vertcat(ses_Matrix.avg_cycle_length) * Samplingrate;
tsc_box_anova(avg_cycle_freq,avg_cycle_freq,1:n_ses);
ylabel('Theta cycle length')

% Panel I
figure
subplot(2,3,1)
longtheta = vertcat(ses_Matrix.longtheta);
tsc_barplot(longtheta(:,1:tSC_num));
ylabel('Ratio of theta cycles')
title('Long theta epochs')

subplot(2,3,2)
shorttheta = vertcat(ses_Matrix.shorttheta);
tsc_barplot(shorttheta(:,1:tSC_num));
title('Short theta epochs')

subplot(2,3,3)
margin_long = vertcat(ses_Matrix.margin_long);
tsc_barplot(margin_long(:,1:tSC_num));
title('Theta transitions')

subplot(2,3,4)
tsc_barplot([longtheta(:,end - 1),1 - longtheta(:,end - 1)]);

subplot(2,3,5)
tsc_barplot([shorttheta(:,end - 1),1 - shorttheta(:,end - 1)]);

subplot(2,3,6)
tsc_barplot([margin_long(:,end - 1),1 - margin_long(:,end - 1)]);

% Panel J
figure
avg_cycle_speed = vertcat(ses_Matrix.avg_cycles_speed);
tsc_box_anova(avg_cycle_speed,avg_cycle_speed,1:size(avg_cycle_speed,2));
ylabel('Speed during theta cycles')

% Panel S1A
tSC_spectrum_plotter(tSC_spect,1:n_ses)

% Panel S1B
figure
avg_allcycle_freq = 1 ./ vertcat(ses_Matrix.avg_allcycle_length) * Samplingrate;
avg_cycle_freq_norm = bsxfun(@rdivide,(bsxfun(@minus, avg_cycle_freq, avg_allcycle_freq)),avg_allcycle_freq);
tsc_box_anova(avg_cycle_freq,avg_cycle_freq_norm,1:n_ses);
ylabel('Theta cycle length')

% Panel S1C
figure
avg_allcycle_speed = vertcat(ses_Matrix.avg_allcycle_speed);
avg_cycle_speed_norm = bsxfun(@rdivide,(bsxfun(@minus, avg_cycle_speed, avg_allcycle_speed)),avg_allcycle_speed);
tsc_box_anova(avg_cycle_speed,avg_cycle_speed_norm,1:size(avg_cycle_speed,2));
ylabel('Speed during theta cycles')

% Panel S1D and S1E
tSC_correlation_fig(ses_Matrix,tSC_num);
%% Figure 2 and S2-S4
% MS single neuron firing is correlated with hippocampal tSCs
%--------------------------------------------------------------------------

% Load firing properties of the neurons during different theta cycles
tSC_num = 5;
frate = vertcat(Matrix.frate);
burstlength = vertcat(Matrix.burstduration);
burstrate = vertcat(Matrix.burstrate);
skipratio = vertcat(Matrix.skipratio);
hang_theta = vertcat(Matrix.hang);
hmvl_theta = vertcat(Matrix.hmvl);
pRayleigh_theta = vertcat(Matrix.pRayleigh);
rythmgroupnum = vertcat(Matrix.rythmgroupnum);
drift = vertcat(Matrix.drift);

% Normalize firing properties 
frate_norm = tsc_firingprop_norm(frate);
burstlength_norm = tsc_firingprop_norm(burstlength);
burstrate_norm = tsc_firingprop_norm(burstrate);
hmvl_theta_norm = tsc_firingprop_norm(hmvl_theta);
skipratio_norm = tsc_firingprop_norm(skipratio);
phase_hist_align_norm = tSC_phasehist_aligner(Matrix,tSC_num);

% Panel A 
tSC_neuron_firingprop(mainpath,'20161695','1824','27',1,9)
tSC_neuron_firingprop(mainpath,'20161750','5054','15',1,17)
tSC_neuron_firingprop(mainpath,'20161869','1314','32',1,3)

% Panel B
filter = frate(:,end) > 3 & (drift == 1); % filter for non-drifting neurons with firing rate above 3 Hz
figure
subplot(2,1,1)
tsc_box_anova(frate,frate_norm,filter);
ylabel('Relative firing rate')
subplot(2,1,2)
scatter(frate_norm(filter,1),frate_norm(filter,4),'k','x')
hold on
scatter(frate_norm([338,309,48],1),frate_norm([342,309,48],4),'r','x')
line([-2,2],[-2,2])
xlim([-1.5,1.5])
ylim([-1.5,1.5])
xlabel('Relative firing rate during tSC1 cycles')
ylabel('Relative firing rate during tSC4 cycles')

% Panel C1 and Fig S2
types = [0,1,13,0];
titles = {'Phase coupled','Pacemaker','Theta follower','tSC1 activated'};
for type = 1:4
    % define filters for the different neuron types
    if type == 1
        filter = pRayleigh_theta(:,end) < 0.01 & frate(:,end) > 3 & (drift == 1);
    elseif type == 4
        [~,MX] = max(frate_norm,[],2);
        filter = pRayleigh_theta(:,end) < 0.01 & frate(:,end) > 3 & (drift == 1) & (MX == 1);
    else
        filter = rythmgroupnum == types(type)  & frate(:,end) > 3 & (drift == 1)  ;
    end
    % plot phase histogram curves
    figure
    for tSC=1:tSC_num + 1
        hold on
        errorshade(1:15,squeeze(mean(phase_hist_align_norm(filter,tSC,:),1))',squeeze(std(phase_hist_align_norm(filter,tSC,:),1)/sqrt(size(phase_hist_align_norm(filter,tSC,:),1)))', 'LineColor',Colors(tSC,:), 'ShadeColor',Colors(tSC,:), 'LineWidth', 3)        
    end
    setmyplot_balazs;
    title(titles(type));
    if type==1 || type==3
        ylabel('Spiking probability')
    end
    set(gca,'xtick',[1,4.5,8,11.5,15]);
    xlabel('Theta phase relative to the preffered phase')
    set(gca,'xticklabel',{'-\pi/2','0','\pi/2','\pi','3\pi/2'});
    xlim([1,15])
end

% Panel C2
filter = pRayleigh_theta(:,end) < 0.01 & (drift == 1) & frate(:,end) > 3;
figure
tsc_box_anova(hmvl_theta,hmvl_theta_norm,filter);
ylabel('Relative strength of theta phase coupling')

% Panel D
figure
Phase_diff = hang_theta(filter,1) - hang_theta(filter,4);
rose(Phase_diff,30)
title('Preferred theta phase differences tSC1-tSC4 cycles')

% Panel E
figure
tsc_box_anova(burstrate,burstrate_norm,filter);
ylabel('Relative strength of theta phase coupling')

% Fig S3 is created with the ms_sync_analysis package by Kocsis et al. 2021
% (https://github.com/hangyabalazs/ms_sync_analysis) 

% Fig S4
figure
maximize_figure
rythm_types = [0,1,13,2];
titles = {'phase coupled','pacemaker','theta rythmic','tonic'};
types_num = length(titles);
for type = 1:types_num
     % define filters for the different neuron types
    if type == 1
        filter = pRayleigh_theta(:,end) < 0.01 & (drift == 1) & frate(:,end) > 3;
    else
        filter = (rythmgroupnum == rythm_types(type)) & (drift == 1)  & frate(:,end) > 3;
    end
    % compare firing during theta cycles with different tSC on box plots
    subplot(5,types_num,0 + type);
    tsc_box_anova(frate,frate_norm,filter);
    set(gca,'xtick',[]);
    if type == 1
        ylabel('Relative firing rate')
    end
    title(titles(type));
    
    subplot(5,types_num,types_num + type);
    tsc_box_anova(hmvl_theta,hmvl_theta_norm,filter);
    set(gca,'xtick',[]);
    if type == 1
        ylabel('Relative strength of theta phase coupling')
    end
    
    if rythm_types(type) ~= 2 %skp burst properties for Tonic neurons
        
        subplot(5,types_num,types_num * 2 + type);
        tsc_box_anova(burstrate,burstrate_norm,filter);
        set(gca,'xtick',[]);
        if type == 1
            ylabel('Relative burstrate')
        end
        
        subplot(5,types_num,types_num * 3 + type);
        tsc_box_anova(burstlength,burstlength_norm,filter);
        set(gca,'xtick',[]);
        if type == 1
            ylabel('Relative burstlength')
        end
        
        subplot(5,types_num,types_num * 4 + type);
        tsc_box_anova(skipratio,skipratio_norm,filter);
        if type == 1
            ylabel('Relative burst-skip ratio')
        end
    end
end

%% Figure 3, S5 and S6
% MS neurons show phase coupling to hippocampal tSCs.
%--------------------------------------------------------------------------

% Panel A
tsc_STA(mainpath, '20161869', '45', '32', 1, 3)

% Find significcant tSC coupling
tSC_num = 4;
p_sig = 0.05;
minspikesnum = 10;
significance_matrix = tSC_coupling_test(Matrix,p_sig,minspikesnum,ones(1,n),tSC_num);

% Panel B
tSC_neuron_coupling(mainpath,'20161989','163164','3',1,0,10)
tSC_neuron_coupling(mainpath,'20161988','281282','12',1,0,2)

% S5
tSC_neuron_coupling(mainpath,'20161869','45','32',1,0,3)
tSC_neuron_coupling(mainpath,'20161865','104105','31',1,0,3)
tSC_neuron_coupling(mainpath,'20161988','281282','12',1,0,1)

% Panel C
tSC_population_psth(Matrix,significance_matrix,tSC_num,160,12);

% Panel D - F
tSC_coupled_neurons_distribution(significance_matrix,tSC_num,vertcat(Matrix.pRayleigh),vertcat(Matrix.rythmgroupnum));

% Fig S6
tSC_coupling_shufflingcontroll(mainpath,p_sig,minspikesnum,tSC_num);
%% Figure 4, S7 and S8
% MS neuron firing predicts tSCs signals.
%--------------------------------------------------------------------------
% F4
significance_matrix = tSC_coupling_test(Matrix,p_sig,minspikesnum,ones(1,n),tSC_num);
tSC_Zshift_plot(Matrix,significance_matrix,tSC_num);

% S7 
significance_matrix = tSC_coupling_test(Matrix,p_sig,minspikesnum,rythmgroupnum == 1,tSC_num);
tSC_Zshift_plot(Matrix,significance_matrix,tSC_num);
tSC_coupled_neurons_distribution(significance_matrix,tSC_num,pRayleigh_theta,rythmgroupnum);

% S8 
significance_matrix = tSC_coupling_test(Matrix,p_sig,minspikesnum,rythmgroupnum == 13,tSC_num);
tSC_Zshift_plot(Matrix,significance_matrix,tSC_num);
tSC_coupled_neurons_distribution(significance_matrix,tSC_num,pRayleigh_theta,rythmgroupnum);

%% Figure 5 and S9-S13
% Optogenetic stimulation of PV-expressing mS neurons evokes tSC-like
% activity patterns in the CA1.
%-------------------------------------------------------------------------

% Figure 5 (except for panel C), S12C and S13 is created in python with the
% plotMeantSCcycle function of the 'Analysis framework for the extraction
% of theta-nested spectral components' package by Lopes-dos-Santos et al.,
% 2018. S9 demonstrates track reconstruction and built up from micropsopy
% images.

% 5C,S11A and S11B is created with Cellbase 
% (https://github.com/hangyabalazs/CellBase) 
choosecb('tSC_stim')
load(getpref('cellbase','fname'));
tsc_interneuron_stim('2391_210112a_1.44');
tsc_interneuron_stim('2394_210202b_1.4');
tsc_interneuron_stim('2391_210112a_1.84');

% S10
mainpath = 'L:\Balint\tsc\STIM_Exp';
animal = '2384';
ch = '86';
time_window = 200;
figure
tmS_STA(mainpath, animal, '2', ch, 2, 2, time_window, 4, 1) % tmS 22 Hz
tmS_STA(mainpath, animal, '2', ch, 3, 3, time_window, 4, 2) % tmS 35 Hz
tmS_STA(mainpath, animal, '2', ch, 4, 4, time_window, 4, 3) % tmS 54 Hz
tmS_STA(mainpath, animal, '3', ch, 1, 4, time_window, 4, 4) % tmS 80 Hz

% S12B 
figure
tmS_STA(mainpath, animal, '2', ch, 1, 1, time_window, 1, 1) % 8 Hz stim

% S11C
figure
nchannels = 276;
samplingrate_ap = 20000;
fr = NaN(1,length(CELLIDLIST));
p2v = NaN(1,length(CELLIDLIST));
for neuron = 1:length(CELLIDLIST)
    [fr(neuron),p2v(neuron)] = get_waveforms(CELLIDLIST(neuron),nchannels,samplingrate_ap,0);    
end
scatter(fr,p2v,'.');
%% Figure 6, S14 -- MS neurons show phase coupling to hippocampal tSCs in
% anesthetized rodents
%--------------------------------------------------------------------------

% Figure 6
mainpath='L:\Balint\tsc\anast_rat\';
tsc_anesthetized_fig(mainpath)

% Figure S14
mainpath='L:\Balint\tsc\anast_mouse\';
tsc_anesthetized_fig(mainpath)