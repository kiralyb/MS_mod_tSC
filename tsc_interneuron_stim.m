tsc_interneuron_stim(cellid)
%TSC_INTERNEURON_STIM Plots for figure S11 presented in the Király
%   et al. manuscript.
%   TSC_INTERNEURON_STIM(CELLIDS) plot burst triggered raster plots and
%   PSTHs, average waveforms, and autocorrelograms for a neuron in
%   Cellbase defined by CELLID.

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   23-Dec-2021

% settings
nChannels = 276;
SR = 20000; % in Hz
acg_timewind = 0.1; % in s

% plot stimulus burst triggered raster plots and PSTHs
figure
viewcell2b(cellid,'TriggerName','PulseOn','SortEvent','BurstOff','ShowEvents',{{'PulseOn','PulseOff','BurstOff'}},...
'eventtype','stim','window',[-0.05 0.25],'dt',0.001,'sigma',0.002,'PSTHstd','on','Partitions','#BurstOn_tSC',...
'EventMarkerWidth',0,'PlotZeroLine','off')

% plot average action potential waveforms
figure
subplot(2,1,1)
get_waveforms(cellid,nChannels,SR,1)

% plot autocorrelogram
subplot(2,1,2)
acg(cellid,acg_timewind)