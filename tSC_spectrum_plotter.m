function tSC_spectrum_plotter(tSC_spect,IDS)
%TSC_SPECTRUM_PLOTTER Plot spectral content of tSC
%   TSC_SPECTRUM_PLOTTER(TSC_SPECT, IDS) plots  
%   average tSC spectras with error shades (SEM) and peak frequencies.

%   Required input arguments:
%       TSC_SPECT: tSC_spectra data (frequencies x tSCs x sessions)
%       IDS: session IDS included from tSC_spect

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   19-Dec-2021

global Colors
freqstart = 200 - size(tSC_spect,1) + 1;
figure
for tSC = 1:size(tSC_spect,2)
    errorshade(freqstart:size(tSC_spect,1)+freqstart - 1,mean(tSC_spect(:,tSC,IDS),3),std(tSC_spect(:,tSC,IDS),0,3) / sqrt(size(tSC_spect,3)),'LineColor',Colors(tSC,:),'ShadeColor',Colors(tSC,:));
    [peak,peakfreq]=max(mean(tSC_spect(:,tSC,IDS),3));
    text(peakfreq + freqstart - 1,peak,[int2str(peakfreq + freqstart - 1),' Hz'])
end
xlabel('Frequency (Hz)')
setmyplot_balazs
