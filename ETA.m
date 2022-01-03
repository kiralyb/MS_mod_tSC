function ETA(triggers,LFP,samplingrate,time_window,columnsnum,column)
%ETA Event triggered average analysis
%   ETA(TRIGGERS,LFP,SAMPLINGRATE,TIME_WINDOW,COLUMNSNUM,COLUMN) performs
%   triggered average analysis on LFP and on the Wavelet power and phase of
%   the LFP. 
%
%   Required input arguments:
%       TRIGGERS: Trigger event times (in data points)
%       LFP1: local field potential data
%       SAMPLINGRATE: sampling rate of the LFP in Hz
%       TIME_WINDOW: % +/- size of the time window around events (ms) 
%       COLUMNSNUM: number of columns plotted
%       COLUMN: number of the current column
%
%   See also EEGWAVELET, ERS and ASTANORM2

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   30-Oct-2021


% Perform wavelet transformation on the LFP data
[pow, phase, freq] = eegwavelet(LFP,samplingrate);

% Calculate spike triggered averages
[sta,~,~,~,~,~,~] = astanorm2(triggers,LFP,samplingrate,time_window * 2);
[pow_sta,~] = ers(triggers,pow,samplingrate,time_window * 2);
[phase_sta,~] = ers(triggers,phase,samplingrate,time_window * 2);
I = log(pow_sta) - log(repmat(mean(pow_sta,2),1,size(pow_sta,2)));

% Calculate spike triggered LFP average
subplot(3,columnsnum + 1,column)
timeVec = -time_window : 1 : time_window; % time vector
plot(timeVec,sta)
hold on
xlim([-time_window,time_window])
line([0,0],get(gca,'YLim'),'Color','red');
ylabel('Average lfp');
xlabel('Time from trigger (ms)');

% Plot spike triggered wavelet power
subplot(3,columnsnum + 1,columnsnum + 1 + column)
[~,uppLim] = min(abs(freq - 100));
fRange = size(I,1);
fValues = freq(uppLim:fRange); % cuts frequency range at 100Hz
fValues = round(fValues,1); % Frequency scale
[~,uppGamma] = min(abs(fValues - 40));
[~,uppTheta] = min(abs(fValues - 12));
[~,uppDelta] = min(abs(fValues - 4));
fPos = [1 uppGamma uppTheta uppDelta length(fValues)];
I2 = I(uppLim:end,:);
imagesc(timeVec,fliplr(fValues),I2)
set(gca,'YTick',fPos)
b_rescaleaxis('Y',round(fValues))
setappdata(gca,'scaley',round(fValues))
b_zoomset_for_wavelet
colormap jet;
cValues = caxis;
set(gca, 'clim', [-min(abs(cValues)), min(abs(cValues))]);
line([0,0],[0,length(fValues)],'Color','red');
ylabel('Frequency (Hz)');
xlabel('Time from trigger (ms)');

% Plot spike triggered wavelet phase
subplot(3,columnsnum + 1,(columnsnum + 1) * 2 + column)
I2 = phase_sta(uppLim:end,:);
imagesc(timeVec,fliplr(fValues),I2)
set(gca,'YTick',fPos)
b_rescaleaxis('Y',round(fValues))
setappdata(gca,'scaley',round(fValues))
b_zoomset_for_wavelet
colormap jet;
cValues = caxis;
set(gca, 'clim', [-min(abs(cValues)), min(abs(cValues))]);
line([0,0],[0,length(fValues)],'Color','red');
ylabel('Frequency (Hz)');
xlabel('Time from trigger (ms)');