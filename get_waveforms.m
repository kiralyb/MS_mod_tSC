function [fr,p2v] = get_waveforms(cellid,nChannels,SR,isplot)
%GET_WAVEFORMS plots average waveform of a cluster's spikes. 
%   GET_WAVEFORMS(CELLID,NCHANNELS,SR,ISPLOT) Calculates average waveform,
%   firing rate, peak to valley time for the CELLID cluster and plots the
%   average waveform on the four channels with the highest amplitude.
%
%   Required input arguments:
%       cellid: cellid from Cellbase
%       NCHANNELS: number of channels
%       SR: sampling rate
%       isplot: logical variable indicating whether average waveforms 
%               should be plotted 
%
%   Bálint Király and Barnabas Kocsis
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   23-Dec-2021

% Load data
txt = split(cellid,'_');
fullpath = fullfile(getpref('cellbase','datapath'), char(txt(1)), char(txt(2)));
ID = split(txt(3),'.');
cluId = str2double(ID(2)); %cluster Id (to analyze, indexing starts from 0!)
datafile = fopen(fullfile(fullpath, 'continuous.dat'), 'r');
clustersIDs = load(fullfile(fullpath,[char(txt(1)), char(txt(2)),'_goodClusterIDs']));
st = readNPY(fullfile(fullpath, 'spike_times.npy')); %spike times
clu = readNPY(fullfile(fullpath,'spike_clusters.npy')); %spike clusters (corresponding to st)
ST = st(clu==clustersIDs.goodCluIds(cluId)); %spike times corresponding to cluId

% Extract waveforms
recPoints = 82;
waveforms = zeros(length(ST), nChannels, recPoints);
cntr = 1;
for buffered = 1:ST(end)/SR
    data = fread(datafile,[nChannels,SR],'int16');
    bufferedST = ST((buffered-1)*SR + 1<ST & ST<buffered*SR);
    for it1 = 1:length(bufferedST)
        tPData = (bufferedST(it1)-SR*(buffered-1)); %timepoint in buffered data
        if tPData + recPoints/2-1 <SR && tPData - recPoints/2-1 > 0 %if sticking out the waveform, renounce it
            waveforms(cntr, :, :) = data(:, tPData-(recPoints/2-1):tPData+recPoints/2);
            cntr = cntr + 1;
        end
    end
end
fclose(datafile);

% calculagte average waveforms
avrgWave = zeros(recPoints, nChannels);
for it1 = 1:nChannels
    avrgWave(:, it1) = mean(squeeze(waveforms(:, it1, :))).'; %average waveform on one channel
    avrgWave(:, it1) = avrgWave(:, it1)-mean(avrgWave(:, it1)); %align around zero
    avrgWave(:, it1) = avrgWave(:, it1) + abs(min(avrgWave(:, it1))); %shift above zero
end

% get firing rate and spike width (peak-to-valley time)
[~,sortedchs] = sort(min(avrgWave)-max(avrgWave));
wv = avrgWave(:,sortedchs(1))*-1;
[~, postvalleytime] = min(wv(recPoints/2+1:end));   % post-valley
p2v = postvalleytime / SR ;
fr = length(ST)/double(st(end)-st(1))*SR;

% plot average waveforms from the four channel with the highest amplitude
if isplot
    plot(avrgWave(:,sortedchs(1:4)) + [0,200,400,600],'r')
end
