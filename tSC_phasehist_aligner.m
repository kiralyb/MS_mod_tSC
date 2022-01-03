function phase_hist_align_norm = tSC_phasehist_aligner(Matrix,tSC_num)
% Align phase histograms to the preferred phase of the neuron and
% normalize with the firing rate

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   01-Jan-2022

n = length(Matrix);
frate = vertcat(Matrix.frate);
phase_hist_align_norm = NaN(n,tSC_num + 1,length(Matrix(1).phase_hist(1,:)));
for neuron = 1:n
    shifter = round((pi / 2 - Matrix(neuron).hang(end)) / 2 / pi * size(Matrix(neuron).phase_hist,2)); 
    for tSC = 1:tSC_num + 1
        phase_hist_align_norm(neuron,tSC,:) = circshift(squeeze(Matrix(neuron).phase_hist(tSC,:)),[0,shifter]) / frate(neuron,end) * mean(frate(:,end));
    end
end