function [burstWind]=burst_detector(TS,maxBurstL,recLength)
%BURST_DETECTOR Detects burst-spikes of a MS unit.
%   [BURSTWIND] = BURST_DETECTOR(TS,MAXBURSTL,RECLENGTH).
%   It creates a RECLENGTH long bimodal state vector BURSTWIND (1 between
%   burst starts and ends, 0 elsewhere).
%   Parameters:
%   TS: vector, spike times of the cell.
%   MAXBURSTL: maximal length between intraburst spikes (unit: ms, e.g. 40).
%   RECLENGTH: total length of the given recording (unit: ms, e.g. 2144000).
%   BURSTWIND: logical vector, 1 inside bursts, 0 outside.
%
%   See also CELL_BURST_PARAMETERS, THETA_SKIPPING.

%   Author: Barnabas Kocsis
%   Institute of Experimental Medicine, MTA
%   Date: 11/02/2019

ISI = diff(TS);
burstInd = ISI<maxBurstL;
burstWindow = diff([0;burstInd]);
burstStart = TS(burstWindow == 1);
burstEnd = TS(burstWindow == -1);
if burstInd(end) == 1 %last spike is in a burst
    burstEnd(end+1) = TS(end)-1;
end

burstWind = zeros(1,recLength);
burstWind(burstStart) = 1;
burstWind(burstEnd+1) = -1;
burstWind = cumsum(burstWind);
end