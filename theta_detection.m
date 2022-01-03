function [theta, thetaTransf, deltaTransf] = theta_detection(fieldPot, thBand, deBand, nsr, tscgwindow, ratioTresh, windowS)
%THETA_DETECTION Detects theta segments in (hippocampal) field potentials.
%   [THETA,THETATRANSF,DELTATRANSF] = THETA_DETECTION(FIELDPOT,
%   THBAND,DEBAND,NSR,TSCGWINDOW,RATIOTRESH,WINDOWS)
%   Additionally, calculates phases of dominant oscillations (hilbert transf.).
%   Filters field data for theta (THBAND). (For this
%   purposes bandpass filter is used.) Centralize, normalize
%   ((feeg-mean(feeg))/standard deviation(eeg)) and hilbert transform
%   (calculating angle and amplitude) them.
%   Definition: if theta's amplitude's ratio>RATIOTRESH, theta 
%   is dominant against delta (ratio is filtered before with moving
%   average).
%   Parameters:
%   FIELDPOT: field potential vector.
%   THBAND: theta band (Hz) (e.g. [3,8]).
%   DEBAND: non-theta (delta) band (Hz) (e.g. [0.5,3]).
%   NSR: sampling rate (e.g. 1000).
%   TSCGWINDOW: window (in sec) for theta/delta amplitude ratio's smoothing
%   (e.g. 5).
%   THRATIOTRESH: theta/delta amplitude ratio threshold for dominant theta 
%   (e.g. 1).
%   WINDOWS: window size for negligibly small theta-delta segments (e.g. 5).
%
%   See also HIPPO_STATE_DETECTION.

%   Author: Barnabas Kocsis
%   Institute of Experimental Medicine, MTA
%   Date: 03/08/2018

%Filter in theta band:
firOrder = 1024;
thetaFilter = fir1(firOrder,thBand/(nsr/2),'bandpass');
thetaFeeg = filtfilt(thetaFilter,1,fieldPot);
thetaSFeeg = (thetaFeeg - mean(thetaFeeg)) ./ std(fieldPot); % standardize feeg
thetaTransf = hilbert(thetaSFeeg); %hilbert transform to frequency domain
thetaAmp = abs(thetaTransf); %theta amplitude

%Filter in delta band:
deltaFilter = fir1(firOrder,deBand/(nsr/2),'bandpass');
deltaFeeg = filtfilt(deltaFilter,1,fieldPot);
deltaDFeeg = (deltaFeeg - mean(deltaFeeg)) ./ std(fieldPot);
deltaTransf = hilbert(deltaDFeeg);
deltaAmp = abs(deltaTransf);

% Theta-delta ratio
ratio_feeg = thetaAmp ./ (deltaAmp);
ratio_feeg(ratio_feeg>10) = 10; %replace too large values

% Smooth ratio with moving average:
wSize = nsr * tscgwindow; %window size
coeff = ones(1, wSize)/wSize;
ratioSFeeg = filtfilt(coeff, 1, ratio_feeg);

% Find transition points:
theta = [0 ratioSFeeg>ratioTresh 0];   % thresholding
[theta, s1, s2] = unifying_and_short_killer(theta, nsr, windowS); %sorting out too short delta and theta segments

timeVec = 1/nsr:1/nsr:length(fieldPot)/nsr; %time vector
standardizedField = (fieldPot-mean(fieldPot))./std(fieldPot); %just for visualization...
hold on
plot(timeVec, standardizedField, 'k')
plot(timeVec, thetaSFeeg, 'Color', [0, 0.4470, 0.7410])
plot(timeVec, deltaDFeeg, 'Color', [0.8500, 0.3250, 0.0980])
domTh = theta*2; %just for visualization...
domTh(s1) = NaN; %just for visualization...
domTh(s2) = NaN; %just for visualization...
plot(timeVec, domTh, 'y')
plot(timeVec, (ratioSFeeg + ones(size(ratioSFeeg))), 'g') %shift the plot to make it visible
legend('standardized field', 'st. and theta filtered field', 'st. and delta filtered field', 'theta', 'averaged ratio');
xlabel('Seconds');
hold off;
end