function [theta, s1, s2] = unifying_and_short_killer(theta, nsr, windowS)
%UNIFYING_AND_SHORT_KILLER Cuts short isolated sections from vector.
%   [THETA,S1,S2] = UNIFYING_AND_SHORT_KILLER(THETA,NSR,WINDOWS) finds and
%   erase too short segments (both theta and delta) in a logical vector.
%   Parameters:
%   THETA: is a logical vector representing two states (0 or 1).
%   NSR: smapling rate
%   WINDOWSS: window size
%   S1: vector, time points of 0 to 1 transitions.
%   S2: vector, time points of 1 to 0 transitions.
%
%   Author: Barnabas Kocsis
%   Institute of Experimental Medicine, MTA
%   Date: 03/08/2018

if nargin < 3
    windowS = 5; %window size (neglect segments shorter than this value)
end

windowS = windowS*nsr;

%unifier (kills short delta segments):
dtheta = diff(theta);
s1 = find(dtheta==1);  % change from non-theta to theta.
s2 = find(dtheta==-1);  % change from theta to non-theta


for it = 1:length(s1)-1
    if s1(it+1)-s2(it) < windowS
        theta(s2(it):s1(it+1)) = 1;
    end
end

%short_killer (kills short theta segments):
theta = [0 theta 0];
dtheta = diff(theta);
s1 = find(dtheta==1);  % change from non-theta to theta.
s2 = find(dtheta==-1);  % change from theta to non-theta

for it = 1:length(s1)-1
    if s2(it)-s1(it) < windowS
        theta(s1(it):s2(it)) = 0;
    end
end

theta = theta(3:end-2);
dtheta = diff(theta);
s1 = find(dtheta==1);
s2 = find(dtheta==-1);
end