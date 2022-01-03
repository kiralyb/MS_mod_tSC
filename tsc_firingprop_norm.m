function y_norm = tsc_firingprop_norm(y)
%TSC_FIRINGPROP Normalize firing propertie (y) during theta cycles with
%different tSC and without tSC as: (property - mean property)/mean property

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   01-Jan-2022

y_norm = bsxfun(@rdivide ,bsxfun(@minus, y(:,1:end - 1), y(:,end)), y(:,end));