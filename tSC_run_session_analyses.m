function tSC_run_session_analyses(mainpath, animal, session, ch, isstim, isrun, exp)
%TSC_RUN_SESSION_ANALYSES analyse tSC -- single unit recording sessions
%
%   TSC_RUN_SESSION_ANALYSES(MAINPATH,ANIMAL,SESSION,CH,ISSTIM,ISRUN,EXP)
%   run all analyses on a session of simultaneous unit and LFP recordings
%   and store results in the Matrix.mat and ses_Matrix.mat files. This code
%   is for further analyses after single unit clusters and tSC were
%   extracted with Kilosort and the tSC extraction package.
%
%   Required input arguments:
%       MAINPATH: the acces route of the folder containing the results
%           matrix and the database: folders with animal IDs containing
%           folders with session IDs, containing the raw data and
%           extracted tSCs in the \raw folder).
%       ANIMAL: ID of the animal
%       SESSION: ID of the session
%       CH: number of the channel used for tSC extraction
%       ISSTIM: logical variable indicating whether there are stimulation
%           periods that need to be exluded from the data.
%       ISRUN: logical variable indicating wheter there is movement data
%           for freely behaving animals.
%       EXP: name of the experiment used for choosing between
%           preset settings for oscillation state detection defined by
%           Kocsis et al, 2020 . Possible inputs: 'awake_mouse',
%           'anesthetized_mouse', 'anesthetized_rat'.
%
%   See also THETA_DETECTION

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   03-Jan-2022

% theta property and tSC presence correlations
tSC_thetaprop(mainpath, animal, session, ch, isstim, isrun);
tSC_ratio_oscstate(mainpath, animal, session, ch, isstim, exp);
% inter-spike-interval histograms
logISI(mainpath, animal, session, isstim)
% Fourier and Wavlet spectras
neuron_spectra(mainpath ,animal, session, ch, isstim);
% firing property and tSC presence correlations
tSC_neuron_firingprop(mainpath,animal, session, ch, isstim);
% tSC copuling of single units
tSC_neuron_coupling(mainpath,animal, session, ch, isstim,1);
% controll analysis with shuffled spikes for tSC coupling
tSC_neuron_coupling(mainpath,animal, session, ch, isstim,0);

