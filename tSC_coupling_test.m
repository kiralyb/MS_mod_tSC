function significance_matrix = tSC_coupling_test(Matrix,p_sig,minspikesnum,filter,tSC_num)
%TSC_COUPLING_TEST Tests tSC coupling of neurons with corrections 
%   SIGNIFICANCE_MATRIX =
%   TSC_COUPLING_TEST[MATRIX,P_SIG,TIME_WINDOW,MINSPIKESNUM,FILTER,TSC_NUM]
%   examines wheter we found significant coupling to tSCs at the P_SIG
%   significance level at any of the temporal shifts pretested with the
%   TSC_NEURON_COUPLING function and stored in MATRIX. Storey false
%   discovery rate method is applied to correct for multiple comparisons.
%   Nerons with less spikes than MINSPIKESNUM during theta cycles with tSC
%   are automatically considered non coupled. Test can be restricted to the
%   neurons defined in the logical arrey FILTER. TSC_NUM defines the number
%   of tSC required to be tested. 
%
%   See also TSC_NEURON_COUPLING

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   06-Noc-2021


% Preallocate memory 
neuron_num = length(Matrix);
significance_matrix = NaN(4,neuron_num);

% Fix random seed for reproducability
rng(7)
warning ('off','all')

% Find significantly tSC coupled neurons
for tsc = 1:tSC_num % tSC loop
    for neuron = 1:neuron_num; % neuron loop
        % Find the number of spikes during theta cycles with tSC, and if it
        % is less then MINSPIKESNUM consider it non coupled 
        Lengths = cellfun(@size,Matrix(neuron).phase_neuron,'uni',false);
        minspikes = min(vertcat(Lengths{1:tSC_num}));
        if minspikes(2) > minspikesnum && filter(neuron) == 1;
            try
                % Apply Storey’s false discovery rate method to adjust p-values
                [~,q] = mafdr((Matrix(neuron).p_rayleigh_corr(tsc,:)));
                % Test significance at any of the temporal shifts.
                significance_matrix(tsc,neuron) = min(q) < p_sig;
            catch
                % If the method can't be applies because of missing data (2
                % neurons in the awake mouse dataset) consider it non
                % coupled.
                significance_matrix(tsc,neuron) = NaN;
            end
        end
    end
end
% Find neurons not coupled to any of the tSCs
significance_matrix(5,:) = sum(significance_matrix) == 0;
significance_matrix(isnan(significance_matrix)) = 0;
significance_matrix = logical(significance_matrix);
