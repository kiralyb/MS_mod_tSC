function [p_tSC_phasepref] = tSC_population_psth(Matrix,significance_matrix,tSC_num,binnum,kernel_size)
%TSC_POPULATION_PSTH plot population PSTHs around tSC troughs.
%   [P_TSC_PHASEPREF] =
%   TSC_POPULATION_PSTH(MATRIX,FILTER,TSC_NUM,BINNUM,KERNEL_SIZE)
%   Plotting tSC trough triggered peri-stimulus phase histograms for tSC
%   coupled neurons.
%
%   Required input arguments:
%       MATRIX: database containing tSC phase values at spike times
%       extracted with the TSC_NEURON_FIRINGPROP for each neuron.
%       SIGNIFICANCE_MATRIX: logic matrix showing whether each neuron is
%       coupled to different tSCs.
%       TSC_NUM: number of tSCs plotted
%       BINNUM: number of bint for the PSTH
%       KERNEL_SIZE: width of the kernel used for smoothing the PSTHs (ms)
%
%   See also  TSC_NEURON_FIRINGPROP and TSC_COUPLING_TEST

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   06-Nov-2021

% Define the phase range we want to plot and smoothing kernel
kernel = normpdf(-kernel_size:kernel_size,0,kernel_size / 2);
phases = linspace(-4 * pi(),4 * pi(),binnum);

% Find preferred tSC phase of the neurons
hang_tsc_Z = NaN(length(Matrix),tSC_num);
for i = 1:length(Matrix)
    for tsc = 1:tSC_num
        ftm = sum(exp(1).^(1i * Matrix(i).phaserel2tsc_Z{tsc})) / length(Matrix(i).phaserel2tsc_Z{tsc});    % first trigonometric moment
        hang_tsc_Z(i,tsc) = angle(ftm) + pi();   % mean angle
    end
end

figure 
for tsc = 1:tSC_num %tSC loop
    inx = find(significance_matrix(tsc,:) == 1); % find neurons coupled to the given tSC
    
    % Calculate Z-scored smoothed PSTHs
    Z_spsth = NaN(length(inx),binnum); %preallocate
    for i=1:length(inx)
        psth = hist(Matrix(inx(i)).phaserel2tsc_A{tsc},binnum);
        spsth = conv(psth,kernel);
        Z_spsth(i,:) = zscore(spsth(kernel_size + 1:end - kernel_size));
    end
    
    % Sort PSTHs based on the preferred phase in two blocks: neurons
    % with firing rate maximum before/after the trigger
    H_temp = hang_tsc_Z(:,tsc) + pi();
    H_temp(H_temp < 0) = abs(H_temp(H_temp < 0));
    [~,ind_temp] = sort(H_temp(inx));
    [~,H] = max(Z_spsth,[],2);
    ind = [ind_temp(H(ind_temp) < binnum / 2 + 1);ind_temp(H(ind_temp) > binnum / 2)];
    
    % Plot population PSTHs as a heatmap
    subplot(2,tSC_num,tsc)
    imagesc(phases(kernel_size / 2 + 1:end - kernel_size / 2),1:length(inx),squeeze(Z_spsth(ind,kernel_size / 2:end - kernel_size / 2)));
    caxis([-2,2]);
    
    % Mark tSc troughs with white lines
    hold on
    line([0,0],[0,length(inx)+1],'color','white','linewidth',1)  
    line([-2 * pi(),-2 * pi()],[0,length(inx) + 1],'color','white','linewidth',1)
    line([2 * pi(),2 * pi()],[0,length(inx) + 1],'color','white','linewidth',1)
    
    % Mark the border of the two blocks
    line([phases(1),phases(end)],[length(ind_temp(H(ind_temp) < binnum / 2 + 1)),length(ind_temp(H(ind_temp) < binnum / 2 + 1))],'color','white','linewidth',0.5)
    %xlim([-7,7])
    
    % Plot preferred tSC phase distribution of the neurons on rose plot
    subplot (2,tSC_num,tsc + tSC_num)  
    rose(hang_tsc_Z(significance_matrix(tsc,:),tsc),18);
    
    % Test the uniformity of the distributions 
    [~,p_tSC_phasepref,~,~] = b_rao(hang_tsc_Z(significance_matrix(tsc,:),tsc)') ; 
end