function [anticipator_ratio,median_lag,average_lag,halfsumpoint,Lags] = tSC_Zshift_plot(Matrix,significance_matrix,tsc_num)
%TSC_ZSHIFT_PLOT Z-shift analysis between neuron population and tSCs 
%   [ANTICIPATOR_RATIO,MEDIAN_LAG,AVERAGE_LAG,HALFSUMPOINT,LAGS] =
%   TSC_ZSHIFT_PLOT[MATRIX,SIGNIFICANCE_MATRIX,TSC_NUM]
%   search for temporal lags between spikes and tSCs leading to the maximal
%   value of Rayleigh's Z for populations of nezurons defined in
%   SIGNIFICANCE_MATRIX based on preprocessed data with the 
%   TSC_NEURON_COUPLING function and stored in MATRIX. TSC_NUM defines the
%   number of tSC required to be tested. 
%
%   See also TSC_NEURON_COUPLING and TSC_COUPLING_TEST

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   26-Dec-2021
global Colors

% load Z-shift data
Zshift_tsc = vertcat(Matrix.Zshift_tsc);
Zshift_tsc = reshape(Zshift_tsc,[size(Matrix(1).Zshift_tsc,1),length(Matrix),size(Matrix(1).Zshift_tsc,2)]);

% plot settings
kernel_size = 3; % in ms 
kernel = normpdf(-kernel_size:kernel_size,0,kernel_size / 2);
time = 1:size(Zshift_tsc,3);
binsize = (size(Zshift_tsc,3)-1) / 15;

% Preallocate space for the results and initialize the figure 
figure
anticipator_ratio = zeros(1,tsc_num);
median_lag = zeros(1,tsc_num);
average_lag = zeros(1,tsc_num);
halfsumpoint = zeros(1,tsc_num);
Lags = cell(1);
for tsc = 1:tsc_num
    % Find maximal Z-values
    [Z_Max,maxind] = max(squeeze(Zshift_tsc(tsc,:,time)),[],2);
    
    % Plot maximal Z histograms
    subplot(3,tsc_num,tsc)
    H = histogram(maxind(significance_matrix(tsc,:)) + time(1) - median(time),time(1) - median(time):binsize:time(end) - median(time),'FaceColor',Colors(tsc,:));
    [~,sorter]=sort(H.Data);
    hold on
    line([0,0],get(gca,'YLim'),'Color','r')
    setmyplot_balazs;
    xlim([time(1) - median(time),time(end) - median(time)])

    % Normalize Z-curves
    filter = find(significance_matrix(tsc,:) == 1);
    zshift_norm = NaN(length(filter),length(time));
    for f = 1:length(filter)
        zshift_norm(f,:) = squeeze(Zshift_tsc(tsc,filter(f),time)) / Z_Max(filter(f));
    end
    
    % Plot all the Z-curves of the population as a heatmap
    subplot(3,tsc_num,4 + tsc)
    imagesc(time-median(time),1:length(filter),zshift_norm(sorter,:));
    line([0,0],get(gca,'YLim'),'Color','r')
    set(gca,'ytick',length(filter)+0.5)
    set(gca,'yticklabel',length(filter))
    
    % Plot the mean Z-curve
    subplot(3,tsc_num,8 + tsc)
    Avg_Z_curve = nanmean(zshift_norm);
    smooth_avg_zcurve = conv(Avg_Z_curve,kernel);
    errorshade(time - median(time),smooth_avg_zcurve(kernel_size + 1:end - kernel_size),nanstd(zshift_norm) / sqrt(length(filter)),'LineColor',Colors(tsc,:),'ShadeColor',Colors(tsc,:))
    hold on
    line([0,0],get(gca,'YLim'),'Color','r')
    setmyplot_balazs;
    
    % Calculate parameters characterising the temporal relationship of the
    % neuron population with the tSCs
    anticipator_ratio(tsc) = sum(H.Data < 0) / length(H.Data);
    median_lag(tsc) = median(H.Data);
    average_lag(tsc) = mean(H.Data);
    [~,halfsumpoint_temp] = min(abs((cumsum(Avg_Z_curve) - (sum(Avg_Z_curve) / 2))));
    halfsumpoint(tsc) = halfsumpoint_temp - median(time) + time(1) - 1;
    Lags{1} = [Lags{1};H.Data];
end