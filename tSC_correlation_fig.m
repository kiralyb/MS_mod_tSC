function tSC_correlation_fig(ses_Matrix,tSC_num)
%TSC_CORRELATION_FIG Wrapper for figure S1D-E presented in the Király et
%al. manuscript.
%   TSC_CORRELATION_FIG(SES_MATRIX) Examining correlation of tSC
%   expression strength with theta frequency and the speed of the animal.
%
%   Required input arguments:
%       SES_MATRIX: Database of theta cycles lengths, tSC expression
%       strengths and the speed of the animal during theta cycles extracted
%       with the TSC_THETAPROP function.
%
%   See also TSCSTRENGTH_CORRELATION and TSC_THETAPROP

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   1-Jan-2022
global Samplingrate

figure
maximize_figure
n_cycles = sum(vertcat(ses_Matrix.numberofallcycles));
n_ses = length(ses_Matrix);
for tsc = 1:tSC_num
    % load data
    lengths = NaN(1,n_cycles);
    strengths = NaN(1,n_cycles);
    length_inx = 1;
    speeds = NaN(1,n_cycles);
    strengths_s = NaN(1,n_cycles);
    lengths_s = NaN(1,n_cycles);
    length_inx_s = 1;
    for session = 1:n_ses
        lengths(length_inx:length_inx+ses_Matrix(session).numberofallcycles-1) = 1 ./ double(ses_Matrix(session).allcyclelength) * Samplingrate;
        strengths(length_inx:length_inx+ses_Matrix(session).numberofallcycles-1) = ses_Matrix(session).tscStrength(:,tsc);
        length_inx = length_inx + ses_Matrix(session).numberofallcycles;
        if ~isempty(ses_Matrix(session).cycles_speed)
            speeds(length_inx_s:length_inx_s+ses_Matrix(session).numberofallcycles-1) = ses_Matrix(session).cycles_speed_all';
            strengths_s(length_inx_s:length_inx_s+ses_Matrix(session).numberofallcycles-1) = ses_Matrix(session).tscStrength(:,tsc);
            lengths_s(length_inx_s:length_inx_s+ses_Matrix(session).numberofallcycles-1) = 1 ./ double(ses_Matrix(session).allcyclelength) * Samplingrate;           
            length_inx_s = length_inx_s + ses_Matrix(session).numberofallcycles;
        end
    end
    strengths_s(isnan(speeds)) = [];
    lengths_s(isnan(speeds)) = [];
    speeds(isnan(speeds)) = [];

    
    % tSC strength - theta frequency correlation    
    subplot(3,tSC_num,tsc)
    R_L = tSCstrength_correlation(lengths',strengths',[4,12]);
    xlabel('Frequency (Hz)')
    title(sprintf('tsc%d',tsc))
      
    % tSC strength - speed correlation
    subplot(3,tSC_num,tsc + tSC_num)
    R_S = tSCstrength_correlation(speeds',strengths_s',[0,25]);
    xlabel('Speed')
    
    % partial correlations
    subplot(3,tSC_num,tsc + tSC_num*2)
    [rho,p] = partialcorr([strengths_s;lengths_s;speeds]');
    rho(1,2) = R_L;
    rho(1,3) = R_S;
    rho(2,3) = corr(lengths_s',speeds','Type','Spearman');

    imagesc(rho);
    caxis([-0.27,0.27]);
    setmyplot_balazs
    
end
