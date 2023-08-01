function burstisi_tSCs(Matrix,inx)
%BURSTISI_tSCs Plotter function for Fig S= 
% Plot the average difference between the inter-spike interval histograms
% during theta cycles with the most and least tSC content
% 
%   Required input arguments:
%       Matrix: Pre-processed data createdwith the
%       TSC_NEURON_FIRINGPROP_ISI function
%       inx: index of the cells considered for averiging

%   Bálint Király
%   Institute of Experimental Medicine, Budapest, Hungary
%   kiraly.balint@koki.hu
%   01-Aug-2023

%initialize histogram difference variable
SN= [];
j=1;
for i = inx
    tt = 1;
    for t = Matrix(i).tSC_pos_1
    %S(j,tt,:) = Matrix(i).HtSC(t,:)-Matrix(i).Hctrl(t,:);
    SN(j,tt,:) = (Matrix(i).HtSC(t,:)-Matrix(i).Hctrl(t,:))/Matrix(i).mfrate(7);
    %H(j,tt,:) = Matrix(i).HtSC(t,:)/sum(Matrix(i).Hctrl(t,:));
    %C(j,tt,:) = Matrix(i).Hctrl(t,:)/sum(Matrix(i).Hctrl(t,:));
    tt = tt+1;
    end
 j=j+1;
end
    
figure
for tSC = 1:4
    
    subplot(1,4,tSC)
    hold on
    
    errorshade(2.5:5:67.5,mean(squeeze(SN(:,tSC,1:14))),std(squeeze(SN(:,tSC,1:14)))/sqrt(j),'LineColor',Colors(tSC,:))
    ylim([-0.4,2.6])
    
    hold  on
    yyaxis right
    plot(2.5:5:67.5,mean(squeeze(SN(:,tSC,1:14))),'Color',Colors(tSC,:));
    set(gca, 'YScale', 'log')
    setmyplot_balazs
    ylim([7*10^-3,3])
    yyaxis left
    xlim([0,70])
end

for i = inx
    tt = 1;
    for t = Matrix(i).tSC_pos_1
    subplot(1,4,tt)
    line([1000/Matrix(i).mainfreqs(t),1000/Matrix(i).mainfreqs(t)],[-0.4,2.6])
    tt = tt+1;
    end
    
    
end