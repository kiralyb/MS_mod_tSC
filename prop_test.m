function [h,p] = prop_test(X,Y)
% Chi-square test

    df=length(X)-1; % samples
    N = sum(X)+sum(Y);
    N_x = sum(X);
    N_y = sum(Y);
    observed = [X Y];
    expected = [(X+Y)/N*N_x (X+Y)/N*N_y];
    
    % Standard Chi-square test
    chi2stat = sum((observed-expected).^2 ./ expected);
    p = chi2pval2(chi2stat,df);
    
    h=0;
    if p<0.05
        h=1;
    end
end