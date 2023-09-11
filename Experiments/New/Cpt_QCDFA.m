function [D_init, D] = Cpt_QCDFA(threshold, quantnum, cp, ci, Beta, Alpha, lambda)
    % False alarm vs threshold:
    % 0.005	0.73
    % 0.01	0.65
    % 0.05	0.48
    % 0.1	0.38
    A = length(ci);
    
    D = zeros(quantnum,1);
    thres_idx = pi2idx(threshold,quantnum);
    for idx = 1: thres_idx
        D(idx) = 1;
    end
    for idx = thres_idx : quantnum
        D(idx) = A;
    end
    D_init = 1;
end