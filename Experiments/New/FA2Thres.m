function opt_thres = FA2Thres(Z,A,T,lambda,quantnum,Beta,Alpha,fa_obj)
    %% QCD FA
    threshold_l = 1/quantnum;
    threshold_r = 1 - 1/quantnum;
    fa_l = sim_fa(threshold_l, quantnum, lambda, Alpha, Beta, Z, A, T);
    fa_r = sim_fa(threshold_r, quantnum, lambda, Alpha, Beta, Z, A, T);
    for idx = 1:100
        threshold_mid = (threshold_r + threshold_l) / 2;
        fa_mid = sim_fa(threshold_mid, quantnum, lambda, Alpha, Beta, Z, A, T);
        if fa_mid > fa_l || fa_mid < fa_r || abs(threshold_l-threshold_r) < 3/quantnum
            opt_thres = threshold_mid;
            break;
        else
            if fa_mid > fa_obj
                threshold_l = threshold_mid;
                fa_l = fa_mid;
            else
                threshold_r = threshold_mid;
                fa_r = fa_mid;
            end
        end
    end
end




