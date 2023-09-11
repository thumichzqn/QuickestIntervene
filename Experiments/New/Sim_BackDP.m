function V_init = Sim_BackDP(T, quantnum, cp, ci, Beta, Alpha, lambda, D, D_init)
    %% Initial
    Z = length(cp);
    A = length(ci);
    rho = 1- 1/T;
    
    %% Back DP
    V=cell(1,Z);
    for zdx = 1:Z
        V{zdx} = zeros(quantnum, T);
        V{zdx}(:,T) = ones(quantnum, 1) * cp(zdx);
    end

    for t = T-1: -1: 1
        for idx = 1 : quantnum
            pi = idx2pi(idx, quantnum);
            pi_hat = pi + lambda * (1 - pi);
            for zdx = 1:Z
                adx = D(idx);
                J = ci(adx) + cp(zdx);
                for zddx=1:Z
                    pi_new = Piupdate(zddx, Beta{adx}, Alpha, lambda, pi);
                    sigma = (1 - pi_hat) * Alpha(zddx) + pi_hat * Beta{adx}(zddx);
                    J = J + sigma * V{zddx}(pi2idx(pi_new,quantnum), t + 1);
                end
                V{zdx}(idx, t) = J;
            end
        end
    end

    t = 0;
    pi = 0;
    pi_hat = pi + lambda * (1 - pi);
    adx = D_init;
    J = ci(adx);
    for zddx=1:Z
        pi_new = Piupdate(zddx, Beta{adx}, Alpha, lambda, pi);
        sigma = (1 - pi_hat) * Alpha(zddx) + pi_hat * Beta{adx}(zddx);
        J = J + sigma * V{zddx}(pi2idx(pi_new,quantnum), t + 1);
    end
    V_init = J;
end