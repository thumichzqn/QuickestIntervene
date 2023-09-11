function [V_init, D_init, D_r] = Cpt_opt(T, quantnum, cp, ci, Beta, Alpha, lambda)
    % Initialization
    Z = length(cp);
    A = length(ci);
    
    V=cell(1,Z);
    D=cell(1,Z);
    for zdx = 1:Z
        V{zdx} = zeros(quantnum, T);
        V{zdx}(:,T) = ones(quantnum, 1) * cp(zdx);
        D{zdx} = zeros(quantnum, T);
    end
    %%
    for t = T-1: -1: 1
        for idx = 1 : quantnum
            pi = idx2pi(idx, quantnum);
            pi_hat = pi + lambda * (1 - pi);
            for zdx = 1:Z
                J = zeros(A,1);
                for adx = 1:A
                    J(adx) = ci(adx) + cp(zdx);
                    for zddx=1:Z
                        pi_new = Piupdate(zddx, Beta{adx}, Alpha, lambda, pi);
                        sigma = (1 - pi_hat) * Alpha(zddx) + pi_hat * Beta{adx}(zddx);
                        J(adx) = J(adx) + sigma * V{zddx}(pi2idx(pi_new,quantnum), t + 1);
                    end
                end
                [value, decision] = min(J);
                V{zdx}(idx, t) = value;
                D{zdx}(idx, t) = decision;
            end
        end
    end
    D_r = D{1};
    %%
    t = 0;
    pi = 0;
    pi_hat = pi + lambda * (1 - pi);
    J = zeros(A,1);
    for adx = 1:A
        J(adx) = ci(adx);
        for zddx=1:Z
            pi_new = Piupdate(zddx, Beta{adx}, Alpha, lambda, pi);
            sigma = (1 - pi_hat) * Alpha(zddx) + pi_hat * Beta{adx}(zddx);
            J(adx) = J(adx) + sigma * V{zddx}(pi2idx(pi_new,quantnum), t + 1);
        end
    end
    [value, decision] = min(J);
    V_init = value;
    D_init = decision;
end