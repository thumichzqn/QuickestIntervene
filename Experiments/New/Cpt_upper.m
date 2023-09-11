function [D_init, D] = Cpt_upper(T, quantnum, cp, ci, Beta, Alpha, lambda)
    %% 
    Z = length(cp);
    A = length(ci);
    rho = 1- 1/T;
    D = zeros(quantnum, 1);
    for a=1:A-1
        pithres(a)=(ci(a+1)-ci(a))/((1-lambda)*rho*(Beta{a}*cp'-Beta{a+1}*cp'))-lambda/(1-lambda);
    end
    pithres(A)=1;
    for a=1:A-1
        if pithres(a)<=0
            pithres(a)=1/quantnum;
        elseif pithres(a)>1
            pithres(a)=1;
        end
    end
    for a=1:A
        if a == 1
            down = 1;
            up = pi2idx(pithres(a),quantnum) - 1;
        else
            down = pi2idx(pithres(a-1),quantnum);
            up = pi2idx(pithres(a),quantnum) - 1;
        end
        if up >= down
            D(down: up) = a;
        end
    end
    D(end) = D(end-1);
    D_init = D(1);
end