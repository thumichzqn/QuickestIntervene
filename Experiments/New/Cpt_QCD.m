function [D_init, D] = Cpt_QCD(T, quantnum, cp, ci, Beta, Alpha, lambda)
    %% Initialization
    Z = length(cp);
    A = length(ci);
    rho = 1 - 1/T;
    
    VA=(ci(A)+rho*Alpha*cp')/(1-rho);
    V=zeros(quantnum,1);
    D=zeros(quantnum,1);

    Vpre=ones(quantnum,1);
    J=0;
    count=0;
    while norm(abs(V-Vpre))>10e-3
        Vpre = V;
        for idx = 1:quantnum
            pi = idx2pi(idx,quantnum);
            pi_hat = pi + lambda * (1 - pi);
            J = ci(1);
            for zddx=1:Z
                pi_new = Piupdate(zddx,Beta{1},Alpha,lambda,pi);
                sigma = (1 - pi_hat) * Alpha(zddx) + pi_hat * Beta{1}(zddx);
                J = J+rho * sigma * (Vpre(pi2idx(pi_new,quantnum))+cp(zddx));
            end
            V(idx)=min(J,VA);
            if V(idx)==J
                D(idx)=1;
            else
                D(idx)=A;
            end
        end
        count=count+1;
        if count==-1
            break;
        end
    end
    D_init = D(1);
end