clear all; close all;

Z = 5;
A = 4;
T = 20;
rho = 1-1/(1*T);
lambdas = logspace(-2,-0.01,5);
quantnum=100;
cp=[0,1,2,3,4];
ci=[0,0.02,0.06,0.2];

Beta=cell(A);

delta=0.02;
Alpha=[0.2,0.2,0.2,0.2,0.2];
Beta{4}=Alpha;
Beta{3}=[0.2-2*delta,   0.2-delta,  0.2,    0.2+delta,  0.2+2*delta];
Beta{2}=[0.2-2*2*delta, 0.2-2*delta,0.2,    0.2+2*delta,0.2+2*2*delta];
Beta{1}=[0.2-2*3*delta, 0.2-3*delta,0.2,    0.2+3*delta,0.2+2*3*delta];


lambda=lambdas(1);
%% Calculate QCD Decision
VA=(ci(A)+rho*Alpha*cp')/(1-rho);
V=zeros(quantnum,1);
D=zeros(quantnum,1);

Vpre=ones(quantnum,1);
J=0;
count=0;
while norm(abs(V-Vpre))>10e-6
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


%% Backward DP
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
idx = 1;
adx = D(idx);
J(adx) = ci(adx);
for zddx=1:Z
    pi_new = Piupdate(zddx, Beta{adx}, Alpha, lambda, pi);
    sigma = (1 - pi_hat) * Alpha(zddx) + pi_hat * Beta{adx}(zddx);
    J(adx) = J(adx) + sigma * V{zddx}(pi2idx(pi_new,quantnum), t + 1);
end
V_init = J;
D_init = adx;
