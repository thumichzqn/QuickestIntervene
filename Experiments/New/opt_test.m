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

%% Optimal Policy
% Ja = cell(A,Z);
% for adx = 1:A
%     for zdx = 1:Z
%         Ja{adx,zdx}=zeros(quantnum,T);
%         Ja{adx,zdx}(:,T)=ones(quantnum,1)*cp(zdx);
%     end
% end
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
