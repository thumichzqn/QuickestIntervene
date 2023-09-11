clear all; close all;

Z = 5;
A = 4;
T = 50;
rho = 1-1/(1*T);
lambda_num = 8;
lambdas = logspace(-2,-0.9,lambda_num);
quantnum = 1000;
cp=[0,1,2,3,4];
ci=[0,0.04,0.12,0.2];

Beta=cell(A);


delta=0.01;
Alpha=[0.2,0.2,0.2,0.2,0.2];
Beta{4}=[delta,   delta,  1-4*delta,    delta,  delta];
Beta{3}=[delta,   delta,  0.5-3*delta,    0.5,  delta];
Beta{2}=[delta, delta, delta, 1-4*delta, delta];
Beta{1}=[delta,   delta,  delta,    delta,  1-4*delta];

fa_obj = [0.005,0.01,0.05,0.1,0.5];

%%
for fa_idx = 1:5
    for lambda_idx = 1: lambda_num
        lambda = lambdas(lambda_idx);
        opt_thres(fa_idx,lambda_idx) = FA2Thres(Z,A,T,lambda,quantnum,Beta,Alpha,fa_obj(fa_idx));
    end
end
%%
sim_fa(opt_thres(3,8), quantnum, lambda, Alpha, Beta, Z, A, T)
%%
save('QCD_thres.mat','opt_thres')
