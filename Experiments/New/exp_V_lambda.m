clear all; close all;

Z = 5;
A = 4;
T = 200;
rho = 1-1/(1*T);
lambda_num = 8;
lambdas = logspace(-2,-0.5,lambda_num);
quantnum = 1000;
unit_cp = 1; 
cp=[0,1,2,3,4] * unit_cp;
ci=[0,0.04,0.12,0.18];

Beta=cell(A);

delta=0.04;
Alpha=[0.2,0.2,0.2,0.2,0.2];
Beta{4}=[0.2 - delta,   0.3 - delta,    0.2 + 2*delta-0.001,  0.2,            0.001];
Beta{3}=[0.2 - delta,   0.2 - delta,    0.2 - delta,    0.4 + 3*delta - 0.001,  0.001];
Beta{2}=[0.2 - delta,   0.2 - delta,    0.2 - delta,    0.2 + 3*delta,  0.2];
Beta{1}=[0.2 - delta,   0.2 - delta,    0.2 - delta,    0.2 - delta,    0.2 + 4*delta];

%% Simulate Traces
trace_num = 10000;

% Observation
traces = cell(trace_num,1);
for trace_idx = 1:trace_num
    traces{trace_idx} = zeros(A,T);
    for adx = 1:A
        traces{trace_idx}(adx,:) = randsample(Z,T,true,Beta{adx});
    end
end

%% Simulate Change point
Ct_opts_v=zeros(1,lambda_num);
Ct_uppers_v=zeros(1,lambda_num);
Ct_QCDs_v=zeros(1,lambda_num);

for lambda_idx = 1:lambda_num
lambda=lambdas(lambda_idx);

geo_pdf = zeros(T+1,1);
for t=1:T
    geo_pdf(t) = geopdf(t-1,lambda);
end
rest = 1- sum(geo_pdf(1:T));
geo_pdf(T+1) = rest;
changepoints = zeros(trace_num,1);
for trace_idx = 1:trace_num
    changepoints(trace_idx) = randsample(T+1,1,true,geo_pdf);
end

%% upper
[D_init, D] = Cpt_upper(T, quantnum, cp, ci, Beta, Alpha, lambda);
V_init = Sim_BackDP(T, quantnum, cp, ci, Beta, Alpha, lambda, D, D_init);
Ct_uppers_v(lambda_idx) = V_init;

%% QCD
[D_init, D] = Cpt_QCD(T, quantnum, cp, ci, Beta, Alpha, lambda);
V_init = Sim_BackDP(T, quantnum, cp, ci, Beta, Alpha, lambda, D, D_init);
Ct_QCDs_v(lambda_idx) = V_init;

%% opt special case: time-varying
[V_init, D_init, D] = Cpt_opt(T, quantnum, cp, ci, Beta, Alpha, lambda);
Ct_opts_v(lambda_idx) = V_init;
end
%%
figure;
hold on ;
plot(lambdas,Ct_opts_v(1,:))
plot(lambdas,Ct_uppers_v(1,:))
plot(lambdas,Ct_QCDs_v(1,:))
set(gca,'Xscale','log')
legend('opt','upper','qcd','opt_v','upper_v','qcd_v')