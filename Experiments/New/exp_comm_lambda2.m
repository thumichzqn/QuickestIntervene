clear all; close all;

clear all; close all;

Z = 5;
A = 4;
T = 50;
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

fa_obj = [0.05, 0.01,0.05,0.1,0.5];

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

%% Compute threshold for QCDFA
load('QCD_thres.mat')

%% Simulation
Ct_opts=zeros(2,lambda_num);        % 1:mean 2:std
Ct_uppers=zeros(2,lambda_num);
Ct_QCDs=zeros(2,lambda_num);
Ct_QCDs1=zeros(2,lambda_num);
Ct_QCDs2=zeros(2,lambda_num);
Ct_QCDs3=zeros(2,lambda_num);
Ct_QCDs4=zeros(2,lambda_num);
Ct_QCDs5=zeros(2,lambda_num);

Ct_opts_v=zeros(1,lambda_num);
Ct_uppers_v=zeros(1,lambda_num);
Ct_QCDs_v=zeros(1,lambda_num);
Ct_QCDs_v1=zeros(1,lambda_num);
Ct_QCDs_v2=zeros(1,lambda_num);
Ct_QCDs_v3=zeros(1,lambda_num);
Ct_QCDs_v4=zeros(1,lambda_num);
Ct_QCDs_v5=zeros(1,lambda_num);


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
Ct = zeros(trace_num,1);
for trace_idx = 1:trace_num
    trace = traces{trace_idx};
    changepoint = changepoints(trace_idx);

    pi = 0;
    at = D_init;
    CPt = 0;
    CIt = ci(at);
    for t = 1:T
        if t >= changepoint
            zt = trace(at, t);
        else
            zt = trace(end, t);
        end
        CPt = CPt + cp(zt);

        if t < T
            pi = Piupdate(zt,Beta{at},Alpha,lambda,pi);
            at = D(pi2idx(pi, quantnum));
            CIt = CIt + ci(at);
        end
    end
    Ct(trace_idx) = CIt + CPt;
end
Ct_uppers(1,lambda_idx)=mean(Ct);
Ct_uppers(2,lambda_idx)=std(Ct)/sqrt(trace_num)*2;

%% QCDFA 0.005
fa_set = opt_thres(1,:);
[D_init, D] = Cpt_QCDFA(fa_set(lambda_idx), quantnum, cp, ci, Beta, Alpha, lambda);
V_init = Sim_BackDP(T, quantnum, cp, ci, Beta, Alpha, lambda, D, D_init);
Ct_QCDs_v1(lambda_idx) = V_init;
Ct = zeros(trace_num,1);
for trace_idx = 1:trace_num
    trace = traces{trace_idx};
    changepoint = changepoints(trace_idx);

    pi = 0;
    at = D_init;
    CPt = 0;
    CIt = ci(at);
    for t = 1:T
        if t >= changepoint
            zt = trace(at, t);
        else
            zt = trace(end, t);
        end
        CPt = CPt + cp(zt);

        if t < T
            pi = Piupdate(zt,Beta{at},Alpha,lambda,pi);
            at = D(pi2idx(pi, quantnum));
            CIt = CIt + ci(at);
        end
    end
    Ct(trace_idx) = CIt + CPt;
end
Ct_QCDs1(1,lambda_idx)=mean(Ct);
Ct_QCDs1(2,lambda_idx)=std(Ct)/sqrt(trace_num)*2;

%% QCDFA 0.01
fa_set = opt_thres(2,:);
[D_init, D] = Cpt_QCDFA(fa_set(lambda_idx), quantnum, cp, ci, Beta, Alpha, lambda);
V_init = Sim_BackDP(T, quantnum, cp, ci, Beta, Alpha, lambda, D, D_init);
Ct_QCDs_v2(lambda_idx) = V_init;
Ct = zeros(trace_num,1);
for trace_idx = 1:trace_num
    trace = traces{trace_idx};
    changepoint = changepoints(trace_idx);

    pi = 0;
    at = D_init;
    CPt = 0;
    CIt = ci(at);
    for t = 1:T
        if t >= changepoint
            zt = trace(at, t);
        else
            zt = trace(end, t);
        end
        CPt = CPt + cp(zt);

        if t < T
            pi = Piupdate(zt,Beta{at},Alpha,lambda,pi);
            at = D(pi2idx(pi, quantnum));
            CIt = CIt + ci(at);
        end
    end
    Ct(trace_idx) = CIt + CPt;
end
Ct_QCDs2(1,lambda_idx)=mean(Ct);
Ct_QCDs2(2,lambda_idx)=std(Ct)/sqrt(trace_num)*2;

%% QCDFA 0.05
fa_set = opt_thres(3,:);
[D_init, D] = Cpt_QCDFA(fa_set(lambda_idx), quantnum, cp, ci, Beta, Alpha, lambda);
V_init = Sim_BackDP(T, quantnum, cp, ci, Beta, Alpha, lambda, D, D_init);
Ct_QCDs_v3(lambda_idx) = V_init;
Ct = zeros(trace_num,1);
for trace_idx = 1:trace_num
    trace = traces{trace_idx};
    changepoint = changepoints(trace_idx);

    pi = 0;
    at = D_init;
    CPt = 0;
    CIt = ci(at);
    for t = 1:T
        if t >= changepoint
            zt = trace(at, t);
        else
            zt = trace(end, t);
        end
        CPt = CPt + cp(zt);

        if t < T
            pi = Piupdate(zt,Beta{at},Alpha,lambda,pi);
            at = D(pi2idx(pi, quantnum));
            CIt = CIt + ci(at);
        end
    end
    Ct(trace_idx) = CIt + CPt;
end
Ct_QCDs3(1,lambda_idx)=mean(Ct);
Ct_QCDs3(2,lambda_idx)=std(Ct)/sqrt(trace_num)*2;
%% QCDFA 0.1
fa_set = opt_thres(4,:);
[D_init, D] = Cpt_QCDFA(fa_set(lambda_idx), quantnum, cp, ci, Beta, Alpha, lambda);
V_init = Sim_BackDP(T, quantnum, cp, ci, Beta, Alpha, lambda, D, D_init);
Ct_QCDs_v4(lambda_idx) = V_init;
Ct = zeros(trace_num,1);
for trace_idx = 1:trace_num
    trace = traces{trace_idx};
    changepoint = changepoints(trace_idx);

    pi = 0;
    at = D_init;
    CPt = 0;
    CIt = ci(at);
    for t = 1:T
        if t >= changepoint
            zt = trace(at, t);
        else
            zt = trace(end, t);
        end
        CPt = CPt + cp(zt);

        if t < T
            pi = Piupdate(zt,Beta{at},Alpha,lambda,pi);
            at = D(pi2idx(pi, quantnum));
            CIt = CIt + ci(at);
        end
    end
    Ct(trace_idx) = CIt + CPt;
end
Ct_QCDs4(1,lambda_idx)=mean(Ct);
Ct_QCDs4(2,lambda_idx)=std(Ct)/sqrt(trace_num)*2;

%% QCDFA 0.5
fa_set = opt_thres(5,:);
[D_init, D] = Cpt_QCDFA(fa_set(lambda_idx), quantnum, cp, ci, Beta, Alpha, lambda);
V_init = Sim_BackDP(T, quantnum, cp, ci, Beta, Alpha, lambda, D, D_init);
Ct_QCDs_v5(lambda_idx) = V_init;
Ct = zeros(trace_num,1);
for trace_idx = 1:trace_num
    trace = traces{trace_idx};
    changepoint = changepoints(trace_idx);

    pi = 0;
    at = D_init;
    CPt = 0;
    CIt = ci(at);
    for t = 1:T
        if t >= changepoint
            zt = trace(at, t);
        else
            zt = trace(end, t);
        end
        CPt = CPt + cp(zt);

        if t < T
            pi = Piupdate(zt,Beta{at},Alpha,lambda,pi);
            at = D(pi2idx(pi, quantnum));
            CIt = CIt + ci(at);
        end
    end
    Ct(trace_idx) = CIt + CPt;
end
Ct_QCDs5(1,lambda_idx)=mean(Ct);
Ct_QCDs5(2,lambda_idx)=std(Ct)/sqrt(trace_num)*2;


%% opt special case: time-varying
[V_init, D_init, D] = Cpt_opt(T, quantnum, cp, ci, Beta, Alpha, lambda);
Ct_opts_v(lambda_idx) = V_init;
Ct = zeros(trace_num,1);
for trace_idx = 1:trace_num
    trace = traces{trace_idx};
    changepoint = changepoints(trace_idx);

    pi = 0;
    at = D_init;
    CPt = 0;
    CIt = ci(at);
    for t = 1:T
        if t >= changepoint
            zt = trace(at, t);
        else
            zt = trace(end, t);
        end
        CPt = CPt + cp(zt);

        if t < T
            pi = Piupdate(zt,Beta{at},Alpha,lambda,pi);
            at = D(pi2idx(pi, quantnum), t);
            CIt = CIt + ci(at);
        end
    end
    Ct(trace_idx) = CIt + CPt;
end
Ct_opts(1,lambda_idx)=mean(Ct);
Ct_opts(2,lambda_idx)=std(Ct)/sqrt(trace_num)*2;
end

%%
figure;
hold on;
grid on;
errorbar(lambdas,Ct_opts(1,:),Ct_opts(2,:));
errorbar(lambdas,Ct_uppers(1,:),Ct_uppers(2,:));
%errorbar(lambdas,Ct_QCDs1(1,:),Ct_QCDs1(2,:));
%errorbar(lambdas,Ct_QCDs2(1,:),Ct_QCDs2(2,:));
%errorbar(lambdas,Ct_QCDs3(1,:),Ct_QCDs3(2,:));
%errorbar(lambdas,Ct_QCDs4(1,:),Ct_QCDs4(2,:));
errorbar(lambdas,Ct_QCDs5(1,:),Ct_QCDs5(2,:));
set(gca,'Xscale','log')
legend('Optimal','Low-complexity',...
    'QCD FA=0.005','QCD FA=0.01','QCD FA=0.05','QCD FA=0.1','QCD FA=0.5')
legend('Location','northwest')
ylabel('Total Cost','Interpreter','Latex');
xlabel('Change-point $\lambda$','Interpreter','Latex');
