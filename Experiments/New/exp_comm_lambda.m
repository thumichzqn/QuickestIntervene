clear all; close all;

Z = 5;
A = 4;
T = 50;
rho = 1-1/(1*T);
lambda_num = 8;
lambdas = logspace(-2.2,-0.9,lambda_num);
quantnum=1000;
cp=[0,1,2,3,4];
ci=[0,0.02,0.06,0.2];

Beta=cell(A);


delta=0.02;
Alpha=[0.2,0.2,0.2,0.2,0.2];
Beta{4}=Alpha;
Beta{3}=[0.2-2*delta,   0.2-delta,  0.2,    0.2+delta,  0.2+2*delta];
Beta{2}=[0.2-2*2*delta, 0.2-2*delta,0.2,    0.2+2*delta,0.2+2*2*delta];
Beta{1}=[0.2-2*3*delta, 0.2-3*delta,0.2,    0.2+3*delta,0.2+2*3*delta];

%% Simulate Traces
trace_num = 20000;

% Observation
traces = cell(trace_num,1);
for trace_idx = 1:trace_num
    traces{trace_idx} = zeros(A,T);
    for adx = 1:A
        traces{trace_idx}(adx,:) = randsample(Z,T,true,Beta{adx});
    end
end

%% Simulate Change point
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
fa_set = [0.632,0.729,0.810,0.874,0.925,0.954,0.975,0.99];
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
fa_set = [0.549,0.652,0.751,0.835,0.891,0.932,0.96,0.98];
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
fa_set = [0.365,0.472,0.578,0.685,0.776,0.848,0.892,0.93];
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
fa_set = [0.28,0.373,0.481,0.587,0.692,0.768,0.830,0.870];
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
fa_set = [0.113,0.161,0.220,0.290,0.360,0.420,0.435,0.440];
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


%% QCD
[D_init, D] = Cpt_QCD(T, quantnum, cp, ci, Beta, Alpha, lambda);
V_init = Sim_BackDP(T, quantnum, cp, ci, Beta, Alpha, lambda, D, D_init);
Ct_QCDs_v(lambda_idx) = V_init;
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
Ct_QCDs(1,lambda_idx)=mean(Ct);
Ct_QCDs(2,lambda_idx)=std(Ct)/sqrt(trace_num)*2;

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
errorbar(lambdas,Ct_opts(1,:),Ct_opts(2,:),'LineWidth',1);
hold on;
grid on;


errorbar(lambdas,Ct_uppers(1,:),Ct_uppers(2,:),'LineWidth',1);
errorbar(lambdas,Ct_QCDs(1,:),Ct_QCDs(2,:),'LineWidth',1);
errorbar(lambdas,Ct_QCDs1(1,:),Ct_QCDs1(2,:),'LineWidth',1);
errorbar(lambdas,Ct_QCDs2(1,:),Ct_QCDs2(2,:),'LineWidth',1);
errorbar(lambdas,Ct_QCDs3(1,:),Ct_QCDs3(2,:),'LineWidth',1);
errorbar(lambdas,Ct_QCDs4(1,:),Ct_QCDs4(2,:),'LineWidth',1);
errorbar(lambdas,Ct_QCDs5(1,:),Ct_QCDs5(2,:),'LineWidth',1);
set(gca,'Xscale','log')
legend('Optimal','Low-complexity','AdaQCD',...
    'QCD FA=0.005','QCD FA=0.01','QCD FA=0.05','QCD FA=0.1','QCD FA=0.5')
legend('Location','northwest')
ylabel('Total Cost','Interpreter','Latex');
xlabel('Change-point $\lambda$','Interpreter','Latex');
magnify
