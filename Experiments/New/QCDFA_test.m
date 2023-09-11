clear all; close all;

Z = 5;
A = 4;
T = 50;
rho = 1-1/(1*T);
lambda_num = 8;
lambdas = logspace(-2,-0.9,lambda_num);
quantnum=1000;
cp=[0,1,2,3,4];
ci=[0,0.04,0.12,0.2];

Beta=cell(A);


delta=0.01;
Alpha=[0.2,0.2,0.2,0.2,0.2];
Beta{4}=[delta,   delta,  1-4*delta,    delta,  delta];
Beta{3}=[delta,   delta,  0.5-3*delta,    0.5,  delta];
Beta{2}=[delta, delta, delta, 1-4*delta, delta];
Beta{1}=[delta,   delta,  delta,    delta,  1-4*delta];


lambda=lambdas(8);
fa_obj = 0.005;
%% QCD FA
threshold_l = 1/quantnum;
threshold_r = 1 - 1/quantnum;
fa_l = sim_fa(threshold_l, quantnum, lambda, Alpha, Beta, Z, A, T);
fa_r = sim_fa(threshold_r, quantnum, lambda, Alpha, Beta, Z, A, T);
for idx = 1:100
    threshold_mid = (threshold_r + threshold_l) / 2;
    fa_mid = sim_fa(threshold_mid, quantnum, lambda, Alpha, Beta, Z, A, T);
    if fa_mid > fa_l || fa_mid < fa_r || abs(threshold_l-threshold_r) < 3/quantnum
        opt_thres = threshold_mid;
        break;
    else
        if fa_mid > fa_obj
            threshold_l = threshold_mid;
            fa_l = fa_mid;
        else
            threshold_r = threshold_mid;
            fa_r = fa_mid;
        end
    end
end
%%
sim_fa(opt_thres,quantnum, lambda, Alpha, Beta, Z, A, T)



%%
function revalue = sim_fa(threshold, quantnum, lambda, Alpha, Beta, Z, A, T)
    D = zeros(quantnum,1);
    thres_idx = pi2idx(threshold,quantnum);
    for idx = 1: thres_idx
        D(idx) = 1;
    end
    for idx = thres_idx : quantnum
        D(idx) = A;
    end

    trace_num = 10000;
    % Observation
    traces = cell(trace_num,1);
    for trace_idx = 1:trace_num
        traces{trace_idx} = zeros(A,T);
        for adx = 1:A
            traces{trace_idx}(adx,:) = randsample(Z,T,true,Beta{adx});
        end
    end

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

    First_times = zeros(trace_num,1);
    for trace_idx = 1:trace_num
        trace = traces{trace_idx};
        changepoint = changepoints(trace_idx);

        pi = 0;
        at = 1;
        first_time = T+1;
        for t = 1:T
            if t >= changepoint
                zt = trace(at, t);
            else
                zt = trace(end, t);
            end
            if t < T
                pi = Piupdate(zt,Beta{at},Alpha,lambda,pi);
                at_pre = at;
                at = D(pi2idx(pi, quantnum));
                if at_pre == 1 && at == A
                    first_time = t;
                end
            end
        end
        First_times(trace_idx) = first_time;
    end
    Falsealarms = (First_times<changepoints);
    revalue = mean(Falsealarms);

end




