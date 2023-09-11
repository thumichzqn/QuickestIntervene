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
