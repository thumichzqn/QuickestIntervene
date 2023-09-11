clear all;
close all;

%% Parameter Setting
rate_set = [1,2,3,4,5];             % Low rate, medium rate, high rate.
action_set = [1,2,3,4];           % doing nothing, Ban high rate, ban medium rate
service_cost_set = [0,1,2,3,4];     % service cost for low, medium and high rate
payment_loss_set = [0,0.02,0.06,0.2];     % payment loss for low, medium and high rate
rate_num = length(rate_set);
action_num = length(action_set);
% rate_dist_viral = [0.4410, 0.4903, 0.0687];     % rate select prob
% rate_dist= [0.5558, 0.3848, 0.0594];            % rate select prob when going viral

quantnum = 1000;
Beta=cell(4,1);

delta=0.02;
Alpha=[0.2,0.2,0.2,0.2,0.2];
Beta{4}=Alpha;
Beta{3}=[0.2-2*delta,   0.2-delta,  0.2,    0.2+delta,  0.2+2*delta];
Beta{2}=[0.2-2*2*delta, 0.2-2*delta,0.2,    0.2+2*delta,0.2+2*2*delta];
Beta{1}=[0.2-2*3*delta, 0.2-3*delta,0.2,    0.2+3*delta,0.2+2*3*delta];

T = 100;
rho = 1-1/(1*T);

%% Trace Creation
trace_num = 10000;
traces = cell(trace_num,1);
for trace_idx = 1:trace_num
    traces{trace_idx} = zeros(action_num,T);
    for action_idx = 1:action_num
        traces{trace_idx}(action_idx,:) = ...
            randsample(rate_num,T,true,Beta{action_idx});
    end
end
pointnum = 10;
viral_prob_set = logspace(-3,-0.1,pointnum);

Service_opt_lambda=zeros(trace_num,pointnum);
Payment_opt_lambda=zeros(trace_num,pointnum);
Total_opt_lambda = zeros(trace_num,pointnum);

Service_Upper_lambda=zeros(trace_num,pointnum);
Payment_Upper_lambda=zeros(trace_num,pointnum);
Total_Upper_lambda = zeros(trace_num,pointnum);

Service_QCD_lambda=zeros(trace_num,pointnum);
Payment_QCD_lambda=zeros(trace_num,pointnum);
Total_QCD_lambda = zeros(trace_num,pointnum);

Service_LB_lambda=zeros(trace_num,pointnum);
Payment_LB_lambda=zeros(trace_num,pointnum);
Total_LB_lambda = zeros(trace_num,pointnum);
%%
for lambdaidx = 1:pointnum
    
    viral_prob = viral_prob_set(lambdaidx);
    geo_pdf = zeros(T,1);
    for t=1:T
        geo_pdf(t) = geopdf(t-1,viral_prob);
    end
    rest = max(1 - sum(geo_pdf),0);
    geo_pdf(T) = geo_pdf(T)+rest;
    viral_points = zeros(trace_num,1);
    for trace_idx = 1:trace_num
        viral_points(trace_idx) = randsample(T,1,true,geo_pdf');
    end

    %% Policy Computation

    [D_opt,D_DQCD,D_Upper,D_approx,V_opt,V_DQCD,V_Upper,V_approx]...
    = Cpt_Decision(service_cost_set,payment_loss_set,rho,viral_prob,...
    Alpha,Beta,quantnum,T);

    %% Optimal Cost By MDP
    service_cost_traces = zeros(trace_num,T);
    payment_loss_traces = zeros(trace_num,T);
    pi_traces = zeros(trace_num,T);
    service_act_traces = zeros(trace_num,T);
    for trace_idx = 1:trace_num
        viral_time = viral_points(trace_idx);
        State_trace = traces{trace_idx};

        user_state = 0;                 % whether the user start to misbehave
        service_act = 1;
        service_cost = 0;
        payment_loss = 0;
        pi = 0;

        for t = 1: T
            %% Viral
            if user_state == 0
                if t >= viral_time
                    user_state = 1;
                end
            end
            %% Rate Selection
            if user_state == 0
                rate_tochoose_prob = Alpha;
                user_rate = State_trace(end,t);
            else
                rate_tochoose_prob = Beta{service_act};
                user_rate = State_trace(service_act,t);
            end

            service_cost = service_cost + service_cost_set(user_rate);
            payment_loss = payment_loss + payment_loss_set(service_act);

            %% Action Selection
            pi = Piupdate(user_rate,rate_tochoose_prob,Alpha,viral_prob,pi);    % posterior update
            action_idx = pi2idx(pi,quantnum);
            service_act = D_opt(action_idx,t);

            %% Recoding
            service_cost_traces(trace_idx,t) = service_cost;
            payment_loss_traces(trace_idx,t) = payment_loss;
            pi_traces(trace_idx,t) = pi;
            service_act_traces(trace_idx,t) = service_act;

        end
    end
    Service_opt_traces = service_cost_traces;
    Payment_opt_traces = payment_loss_traces;
    Pi_opt_traces = pi_traces;
    Act_opt_traces = service_act_traces;

    %% Cost By LowComplex
    service_cost_traces = zeros(trace_num,T);
    payment_loss_traces = zeros(trace_num,T);
    pi_traces = zeros(trace_num,T);
    service_act_traces = zeros(trace_num,T);
    for trace_idx = 1:trace_num
        viral_time = viral_points(trace_idx);
        State_trace = traces{trace_idx};

        user_state = 0;                 % whether the user start to misbehave
        service_act = 1;
        service_cost = 0;
        payment_loss = 0;
        pi = 0;

        for t = 1: T
            %% Viral
            if user_state == 0
                if t >= viral_time
                    user_state = 1;
                end
            end
            %% Rate Selection
            if user_state == 0
                rate_tochoose_prob = Alpha;
                user_rate = State_trace(end,t);
            else
                rate_tochoose_prob = Beta{service_act};
                user_rate = State_trace(service_act,t);
            end

            service_cost = service_cost + service_cost_set(user_rate);
            payment_loss = payment_loss + payment_loss_set(service_act);

            %% Action Selection
            pi = Piupdate(user_rate,rate_tochoose_prob,Alpha,viral_prob,pi);    % posterior update
            action_idx = pi2idx(pi,quantnum);
            service_act = D_Upper(action_idx);

            %% Recoding
            service_cost_traces(trace_idx,t) = service_cost;
            payment_loss_traces(trace_idx,t) = payment_loss;
            pi_traces(trace_idx,t) = pi;
            service_act_traces(trace_idx,t) = service_act;

        end
    end
    Service_Upper_traces = service_cost_traces;
    Payment_Upper_traces = payment_loss_traces;
    Pi_Upper_traces = pi_traces;
    Act_Upper_traces = service_act_traces;


    %% Cost By QCD
    service_cost_traces = zeros(trace_num,T);
    payment_loss_traces = zeros(trace_num,T);
    pi_traces = zeros(trace_num,T);
    service_act_traces = zeros(trace_num,T);
    for trace_idx = 1:trace_num
        viral_time = viral_points(trace_idx);
        State_trace = traces{trace_idx};

        user_state = 0;                 % whether the user start to misbehave
        service_act = 1;
        service_cost = 0;
        payment_loss = 0;
        pi = 0;

        for t = 1: T
            %% Viral
            if user_state == 0
                if t >= viral_time
                    user_state = 1;
                end
            end
            %% Rate Selection
            if user_state == 0
                rate_tochoose_prob = Alpha;
                user_rate = State_trace(end,t);
            else
                rate_tochoose_prob = Beta{service_act};
                user_rate = State_trace(service_act,t);
            end

            service_cost = service_cost + service_cost_set(user_rate);
            payment_loss = payment_loss + payment_loss_set(service_act);

            %% Action Selection
            pi = Piupdate(user_rate,rate_tochoose_prob,Alpha,viral_prob,pi);    % posterior update
            action_idx = pi2idx(pi,quantnum);
            service_act = D_DQCD(action_idx);

            %% Recoding
            service_cost_traces(trace_idx,t) = service_cost;
            payment_loss_traces(trace_idx,t) = payment_loss;
            pi_traces(trace_idx,t) = pi;
            service_act_traces(trace_idx,t) = service_act;

        end
    end
    Service_QCD_traces = service_cost_traces;
    Payment_QCD_traces = payment_loss_traces;
    Pi_QCD_traces = pi_traces;
    Act_QCD_traces = service_act_traces;



    %% Lower_Bound
    service_cost_traces = zeros(trace_num,T);
    payment_loss_traces = zeros(trace_num,T);
    for trace_idx = 1:trace_num
        viral_time = viral_points(trace_idx);
        State_trace = traces{trace_idx};
        service_cost = 0;
        payment_loss = 0;
        for t=1:T
            user_rate = State_trace(end,t);
            if t<viral_time
                service_act = 1;
                service_cost = service_cost + service_cost_set(user_rate);
                payment_loss = payment_loss + payment_loss_set(service_act);
            else
                service_act = action_num;
                service_cost = service_cost + service_cost_set(user_rate);
                payment_loss = payment_loss + payment_loss_set(service_act);
            end
            service_cost_traces(trace_idx,t) = service_cost;
            payment_loss_traces(trace_idx,t) = payment_loss;
        end
    end
    Service_LB_traces = service_cost_traces;
    Payment_LB_traces = payment_loss_traces;

    %%
    Service_opt_lambda(:,lambdaidx) = Service_opt_traces(:,T);
    Payment_opt_lambda(:,lambdaidx) = Payment_opt_traces(:,T);
    Total_opt_lambda(:,lambdaidx) = Service_opt_traces(:,T) + Payment_opt_traces(:,T);

    Service_Upper_lambda(:,lambdaidx) = Service_Upper_traces(:,T);
    Payment_Upper_lambda(:,lambdaidx) = Payment_Upper_traces(:,T);
    Total_Upper_lambda(:,lambdaidx) = Service_Upper_traces(:,T) + Payment_Upper_traces(:,T);

    Service_QCD_lambda(:,lambdaidx) = Service_QCD_traces(:,T);
    Payment_QCD_lambda(:,lambdaidx) = Payment_QCD_traces(:,T);
    Total_QCD_lambda(:,lambdaidx) = Service_QCD_traces(:,T) + Payment_QCD_traces(:,T);

    Service_LB_lambda(:,lambdaidx) = Service_LB_traces(:,T);
    Payment_LB_lambda(:,lambdaidx) = Payment_LB_traces(:,T);
    Total_LB_lambda(:,lambdaidx) = Service_LB_traces(:,T) + Payment_LB_traces(:,T);
    
end
save('data.mat','viral_prob_set',...
    'Service_opt_lambda','Payment_opt_lambda','Total_opt_lambda',...
    'Service_Upper_lambda','Payment_Upper_lambda','Total_Upper_lambda',...
    'Service_QCD_lambda','Payment_QCD_lambda','Total_QCD_lambda',...
    'Service_LB_lambda','Payment_LB_lambda','Total_LB_lambda')

%% functions
function Pi=Piupdate(zd,Beta,Alpha,lambda,pi)
    Pi=(pi+(1-pi)*lambda)*Beta(zd)/((pi+(1-pi)*lambda)*Beta(zd)+(1-pi)*(1-lambda)*Alpha(zd));
end
function pi=idx2pi(idx,quantnum)
    pi=(idx-0.5)/quantnum;
end
function idx=pi2idx(pi,quantnum)
    idx=ceil(pi*quantnum);
    if idx == 0
        idx=1;
    end
end
function dist=KL(Alpha,Beta)
    L=length(Alpha);
    dist=0;
    for idx=1:L
        dist=dist+Alpha(idx)*log(Alpha(idx)/Beta(idx));
    end
end



