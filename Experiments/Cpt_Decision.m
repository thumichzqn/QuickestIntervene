function [D_opt,D_DQCD,D_Upper,D_approx,V_opt,V_DQCD,V_Upper,V_approx]...
    = Cpt_Decision(cp,ci,rho,lambda,Alpha,Beta,quantnum,T)
%%
Z = length(cp);
A = length(ci);

%% Optimal Policy
V=zeros(quantnum,T+1);
Ja=cell(A,1);
for adx=1:A
    Ja{adx}=zeros(quantnum,T+1);
end
D=zeros(quantnum,T+1);
Vpre=ones(quantnum,T+1);
count=0;
%%
while norm(abs(V-Vpre))>10e-6
    Vpre=V;
    for t = T:-1:1
        for idx=1:quantnum
            pi=idx2pi(idx,quantnum);
            J=zeros(A,1);
            for adx=1:A
                J(adx)=ci(adx);
                for zddx=1:Z
                    pi_new=Piupdate(zddx,Beta{adx},Alpha,lambda,pi);
                    J(adx)=J(adx)+((1-pi)*(1-lambda)*Alpha(zddx)+(pi+(1-pi)*lambda)*Beta{adx}(zddx))...
                        *(Vpre(pi2idx(pi_new,quantnum),t+1)+cp(zddx));
                end
                Ja{adx}(idx,t)=J(adx);
            end 
            [V(idx,t),D(idx,t)]=min(J);
        end
    end
    count=count+1;
    if count==-1
        break;
    end
end
V_opt=V;
D_opt=D;

%% Direct QCD
VA=(ci(A)+rho*Alpha*cp')/(1-rho);
V=zeros(quantnum,1);
Ja=zeros(quantnum,1);
D=zeros(quantnum,1);

Vpre=ones(quantnum,1);
J=zeros(1,1);
count=0;
while norm(abs(V-Vpre))>10e-7
    Vpre=V;
    for idx=1:quantnum
        pi=idx2pi(idx,quantnum);
        J=ci(1);
        for zddx=1:Z
            pi_new=Piupdate(zddx,Beta{1},Alpha,lambda,pi);
            J=J+rho*((1-pi)*(1-lambda)*Alpha(zddx)+(pi+(1-pi)*lambda)*Beta{1}(zddx))...
                *(Vpre(pi2idx(pi_new,quantnum))+cp(zddx));
        end
        Ja(idx)=J;
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
V_DQCD=V;
D_DQCD=D;

%% Upperbound Policy
V=zeros(quantnum,1);
Ja=cell(A,1);
for adx=1:A
    Ja{adx}=zeros(quantnum,1);
end
D=zeros(quantnum,1);
for a=1:A-1
    pithres(a)=(ci(a+1)-ci(a))/((1-lambda)*rho*(Beta{a}*cp'-Beta{a+1}*cp'))-lambda/(1-lambda);
end
pithres(A)=1;
for a=1:A-1
    if pithres(a)<=0
        pithres(a)=1/quantnum;
    elseif pithres(a)>1
        pithres(a)=1;
    end
end
for a=1:A
    if a == 1
        D(1:pi2idx(pithres(a),quantnum))=a;
    else
        D(pi2idx(pithres(a-1),quantnum):pi2idx(pithres(a),quantnum))=a;
    end
end
Vpre=ones(quantnum,1);
count=0;
while norm(abs(V-Vpre))>10e-7
    Vpre=V;
    for idx=1:quantnum
        pi=idx2pi(idx,quantnum);
        J=zeros(A,1);
        for adx=1:A
            J(adx)=ci(adx);
            for zddx=1:Z
                pi_new=Piupdate(zddx,Beta{adx},Alpha,lambda,pi);
                J(adx)=J(adx)+rho*((1-pi)*(1-lambda)*Alpha(zddx)+(pi+(1-pi)*lambda)*Beta{adx}(zddx))...
                    *(Vpre(pi2idx(pi_new,quantnum))+cp(zddx));
            end
            Ja{adx}(idx)=J(adx);
        end        
        V(idx)=J(D(idx));
    end
    count=count+1;
    if count==-1
        break;
    end
end
V_Upper=V;
D_Upper=D;

%% Approx
V=zeros(quantnum,1);
Ja=cell(A,1);
for adx=1:A
    Ja{adx}=zeros(quantnum,1);
end
D=zeros(quantnum,1);

Vpre=ones(quantnum,1);
count=0;
while norm(abs(V-Vpre))>10e-6
    Vpre=V;
    for idx=1:quantnum
        pi=idx2pi(idx,quantnum);
        J=zeros(A,1);
        for adx=1:A
            J(adx)=ci(adx)+rho*Alpha*cp'+rho*(pi+(1-pi)*lambda)*(Beta{adx}-Alpha)*cp'...
                +rho*Vpre(pi2idx(pi+(1-pi)*lambda,quantnum));
            Ja{adx}(idx)=J(adx);
        end 
        [V(idx),D(idx)]=min(J);

    end
    count=count+1;
    if count==-1
        break;
    end
end
V_approx=V;
D_approx=D;
end



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