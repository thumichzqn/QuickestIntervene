clear all; close all;

Z = 5;
A = 4;
rho = 0.95;
lambda=0.1;
quantnum=1000;
cp=[0,1,2,3,4];
deltas=[0.03,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001];

%%
% KL(Alpha,Beta{1})


for idxdelta = 1:length(deltas)
    delta=deltas(idxdelta);
    
    ci=[0,0.02,0.06,0.2]/(0.02)*delta;
    Beta=cell(A);
    Alpha=[0.2,0.2,0.2,0.2,0.2];
    Beta{4}=Alpha;
    Beta{3}=[0.2-2*delta,   0.2-delta,  0.2,    0.2+delta,  0.2+2*delta];
    Beta{2}=[0.2-2*2*delta, 0.2-2*delta,0.2,    0.2+2*delta,0.2+2*2*delta];
    Beta{1}=[0.2-2*3*delta, 0.2-3*delta,0.2,    0.2+3*delta,0.2+2*3*delta];
    
    %% Optimal Policy
    V=zeros(quantnum,A);
    Ja=cell(A,1);
    for adx=1:A
        Ja{adx}=zeros(quantnum,1);
    end
    D=zeros(quantnum,A);

    Vpre=ones(quantnum,A);
    J=zeros(A,1);
    count=0;
    while norm(abs(V-Vpre))>10e-6
        Vpre=V;
        for idx=1:quantnum
            pi=idx2pi(idx,quantnum);
            J=zeros(A,1);
            for adx=1:A
                J(adx)=ci(adx);
                for zddx=1:Z
                    pi_new=Piupdate(zddx,Beta{adx},Alpha,lambda,pi);
                    J(adx)=J(adx)+rho*((1-pi)*(1-lambda)*Alpha(zddx)+(pi+(1-pi)*lambda)*Beta{adx}(zddx))...
                        *(Vpre(pi2idx(pi_new,quantnum),adx)+cp(zddx));
                end
                Ja{adx}(idx)=J(adx);
            end 
            for adx=1:A
                if adx==A
                    [V(idx,A),D(idx,A)]=min(J(A:end));
                    D(idx,A)=D(idx,A)+A-1;
                else
                    [V(idx,adx),D(idx,adx)]=min(J(adx:adx+1));
                    D(idx,adx)=D(idx,adx)+adx-1;
                end
            end

        end
        count=count+1;
        if count==-1
            break;
        end
    end
    V_prop=V;
    D_prop=D;


    %% QCD Policy
    V=zeros(quantnum,A);
    Ja=cell(A,1);
    for adx=1:A
        Ja{adx}=zeros(quantnum,1);
    end
    D=zeros(quantnum,A);

    Vpre=ones(quantnum,A);
    J=zeros(A,1);
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
                        *(Vpre(pi2idx(pi_new,quantnum),adx)+cp(zddx));
                end
                Ja{adx}(idx)=J(adx);
            end

            for adx=1:A
                if adx==A
                    [V(idx,A),D(idx,A)]=min(J(A:end));
                    D(idx,A)=D(idx,A)+A-1;
                elseif adx==1
                    [V(idx,adx),D(idx,adx)]=min(J(adx:adx+1));
                    D(idx,adx)=D(idx,adx)+adx-1;
                else
                    V(idx,adx)=J(adx+1);
                    D(idx,adx)=adx+1;
                end
            end

        end
        count=count+1;
        if count==-1
            break;
        end
    end
    V_QCD=V;
    D_QCD=D;
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
    V=zeros(quantnum,A);
    Ja=cell(A,1);
    for adx=1:A
        Ja{adx}=zeros(quantnum,1);
    end
    D=zeros(quantnum,A);
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
    for a=1:A-1
        D(1:pi2idx(pithres(a),quantnum)-1,a)=a;
        D(pi2idx(pithres(a),quantnum):end,a)=a+1;
    end
    D(:,A)=A;
    Vpre=ones(quantnum,A);
    J=zeros(A,1);
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
                        *(Vpre(pi2idx(pi_new,quantnum),adx)+cp(zddx));
                end
                Ja{adx}(idx)=J(adx);
            end           
            for adx=1:A
                V(idx,adx)=J(D(idx,adx));
            end
        end
        count=count+1;
        if count==-1
            break;
        end
    end
    V_Upper=V;
    D_Upper=D;
    

    %% Approx
    V=zeros(quantnum,A);
    Ja=cell(A,1);
    for adx=1:A
        Ja{adx}=zeros(quantnum,1);
    end
    D=zeros(quantnum,A);

    Vpre=ones(quantnum,A);
    J=zeros(A,1);
    count=0;
    while norm(abs(V-Vpre))>10e-9
        Vpre=V;
        for idx=1:quantnum
            pi=idx2pi(idx,quantnum);
            J=zeros(A,1);
            for adx=1:A
                J(adx)=ci(adx)+rho*Alpha*cp'+rho*(pi+(1-pi)*lambda)*(Beta{adx}-Alpha)*cp'...
                    +rho*Vpre(pi2idx(pi+(1-pi)*lambda,quantnum),adx);
                Ja{adx}(idx)=J(adx);
            end 
            for adx=1:A
                if adx==A
                    [V(idx,A),D(idx,A)]=min(J(A:end));
                    D(idx,A)=D(idx,A)+A-1;
                else
                    [V(idx,adx),D(idx,adx)]=min(J(adx:adx+1));
                    D(idx,adx)=D(idx,adx)+adx-1;
                end
            end

        end
        count=count+1;
        if count==-1
            break;
        end
    end
    V_approx=V;
    D_approx=D;

    thres_approx_idx=zeros(A-1,1);
    thres_approx=zeros(A-1,1);
    t_a=zeros(A-1,1);
    D_i=zeros(A-1,1);
    D_p=zeros(A-1,1);
    for a=1:A-1
        thres_approx_idx(a)=find(D_approx(:,a)==a+1,1);
        if isempty(thres_approx_idx(a))
            thres_approx_idx(a)=quantnum;
        end
        thres_approx(a)=idx2pi(thres_approx_idx(a),quantnum);
        t_a(a)=ceil(log(1-thres_approx(a))/log(1-lambda));
    end
    for a=2:A
        D_i(a-1)=ci(a)-ci(a-1);
        D_p(a-1)=(Beta{a}-Beta{a-1})*cp';
    end
    V_approxclose=1/(1-rho)*(rho.^t_a)'*D_i+rho/(1-rho)*Alpha*cp'...
        + lambda*rho/(1-rho)/(1-rho*(1-lambda))*((Beta{1}-Alpha)*cp')...
        +rho/(1-rho)*(rho.^t_a)'*D_p ...
        - rho*(1-lambda)/(1-rho*(1-lambda))*((rho*(1-lambda)).^t_a)'*D_p;
%%

    V_Oracle=rho*(Alpha*cp')/(1-rho)+((Beta{A}-Alpha)*cp'+ci(A))*(rho/(1-rho)+(2*lambda-1)*rho/(1-rho*(1-lambda)));
    Reg_approx(idxdelta)=V_approx(1,1)-V_Oracle;
    KLs(idxdelta)=KL(Alpha,Beta{1});
    Reg_prop(idxdelta)=V_prop(1,1)-V_Oracle;
    Reg_Upper(idxdelta)=V_Upper(1,1)-V_Oracle;
    Reg_QCD(idxdelta)=V_QCD(1,1)-V_Oracle;
    Reg_DQCD(idxdelta)=V_DQCD(1,1)-V_Oracle;
    RegProportion1(idxdelta)=(Reg_QCD(idxdelta)-Reg_prop(idxdelta))/Reg_QCD(idxdelta);
    RegProportion1Up(idxdelta)=(Reg_QCD(idxdelta)-Reg_Upper(idxdelta))/Reg_QCD(idxdelta);
    RegProportion2(idxdelta)=(Reg_DQCD(idxdelta)-Reg_prop(idxdelta))/Reg_DQCD(idxdelta);
    RegProportion2Up(idxdelta)=(Reg_DQCD(idxdelta)-Reg_Upper(idxdelta))/Reg_DQCD(idxdelta);
end
%%
figure;
semilogx(KLs,Reg_prop,'Linewidth',1.4,'Marker','o');
hold on;
grid on;
semilogx(KLs,Reg_Upper,'Linewidth',1.4,'Marker','^');
semilogx(KLs,Reg_QCD,'Linewidth',1.4,'Marker','*');
semilogx(KLs,Reg_DQCD,'Linewidth',1.4,'Marker','s');
semilogx(KLs,Reg_approx,'Linewidth',1.4,'Marker','d');
%set(gca,'YScale','log');
legend('Optimal','Low-complexity','QCD','DQCD','Approximated','Interpreter','Latex')
legend('Location','northwest')
xlabel('KL divergence $D(\alpha,\beta_0)$','Fontname','Times New Roman','Interpreter','Latex'),
ylabel('Regret','Interpreter','Latex'),
set(gca,'Ylim',[0,1.5])
%%

% figure;
% semilogx(KLs,RegProportion1,'Marker','o');
% hold on
% semilogx(KLs,RegProportion1Up,'Marker','*');
% legend('prop on QCD','Upper on QCD')
% xlabel('KL divergence'),
% ylabel('Regret Improvement'),
% set(gca,'Ylim',[0.21,0.23])
% 
% figure;
% semilogx(KLs,RegProportion2,'Marker','o');
% hold on
% semilogx(KLs,RegProportion2Up,'Marker','*');
% legend('prop on DQCD','Upper on DQCD')
% xlabel('KL divergence'),
% ylabel('Regret Improvement'),
% set(gca,'Ylim',[0.24,0.32])


%%
function Trans=Normalize(Trans)
    [m,n]=size(Trans);
    for idx1=1:m
        Trans(idx1,:)=Trans(idx1,:)/(sum(Trans(idx1,:)));
    end
end
function Pi=Piupdate(zd,Beta,Alpha,lambda,pi)
    Pi=(pi+(1-pi)*lambda)*Beta(zd)/((pi+(1-pi)*lambda)*Beta(zd)+(1-pi)*(1-lambda)*Alpha(zd));
end
function pi=idx2pi(idx,quantnum)
    pi=(idx-0.5)/quantnum;
end
function idx=pi2idx(pi,quantnum)
    idx=ceil(pi*quantnum);
end
function dist=KL(Alpha,Beta)
    L=length(Alpha);
    dist=0;
    for idx=1:L
        dist=dist+Alpha(idx)*log(Alpha(idx)/Beta(idx));
    end
end