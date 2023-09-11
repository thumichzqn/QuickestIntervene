clear all; close all;

Z = 2;
A = 5;
rho = 0.95;
cp=[0,1];
%ci=[0.05,0.1,0.15,0.2,0.25];
ci=[1:A]*0.005;
lambda=0.1;
quantnum=1000;

%% Transition Matrix
% Beta=cell(3);
% Alpha=[0.5, 0.5, 0.5 ; 0.5, 0.5, 0.5; 0.5, 0.5, 0.5];
% Beta{1}=[0.5, 0.5, 0.5; 0.3, 0.5, 0.7; 0.1, 0.5, 0.9];
% Beta{2}=[0.5, 0.5, 0.5; 0.4, 0.5, 0.6; 0.3, 0.5, 0.7];
% Beta{3}=Alpha;
% Alpha=Normalize(Alpha);
% for idx=1:A
%     Beta{idx}=Normalize(Beta{idx});
% end

Beta=cell(A);
Alpha=[0.5, 0.5];
Beta{1}=[0.46, 0.54];
Beta{2}=[0.47, 0.53];
Beta{3}=[0.48, 0.52];
Beta{4}=[0.49, 0.51];
Beta{5}=Alpha;
Alpha=Normalize(Alpha);
for idx=1:A
    Beta{idx}=Normalize(Beta{idx});
end
%%

V=zeros(quantnum,A);
Ja=cell(A,1);
for adx=1:A
    Ja{adx}=zeros(quantnum,1);
end
D=zeros(quantnum,A);

Vpre=ones(quantnum,A);
J=zeros(A,1);
count=0;
while norm(V-Vpre)>1e-3
    Vpre=V;
    for idx=1:quantnum
        pi=idx2pi(idx,quantnum);
        for adx=1:A
            J(adx)=ci(adx);
            for zddx=1:Z
                pi_new=Piupdate(zddx,Beta{adx},Alpha,lambda,pi);
                J(adx)=J(adx)+(pi*Beta{adx}(zddx)+(1-pi)*Alpha(zddx))*cp(zddx)...
                    +rho*((1-pi)*(1-lambda)*Alpha(zddx)+(pi+(1-pi)*lambda)*Beta{adx}(zddx))...
                    *Vpre(pi2idx(pi_new,quantnum),adx);
            end
            Ja{adx}(idx)=J(adx);
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
%%
figure;
imagesc(D);
colorbar;

figure;
for a=1:A-1
    subplot(1,A-1,a)
    hold on
    plot([0:1/quantnum:1-1/quantnum],Ja{a});
    plot([0:1/quantnum:1-1/quantnum],Ja{a+1});
    legend([num2str(a)],[num2str(a+1)])
end
%%
figure;
plot(Ja{4}-Ja{3})
%%
% a=2;
% for idx=1:quantnum
%     idxd1=pi2idx(Piupdate(1,Beta{a},Alpha,lambda,idx2pi(idx,quantnum)),quantnum);
%     idxd2=pi2idx(Piupdate(2,Beta{a},Alpha,lambda,idx2pi(idx,quantnum)),quantnum);
%     idxd11=pi2idx(Piupdate(1,Beta{a+1},Alpha,lambda,idx2pi(idx,quantnum)),quantnum);
%     idxd12=pi2idx(Piupdate(2,Beta{a+1},Alpha,lambda,idx2pi(idx,quantnum)),quantnum);
%     Difmin(idx)=min(Ja{a}(idxd2),Ja{a+1}(idxd2))-min(Ja{a}(idxd1),Ja{a+1}(idxd1))...
%         -min(Ja{a+1}(idxd12),Ja{a+2}(idxd12))+min(Ja{a+1}(idxd11),Ja{a+2}(idxd11));
% end
% figure;
% plot(Difmin)

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