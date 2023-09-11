clear all; close all;

Z = 2;
A = 2;
rho = 0.95;
cp=[0,1];
ci=[0,0.005];
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

Beta=cell(2);
Alpha=[0.5, 0.5];
Beta{1}=[0.49, 0.51];
Beta{2}=Alpha;
Alpha=Normalize(Alpha);
for idx=1:A
    Beta{idx}=Normalize(Beta{idx});
end
%%

V=zeros(quantnum,1);
Ja=cell(A);
for adx=1:A
    Ja{adx}=zeros(quantnum,1);
end
D=zeros(quantnum,1);

Vpre=ones(quantnum,1);
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
                    *Vpre(pi2idx(pi_new,quantnum));
            end
            Ja{adx}(idx)=J(adx);
        end           
        [V(idx),D(idx)]=min(J);
    end
    count=count+1;
    if count==-1
        break;
    end
end
%%
figure;
plot([0:1/quantnum:1-1/quantnum],V);
hold on;
a=1;

%%
figure;
hold on;
plot([0:1/quantnum:1-1/quantnum],Ja{1})
plot([0:1/quantnum:1-1/quantnum],Ja{2})
legend('1','2');
%%
figure;
hold on
plot([0:1/quantnum:1-1/quantnum],Ja{2}-Ja{1})
C=zeros(quantnum,1);
a=1;
for idx=1:quantnum
    pi=idx2pi(idx,quantnum);
    for z=1:Z
        C(idx)=C(idx)+pi*(Beta{a+1}(z)-Beta{a}(z))*cp(z); 
    end
end
%plot([0:1/quantnum:1-1/quantnum],C)
%%
figure;
imagesc(D);
colorbar;
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