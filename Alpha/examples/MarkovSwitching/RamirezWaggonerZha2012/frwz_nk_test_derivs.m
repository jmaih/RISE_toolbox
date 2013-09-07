clc
syms PAIp1 Rp1 Yp1 PAI R Y PAIm1 Rm1 Ym1 EPS_R betta eta kappa rhor sigr mu psi_ PAIss Rss

f1=sym('1-betta*(1-.5*kappa*(PAI-1)^2)*Y*R/((1-.5*kappa*(PAIp1-1)^2)*Yp1*exp(mu)*PAIp1)');

f2=sym('1-eta+eta*(1-.5*kappa*(PAI-1)^2)*Y+betta*kappa*(1-.5*kappa*(PAI-1)^2)*(PAIp1-1)*PAIp1/(1-.5*kappa*(PAIp1-1)^2)-kappa*(PAI-1)*PAI');

f3=('(Rm1/Rss)^rhor*(PAI/PAIss)^((1-rhor)*psi_)*exp(sigr*EPS_R)-R/Rss');

f=[f1;f2;f3];
%%
wrt=[PAIp1 Rp1 Yp1 PAI R Y PAIm1 Rm1 Ym1 EPS_R, mu psi_];
J=jacobian(f,wrt);
J=J.';

H=jacobian(J(:),wrt);
H=H.';

%%
Jfunc=matlabFunction(J);
Hfunc=matlabFunction(H);

%%
EPS_R=0;PAIss=1;PAI=PAIss;PAIp1=PAIss;Rss=1.030506;R=Rss;Rm1=Rss;
Yss=0.9;Y=Yss;Yp1=Yss;betta=0.99;eta=10;kappa=161;mu=[0.03 0.01];rhor=0.8;
sigr=0.0025;psi_=[3.1 0.9];
JJ=Jfunc(EPS_R,PAI,PAIp1,PAIss,R,Rm1,Rss,Y,Yp1,betta,eta,kappa,mu(1),psi_(1),rhor,sigr);
JJ(:,:,2)=Jfunc(EPS_R,PAI,PAIp1,PAIss,R,Rm1,Rss,Y,Yp1,betta,eta,kappa,mu(2),psi_(2),rhor,sigr);
% J2=reshape(J2,3,[]);
H1=Hfunc(EPS_R,PAI,PAIp1,PAIss,R,Rm1,Rss,Y,Yp1,betta,eta,kappa,mu(1),psi_(1),rhor,sigr);
HH=reshape(H1,[],3);
H2=Hfunc(EPS_R,PAI,PAIp1,PAIss,R,Rm1,Rss,Y,Yp1,betta,eta,kappa,mu(2),psi_(2),rhor,sigr);
HH(:,:,2)=reshape(H2,[],3);
%%
np=3;nc=3;nm=3;ne=1;nt=2;
n=[np,nc,nm,ne,nt];
labels={'p','c','m','e','t'};
longLabels='';
for ilab=1:numel(labels)
    longLabels=[longLabels,labels{ilab}*ones(1,n(ilab))];
end
%%
wrt_list={'PAI{+1}','R{+1}','Y{+1}','PAI','Y','R','PAI{-1}','R{-1}','Y{-1}','EPS_R','mu','psi'};
H_list={};
iter=0;
partitions=struct();
for irow=1:numel(wrt_list)
    first=longLabels(irow);
    iter=iter+1;
    if ~isfield(partitions,first)
        partitions.(first)=[];
    end
    partitions.(first)=[partitions.(first),iter];
end
iter=0;
for irow=1:numel(wrt_list)
    first=longLabels(irow);
    for icol=1:numel(wrt_list)
        iter=iter+1;
        H_list=[H_list,[wrt_list{irow},'&',wrt_list{icol}]];
        second=longLabels(icol);
        first_second=[first,second];
        if ~isfield(partitions,first_second)
            partitions.(first_second)=[];
        end
        partitions.(first_second)=[partitions.(first_second),iter];
    end
end
%%
[wrt_list',num2cell(JJ(:,:,1))]
[wrt_list',num2cell(JJ(:,:,2))]
[H_list',num2cell(HH(:,:,1))]
[H_list',num2cell(HH(:,:,2))]
%%
frwz=rise('frwz_nk');
frwzs=solve(frwz,'solve_order',2);
syst=get(frwzs,'structure');

fields=fieldnames(partitions);
%%
for ii=1:2
    index=ii;
    for ifield=1:numel(fields)
        disp([' ------------- ',fields{ifield},'(',int2str(index),') ------------- '])
        if size(fields{ifield},2)==1
            first=JJ(partitions.(fields{ifield}),:,index)';
        else
            first=HH(partitions.(fields{ifield}),:,index)';
        end
        second=full(syst.(['G',fields{ifield}]){index,index})/frwzs.solution.Q(index,index);
        disp(max(abs(first(:)-second(:))))
    end
end

