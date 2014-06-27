function [T,R,SS,retcode]=dsge_solve_aim(Aplus,A0,Aminus,B,C,T_only)

if nargin<6
    T_only=false;
    if nargin<5
        CC=[];
        if nargin<4
            B=[];
        end
    end
end
retcode=1;
T=[];SS=[];R=[];

lags=1; % no of lags and leads
leads=1;
nn= size(A0,1); % no of equations
theAIM_H=[Aminus,A0,Aplus];

condn  = 1.e-10;%SPAmalg uses this in zero tests
uprbnd = 1 + 1.e-6;%allow unit roots

% [T,rts,ia,nexact,nnumeric,lgroots,aimcode] =...
[T,~,~,~,~,~,aimcode] =SPAmalg(theAIM_H,nn,lags,leads,condn,uprbnd);

retcode=22;
if(aimcode==1)
    e='Aim: unique solution.'; %#ok<*NASGU>
retcode=0;
elseif(aimcode==2)
    e='Aim: roots not correctly computed by real_schur.';
elseif(aimcode==3)
    e='Aim: too many big roots.';
elseif(aimcode==35)
    e='Aim: too many big roots, and q(:,right) is singular.';
elseif(aimcode==4)
    e='Aim: too few big roots.';
    retcode=21;
elseif(aimcode==45)
    e='Aim: too few big roots, and q(:,right) is singular.';
    retcode=21;
elseif(aimcode==5)
    e='Aim: q(:,right) is singular.';
elseif(aimcode==61)
    e='Aim: too many exact shiftrights.';
elseif(aimcode==62)
    e='Aim: too many numeric shiftrights.';
else
    e='Aimerr: return code not properly specified';
end

if retcode
    return
end
if T_only
    return
end
theAIM_Psi= - B;%

scof = SPObstruct(theAIM_H,T,nn,lags,leads);
scof1= scof(:,lags*nn+1:end);
R =scof1\theAIM_Psi;
SS=(eye(nn)-T)\C;

