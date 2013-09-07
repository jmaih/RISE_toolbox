function [a,P,Tt,Rt,Record]=kalman_prediction(T,R,att,Ptt,MUt,OMGt,DPHI,DT,Record,ExpandedFlag)
if nargin==4
    bt=0;
    Tt=T;
    RR=R;
    Record=[];
    Rt=[];
elseif nargin==10
    [Tt,Rt,bt,~,Record]=ConditionalStateMatrices(T,R,MUt,OMGt,DPHI,DT,Record,ExpandedFlag);
    RR=Rt*Rt';
else
    error([mfilename,':: number of input arguments must be 4 or 10'])
end
% only compute the places where there is some action
test=true;
if test
    a=bt+Tt*att;
    P=Tt*Ptt*transpose(Tt)+RR;
else
    kk=any(Tt);
    a=bt+Tt(:,kk)*att(kk,:);
    P=Tt(:,kk)*Ptt(kk,kk)*transpose(Tt(:,kk))+RR;
end

% Make sure we remain symmetric
P=symmetrize(P);
end
