function [R,HHRestr,NumberOfCondShocksPeriods]=BuildShockRestrictions(H,G,yrest_id,xrest_id,ncp,NumberOfAnticipatedPeriods,Hypothesis)

Hypothesis=upper(Hypothesis);
switch Hypothesis
    case {'NCP','IWB'}
        NumberOfCondShocksPeriods=ncp;
    case {'NAS','IRIS'}
        NumberOfCondShocksPeriods=NumberOfAnticipatedPeriods;
        if ~isequal(ncp,NumberOfAnticipatedPeriods)
            error([mfilename,':: for the NAS or IRIS assumption, you need # anticipated steps = # conditioning periods'])
        end
    case {0,'JMA'}
        NumberOfCondShocksPeriods=NumberOfAnticipatedPeriods+ncp;
    otherwise
        error([mfilename,':: Unknown option for the anticipation hypothesis'])
end

sizg=size(G);
endo_nbr=sizg(1);
exo_nbr=sizg(2);
RestEndo_nbr=numel(yrest_id);
cutoff=RestEndo_nbr*ncp;
R=zeros(cutoff,exo_nbr*NumberOfCondShocksPeriods);
RestShocks_nbr=numel(xrest_id);
S=zeros(RestShocks_nbr*ncp,exo_nbr*NumberOfCondShocksPeriods);

Hj=eye(endo_nbr);
HH=zeros(endo_nbr,endo_nbr,ncp);
HHRestr=nan(cutoff,endo_nbr);
for j=1:ncp
    Hj=Hj*H;
    HH(:,:,j)=Hj;
    HHRestr((j-1)*RestEndo_nbr+1:j*RestEndo_nbr,:)=HH(yrest_id,:,j);
end

for k=1:ncp
    tmp=zeros(endo_nbr,exo_nbr*NumberOfCondShocksPeriods);
    for j=1:NumberOfCondShocksPeriods
        PHI_k_j=zeros(endo_nbr,exo_nbr);
        for i=1:j
            sss=k-i;
            if sss>0
                Hs=HH(:,:,sss);
            elseif sss==0
                Hs=1;
            else
                Hs=0;
            end
            ttt=j-i+1;
            test=ttt<=NumberOfAnticipatedPeriods && sss>=0;
            if  test
                PHI_k_j=PHI_k_j+Hs*G(:,:,ttt);
            end
        end
        tmp(:,(j-1)*exo_nbr+1:j*exo_nbr)=PHI_k_j;
    end
    R((k-1)*RestEndo_nbr+(1:RestEndo_nbr),:)=tmp(yrest_id,:);
    S((k-1)*RestShocks_nbr+(1:RestShocks_nbr),(k-1)*exo_nbr+xrest_id)=eye(RestShocks_nbr);
end
clear tmp

R=[R;S];


end
