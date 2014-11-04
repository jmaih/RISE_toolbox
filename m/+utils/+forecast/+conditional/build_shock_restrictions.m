function [R,HHRestr,ncsp]=build_shock_restrictions(H,G,yrest_id,xrest_id,ncp,nap,hypo)
% build_shock_restrictions - build restriction conditions on the shocks
%
% Syntax
% -------
% ::
%
%   [R,HHRestr,ncsp]=build_shock_restrictions(H,G,yrest_id,xrest_id,ncp)
%   [R,HHRestr,ncsp]=build_shock_restrictions(H,G,yrest_id,xrest_id,ncp,nap)
%   [R,HHRestr,ncsp]=build_shock_restrictions(H,G,yrest_id,xrest_id,ncp,nap,hypo)
%
% Inputs
% -------
%
% - **H** [matrix]: Autoregressive part of state matrix with format ny x ny,
%   where ny is the number of endogenous variables. The assumed model is of
%   the form y_t=H*y_{t-1}+G*e_t  
%
% - **G** [matrix]: Shock impact, with format ny x nx x nh, where nx is the
%   number of exogenous variables and nh the horizon of anticipation +1. If
%   nh=1, there is no anticipation of future events. So, the number of
%   anticipated events is k=nh-1. The assumed model is of the form 
%   y_t=H*y_{t-1}+G*e_t 
%
% - **yrest_id** [vector|scalar]: location of the restricted endogenous in
%   the rows of H and G
%
% - **xrest_id** [vector|scalar]: location of the restricted shocks in
%   the columns of G
%
% - **ncp** [integer]: number of conditioning periods
%
% - **nap** [integer|size(G,3)]: number of anticipated periods
%
% - **hypo** [NCP|NAS|{JMA}]: forecasting hypothesis, determining the
%   number of periods over which future shocks will be drawn or assumed
%   known to the agents.
%   - **NCP** : under this scheme, the number of periods of future shocks 
%       is equal to the number of periods over which conditioning
%       information is available. Beyond the conditioning period, the
%       shocks are 0.
%   - **NAS** : under this scheme, the number of periods of future shocks 
%       is equal to the number of anticipated steps. Beyond that number,
%       the shocks are expected to be 0.
%   - **JMA** : under this scheme, the number of periods of future shocks 
%       is equal to the number of anticipated steps plus the number of
%       conditioning periods. It is assumed that every condition is
%       determined by the same number of future shocks
%
% Outputs
% --------
%
% - **R** [matrix]: Convolutions of H and G for the restricte endogenous
%   and convolution of future shocks
%
% - **HHRestr** [matrix]: convolutions of powers of H for the restricted
%   endogenous  
%
% - **ncsp** [integer]: number of periods over which future shocks are
%   computed 
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

narginchk(5,7)

if nargin<7
    hypo='JMA';
    if nargin<6
        nap=size(G,3);
    end
end

hypo=upper(hypo);
switch hypo
    case {'NCP'}
        ncsp=ncp;
    case {'NAS'}
        ncsp=nap;
        if ~isequal(ncp,nap)
            error([mfilename,':: for the NAS assumption, you need # anticipated steps = # conditioning periods'])
        end
    case {0,'JMA'}
        ncsp=nap+ncp;
    otherwise
        error([mfilename,':: Unknown option for the anticipation hypothesis'])
end

sizg=size(G);
endo_nbr=sizg(1);
exo_nbr=sizg(2);
RestEndo_nbr=numel(yrest_id);
cutoff=RestEndo_nbr*ncp;
R=zeros(cutoff,exo_nbr*ncsp);
RestShocks_nbr=numel(xrest_id);
S=zeros(RestShocks_nbr*ncp,exo_nbr*ncsp);

Hj=eye(endo_nbr);
HH=zeros(endo_nbr,endo_nbr,ncp);
HHRestr=nan(cutoff,endo_nbr);
for jj=1:ncp
    Hj=Hj*H;
    HH(:,:,jj)=Hj;
    HHRestr((jj-1)*RestEndo_nbr+1:jj*RestEndo_nbr,:)=HH(yrest_id,:,jj);
end

for k=1:ncp
    tmp=zeros(endo_nbr,exo_nbr*ncsp);
    for jj=1:ncsp
        PHI_k_j=zeros(endo_nbr,exo_nbr);
        for ii=1:jj
            sss=k-ii;
            if sss>0
                Hs=HH(:,:,sss);
            elseif sss==0
                Hs=1;
            else
                Hs=0;
            end
            ttt=jj-ii+1;
            test=ttt<=nap && sss>=0;
            if  test
                PHI_k_j=PHI_k_j+Hs*G(:,:,ttt);
            end
        end
        tmp(:,(jj-1)*exo_nbr+1:jj*exo_nbr)=PHI_k_j;
    end
    R((k-1)*RestEndo_nbr+(1:RestEndo_nbr),:)=tmp(yrest_id,:);
    S((k-1)*RestShocks_nbr+(1:RestShocks_nbr),(k-1)*exo_nbr+xrest_id)=eye(RestShocks_nbr);
end
clear tmp

R=[R;S];

end
