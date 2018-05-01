function [a,P,Tt,Rt,Record]=prediction_step(T,R,att,Ptt,MUt,OMGt,DPHI,DT,Record,ExpandedFlag)
% prediction_step - prediction step for Kalman filter
%
% ::
%
%
%   [a,P,Tt,Rt,Record]=prediction_step(T,R,att,Ptt)
%   [a,P,Tt,Rt,Record]=prediction_step(T,R,att,Ptt,MUt,OMGt,DPHI,DT,Record,ExpandedFlag)
%
% Args:
%
%    - **T** [matrix] : m x m state matrix (autoregressive part)
%
%    - **R** [matrix] : m x n state matrix (shock impacts)
%
%    - **att** [vector] : m x 1 state update
%
%    - **Ptt** [matrix] : m x m covariance matrix of state update
%
%    - **MUt** [matrix] : k x l matrix with forward information to match
%
%    - **OMGt** [matrix|{[]}] : kl x kl covariance matrix of forward
%      information
%
%    - **DPHI** [matrix] : Matrix representing the restrictions on future
%      shocks
%
%    - **DT** [matrix] : Convoluted impact of initial conditions
%
%    - **Record** [] : Holder of invariant information
%
%    - **ExpandedFlag** [true|false] : if true, returns the expanded state
%      vector including the future shocks. If false, returns only the
%      endogenous variables
%
% Returns:
%    :
%
%    - **a** [vector] : m x 1 or mm x 1 vector of predictions
%
%    - **P** [matrix] : m x m or mm x mm covariance of predictions
%
%    - **Tt** [matrix] : m x m or mm x mm time-varying matrix, modifying the
%      impact of lagged endogenous (T)
%
%    - **Rt** [] : m x n or mm x nn matrix of impact of shocks
%
%    - **Record** [] : Holder of invariant information
%
% Note:
%
% Example:
%
%    See also:

narginchk(4,10);
reduced=nargin==4;
bt=0;
Tt=T;
Rt=R;
if ~reduced
    reduced=reduced||(~ExpandedFlag && all(isnan(MUt(:))));
    if ~reduced
        [Tt,Rt,bt,~,Record]=utils.forecast.rscond.state_matrices(T,R,MUt,OMGt,DPHI,DT,Record,ExpandedFlag);
    end
end

if reduced
    Record=[];
end

RR=Rt*Rt';
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
P=utils.cov.symmetrize(P);
end