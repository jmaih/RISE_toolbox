%  INTERNAL FUNCTION: Prediction step for Kalman filter
% 
%  ::
% 
%    [a,P,Tt,Rt,Record]=prediction_step(T,R,att,Ptt)
%    [a,P,Tt,Rt,Record]=prediction_step(T,R,att,Ptt,MUt,OMGt,DPHI,DT,Record,ExpandedFlag)
% 
%  Args:
% 
%     - **T** [matrix] : m x m state matrix (autoregressive part)
%     - **R** [matrix] : m x n state matrix (shock impacts)
%     - **att** [vector] : m x 1 state update
%     - **Ptt** [matrix] : m x m covariance matrix of state update
%     - **MUt** [matrix] : k x l matrix with forward information to match
%     - **OMGt** [matrix|{[]}] : kl x kl covariance matrix of forward
%       information
%     - **DPHI** [matrix] : Matrix representing the restrictions on future
%       shocks
%     - **DT** [matrix] : Convoluted impact of initial conditions
%     - **Record** [] : Holder of invariant information
%     - **ExpandedFlag** [true|false] : if true, returns the expanded state
%       vector including the future shocks. If false, returns only the
%       endogenous variables
% 
%  Returns:
%     :
% 
%     - **a** [vector] : m x 1 or mm x 1 vector of predictions
%     - **P** [matrix] : m x m or mm x mm covariance of predictions
%     - **Tt** [matrix] : m x m or mm x mm time-varying matrix, modifying the
%       impact of lagged endogenous (T)
%     - **Rt** [] : m x n or mm x nn matrix of impact of shocks
%     - **Record** [] : Holder of invariant information
% 
%