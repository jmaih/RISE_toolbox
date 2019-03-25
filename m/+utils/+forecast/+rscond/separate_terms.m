%  INTERNAL FUNCTION: separate terms for conditional forecasting
% 
%  ::
% 
%    [Ty,Te,Tsig,C]=SEPARATE_TERMS(T,sstate,ny,nx,nshocks,k,h)
% 
%  Args:
% 
%     - **T** [cell array]: solution impact
%     - **sstate** [cell array]: steady state
%     - **state_cols** [vector]: location of state variables excluding shocks
%     - **k** [numeric]: number of forward steps (anticipation)
%     - **nshocks** [numeric]: number of shocks
% 
%  Returns:
%     :
% 
%     - **Ty** [cell array]: includes square matrices for the impact of
%       endogenous variables
%     - **Te** [cell array]: includes matrices for the impact of shocks
%     - **Tsig** [cell array]: includes vectors summing the impact of
%       uncertainty and the trend
%     - **C** [cell array]: includes vectors summing the impact of
%       uncertainty, the trend and the steady state
% 
%