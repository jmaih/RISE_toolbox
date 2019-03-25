%  INTERNAL FUNCTION: initialize the storage of filtering
% 
%  ::
% 
%    [Filters,K_store,iF_store,v_store]=initialize_storage(a,P,PAI,Q,p0,...
%    exo_nbr,horizon,m,nsteps,smpl,h,store_filters)
% 
%  Args:
% 
%     - **a** [cell]: initial conditions of the filter for each regime
%     - **P** [cell]: initial covariance of the filter for each regime
%     - **PAI** [vector]: initial probability distribution of regimes
%     - **p0** [scalar]: number of observables
%     - **exo_nbr** [scalar]: number of exogenous
%     - **horizon** [{1}|scalar]: number of anticipated steps + 1
%     - **nsteps** [scalar]: number of forecast steps
%     - **smpl** [scalar]: number of observations
%     - **store_filters** [0|1|2|3]: 0 (no storage), 1(predicted only),
%       2(predicted and updated), 3(predicted, updated and smoothed)
% 
%  Returns:
%     :
% 
%     - **Filters** [struct]: structure with different fields
%     - **K_store** [cell]: Place holder for Kalman gains
%     - **iF_store** [cell]: Place holder for inverses of covariance matrices
%       of forecast errors
%     - **v_store** [cell]: Place holder for forecast errors
% 
%