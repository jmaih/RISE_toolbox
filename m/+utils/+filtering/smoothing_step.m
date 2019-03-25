%  INTERNAL FUNCTION: smoothing step
% 
%  ::
% 
%    [atT,etat,rlag]=smoothing_step(a,r,K,P,T,R,Z,iF,v)
% 
%  Args:
% 
%     - **a** [vector] : m x 1 filtered state
%     - **r** [vector] : m x 1
%     - **K** [matrix] : m x p kalman gain
%     - **P** [matrix] : m x m covariance of filtered state
%     - **T** [matrix] : m x m state matrix of all endogenous
%     - **R** [matrix] : m x n shock impact matrix
%     - **Z** [square matrix] : selection matrix of observed into grand state
%     - **iF** [matrix] : p x p inverse of covariance matrix of forecast errors
%     - **v** [vector] : p x 1  forecast errors
% 
%  Returns:
%     :
% 
%     - **atT** [vector] : smoothed state
%     - **etat** [square matrix] : smoothed errors
%     - **rlag** [square matrix]
% 
%