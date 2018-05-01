function [atT,etat,rlag]=smoothing_step(a,r,K,P,T,R,Z,iF,v)
% smoothing_step - smoothing step
%
% ::
%
%
%   [atT,etat,rlag]=smoothing_step(a,r,K,P,T,R,Z,iF,v)
%
% Args:
%
%    - **a** [vector] : m x 1 filtered state
%
%    - **r** [vector] : m x 1
%
%    - **K** [matrix] : m x p kalman gain
%
%    - **P** [matrix] : m x m covariance of filtered state
%
%    - **T** [matrix] : m x m state matrix of all endogenous
%
%    - **R** [matrix] : m x n shock impact matrix
%
%    - **Z** [square matrix] : selection matrix of observed into grand state
%
%    - **iF** [matrix] : p x p inverse of covariance matrix of forecast errors
%
%    - **v** [vector] : p x 1  forecast errors
%
% Returns:
%    :
%
%    - **atT** [vector] : smoothed state
%
%    - **etat** [square matrix] : smoothed errors
%
%    - **rlag** [square matrix]
%
% Note:
%
% Example:
%
%    See also:

L=T-T*K*Z;
% Note that Durbin and Koopman define K=T*P*Z'*iF, while here it is defined
% as K=P*Z'*iF. Hence, in the definition of Lt, I have to premultiply K by
% T
rlag=Z'*iF*v+L'*r;
atT=a+P*rlag;
% Q=eye(exo_nbr) in this algorithm...
etat=R'*rlag; % <--- etat=Rt'*rt;
% The state equation in Durbin and Koopman is
% a_{t+1}=T_{t}*a_{t}+R_{t}*eta_{t}, whereas the state equation in our
% models is a_{t}=T_{t}*a_{t-1}+R_{t}*eta_{t}, and this explains the change
% I made to the shock smoothing equation. With this change, if we define an
% auxiliary variable to be equal to a shock, we should retrieve the same
% values for both in the smoothing. Without the change, there will be a
% "delay"
end