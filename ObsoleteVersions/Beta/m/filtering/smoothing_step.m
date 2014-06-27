function [alphat,etat,rlag]=smoothing_step(at,rt,Kt,Pt,Tt,Rt,Zt,iFt,vt)
Lt=Tt-Tt*Kt*Zt;
% Note that Durbin and Koopman define K=T*P*Z'*iF, while here it is defined
% as K=P*Z'*iF. Hence, in the definition of Lt, I have to premultiply K by
% T
rlag=Zt'*iFt*vt+Lt'*rt;
alphat=at+Pt*rlag;
% Q=eye(exo_nbr) in this algorithm...
etat=Rt'*rlag; % <--- etat=Rt'*rt;
% The state equation in Durbin and Koopman is
% a_{t+1}=T_{t}*a_{t}+R_{t}*eta_{t}, whereas the state equation in our
% models is a_{t}=T_{t}*a_{t-1}+R_{t}*eta_{t}, and this explains the change
% I made to the shock smoothing equation. With this change, if we define an
% auxiliary variable to be equal to a shock, we should retrieve the same
% values for both in the smoothing. Without the change, there will be a
% "delay"
end
