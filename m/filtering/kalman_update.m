function [att,Ptt,K]=kalman_update(a,P,iF,v,obs_id)
%         K=P(:,obs_id);
%         att=a+K*iF*v;
%         Ptt=P-K*iF*transpose(K);
K=P(:,obs_id)*iF;
att=a+K*v;
Ptt=P-K*transpose(P(:,obs_id));
% Ptt=symmetrize(Ptt);
end
