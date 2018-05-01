function [alpha_t]=smoothing_step_classical(T,att,Ptt,Ptplus1,alpha_plus,a_plus)

% smoothing_step - smoothing step
%
% ::
%
%
%   [atT,etat,rlag]=smoothing_step(a,r,K,P,T,R,Z,iF,v)
%
% Args:
%
%    - **att** [vector] : m x 1 updated state
%
%    - **Ptt** [matrix] : m x m covariance matrix of updated state
%
%    - **Ptplus1** [matrix] : m x m covariance matrix of filtered state next
%    period
%
%    - **alpha_plus** [vector] : m x 1 smoothed state next period
%
%    - **a_plus** [vector] : m x 1 filtered state next period
%
% Returns:
%    :
%
%    - **alpha_t** [vector] : m x 1 smoothed state
%
% Note:
%
% Example:
%
%    See also:

alpha_t=att;

if ~isempty(Ptt)
    alpha_t=alpha_t+Ptt*T.'*pinv(Ptplus1)*(alpha_plus-a_plus);
%     alpha_t=alpha_t+Ptt*T.'*Ptplus1\(alpha_plus-a_plus);
end

end