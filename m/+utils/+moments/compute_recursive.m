%  INTERNAL FUNCTION: Recursively computes moments
% 
%  ::
% 
%     [m,V]=compute_recursive(X)
% 
%  Args:
%     X (matrix of doubles): Data points of the parameters [Xn_{ij}] where
%       i denote parameters and j denote data sample
% 
%  Returns:
%     :
% 
%         - m (vector of doubles): mean
%         - V (matrix of doubles): variance
% 
%  Example:
% 
%     ::
% 
%        npar=11;
%        ndraws=3000;
%        Params=rand(npar,ndraws);
% 
%        % direct calculation
%        %-------------------
%        [m00,V00]=utils.moments.recursive(m,V,Params,1);
% 
%        % one at a time
%        %--------------
%        [m,V]=utils.moments.compute_recursive(Params)
% 
%        max(max(abs(V-cov(Params',1))))
%        max(max(abs(V00-cov(Params',1))))
% 
%  See also : utils.moments.recursive
%