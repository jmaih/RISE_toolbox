%  INTERNAL FUNCTION: Recursively update moments
% 
%  ::
% 
%     [m1,V1,n] = recursive(m0,V0,Xn,n);
%     [m1,V1,n] = recursive(m0,V0,Xn,n,n0);
% 
%  Args:
%     m0 (vector of doubles): Mean of the distribution so far
%     V0 (matrix of doubles): Variance of the distribution so far
%     Xn (matrix of doubles): Data points of the parameters [Xn_{ij}] where
%     i denote parameters and j denote data sample
%     n (integer): Weight of the mean/variance approximation so far, i.e.,
%     the number of sample points used in computing the mean/variance prior
%     to the function call
%     n0 (integer): past n, defaults to n-1 if not provided.
% 
%  Returns:
%     :
% 
%         - m1 (vector of doubles): updated mean
%         - V1 (matrix of doubles): updated variance
%         - n (integer): updated weight of the distribution
% 
%  Example:
% 
%     ::
% 
%        npar=11;
%        ndraws=3000;
%        Params=rand(npar,ndraws);
%        m = 0;
%        V = 0;
% 
%        % direct calculation
%        %-------------------
%        [m00,V00]=utils.moments.recursive(m,V,Params,1);
% 
%        % one at a time
%        %--------------
%        for ii=1:ndraws
%            [m,V]=utils.moments.recursive(m,V,Params(:,ii),ii);
%        end
%        max(max(abs(V-cov(Params',1))))
%        max(max(abs(V00-cov(Params',1))))
% 
%