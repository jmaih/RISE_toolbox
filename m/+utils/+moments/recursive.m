function [m1,V1,n]=recursive(m0,V0,Xn,n)
% INTERNAL FUNCTION: Recursively update moments
%
% ::
%
%    [m1,V1,n] = recursive(m0,V0,Xn,n);
%
% Args:
%    m0 (vector of doubles): Mean of the distribution so far
%    V0 (matrix of doubles): Variance of the distribution so far
%    Xn (matrix of doubles): Data points of the parameters [Xn_{ij}] where i denote parameters and j denote data sample
%    n (integer): Weight of the mean/variance approximation so far, i.e., the number of sample points used in computing the mean/variance prior to the function call
%
% Returns:
%    :
%
%        - m1 (vector of doubles): updated mean
%        - V1 (matrix of doubles): updated variance
%        - n (integer): updated weight of the distribution
%
% Example:
%
%    ::
%
%       npar=11;
%       ndraws=3000;
%       Params=rand(npar,ndraws);
%       m = 0;
%       V = 0;
%
%       % direct calculation
%       %-------------------
%       [m00,V00]=utils.moments.recursive(m,V,Params,1);
%
%       % one at a time
%       %--------------
%       for ii=1:ndraws
%           [m,V]=utils.moments.recursive(m,V,Params(:,ii),ii);
%       end
%       max(max(abs(V-cov(Params',1))))
%       max(max(abs(V00-cov(Params',1))))
%

if isempty(n)
    n=1;
end

m1=m0;
V1=V0;
covariance_also=nargout>1;
ncols=size(Xn,2);
for iup=1:ncols
    Xi=Xn(:,iup);
    do_one_update(Xi,n)
    if iup<ncols
        n=n+1;
        m0=m1;
        if covariance_also
            V0=V1;
        end
    end
end
    function do_one_update(Xi,n)
        m1=1/n*(Xi+(n-1)*m0);
        if covariance_also
            V1=(n-1)/n*(V0+m0*m0')+1/n*(Xi*Xi')-m1*m1';
        end
    end
end