function [mu,SIG]=update_moments(mu,SIG,x,n,xi2,c2,e)
% update_moments -- updates the moments for mcmc algorithms
%
% ::
%
%
%   [mu,SIG]=update_moments(mu,SIG,x,n,xi2)
%
%   [mu,SIG]=update_moments(mu,SIG,x,n,xi2,c2)
%
%   [mu,SIG]=update_moments(mu,SIG,x,n,xi2,c2,e)
%
% Args:
%
%    - **mu** [vector]: mean
%
%    - **SIG** [matrix]: covariance
%
%    - **x** [vector|matrix]: new parameter vectors
%
%    - **n** [integer]: iteration
%
%    - **xi2** [numeric]: appears in gam2=c2*(n+1)^(-xi2). Must be in (0.5,1)
%
%    - **c2** [numeric|{1}]: appears in gam2=c2*(n+1)^(-xi2)
%
%    - **e** [numeric]: minimum eigenvalue of the updated covariance matrix
%
% Returns:
%    :
%
%    - **mu** [vector]: updated mean
%
%    - **SIG** [matrix]: updated covariance matrix
%
% Note:
%
% Example:
%
%    See also:

% Applies formulae 23 and 24 in Blazej Miasojedow, Eric Moulines and Matti
% Vihola (2012): "Adaptive Parallel Tempering Algorithm"
if nargin<7
    
    e=sqrt(eps);
    
    if nargin<6
        
        c2=1;
        
    end
    
end

if c2<=0||c2>1
    
    error('c2 must be in (0,1]')
    
end

if xi2<=0.5||xi2>=1
    
    error('xi2 must be in (0.5,1)')
    
end

L=size(x,2);

L_xbar=sum(x,2);

L_C=0;

for icol=1:L
    
    xm=x(:,icol)-mu;
    
   L_C=L_C+(xm*xm.'); 
   
end

gam2=c2*(n+1)^(-xi2);

SIG=(1-gam2)*SIG+gam2/L*L_C;

SIG = utils.cov.project(SIG,e);

mu=(1-gam2)*mu+gam2/L*L_xbar;

end
