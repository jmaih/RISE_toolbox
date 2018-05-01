function log_c=update_scaling(log_c,accept_ratio,alpha_range,fixed_scaling,...
n,xi3,c3,c_range)
% update_scaling -- updates the log of the scaling parameter for mcmc
% algorithms
%
% ::
%
%
%   log_c=update_scaling(log_c,accept_ratio,alpha_range,fixed_scaling,n,xi3)
%
%   log_c=update_scaling(log_c,accept_ratio,alpha_range,fixed_scaling,n,xi3,c3)
%
%   log_c=update_scaling(log_c,accept_ratio,alpha_range,fixed_scaling,n,xi3,c3,c_range)
%
% Args:
%
%    - **log_c** [numeric]: initial scaling
%
%    - **accept_ratio** [numeric]: acceptance rate
%
%    - **alpha_range** [interval|{[]}]:
%
%    - **fixed_scaling** [true|false]:
%
%    - **n** [integer]: iteration
%
%    - **xi3** [numeric]: appears in gam3 = c3*(n+1)^(-xi3). must lie in
%    (0.5,1)
%
%    - **c3** [numeric|{1}]: appears in gam3 = c3*(n+1)^(-xi3);
%
%    - **c_range** [interval|{[sqrt(eps),100]}]: range of variation of the
%    scaling parameter (not its log!!!)
%
% Returns:
%    :
%
%    - **log_c** [numeric]: updated scaling
%
% Note:
%
% Example:
%
%    See also:

% Adapts formulae 20 in Blazej Miasojedow, Eric Moulines and Matti
% Vihola (2012): "Adaptive Parallel Tempering Algorithm"
if nargin<8
    
    c_range=[sqrt(eps),100];
    
    if nargin<7
        
        c3=[];
        
    end
    
end

if isempty(c3)
    
    c3=1;
    
end

if xi3<=0.5||xi3>=1
    
    error('xi3 must be in (0.5,1)')
    
end

gam3 = c3*(n+1)^(-xi3);

alpha_1=min(alpha_range);

alpha_2=max(alpha_range);

alpha=.5*(alpha_1+alpha_2);

alpha_diff = (~fixed_scaling).*(accept_ratio-alpha);

log_c = log_c + gam3*alpha_diff;

projection_facility()

    function projection_facility()
        
        if ~fixed_scaling
            
            log_c_range=log(c_range);
            
            if log_c<log_c_range(1)
                
                log_c=log_c_range(1);
                
            end
            
            if log_c>log_c_range(2)
                
                log_c=log_c_range(2);
                
            end
            
        end
        
    end

end
