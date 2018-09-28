function [T,itercode,retcode]=fix_point_iterator(iterate_func,T0,options,varargin)
% fix_point_iterator solves the fix point of a function
%
% ::
%
%   [T,itercode,retcode]=fix_point_iterator(iterate_func,T0,options,varargin)
%
% Args:
%    - iterate_func : [func_handle]: a function handle which returns 2
%      elements: The improved value of T and F0, the value of the function
%      whose discrepancy is to be minmized
%    - T0 : is the initial guess for the solution
%    - options : [struct] with fields
%      - fix_point_explosion_limit : [positive scalar |{1e+6}] : maximum
%        divergence of F0
%      - fix_point_TolFun : [positive scalar |{sqrt(eps)}] : tolerance
%        criterion
%      - fix_point_maxiter : [positive scalar |{1000}] : maximum number of
%        iterations
%      - fix_point_verbose : [true|{false}] : show iterations or not
%    - varargin
%
% Returns:
%    :
%    - T : final solution
%    - itercode : final number of iterations
%    - retcode : return code
%      - 0 : successful
%      - 21 : maximum number of iterations reached
%      - 22 : nan or inf in F0 or T
%      - 23 : divergence
%
% Note:
%
% Example:
%
%    See also:

mydefaults=the_defaults();
    
if nargin==0
    
    if nargout
        
        T=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

default_solve=disp_defaults(mydefaults);

options=utils.miscellaneous.setfield(default_solve,options);

% solving for T
F0=0.1*options.fix_point_explosion_limit;

iter=0;

conv_T=inf;

conv_F=F0;

T00=T0;

fix_point_verbose=options.fix_point_verbose;

while max(conv_T,conv_F)>options.fix_point_TolFun && ...
        iter<options.fix_point_maxiter && ...
        conv_F<options.fix_point_explosion_limit && ...
        utils.error.valid(F0)
    
    iter=iter+1;
    
    [T,F0] = iterate_func(T0,varargin{:});
    
    conv_T=max(abs(T(:)-T0(:)));
    
    conv_F=max(abs(F0(:)));
    
    if fix_point_verbose
        
        fprintf(1,'iter # %0.0f : conv(x)=%0.8f, conv(F(x))=%0.8f\n',iter,full(conv_T),conv_F); 
        
    end
    
    T0=T;
    
    if ~utils.error.valid(T)
        
        break
        
    end
    
end

retcode=0;

itercode=iter;

if iter>=options.fix_point_maxiter
    
    retcode=21;
    
elseif ~utils.error.valid(F0)||~utils.error.valid(T)
    
    retcode=22;
    
elseif conv_F>=options.fix_point_explosion_limit
    
    retcode=23;
    
end

if retcode==0 && ~utils.error.valid(T)
    
    assignin('base','T0',T00)
    
    assignin('base','iterate_func',iterate_func)
    
    assignin('base','options_',options)
    
    assignin('base','varargin_',varargin)
    
    disp(iter)
    
    keyboard
    
end

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;

d={
    'fix_point_TolFun(r)',sqrt(eps),@(x)num_fin(x)&& x>0,...
    'fix_point_TolFun must be a finite and positive scalar'
    
    'fix_point_maxiter',1000,@(x)num_fin_int(x),...
    'fix_point_maxiter must be a finite and positive integer'
    
    'fix_point_verbose(r)',false,@(x)islogical(x),...
    'fix_point_verbose must be a logical'
    
    'fix_point_explosion_limit(r)',1e+12,@(x)num_fin_int(x),...
    'fix_point_explosion_limit must be a finite and positive integer'
    };

end


