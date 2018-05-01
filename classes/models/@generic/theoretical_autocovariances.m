function [V,retcode]=theoretical_autocovariances(obj,varargin)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
%    - The variances of the nonstationary variables are automatically set to
%      1. The nonstationary variables are detected by the program when option
%      "solve_bgp" is set to true.
%
% Example:
%
%    See also:

if isempty(obj)
    
    mydefaults=the_defaults();
    
    if nargout
        
        V=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

nobj=numel(obj);

if nobj>1
    
    V=cell(1,nobj);
    
    retcode=cell(1,nobj);
    
    for iobj=1:nobj
        
        [V{iobj},retcode{iobj}]=theoretical_autocovariances(obj(iobj),varargin{:});
   
    end
    
    return
    
end

obj=set(obj,varargin{:});

autocov_ar=obj.options.autocov_ar;

retcode=0;

V=[];

if obj.options.autocov_model_resolve
    
    [obj,retcode]=solve(obj);
    
end

if retcode
    
    return
    
end

n=obj.endogenous.number;

% 0- get the solution in alphabetical order
%------------------------------------------
[Tz,Re]=set_solution_to_companion(obj);

% 1- aggregate
%--------------
[T,R]=utils.miscellaneous.integrate_regimes(obj.solution.transition_matrices.Q,Tz,Re);

isgood=stationary_index(obj);
    
% 3- compute covariance of the state variables first
%----------------------------------------------------
[Vx,retcode]=lyapunov_equation(T,R*R',obj.options,isgood);

if ~retcode
    % recompute covariance of the whole vector conditional on the state
    % covariance. Should/could be used also during estimation. Need a way of
    % separating out the stationary and the nonstationary
    %-----------------------------------------------------------------------
    V=zeros(n,n,autocov_ar+1);
    
    for ii=1:autocov_ar+1
        
        V(:,:,ii)=Vx(1:n,1:n);
        
        Vx=T*Vx;
        
    end
    
end

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;

d={
    'autocov_ar',5,@(x)num_fin_int(x),...
    'autocov_ar must be a finite and positive integer'
    
    'autocov_model_resolve',true,@(x)islogical(x),...
    'autocov_model_resolve should be true or false'
    };

end

