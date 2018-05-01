function [V,retcode]=lyapunov_equation(T,Q,options,stationary)
% lyapunov_equation solves the equation V=T*V*T'+Q
%
% ::
%
%   [V,retcode]=lyapunov_equation(T,Q)
%   [V,retcode]=lyapunov_equation(T,Q,options)
%
% Args:
%    - T :
%    - Q :
%    - options :
%
% Returns:
%    :
%    - V :
%    - retcode :
%
% Note:
%
% Example:
%
%    See also:

% example
% we would like to compute the variance of the system
% y_t=ay*y_{t-1}+ax*x_{t-1}+ae*e_t
% x_t=by*y_{t-1}+bx*x_{t-1}+be*e_t
% the matrix representation is
% X_t=[ay ax;by bx]*X_{t-1}+[ae be]'*e_t;

mydefaults=the_defaults();

if nargin==0
    
    if nargout
        
        V=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

defaults=disp_defaults(mydefaults);

if nargin<4
    
    stationary=[];
    
    if nargin<3
        
        options=struct();
        
    end
    
end

n=size(T,1);

if isempty(stationary)
    
    stationary=true(1,n);
    
end

fnames=fieldnames(defaults);

for ii=1:numel(fnames)
    
    if ~isfield(options,fnames{ii})
        
        options.(fnames{ii})=defaults.(fnames{ii});
        
    end
    
end

algo=options.lyapunov_algo;

% variance of nonstationary terms
%--------------------------------
V=diag(options.lyapunov_diffuse_factor*ones(n,1));

V(stationary,stationary)=lyapunov_engine(T(stationary,stationary),Q(stationary,stationary));

    function v=lyapunov_engine(T,Q)
        
        % locate_state_variables
        %-----------------------
        t_pb=any(T,1);
        
        switch algo
            
            case 'schur'
                
                error('Schur algorithms are under revision')
                
            case 'doubling'
                
                [v,retcode]=doubling_solve(T(t_pb,t_pb),[],Q(t_pb,t_pb),options);
                
            otherwise
                
                error([mfilename,':: unknown lyapunov algorithm option ',algo])
                
        end
        
        v=T(:,t_pb)*v*T(:,t_pb).'+Q;
        
    end

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

d={
    'lyapunov_algo','doubling',@(x)ismember(x,{'doubling','schur'}),...
    'lyapunov_algo must be doubling or schur'
    
    'lyapunov_diffuse_factor',1,@(x)num_fin(x) && x>=0,...
    'lyapunov_diffuse_factor must be >=0'
    
    'fix_point_verbose(r)',false,@(x)islogical(x),'fix_point_verbose must be a logical'
    };

end

