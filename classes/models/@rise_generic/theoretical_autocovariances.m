function [V,retcode]=theoretical_autocovariances(obj,varargin)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    V=struct('autocov_ar',5);
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

[obj,retcode]=solve(obj);
V=[];
if retcode
    return
end

n=obj.endogenous.number(end);

V=zeros(n,n,autocov_ar+1);

% 0- get the solution in alphabetical order
%------------------------------------------
[Tz,Re]=set_solution_to_companion(obj);

% 1- aggregate
%--------------
[T,R]=utils.miscellaneous.integrate_regimes(obj.solution.transition_matrices.Q,Tz,Re);

% 2- locate state variables
%--------------------------
t_pb=any(T,1);
T=T(:,t_pb);

% 3- compute covariance of the state variables first
%----------------------------------------------------
[Vx,retcode]=lyapunov_equation(T(t_pb,:),R(t_pb,:)*R(t_pb,:)',obj.options);

if ~retcode
    % recompute covariance of the whole vector conditional on the state
    % covariance. Should/could be used also during estimation. Need a way of
    % separating out the stationary and the nonstationary
    %-----------------------------------------------------------------------
    for ii=1:autocov_ar+1
        if ii==1
            V0=T*Vx*T'+R*R';
        else
            V0=T*V0(t_pb,:);
        end
        V(:,:,ii)=V0(1:n,1:n);
    end
end
end

