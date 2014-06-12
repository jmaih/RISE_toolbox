function [V,retcode]=theoretical_autocovariances(obj,varargin)
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

reg_nbr=obj.markov_chains.regimes_number;
shock_horizon=max(obj.exogenous.shock_horizon);
if shock_horizon>1
    error([mfilename,':: autocovariances for expansion autocov_ar greater than 1 not ready'])
end
% 0- get the solution and put it in the order_var form
%-----------------------------------------------------
ov=obj.order_var.after_solve;
z_pb=obj.locations.after_solve.z.pb;
t_pb=obj.locations.after_solve.t.pb;
% t_sf=[obj.locations.after_solve.t.s,obj.locations.after_solve.t.f];
e_0=obj.locations.after_solve.z.e_0;
Re=cell(1,reg_nbr);
Tz=Re;
for isol=1:reg_nbr
    tmp=obj.solution.Tz{isol}(ov,:);
    % separate autoregressive part from shocks
    %----------------------------------------
    Tz{isol}=tmp(:,z_pb);
    Re{isol}=tmp(:,e_0);
end

% aggregate
%----------
Q=obj.solution.transition_matrices.Q;
pai=[eye(reg_nbr)-Q'
    ones(1,reg_nbr)]\[zeros(reg_nbr,1)
    1];
T=0;
R=0;
for ireg=1:reg_nbr
    T=T+pai(ireg)*Tz{ireg};
    R=R+pai(ireg)*Re{ireg};
end

% 1-compute covariance of the state variables first
%--------------------------------------------------
n=obj.endogenous.number(end);
V=zeros(n,n,autocov_ar+1);
[Vx,retcode]=lyapunov_equation(T(t_pb,:),R(t_pb,:)*R(t_pb,:)',obj.options);
if ~retcode
    % 2-recompute covariance of the whole vector conditional on the state
    % covariance. Should/could be used also during estimation. Need a way of
    % separating out the stationary and the nonstationary
    %-----------------------------------------------------------------------
    V(:,:,1)=T*Vx*T'+R*R';
    for ii=2:autocov_ar+1
        V(:,:,ii)=T*V(t_pb,:,ii-1);
    end
    % re-order alphabetically
    %------------------------
    inv_order_var=obj.inv_order_var.after_solve;
    V=V(inv_order_var,inv_order_var,:);
end

end

