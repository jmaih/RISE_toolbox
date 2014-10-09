function [V,retcode]=lyapunov_equation(T,Q,options)
% lyapunov_equation solves the equation V=T*V*T'+Q
%
% Syntax
% -------
% ::
%   [V,retcode]=lyapunov_equation(T,Q)
%   [V,retcode]=lyapunov_equation(T,Q,options)
%
% Inputs
% -------
% - T :
% - Q :
% - options :
%
% Outputs
% --------
% - V :
% - retcode :
% 
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

% example
% we would like to compute the variance of the system
% y_t=ay*y_{t-1}+ax*x_{t-1}+ae*e_t
% x_t=by*y_{t-1}+bx*x_{t-1}+be*e_t
% the matrix representation is 
% X_t=[ay ax;by bx]*X_{t-1}+[ae be]'*e_t;

% 
defaults=struct('lyapunov_algo','doubling',...
    'lyapunov_diffuse_factor',0,...
    'lyapunov_fast_doubling',true);
if nargin==0
	if nargout>1
		error([mfilename,':: number of output arguments cannot exceed 1 if there are no inputs'])
	end
	V=defaults;
	return
end

if nargin<3
    options=[];
end
if isfield(options,'lyapunov_fast_doubling')
    lyapunov_fast_doubling=options.lyapunov_fast_doubling;
else
    lyapunov_fast_doubling=defaults.lyapunov_fast_doubling;
end
if isfield(options,'lyapunov_algo')
    algo=options.lyapunov_algo;
else
    algo=defaults.lyapunov_algo;
end
if isfield(options,'lyapunov_diffuse_factor')
    lyapunov_diffuse_factor=options.lyapunov_diffuse_factor;
else
    lyapunov_diffuse_factor=defaults.lyapunov_diffuse_factor;
end

switch algo
    case 'schur'
        [Ts,Qs,Us,stable]=lower_triangularize(T,Q,options);
        [Vs,retcode]=solve_triangular_lyapunov(Ts(stable,stable),Qs(stable,stable));
        V=[];
        if ~retcode
            V=zeros(size(T));
            V(stable,stable)=Vs;
            V(~stable,~stable)=diag(lyapunov_diffuse_factor*ones(sum(~stable),1));
            % complex schur was taken earlier, reverse it by taking real.
            % But we only take the stable guys...
            V=real(Us'*V*Us);
%             V(stable,stable)=real(Us(stable,stable)'*V(stable,stable)*Us(stable,stable));
        end
    case 'doubling'
        if lyapunov_fast_doubling
            tmpT=T;
            tmpQ=Q;
            state_vars=any(T);
            T=T(state_vars,state_vars);
            Q=Q(state_vars,state_vars);
        end
        V0=Q;
        G0=T;
        [V,~,retcode]=fix_point_iterator(@iterator,V0,options);
        if retcode
            retcode=280+retcode;
        else
            if lyapunov_fast_doubling
                V=tmpT(:,state_vars)*V*tmpT(:,state_vars)'+tmpQ;
            end
        end
    otherwise
        error([mfilename,':: unknown lyapunov algorithm option ',algo])
end
    function [V,F0]=iterator(V0)
        V=V0+G0*V0*G0';
        G0=G0*G0;
        F0=V-V0;
   end
end

function [V,retcode]=solve_triangular_lyapunov(T,Q)
% check that T is triangular
if ~is_lower_triangular(T)
    error([mfilename,':: T must be upper triangular'])
end
R=transpose(T);
n=size(T,1);
I=eye(n);
V=nan(n);
TV=nan(n);
retcode=0;
for k=1:n
    tmp=(I-R(k,k)*T);
    if rcond(tmp)<1e-12
        V=nan(n);
        retcode=302;
        return
    end
    V(:,k)=tmp\(TV(:,1:k-1)*R(1:k-1,k)+Q(:,k));
    if k<n
        % in order to avoid recomputing TV all the time.
        TV(:,k)=T*V(:,k);
    end
end
end

function [Ts,Qs,Us,stable]=lower_triangularize(T,Q,options)
if nargin<3
    Q=[];
    if nargin<2
        options=[];
    end
end
if isempty(options)
    options=1e-9;
end
switch class(options)
    case 'double'
        tol=options;
    case 'struct'
        if isfield(options,'fix_point_TolFun')
            tol=options.fix_point_TolFun;
        else
            tol=1e-9;
        end
    otherwise
        error([mfilename,':: class ''',class(options),''' not supported for input argument options'])
end
% check whether T is quasi-triangular
% Ts=T;Qs=Q;Us=1;stable=true(1,size(T,1));
% if ~is_lower_triangular(T)
    [Up_s,Tp_s]=schur(transpose(T),'complex'); 
    % go complex to ensure the matrix is triangular. Otherwise the matrix
    % may be quasi-triangular and we might have to solve for blocks, which
    % I don't feel like doing today.
    % take the eigenvalues which are now on the diagonal of Tp_s
    E=diag(Tp_s); % <--- E=ordeig(Tp_s);
    % separate stable from explosive and unit
    unit=abs(E)>1-tol;
    stable=~unit;
    if any(unit)
        clusters=zeros(size(E));
        clusters(unit)=1;
        % things are sorted in descending order of the elements in clusters
        % and so the explosive and unit roots will come first.
        [Up_s,Tp_s]=ordschur(Up_s,Tp_s,clusters);
        stable=[false(1,sum(unit)),true(1,sum(stable))];
    end
    % Now, Tp_s is upper triangular and we have transpose(T)=Up_s*Tp_s*Up_s' or
    % alternatively T=Up_s*Tp_s'*Up_s'
    Ts=Tp_s';
    Us=Up_s';
    Qs=[];
    if ~isempty(Q)
        Qs=Us'*Q*Us;
    end
    % we have Us'*Ts*Us-T=0
% end
end

function flag=is_lower_triangular(T)
n=size(T,1);
flag=true;
for ii=1:n
    if any(T(1:ii-1,ii))
        flag=false;
        break
    end
end
end
