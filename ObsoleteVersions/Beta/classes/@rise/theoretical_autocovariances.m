function [A,retcode]=theoretical_autocovariances(obj,varargin)%,resolve_flag
default_ar=5;
if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    A=struct('ar',default_ar);
    return
end

nobj=numel(obj);
if nobj>1
    A=cell(1,nobj);
    retcode=cell(1,nobj);
    for iobj=1:nobj
        [A{iobj},retcode{iobj}]=theoretical_autocovariances(obj(iobj),varargin{:});
    end
    return
end
obj=set_options(obj,varargin{:});
ar=obj.options.ar;

[obj,retcode]=solve(obj);
A=[];
if retcode
    %     error([mfilename,':: no solution: theoretical autocovariances cannot be computed'])
    return
end

T0=obj.T;
R0=obj.R;
Q=obj.Q;
[n,~,expansion_order,reg_nbr]=size(R0);
if expansion_order>1
    error([mfilename,':: autocovariances for expansion ar greater than 1 not ready'])
end
pai=[eye(reg_nbr)-Q'
    ones(1,reg_nbr)]\[zeros(reg_nbr,1)
    1];
T=0;
R=0;
for ireg=1:reg_nbr
    T=T+pai(ireg)*T0(:,:,ireg);
    R=R+pai(ireg)*R0(:,:,ireg);
end
A=nan(n,n,ar+1);
% if reg_nbr==1
[P,retcode]=lyapunov_equation(T,R*R',obj.options);
if ~retcode
    A(:,:,1)=P;
    for ii=2:ar+1
        A(:,:,ii)=T*A(:,:,ii-1);
    end
end
% else
%     % We set the first flag (measurement errors) to false. The model has
%     % already been solved above and we don't need to resolve it. And so,
%     % the second flag (resolve) is set to false as well.
%     N=1500;
%     [obj,State] = simulate(obj,N,false,false);
%     data=vertcat(obj.varendo.value);
%     steady_state=vertcat(obj.varendo.det_steady_state);
%     % use the steady steady state you were in to deflate the data
%     for ii=1:reg_nbr
%         reg_loc=State==ii;
%         if any(reg_loc)
%             data(:,reg_loc)=bsxfun(@minus,data(:,reg_loc),steady_state(:,ii));
%         end
%     end
%     % Now, in theory at least, the simulated data have zero mean
%     for ii=1:ar+1
%         A(:,:,ii)=data(:,ar+1:end)*data(:,(ar+1:end)-ii+1)';
%     end
%     A=A/(N-ar);
% end

