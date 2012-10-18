function [A,retcode]=theoretical_autocovariances(obj,ar)%,resolve_flag
default_ar=5;
if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    A=struct('ar',default_ar);
    return
end
if nargin<2
    ar=[];
end
if isempty(ar)
    ar=obj.options.ar;
end

[obj,retcode]=solve(obj);
A=[];
if retcode
%     error([mfilename,':: no solution: theoretical autocovariances cannot be computed'])
    return
end
    
T=obj.T;
R=obj.R;
[n,junk,expansion_order,reg_nbr]=size(R);
if expansion_order>1
    error([mfilename,':: autocovariances for expansion ar greater than 1 not ready'])
end
A=nan(n,n,ar+1);
if reg_nbr==1
    [P,retcode]=lyapunov_equation(T,R*R',obj.options);
    if ~retcode
        A(:,:,1)=P;
        for ii=2:ar+1
            A(:,:,ii)=T*A(:,:,ii-1);
        end
    end
else
    % We set the first flag (measurement errors) to false. The model has
    % already been solved above and we don't need to resolve it. And so,
    % the second flag (resolve) is set to false as well.
    N=1500;
    [obj,State] = simulate(obj,N,false,false);
    data=vertcat(obj.varendo.value);
    steady_state=vertcat(obj.varendo.det_steady_state);
    % use the steady steady state you were in to deflate the data
    for ii=1:reg_nbr
        reg_loc=State==ii;
        if any(reg_loc)
            data(:,reg_loc)=bsxfun(@minus,data(:,reg_loc),steady_state(:,ii));
        end
    end
    % Now, in theory at least, the simulated data have zero mean
    for ii=1:ar+1
        A(:,:,ii)=data(:,ar+1:end)*data(:,(ar+1:end)-ii+1)';
    end
    A=A/(N-ar);
end

