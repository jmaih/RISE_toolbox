function J=first_order_derivatives(func,x,varargin)
x=x(:);
nx=numel(x);
vectorized=nx<300;
if vectorized
    x0=automatic(x,speye(nx));
    out=func(x0,varargin{:});
    J=vertcat(out.dx);
else
    x0=automatic(x,zeros(nx,1));
    ftest=func(x,varargin{:});
    nr=numel(ftest);
    J=nan(nr,nx);
    for ii=1:nx
        x0_i=x0;
        % activate variable ii
        x0_i(ii).dx=1;
        out=func(x0_i,varargin{:});
        J(:,ii)=vertcat(out.dx);
    end
end
