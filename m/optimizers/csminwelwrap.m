function [x1,f1,exitflag]=csminwelwrap(fh_handle,x0,lb,ub,options,varargin)
% translate the inputs from RISE into inputs for csminwel
%--------------------------------------------------------
nx=numel(x0);
H0 = 1e-4*eye(nx);
crit = 1e-7;
nit = options.MaxIter; % using Matlab's MaxIter, which RISE will send
analytic_grad=[];

[fh,xh,gh,H,itct,funcCount,retcodeh] = csminwel(...
    @myobjective,x0,H0,analytic_grad,crit,nit,varargin{:});

% translate the outputs from csminwel into outputs for RISE
%----------------------------------------------------------
exitflag=1;
x1=xh;
f1=fh;

    function fval=myobjective(x)
        % RISE expect your optimizer to correctly handle parameter bounds
        % restrictions. If your optimizer does not, you have to do
        % something about it yourself as I do in this subfunction.
        
        if any(x<lb)||any(x>ub)
            
            fval=1e+12;
            
        else
            
            fval=fh_handle(x,varargin{:});
            
        end
        
    end

end
