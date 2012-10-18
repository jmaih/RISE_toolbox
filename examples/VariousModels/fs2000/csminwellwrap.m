function [xh,fh,exitflag]=csminwellwrap(fh_handle,x0,lb,ub,options,varargin)

nx=numel(x0);
H0 = 1e-4*eye(nx);
crit = 1e-7;
nit = 1000;
analytic_grad=[];
exitflag=1;

[fh,xh,gh,H,itct,fcount,retcodeh] = csminwel(...
    @myobjective,x0,H0,analytic_grad,crit,nit,varargin{:});

    function fval=myobjective(x)
        if any(x<lb)||any(x>ub)
            fval=1e+12;
        else
            fval=fh_handle(x,varargin{:});
        end
    end

end
