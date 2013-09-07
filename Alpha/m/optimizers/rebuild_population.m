function [xx,ff,funevals]=rebuild_population(xx,ff,Objective,lb,ub,tol,varargin)
[npar,MaxNodes]=size(xx);
ul=ub-lb;
funevals=0;
newstyle=true;
for ii=1:MaxNodes
    xi=(xx(:,ii)-lb)./ul;
    if isnan(ff(ii))
        ff(ii)=Objective(xx(:,ii),varargin{:});
        funevals=funevals+1;
    end
    for jj=ii+1:MaxNodes
        xj=(xx(:,jj)-lb)./ul;
        if distance(xi,xj)<tol
            if ff(ii)<ff(jj)
                change=jj;
            else
                change=ii;
            end
            if newstyle
                xx(:,change)=lb+ul.*rand(npar,1);
            else
                u=rand(npar,1);
                shift=.5*u+.25*u.*randn(npar,1);
                xx(:,change)=recenter(xx(:,change)+shift,lb,ub);
            end
            ff(change)=Objective(xx(:,change),varargin{:});
            funevals=funevals+1;
        elseif isnan(ff(jj))
            ff(jj)=Objective(xx(:,jj),varargin{:});
            funevals=funevals+1;
        end
    end
end