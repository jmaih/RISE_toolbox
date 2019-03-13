function [x,fval,exitflag] = mynewton(f,x,opts,varargin)

if nargin<3
    
    opts=[];
    
end

fcount=0;

nv=length(x);

objfun=@objective;

opts=nonlinsyst.set_defaults(opts);

% differencing interval for numerical gradient
%---------------------------------------------
delta = opts.TolX;

tvec=delta*eye(nv);

is_verbose=strcmp(opts.Display,'iter');

exitflag=1;

line_search_maxsteps=25;

iter=0;

xbest=struct('x',x,'fval',nan,'fnorm',inf);

while iter < opts.MaxIter
    
    iter=iter+1;
    
    if iter==1 && is_verbose
        
        % Print header.
        fprintf('\n Iteration  Func-count     norm f(x)      backstep\n');
        formatstr = '%5.0f      %5.0f   %15.6g    %8.0f\n';
        
    end
    
    [fval,fjac] = objfun(x);
    
    fnorm = nonlinsyst.norm(fval);
    
    if isnan(fnorm)||isinf(fnorm)
        
        exitflag=-2;
        
        break
        
    elseif fnorm<1e-4*opts.TolFun
        
        if is_verbose
            
            fprintf(formatstr,iter,fcount,fnorm,0);
            
        end
        
        x=xbest.x;
        
        fval=xbest.fval;
        
        return
        
    end
    
    dx = -(fjac\fval);
    
    if any(isnan(dx)|isinf(dx))
        
        exitflag=-2;
        
        break
        
    end
    
    % perform linesearch
    [dx,~,~,flag,backstep]= nonlinsyst.linesearch(objfun,x,dx,fnorm,line_search_maxsteps);
    
%     if flag
%         
%         dx=norm(x)./randn(size(x));
%         
%     end
    
    x = x+dx;
    
    if is_verbose
        
        fprintf(formatstr,iter,fcount,xbest.fnorm,backstep);
        
    end
    
end

if iter>=opts.MaxIter
    
    exitflag=0;
    
end

x=xbest.x;

fval=xbest.fval;

    function [o,g]=objective(x)
        
        o=f(x,varargin{:});
        
        fcount=fcount+1;
        
        nn=nonlinsyst.norm(o);
        
        if nn<xbest.fnorm
            
            xbest.x=x;
            
            xbest.fval=o;
            
            xbest.fnorm=nn;
            
        end
        
        if nargout>1
            
            g = o(:,ones(nv,1));
            
            f0=o;
            
            for i=1:nv
                
                f1=objective(x+tvec(:,i));
                
                g(:,i) = (f1-f0)/delta;
                
            end
            
        end
        
    end

end
