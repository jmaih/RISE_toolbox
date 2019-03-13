function [x,fval,exitflag] = mycsolve(FUN,x0,opts,varargin)
% Inspired by Chris Sims' csolve

if nargin<3
    
    opts=[];
    
end

%------------ alpha ------------------
% tolerance on rate of descent
alpha=1e-3;
%---------------------------------------

opts=set_defaults(opts);

fcount=0;

nv=length(x0);

% differencing interval for numerical gradient
%---------------------------------------------
delta = opts.TolX;

tvec=delta*eye(nv);

done=false;

x=x0;

f0=objective(x);

af0=sum(abs(f0));

af00=af0;

iter=0;

is_verbose=strcmp(opts.Display,'iter');

if is_verbose
    
    % Print header.
    fprintf(['\n                                          Sum of       Step               \n',...
                ' Iteration  Func-count     f(x)          abs(f)       Length    retcode\n']);
    formatstr = ' %5.0f      %5.0f   %13.6g  %13.6g   %8.3g    %6.0f\n';
    
end

while ~done
    
    if iter>3 && af00-af0<opts.TolFun*max(1,af0) && rem(iter,2)==1
        
        randomize=1;
        
    else
        
        [~,grad]=objective(x,f0);
        
        if isreal(grad)
            
            dx0=-grad\f0;
            
            randomize=0;
            
        else
            
            if is_verbose
                
                disp('gradient imaginary')
                
            end
            
            randomize=true;
            
        end
        
    end
    
    if randomize
        
        if is_verbose
            
            fprintf(1,' Random Search\n')
            
        end
        
        dx0=norm(x)./randn(size(x));
        
    end
    
    [lambdamin,xmin,fmin,afmin,rc]=sub_problem();
    
    iter=iter+1;
    
    fval=max(abs(fmin));
    
    if is_verbose
        
      fprintf(formatstr,iter,fcount,fval,afmin,lambdamin,rc);
      
    end
    
    x=xmin;
    
    f0=fmin;
    
    af00=af0;
    
    af0=afmin;
    
    if iter >= opts.MaxIter
        
        done=true;
        
        rc=4;
        
    elseif af0<opts.TolFun
        
        done=true;
        
        rc=0;
        
    end
    
end

switch rc
    
    case 0
        
        exitflag=1;
        
    case {1,2,3}
        % 1 and 3 mean no solution despite extremely fine adjustments in
        % step length (very likely a numerical problem, or a discontinuity) 
        
        exitflag=-2;
        
    case 4
        %4 means opts.MaxIter
        
        exitflag=0;
        
end

    function [fx,g]=objective(x,fx)
        
        if nargin<2
            
            fx=FUN(x,varargin{:});
            
            fcount=fcount+1;
            
        end
        
        if nargout>1
                    
            g = fx(:,ones(nv,1));
                        
            for ii=1:nv
                
                f1=objective(x+tvec(:,ii));
                
                g(:,ii) = (f1-fx)/delta;
                
            end
            
        end

    end

    function [lambdamin,xmin,fmin,afmin,rc]=sub_problem()
        
        lambda=1;
        
        lambdamin=lambda;
        
        fmin=f0;
        
        xmin=x;
        
        afmin=af0;
        
        dxSize=norm(dx0);
        
        factor=.6;
        
        shrink=1;
        
        subDone=false;
        
        while ~subDone
            
            dx=lambda*dx0;
            
            f=objective(x+dx);
            
            af=sum(abs(f));
            
            if af<afmin
                
                afmin=af;
                
                fmin=f;
                
                lambdamin=lambda;
                
                xmin=x+dx;
                
            end
            
            if (lambda >0 && af0-af < alpha*lambda*af0) ||...
                    (lambda<0 && af0-af < 0)
                
                if ~shrink
                    
                    factor=factor^.6;
                    
                    shrink=true;
                    
                end
                
                if abs(lambda*(1-factor))*dxSize > .1*delta
                    
                    lambda = factor*lambda;
                    
                elseif lambda > 0 && factor==.6 %i.e., we've only been shrinking
                    
                    lambda=-.3;
                    
                else
                    
                    subDone=true;
                    
                    if lambda > 0
                        
                        if factor==.6
                            
                            rc = 2;
                            
                        else
                            
                            rc = 1;
                            
                        end
                        
                    else
                        
                        rc=3;
                        
                    end
                    
                end
                
            elseif lambda >0 && (af-af0 > (1-alpha)*lambda*af0)
                
                if shrink
                    
                    factor=factor^.6;
                    
                    shrink=false;
                    
                end
                
                lambda=lambda/factor;
                
            else % good value found
                subDone=true;
                
                rc=0;
                
            end
            
        end % while ~subDone
        
    end

end

function opts=set_defaults(opts)

options=struct('Display','off',...% if this is set to zero, all screen output is suppressed
    'TolFun',sqrt(eps),...
    'MaxIter',1000);

if isempty(opts)
    
    opts=struct();
    
elseif ~isstruct(opts)
    
    error('opts must be empty or a structure')
    
end

fopts=fieldnames(options);

for iopt=1:numel(fopts)
    
    if ~isfield(opts,fopts{iopt})
        
        opts.(fopts{iopt})=options.(fopts{iopt});
        
    end
    
end

end