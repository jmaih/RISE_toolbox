function [x1,f1,H,x0,f0,viol,funevals,issue,obj]=find_posterior_mode(obj,x0,lb,ub)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

nobj=numel(obj);

if nobj==0
    
    if nargout>1
        
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    
    end
    
    x1=struct();
    
    return

end

% initialize the number of function calls
%----------------------------------------
funevals=0;

Nsim=max(1,obj(1).options.estim_parallel);
% find the linear restrictions and shorten everything
%----------------------------------------------------
npar=size(x0,1);

x0=[x0,nan(npar,Nsim-1)];

f0=nan(1,Nsim);

% all objects are evaluated at the same point. Without a second argument,
% this is exactly what will happen.
[x0(:,1),f0(1),~,retcode0,viol,obj,funevals]=estimation_wrapper(obj,[],x0(:,1),lb,ub,funevals);

if retcode0||any(viol>0)
    % first check constraint violations
    f0(1)=obj(1).options.estim_penalty;

end
% the objects are automatically updated and potentially contain crucial
% information going forward. In particular, they contain information
% about whether the models are stationary or not.

if f0(1)<obj(1).options.estim_penalty
    
    beg=2;

else
    
    beg=1;

end

%% disable those elements
warning('off','MATLAB:nearlySingularMatrix')

warning('off','MATLAB:illConditionedMatrix')

fprintf(1,'%s\n','Looking for good enough start values. Please wait...');

gen_start=@()estimation_wrapper(obj,'draw',[],lb,ub,funevals);

for ii=beg:Nsim
    
    NotDone=true;
    
    iter=0;
    
    while NotDone
        
        [xtest,ftest,retcode]=...
            utils.estim.generate_starting_point(gen_start);
        
        if ftest<obj(1).options.estim_penalty
            
            NotDone=false;
            
            f0(ii)=ftest;
            
            x0(:,ii)=xtest;
        
        end
        
        iter=iter+1;
        
        if iter>=obj(1).options.estim_max_trials
            
            error([mfilename,':: No admissible starting value after ',...
                int2str(obj(1).options.estim_max_trials),' trials'])
        
        else
            
            fprintf(1,'%3.0d :: %s\n',iter,utils.error.decipher(retcode));
        
        end
        
    end
    
    disp(['Starting value # ',int2str(ii),' found after ',int2str(iter),' iterations'])
    
    ratio=ii/Nsim;
    
    fprintf(1,'%s\n',['...', num2str(100*ratio),'% done']);
    
end

[x1,f1,H,issue,viol,obj,funevals]=estimation_wrapper(obj,'estimate',x0,lb,ub,funevals);

viol=viol(viol>0);

warning('on','MATLAB:nearlySingularMatrix')

warning('on','MATLAB:illConditionedMatrix')

end
