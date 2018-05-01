function [x,f,viol,funevals]=generate_candidates(objective,lb,ub,n,...
max_trials,restrictions,opt,penalty,varargin)
% generate_candidates -- generate candidates for optimization
%
% ::
%
%
%   [x,f,viol,funevals]=generate_candidates(objective,lb,ub,n,...
%       restrictions,penalty,varargin)
%
% Args:
%
%    - **objective** [function_handle]: objective to minimize
%
%    - **lb** [vector]: lower bound of the search space
%
%    - **ub** [vector]: upper bound of the search space
%
%    - **n** [integer]: number of candidates to generate
%
%    - **max_trials** [integer]: number of trials after which the procedure
%    crashes
%
%    - **restrictions** [empty|function_handle]: function evaluating the
%    violations
%
%    - **opt** [struct]: structure with fields
%      - **restrictions_in_objective** [true|false]:
%      - **returns_retcode** [true|false]:
%      - **restrictions_same_weights** [true|{false}]:
%      - **allow_restrictions_violations** [true|{false}]:
%
%    - **penalty** [numeric]: value functions exceeding this value in
%    absolute value are assigned this value
%
%    - **varargin** []: additional input arguments for the objective function
%
% Returns:
%    :
%
%    - **x** [d x n matrix]: parameter vectors
%
%    - **f** [row vector]: value function for each parameter vector
%
%    - **viol** [row vector]: violation penalty for each parameter vector
%
%    - **funevals** [integer]: number of function evaluations
%
% Note:
%
% Example:
%
%    See also:

if isempty(restrictions)
    
    restrictions=@(z)[];
    
end

if isempty(max_trials)
    
max_trials=50;

end

MaxIter=50;

npar=size(lb,1);

x=nan(npar,n);

f=nan(1,n);

fcount=nan(1,n);

viol=zeros(1,n);

funevals=0;

% can we get 2 outputs?

if isempty(opt)
    
    opt=struct('restrictions_in_objective',false,...
    'returns_retcode',false,...
    'restrictions_same_weights',false,...
    'allow_restrictions_violations',false);

end

if opt.restrictions_in_objective
    
    opt.returns_retcode=true;
    
else
    
    try %#ok<TRYNC>
        
        xtest=lb+(ub-lb).*rand(npar,1);
        
        [~,rcode]=objective(xtest,varargin{:});
        
        try %#ok<TRYNC>
            
            if isnumeric(rcode) && isscalar(rcode)
                
                msg=decipher(rcode);
                
                opt.returns_retcode=true;
                
            end
            
        end
        
    end

end

% success=nargout(objective)>=2;
msg='';

myinvalid=@(v,mess,f)(~opt.allow_restrictions_violations && m>0)...
                ||~isempty(mess)||f>=penalty;

the_loop=@loop_body;

for ii=1:n
    
    [x(:,ii),f(ii),viol(ii),fcount(ii)]=the_loop(ii);
    
end

funevals=funevals+sum(fcount);

    function [xi,fi,violi,fcount]=loop_body(ii)
        
        invalid=true;
        
        iter=0;
        
        fcount=0;
        
        while invalid
            
            if iter>=MaxIter
                
                error([mfilename,':: could not generate a valid candidate after ',...
                    int2str(MaxIter*max_trials),' attempts. Exiting with the best vector'])
                
            end
            
            [xi,fi,violi]=draw_and_evaluate_vector();
            
            fcount=fcount+1;
            
            iter2=0;
            
            invalid=myinvalid(violi,msg,fi);
            
            while iter2<max_trials && invalid
                
                iter2=iter2+1;
                
                [xi,fi,violi]=draw_and_evaluate_vector();
                
                fcount=fcount+1;
                
                invalid=myinvalid(violi,msg,fi);
                
            end
            
            if invalid
                
                fprintf(1,'%5s %3.0d/%3.0d %8s %5.0d %8s %5.0d %8s %s \n',...
                    'warrior #',ii,n,'iter',iter,'interm runs',iter2,'pb',msg);
                
            end
            
            iter=iter+1;
            
        end
        
    end

    function [c,fc,viol]=draw_and_evaluate_vector()
        
        c=lb+(ub-lb).*rand(npar,1);
        
        [fc,viol,msg]=utils.estim.eval_objective_and_restrictions(c,...
            objective,restrictions,opt,varargin{:});
        
    end


end