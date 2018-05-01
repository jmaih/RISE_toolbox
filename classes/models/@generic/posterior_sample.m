function [result,time_it_took]=posterior_sample(m,pop,dowhat,howmany,ouf,varargin)
% POSTERIOR_SAMPLE -- computes a sample of any quantity of interest using
% parameter draws from a population e.g. a posterior simulation
%
% ::
%
%
%   [result]=POSTERIOR_SAMPLE(m,pop,dowhat)
%
%   [result]=POSTERIOR_SAMPLE(m,pop,dowhat,howmany)
%
%   [result]=POSTERIOR_SAMPLE(m,pop,dowhat,howmany,ouf)
%
%   [result]=POSTERIOR_SAMPLE(m,pop,dowhat,howmany,ouf,varargin)
%
%   [result,is_failed,time_it_took]=POSTERIOR_SAMPLE(...)
%
% Args:
%
%    - **m** [rise|dsge|svar|rfvar|valid rise object]: model object
%
%    - **pop** [m x n struct]: parameter draws, with "m" the number of chains
%    and "n" the number of draws in each chain. Each element of "pop" is a
%    structure with fields "f" (not used), the value of the posterior and "x"
%    the parameter vector
%
%    - **dowhat** [fhandle]: function (handle) to apply to each parameterized
%    model object. e.g. dowhat=@irf, dowhat=@simulate, dowhat=@forecast, etc.
%    "dowhat" need not be a method of "m": it represents the quantity of
%    interest.
%
%    - **howmany** [integer|{m x n}]: number of draws to use in the
%    calculation
%
%    - **ouf** [fhandle|{[]}]: output update function. Function that updates
%    the output before storing it. e.g. if dowhat=@filter, one may be
%    interested in the the filters only and in that case ouf=@(x)x.filtering.
%
%    - **varargin** [pairwise args]: valid pairwise arguments for the model
%    object
%
% Returns:
%    :
%
%    - **result** [1 x howmany cell]: container of the various applications of
%    the "dowhat" handle
%
%    - **time_it_took** [numeric]: number of seconds needed to run all the
%    simulations.
%
% Note:
%
%    - the function will exploit parallel computation if there are workers
%    idle.
%
%    - Because the solving of the model is sometimes iterative, a change of
%    solver or of the settings of the solver can result in the model not
%    being solved or more generally simulations failures. The algorithm will
%    loop until the requested number of simulations is obtained. But it will
%    not point to the parameter vectors that fail.
%
% Example:
%
%    See also:
    
% hard coded
%-----------
max_while_loops=20;

if isempty(m)
    
    result=cell(0,4);
    
else
    
    if nargin<5
        
        ouf=[];
        
        if nargin<4
            
            howmany=[];
            
        end
        
    end
    
    if isempty(ouf)
        
        ouf=@(x)x;
        
    end
    
    m=set(m,varargin{:});
    
    m=setup_restrictions(m);
    
    pop=pop(:).';
    
    npop=numel(pop);
    
    if isempty(howmany)
        
        howmany=npop;
                
    end
    
    result=cell(1,howmany);
    
    objective=@do_it;
    
    NumWorkers=utils.parallel.get_number_of_workers();
    
    tic
    
    offset=0;
    
    howmany0=howmany;
    
    iter=0;
    
    while howmany0 && iter<max_while_loops
        
        iter=iter+1;
        
        small_result=do_one_pass(howmany0);
        
        ngood=numel(small_result);
        
        result(offset+(1:ngood))=small_result;
        
        offset=offset+ngood;
        
        howmany0=howmany-offset;
        
    end
    
    if iter==max_while_loops
        
        error('Too many vectors fail. This can occur if e.g. the solver or the solver''s settings are changed')
        
    end
    
    time_it_took=toc;
    
end

    function small_result=do_one_pass(hm)
        % update the number in case population is shrunk
        %-----------------------------------------------
        npop=numel(pop);
        
        if hm>npop
            
            error('howmany exceeds the number of elements in the population... possibly after chop off')
            
        end
        % reorder
        %-------------
        order=randperm(npop);
        
        pop=pop(order);
        
        small_result=cell(1,hm);
        
        is_failed_this_pass=true(1,hm);
        
        parfor(ii=1:hm,NumWorkers)
            
            if is_failed_this_pass(ii)
                
                x=pop(ii).x;
                
                try
                    
                    small_result{ii}=objective(x); %#ok<PFOUS>
                    
                    is_failed_this_pass(ii)=false;
                    
                catch me
                    
                    disp(ii),disp(me.message),disp(' ')
                    
                end
                
            end
            
        end
        
        % bump the failed
        %----------------
        failed=find(is_failed_this_pass);
        
        small_result(is_failed_this_pass)=[];
        
        if ~isempty(failed)
            
            warning(sprintf('%s:: %0.0f failed simulations',mfilename,numel(failed))) %#ok<SPWRN>
            
        end
        
        pop=pop(hm+1:end);
        
    end

    function out=do_it(x)
        
        obj=assign_estimates(m,x);
        
        out=dowhat(obj);
        
        out=ouf(out);
        
    end

end