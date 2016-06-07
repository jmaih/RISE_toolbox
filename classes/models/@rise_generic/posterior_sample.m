function [result,is_failed,time_it_took]=posterior_sample(m,pop,dowhat,howmany,ouf,varargin)
% POSTERIOR_SAMPLE -- computes a sample of any quantity of interest using
% parameter draws from a population e.g. a posterior simulation
%
% Syntax
% -------
% ::
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
% Inputs
% -------
%
% - **m** [rise|dsge|svar|rfvar|valid rise object]: model object
%
% - **pop** [m x n struct]: parameter draws, with "m" the number of chains
% and "n" the number of draws in each chain. Each element of "pop" is a
% structure with fields "f" (not used), the value of the posterior and "x"
% the parameter vector
%
% - **dowhat** [fhandle]: function (handle) to apply to each parameterized
% model object. e.g. dowhat=@irf, dowhat=@simulate, dowhat=@forecast, etc.
% "dowhat" need not be a method of "m": it represents the quantity of
% interest.
%
% - **howmany** [integer|{m x n}]: number of draws to use in the
% calculation
%
% - **ouf** [fhandle|{[]}]: output update function. Function that updates
% the output before storing it. e.g. if dowhat=@filter, one may be
% interested in the the filters only and in that case ouf=@(x)x.filtering.
%
% - **varargin** [pairwise args]: valid pairwise arguments for the model
% object
%
% Outputs
% --------
%
% - **result** [1 x howmany cell]: container of the various applications of
% the "dowhat" handle
%
% - **is_failed** [1 x howmany logical]: detector of failed simulations
%
% - **time_it_took** [numeric]: number of seconds needed to run all the
% simulations.
%
% More About
% ------------
%
% - the function will exploit parallel computation if there are workers
% idle.
%
% Examples
% ---------
%
% See also:

if isempty(m)
    
    result=struct();
    
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
        
    else
        
        if howmany>npop
            
            error('howmany cannot exceed the number of elements in the population')
            
        end
        
    end
    
    order=randperm(npop);
    
    order=order(1:howmany);
    
    result=cell(1,howmany);
    
    % cut the crap
    %-------------
    pop=pop(order);
    
    objective=@do_it;
    
    NumWorkers=utils.parallel.get_number_of_workers();
    
    is_failed=false(1,howmany);
    
    tic
    
    parfor(ii=1:howmany,NumWorkers)
        
        x=pop(ii).x;
        
        try
            
            result{ii}=objective(x);
            
        catch me
            
            disp(ii),disp(me.message),disp(' ')
            
            is_failed(ii)=true;
            
        end
        
    end
    
    time_it_took=toc;
    
end

    function out=do_it(x)
        
        obj=assign_estimates(m,x);
        
        out=dowhat(obj);
        
        out=ouf(out);
        
    end

end