function mycontrib=historical_decomposition(m,sim_engine,nsim,varargin)
% historical_decomposition Computes historical decomposition of a nonlinear
% DSGE model 
%
% Syntax
% -------
% ::
%
%   mycontrib=historical_decomposition(m)
%   mycontrib=historical_decomposition(m,sim_engine)
%   mycontrib=historical_decomposition(m,sim_engine,nsim)
%   mycontrib=historical_decomposition(m,sim_engine,nsim,varargin)
%
% Inputs
% -------
%
% - m : [rise|dsge] model(s) for which to compute the historical 
%   decomposition. m could be a vector of models. The computation of the
%   historical decomposition requires a computation of smoothed history.
%   Hence either m should contain all the information needed for that
%   purpose or the information should be passed along through varargin below
%
% - sim_engine : [empty|function handle] function computing forecasts
%
% - nsim : [empty|numeric|{100}] Number of simulation to consider for the
%   integration exercise. nsim is automatically set to 1 if the model is
%   detected to be solved at order 1 and not contain regime switches.
%
% - varargin : additional information needed for the computation of the
%   smoothed quantities through filtration.
%
% Outputs
% --------
%
% - mycontrib : [struct|cell array] structure or cell array of structures
%   with the contributions in each model and each solution. The
%   decompositions are given in terms of:
%   - the exogenous variables
%   - **init** : the effect of initial conditions, which includes the
%     steady state!!!
%
% Remarks
% --------
%
% - **N.B** : For a nonlinear model (e.g linear but switching), the type of
%   decomposition/historical_decomposition we do for linear/linearized
%   constant-parameter models is not feasible. RISE uses a monte carlo
%   integration to provide an approximation to the
%   decomposition/historical_decomposition
%
% - **N.B** : if m is a vector of models, then each model should return a
%   unique solution (implying a unique filtration), else the concatenation
%   of decompositions will fail. In that case it is better to run one model
%   at a time.
%
% - **N.B** : When a model has multiple solutions, there is no guarantee
%   that every solution will have succeed at filtration. For this reason,
%   the result of each decomposition is stored in a separate cell. Empty
%   cells represent the solutions for which the filtration could not be
%   obtained.
%
% Examples
% ---------
%
% See also: dsge/historical_decomposition

if isempty(m)
    
    mydefaults=cell(0,4);
    
    if nargout
        
        mycontrib=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

if nargin<3
    
    nsim=[];
    
    if nargin<2
        
        sim_engine=[];
        
    end
    
end

nobj=numel(m);

if nobj>1
    
    tmp=cell(1,nobj);
    
    for ii=1:nobj
        
        tmp{ii}=nonlinear_historical_decomposition(m(ii),sim_engine,nsim,...
            varargin{:});
        
        if ii==1
            
            vnames=fieldnames(tmp{ii});
            
        else
            
            vnames=intersect(vnames,fieldnames(tmp{ii}));
            
        end
        
    end
    
    mycontrib=struct();
    
    for iname=1:numel(vnames)
        
        v=vnames{iname};
        
        db=tmp{1}.(v);
        
        for ii=2:nobj
            
            db=cat(3,db,tmp{ii}.(v));
            
        end
        
        mycontrib.(v)=db;
        
    end
    
    return
    
end

m=set(m,varargin{:});

mydefaults=the_defaults();

allshks=get(m(1),'exo_list');

recheck_inputs()

[filtration,~,~,retcode,m]=filter(m);

if all(retcode)

	error(decipher(retcode))

end

[T,~,Qfunc,ss,~,is_state]=load_solution(m,'ov');

nsols=m.nsols;

if nsols>1
    
    mycontrib=cell(1,nsols);
    
    for isol=1:nsols
        
        if retcode(isol)
            
            continue
            
        end
        
        mycontrib{isol}=do_one_contribution(filtration{isol},T(:,:,isol));
        
    end
    
else
    
     mycontrib=do_one_contribution(filtration,T);
    
end

    function mycontrib=do_one_contribution(filtration,T)
        
        [smoothed,shocks,pai0,ordered_endos,start_date]=...
            utils.filtering.smoothed_quantities(m,filtration);
        
        contrib=utils.filtering.nonlinear_shock_decomp(sim_engine,smoothed,shocks,nsim,...
            ss,T,Qfunc,is_state,pai0);
        
        cnames=[allshks,'init'];
        
        tmplt=ts(start_date,permute(contrib(1,:,:),[2,3,1]),cnames);
        
        mycontrib=struct();
        
        for iii=1:numel(ordered_endos)
            
            v=ordered_endos{iii};
            
            mycontrib.(v)=set(tmplt,'data',permute(contrib(iii,:,:),[2,3,1]));
            
        end
        
    end


    function recheck_inputs()
        
        myopts=struct('histdec_nsim',nsim,...
            'histdec_sim_engine',sim_engine);
        
        opts=disp_defaults(mydefaults);
        
        fields=fieldnames(opts);
        
        for iii=1:numel(fields)
            
            fname=fields{iii};
            
            if ~isempty(myopts.(fname))
                
                opts.(fname)=myopts.(fname);
                
            end
            
        end
        
        opts=parse_arguments(mydefaults,opts);
        
        nsim=opts.histdec_nsim;
        
        if ~any(m.parameters.is_switching) && m.options.solve_order==1
            
            nsim=1;
            
        end
        
        sim_engine=opts.histdec_sim_engine;
        
    end

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=1;

d={
    'histdec_sim_engine',@utils.filtering.sim_engine,@(x)isempty(x)||isa(x,'function_handle'),...
    'sim_engine must be a function handle'
    'histdec_nsim',100,@(x)isempty(x)||num_fin_int(x),...
    'nsim must be a finite and positive integer'
    };
end
