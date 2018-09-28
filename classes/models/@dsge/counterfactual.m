function counterf=counterfactual(m,sim_engine,nsim,shock_names,...
    varargin)
% counterfactual Computes counterfactual history of a nonlinear DSGE model
%
% Syntax
% -------
% ::
%
%   counterf=counterfactual(m)
%   counterf=counterfactual(m,sim_engine)
%   counterf=counterfactual(m,sim_engine,nsim)
%   counterf=counterfactual(m,sim_engine,nsim,shock_names)
%   counterf=counterfactual(m,sim_engine,nsim,shock_names,varargin)
%
% Inputs
% -------
%
% - m : [rise|dsge] model(s) for which to compute the
%   counterfactual. m could be a vector of models. The computation of the
%   counterfactual requires a computation of smoothed history. Hence either
%   m should contain all the information needed for that purpose or the
%   information should be passed along through varargin below
%
% - sim_engine : [empty|function handle] function computing forecasts
%
% - nsim : [empty|numeric|{100}] Number of simulation to consider for the
%   integration exercise. nsim is automatically set to 1 if the model is
%   detected to be solved at order 1 and not contain regime switches.
%
% - shock_names : [empty|char|cellstr] list of shocks to consider in the
%   computation of the counterfactual
%
% - varargin : additional information needed for the computation of the
%   smoothed quantities through filtration.
%
% Outputs
% --------
%
% - counterf : [struct|cell array] structure or cell array of structures
%   with the counterfactuals in each model and each solution.
%
% Remarks
% --------
%
% - **N.B** : For a nonlinear model (e.g linear but switching), the type of
%   decomposition/counterfactual we do for linear/linearized
%   constant-parameter models is not feasible. RISE uses a monte carlo
%   integration to provide an approximation to the
%   decomposition/counterfactual
%
% - **N.B** : if m is a vector of models, then each model should return a
%   unique solution (implying a unique filtration), else the concatenation
%   of counterfactuals will fail. In that case it is better to run one
%   model at a time.
%
% - **N.B** : When a model has multiple solutions, there is no guarantee
%   that every solution will have succeed at filtration. For this reason,
%   the result of each counterfactual is stored in a separate cell. Empty
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
        
        counterf=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

if nargin<4
    
    shock_names=[];
    
    if nargin<3
        
        nsim=[];
        
        if nargin<2
            
            sim_engine=[];
            
        end
        
    end
    
end

nobj=numel(m);

if nobj>1
    
    tmp=cell(1,nobj);
    
    for isol=1:nobj
        
        tmp{isol}=counterfactual(m(isol),sim_engine,nsim,shock_names,...
            varargin{:});
        
        if isol==1
            
            vnames=fieldnames(tmp{isol});
            
        else
            
            vnames=intersect(vnames,fieldnames(tmp{isol}));
            
        end
        
    end
    
    counterf=struct();
    
    for iname=1:numel(vnames)
        
        v=vnames{iname};
        
        db=tmp{1}.(v);
        
        for isol=2:nobj
            
            db=[db,tmp{isol}.(v)]; %#ok<AGROW>
            
        end
        
        counterf.(v)=db;
        
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

active=ismember(allshks,shock_names);

[T,~,Qfunc,ss,~,is_state]=load_solution(m,'ov');

nsols=m.nsols;

if nsols>1
    
    counterf=cell(1,nsols);
    
    for isol=1:nsols
        
        if retcode(isol)
            
            continue
            
        end
        
        counterf{isol}=do_one_counterfactual(filtration{isol},T(:,:,isol));
        
    end
    
else
    
    counterf=do_one_counterfactual(filtration,T);
    
end

    function counterf=do_one_counterfactual(filtration,T)
        
        [smoothed,shocks,pai0,ordered_endos,start_date]=...
            utils.filtering.smoothed_quantities(m,filtration);
        
        % the initial condition, in last position, is always active
        %----------------------------------------------------------
        active=[active(:).',true];
        
        contrib=utils.filtering.nonlinear_shock_decomp(sim_engine,smoothed,shocks,nsim,...
            ss,T,Qfunc,is_state,pai0);
        
        contrib(:,:,~active)=0;
        
        CF=sum(contrib,3);
        
        tmplt=ts(start_date,CF(1,:).');
        
        counterf=struct();
        
        for iii=1:numel(ordered_endos)
            
            v=ordered_endos{iii};
            
            counterf.(v)=set(tmplt,'data',CF(iii,:).');
            
        end
        
    end

    function recheck_inputs()
        
        if ~isempty(shock_names)
            
            if ischar(shock_names)
                
                shock_names=cellstr(shock_names);
                
            end
            
            good=ismember(shock_names,allshks);
            
            if any(~good)
                
                disp(shock_names(~good))
                
                error('The shocks above are not declared in the model')
                
            end
            
        end
        
        myopts=struct('ctrf_nsim',nsim,...
            'ctrf_sim_engine',sim_engine);
        
        myopts.ctrf_shock_list=shock_names;
        
        opts=disp_defaults(mydefaults);
        
        fields=fieldnames(opts);
        
        for iii=1:numel(fields)
            
            fname=fields{iii};
            
            if ~isempty(myopts.(fname))
                
                opts.(fname)=myopts.(fname);
                
            end
            
        end
        
        opts=parse_arguments(mydefaults,opts);
        
        nsim=opts.ctrf_nsim;
        
        if ~any(m.parameters.is_switching) && m.options.solve_order==1
            
            nsim=1;
            
        end
        
        sim_engine=opts.ctrf_sim_engine;
        
        shock_names=opts.ctrf_shock_list;
        
    end

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=1;

d={
    'ctrf_sim_engine',@utils.filtering.sim_engine,@(x)isempty(x)||isa(x,'function_handle'),...
    'sim_engine must be a function handle'
    'ctrf_nsim',100,@(x)isempty(x)||num_fin_int(x),...
    'nsim must be a finite and positive integer'
    'ctrf_shock_list','',@(x)isempty(x)||ischar(x)||iscellstr(x),...
    'shock_list must be char or cellstr'
    }; %#ok<ISCLSTR>

end
