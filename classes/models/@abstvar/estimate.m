function self=estimate(self,varargin)
% Estimates VAR parameters using frequentist/Bayesian interpretation
%
% ::
%
%    var = estimate(var);
%    var = estimate(var, data);
%    var = estimate(var, data, data_range);
%    var = estimate(var, data, data_range, prior);
%    var = estimate(var, data, data_range, prior, restrictions);
%    var = estimate(var, data, data_range, prior, restrictions,optimizer);
%    var = estimate(var, data, data_range, prior, restrictions,optimizer,is_fixed_regime);
%
% Args:
%
%    var (var object): var object
%
%    data (struct or ts obect): data for estimation
%
%    data_range (serial): (optional) date_range
%
%    prior (string): priors for parameters
%
%       - If no prior is given, maximum likehood estimators (frequentist)
%       are set for parameters 
%       - With priors, posterior mode (Bayesian) values are set for parameters
%
%    restrictions : restrictions. Refer to XXXXXXXX for analysis
%
%    optimizer (char|function_handle|cell|{fmincon}) : optimization
%       procedure. Used with optimization is required. e.g. markov
%       switching. This can be the name of a standard matlab optimizer or
%       RISE optimization routine or a user-defined optimization procedure
%       available of the matlab search path. If the optimzer is provided as
%       a cell, then the first element of the cell is the name of the
%       optimizer or its handle and the remaining entries in the cell are
%       additional input arguments to the user-defined optimization
%       routine. A user-defined optimization function should have the
%       following syntax :: 
%
%            [xfinal,ffinal,exitflag,H]=optimizer(fh,x0,lb,ub,options,varargin);
%
%         That is, it accepts as inputs:
%
%             - **fh**: the function to optimize
%             - **x0**: a vector column of initial values of the parameters
%             - **lb**: a vector column of lower bounds
%             - **ub**: a vector column of upper bounds
%             - **options**: a structure of options whose fields will be similar
%               to matlab's optimset
%             - **varargin**: additional arguments to the user-defined
%               optimization procedure
%
%         That is, it provides as outputs:
%
%             - **xfinal**: the vector of final values
%             - **ffinal**: the value of **fh** at **xfinal**
%             - **exitflag**: a flag similar to the ones provided by matlab's
%               optimization functions.
%             - **H**: an estimate of the Hessian
%
%    is_fixed_regime : (true|{false}): if true, the regimes are known in
%       advance. In that case the data should contain time series for a
%       variable called "hist_regimes"
%
% Returns:
%    : var object with parameters estimated based on data
%
% See also:
%    - :func:`posterior_mode <var.posterior_mode>`
%    - :func:`identification <var.identification>`
%

data=[];date_range=[];prior=[];restrictions=[];optimizer=[]; is_fixed_regime=[];

n=length(varargin);

if n
    
    data=varargin{1};
    
    if n>1
        
        date_range=varargin{2};
        
        if n>2
            
            prior=varargin{3};
            
            if n>3
                
                restrictions=varargin{4};
                
                if n>4
                    
                    optimizer=varargin{5};
                    
                    if n>5
                        
                        is_fixed_regime=varargin{6};
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

nobj=numel(self);

if nobj>1
    
    for ii=1:nobj
        
        self(ii)=estimate(self(ii),data,date_range,prior,restrictions,optimizer,is_fixed_regime);
        
    end
    
    return
    
end

% set the prime-time restrictions that are associated with the
% construction of the model
self=prime_time(self);

self=set_inputs(self,'linear_restrictions',restrictions,...
    'data',data,'prior',prior);

[y,x,self.estim_.date_range,fixed_regimes]=collect_data(self,date_range,is_fixed_regime);

% N.B: yx does not include the constant while x may include it.
self = abstvar.embed(self,y,x);

% self.mapping
links=struct();

[links.estimList,links.start,links.elb,links.eub,links.theMap,...
    links.transProbs]=abstvar.map_estimation(self.markov_chains,...
    self.mapping.regimes,self.param_guide);

self.estim_.links=links;

% process the restrictions above + additional ones set by the user
self=process_linear_restrictions(self);

if self.is_switching||self.optimize
    
    self=switching_parameter_estimation(self,y,optimizer,fixed_regimes);
    
else
    
    self=constant_parameter_estimation(self,y);
    
end

end

function self=switching_parameter_estimation(self,y,optimizer,fixed_regimes)

prior_trunc=1e-10;

penalty=1000;

if ~isempty(fixed_regimes)
    % lose nlags
    fixed_regimes=fixed_regimes(self.nlags+1:end);
    
end

[Dummies,epdata,endopri]=load_priors();

is_endo_priors=~isempty(endopri);

is_filt_required=false;

filt=[];

if is_endo_priors
    
    is_filt_required=endopri.is_filt_required;
    
end

XX0=self.estim_.X;

YY0=self.estim_.Y;

T0=Dummies.T;

if T0
    
    XX0=[XX0,Dummies.X];
    
    YY0=[YY0,Dummies.Y];
    
    if ~isempty(fixed_regimes)
        
        fixed_regimes=[fixed_regimes(:).',fixed_regimes(end)*ones(1,T0)];
        
    end
    
end

is_time_varying_trans_prob=self.is_time_varying_trans_prob;

nonlinres=self.estim_.nonlinres;

mapping=self.mapping;

markov_chains=self.markov_chains;

% linres=self.estim_.linres;

theMap=self.estim_.links.theMap;

objfun=@engine;

% using the linres, transform x0, lb and ub

x0=transform(self,self.estim_.links.start);

% I am not sure about the transformations on the bounds but in case there
% is an issue, one can always check the bounds after the transformations,
% possibly penalizing the objective function in case of a violation.

lb=transform(self,self.estim_.links.elb);

ub=transform(self,self.estim_.links.eub);

estim_blocks=[];

options=struct();
options.Display='iter';
options.MaxFunEvals=50000*6;
options.MaxIter=50000;

PROBLEM_=struct('objective',objfun,...
    'x0',x0,...
    'lb',lb,...
    'ub',ub,...
    'nonlcon',[],... % the nonlinear constraints restrictions take the same inputs as fh_wrapper
    'options',options,...obj(1).options.optimset
    'solver',optimizer);

wd=utils.estim.warnings_disable();

[x1,f1,H]=optimization.estimation_engine(PROBLEM_,estim_blocks); %#ok<ASGLU>

if any(isnan(H(:)))
    
    H=utils.numdiff.hessian(@engine,x1);
    
end

self.estim_.vcov=H\eye(numel(x1));

self.estim_.estim_param=untransform(self,x1);

self.estim_.objfun=objfun;

utils.estim.warnings_enable(wd)

    function varargout=engine(params0)
        
        Lpost=uminus(1e+8); Lprior=[0,0]; Incr=[];
        
        params1=untransform(self,params0);
        
        % Evaluate prior first and only evaluate the likelihood if prior
        % does not fail.
        [Lprior(1),retcode]=utils.estim.prior_evaluation_engine(epdata,...
            params1,Lprior(1));
        
        if ~retcode
            
            M=vartools.estim2states(params1,theMap,mapping.nparams,...
                mapping.nregimes);
            
            pen=utils.estim.penalize_violations2(M,nonlinres,penalty);
            
            if is_filt_required
                
                [LogLik,Incr,retcode,filt]=vartools.likelihood(M,mapping,...
                    YY0,XX0,is_time_varying_trans_prob,...
                    markov_chains,fixed_regimes);
                
            else
                
                [LogLik,Incr,retcode]=vartools.likelihood(M,mapping,...
                    YY0,XX0,is_time_varying_trans_prob,...
                    markov_chains,fixed_regimes);
                
            end
            
            if ~retcode
                
                if T0
                    
                    LogLik0=sum(Incr(1:end-T0));
                    
                    Lprior(1)=Lprior(1)+(LogLik-LogLik0);
                    
                    LogLik=LogLik0;
                    
                end
                
                if is_endo_priors
                    % do the uncorrelated guys
                    %--------------------------
                    pp=endopri.estim_endogenous_priors(self,filt);
                    
                    [Lprior(2),retcode]=utils.estim.prior.evaluate_uncorrelated(...
                        endopri.estim_endogenous_priors_data,pp,Lprior(2));
                    
                    if retcode
                        
                        retcode=309;
                        
                    end
                    
                end
                
                if ~retcode
                    
                    Lpost=LogLik+sum(Lprior)+pen;
                    
                end
                
            end
            
        end
        
        vout={Lpost,Incr,retcode};
        
        varargout=vout(1:nargout);
        
        % negative of likelihood
        varargout{1}=-varargout{1};
        
    end

    function [Dummies,indPriors,endopri]=load_priors()
        
        Dummies=do_var_prior();
        
        indPriors=do_nonvar_priors();
        
        endopri=do_endogenous_priors();
        
        function endopri=do_endogenous_priors()
            
            endopri=[];
            
            if ~isfield(self.estim_.prior,'endogenous')||...
                    isempty(self.estim_.prior.endogenous)
                
                return
                
            end
            
            endopri=struct();
            
            fh=self.estim_.prior.endogenous;
            
            test=fh();
            
            endopri.is_filt_required=test.kf_filtering_level>0;
            
            [endopri.estim_endogenous_priors_data,...
                endopri.endogenous_priors,endopri.estim_endogenous_priors]...
                =utils.prior.setup_endogenous_priors_engine(prior_trunc,fh);
            
        end
        
        function newpri=do_nonvar_priors()
            
            estimList=self.estim_.links.estimList;
            
            nest=numel(estimList);
            
            processed=true(1,nest);
            
            nonvar_id=locate_variables(self.nonvar_parameters,estimList);
            
            processed(nonvar_id)=false;
            
            newpri=[];
            
            pri=[];
            
            d=[]; shortcut_d=[];
            
            if ~isfield(self.estim_.prior,'nonvar')||...
                    isempty(self.estim_.prior.nonvar)
                
                return
                
            end
            
            nonvar=self.estim_.prior.nonvar;
            
            estnames=fieldnames(nonvar);
            
            is_dirichlet=strncmp(estnames,'dirichlet',9);
            
            for est_id=1:numel(estnames)
                
                fildname=estnames{est_id};
                
                position=find(strcmp(fildname,estimList));
                
                tmp=nonvar.(fildname);
                
                if isempty(position)
                    
                    if ~is_dirichlet(est_id)
                        
                        error(['unknown parameter "',fildname,'"'])
                        
                    end
                    
                    [d,shortcut_d]=utils.estim.format_dirichlet(d,tmp,...
                        estimList,shortcut_d);
                    
                    lastlocs=d(end).location;
                    
                    if any(processed(lastlocs))
                        
                        disp(estimList(lastlocs))
                        
                        error('the parameters above appear to be set multiple times')
                        
                    end
                    
                    self.estim_.links.elb(lastlocs)=0;
                    
                    self.estim_.links.start(lastlocs)=...
                        d(end).moments.mean(d(end).pointers);
                    
                    self.estim_.links.eub(lastlocs)=1;
                    
                    processed(lastlocs)=true;
                    
                else
                    
                    if processed(position)
                        % this cannot happen since it is a field name but
                        % if the user sets a VAR parameter, this is where
                        % he will be caught
                        
                        error(['parameter "',fildname,'" is duplicated'])
                        
                    end
                    
                    block=utils.prior.cell2block(tmp,fildname,'const',1,'');
                    
                    pri=utils.prior.prior_setting_engine(pri,block,est_id,...
                        prior_trunc);
                    
                    self.estim_.links.elb(position)=pri(end).lower_bound;
                    
                    self.estim_.links.start(position)=pri(end).start;
                    
                    self.estim_.links.eub(position)=pri(end).upper_bound;
                    
                    processed(position)=true;
                    
                end
                
            end
            
            nonprocessed_params=estimList(~processed);
            
            if ~isempty(nonprocessed_params)
                
                disp(nonprocessed_params)
                
                error('Missing priors for the parameters above')
                
            end
            % Not all parameters are estimated and so the location may be
            % off, it has to be corrected
            
            if ~isempty(pri)
                
                prilocs=locate_variables({pri.name},self.estim_.links.estimList);
                
                newpri=utils.prior.load_priors(struct(),pri,shortcut_d,prilocs);
                
            end
            
        end
        
        function D=do_var_prior()
            
            D=struct('T',0);
            
            if ~isfield(self.estim_.prior,'var')||...
                    isempty(self.estim_.prior.var)
                
                return
                
            end
            
            var_prior=self.estim_.prior.var;
            
            switch var_prior.type
                
                case {'sims-zha','sz'}
                    
                    sig=std(y,[],2);
                    
                    ybar=mean(y,2);
                    
                    [D.Y,D.X]=vartools.sims_zha_dummies(var_prior,...
                        self.nvars*self.ng,self.nx*self.ng,self.nlags,...
                        sig,ybar);
                    
                    D.T=size(D.Y,2);
                    
                case 'minnesota'
                    
                    error(['prior "',var_prior.type,'" not ready']);
                    
                case {'normal-wishart','nw'}
                    
                    error(['prior "',var_prior.type,'" not ready']);
                    
                case {'indep-normal-wishart','inw'}
                    
                    error(['prior "',var_prior.type,'" not ready']);
                    
                case 'jeffrey'
                    
                    error(['prior "',var_prior.type,'" not ready']);
                    
                otherwise
                    
                    error(['unknown prior type ',var_prior.type])
                    
            end
            
        end
        
    end

end

function self=constant_parameter_estimation(self,y)

warmup= vartools.ols(self);

if ~isfield(self.estim_.prior,'var')||isempty(self.estim_.prior.var)
    
    B = warmup.B;
    
    SIG = warmup.Sigma;
    
else
    
    var_prior=self.estim_.prior.var;
    
    switch var_prior.type
        
        case 'minnesota'
            
            [abar,SIG,self.estim_.sampler]=vartools.minnesota_prior(...
                self,y,warmup.Sigma,var_prior);
            
            B = reshape(abar,self.nvars,[]);
            
        case {'normal-wishart','nw'}
            
            [abar,SIG,self.estim_.sampler]=vartools.normal_wishart_prior(...
                self,y,warmup.Sigma,var_prior);
            
            B = reshape(abar,self.nvars,[]);
            
        case {'indep-normal-wishart','inw'}
            
            [abar,SIG,self.estim_.sampler]=vartools.independent_normal_wishart_prior(...
                self,y,warmup.Sigma,var_prior);
            
            B = reshape(abar,self.nvars,[]);
            
        case {'sims-zha','sz'}
            
            [abar,SIG,self.estim_.sampler]=vartools.sims_zha_prior(...
                self,y,warmup.Sigma,var_prior);
            
            B = reshape(abar,self.nvars,[]);
            
        case 'jeffrey'
            
            error('jeffrey priors not yet implemented')
            
        otherwise
            
            error(['unknown prior type ',var_prior.type])
            
    end
    
end

self.estim_.estim_param=[B(:);vech(SIG)];

end
