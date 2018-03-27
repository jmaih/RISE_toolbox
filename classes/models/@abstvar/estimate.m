function self=estimate(self,varargin)

% data,date_range,prior,restrictions,
data=[];date_range=[];prior=[];restrictions=[];

n=length(varargin);

if n
    
    data=varargin{1};
    
    if n>1
        
        date_range=varargin{2};
        
        if n>2
            
            prior=varargin{3};
            
            if n>3
                
                restrictions=varargin{4};
                
            end
            
        end
        
    end
    
end

nobj=numel(self);

if nobj>1
    
    for ii=1:nobj
        
        self(ii)=estimate(self(ii),data,date_range,prior,restrictions);
        
    end
    
    return
    
end

% set the prime-time restrictions that are associated with the
% construction of the model
self=prime_time(self);

self=set_inputs(self,'linear_restrictions',restrictions,...
    'data',data,'prior',prior);

[y,x,self.estim_.date_range]=collect_data(self,date_range);

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
    
    self=switching_parameter_estimation(self,y);
    
else
    
    self=constant_parameter_estimation(self,y);

end

end

function self=switching_parameter_estimation(self,y)

prior_trunc=1e-10;

penalty=1000;

[Dummies,epdata]=load_priors();

XX0=self.estim_.X;

YY0=self.estim_.Y;

T0=Dummies.T;

if T0
    
    XX0=[XX0,Dummies.X];
    
    YY0=[YY0,Dummies.Y];
    
end

is_time_varying_trans_prob=self.is_time_varying_trans_prob;

nonlinres=self.estim_.nonlinres;

mapping=self.mapping;

markov_chains=self.markov_chains;
                
linres=self.estim_.linres;

theMap=self.estim_.links.theMap;

objfun=@engine;

% using the linres, transform x0, lb and ub

x0=transform(self.estim_.links.start);

% I am not sure about the transformations on the bounds but in case there
% is an issue, one can always check the bounds after the transformations,
% possibly penalizing the objective function in case of a violation.

lb=transform(self.estim_.links.elb);

ub=transform(self.estim_.links.eub);

options=struct();
options.Display='iter';
options.MaxFunEvals=50000*6;
options.MaxIter=50000;

out=fmincon(objfun,x0,[],[],[],[],lb,ub,[],options);

self.estim_.estim_param=untransform(out);

    function varargout=engine(params0)
        
        Lpost=uminus(1e+8); Lprior=0; Incr=[];
        
        params1=untransform(params0);
        
        % Evaluate prior first and only evaluate the likelihood if prior
        % does not fail.
        [Lprior,retcode]=utils.estim.prior_evaluation_engine(epdata,...
            params1,Lprior);
        
        if ~retcode
            
            M=vartools.estim2states(params1,theMap,mapping.nparams,...
                mapping.nregimes);
            
            pen=utils.estim.penalize_violations2(M,nonlinres,penalty);
            
            [LogLik,Incr,retcode]=vartools.likelihood(M,mapping,...
                YY0,XX0,is_time_varying_trans_prob,...
                markov_chains);
            
            if ~retcode
                
                if T0
                    
                    LogLik0=sum(Incr(1:end-T0));
                    
                    Lprior=Lprior+(LogLik-LogLik0);
                    
                    LogLik=LogLik0;
                    
                end
                
                if ~retcode
                    
                    Lpost=LogLik+Lprior+pen;
                    
                end
                
            end
            
        end
        
        vout={Lpost,Incr,retcode};
        
        varargout=vout(1:nargout);
        
        % negative of likelihood
        varargout{1}=-varargout{1};
        
    end

    function p=transform(p)
        
        if ~isempty(linres)
            
            p=linres.a_to_a2tilde(p);
            
        end
        
    end

    function p=untransform(p)
        % using the linres, untransform params before evaluating the
        % likelihood and later on the prior as well.
        
        if ~isempty(linres)
            
            p=linres.a2tilde_to_a(p);
            
        end
        
    end

    function [Dummies,indPriors]=load_priors()
        
        Dummies=do_var_prior();
        
        indPriors=do_nonvar_priors();
        
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
