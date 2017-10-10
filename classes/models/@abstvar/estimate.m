function self=estimate(self,date_range,varargin)

if nargin < 2
    
    date_range=[];
    
end

self=set(self,varargin{:});%self=setOptions(self,varargin{:});

% set the prime-time restrictions that are associated with the
% construction of the model
self=prime_time(self);

[y,x,self.date_range]=collect_data(self,date_range);

self = abstvar.embed(self,y,x);

% process the restrictions above + additional ones set by the user
self=process_linear_restrictions(self);

if self.is_switching||self.optimize
    
    self=switching_parameter_estimation(self,y);
    
else
    
    self=constant_parameter_estimation(self,y);

end

end

function self=switching_parameter_estimation(self,y)

[kdata,Tprior]=load_priors();

linres=self.linres;

objfun=@engine;

% using the linres, transform x0, lb and ub

x0=transform(self.mapping.start);

% I am not sure about the transformations on the bounds but in case there
% is an issue, one can always check the bounds after the transformations,
% possibly penalizing the objective function in case of a violation.

lb=transform(self.mapping.elb);

ub=transform(self.mapping.eub);

options=struct();
options.Display='iter';
options.MaxFunEvals=50000;

out=fmincon(objfun,x0,[],[],[],[],lb,ub,[],options);

self.estim_param=untransform(out);

    function varargout=engine(params0)
        
        params1=untransform(params0);
        
        [varargout{1:nargout}]=vartools.likelihood(params1,...
            kdata.mapping,kdata,kdata.markov_chains);
        
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

    function [kdata,Tprior]=load_priors()
        
        kdata=self;
        
        if isempty(self.prior)
            
            Tprior=0;
            
            return
            
        end
        
        switch self.prior.type
            
            case 'minnesota'
                
                error(['prior "',self.prior.type,'" not ready']);
                
            case {'normal-wishart','nw'}
                
                error(['prior "',self.prior.type,'" not ready']);
                
            case {'indep-normal-wishart','inw'}
                
                error(['prior "',self.prior.type,'" not ready']);
                
            case {'sims-zha','sz'}
                
                sig=std(y,[],2);
                
                ybar=mean(y,2);
                
                [Y,X]=vartools.sims_zha_dummies(kdata.prior,...
                    kdata.nvars*kdata.ng,kdata.nx*kdata.ng,kdata.nlags,...
                    sig,ybar);
                
                Tprior=size(Y,2);
                
                kdata.X=[kdata.X,X];
                
                kdata.Y=[kdata.Y,Y];
                
                kdata.T=size(kdata.Y,2);
                
            case 'jeffrey'
                
                error(['prior "',self.prior.type,'" not ready']);
                
            otherwise
                
                error(['unknown prior type ',self.prior.type])
                
        end
        
    end

end

function self=constant_parameter_estimation(self,y)

warmup= vartools.ols(self);

if isempty(self.prior)
    
    B = warmup.B;
    
    SIG = warmup.Sigma;
    
else
    
    switch self.prior.type
        
        case 'minnesota'
            
            [abar,SIG,self.sampler]=vartools.minnesota_prior(...
                self,y,warmup.Sigma,self.prior);
            
            B = reshape(abar,self.nvars,[]);
            
        case {'normal-wishart','nw'}
            
            [abar,SIG,self.sampler]=vartools.normal_wishart_prior(...
                self,y,warmup.Sigma,self.prior);
            
            B = reshape(abar,self.nvars,[]);
            
        case {'indep-normal-wishart','inw'}
            
            [abar,SIG,self.sampler]=vartools.independent_normal_wishart_prior(...
                self,y,warmup.Sigma,self.prior);
            
            B = reshape(abar,self.nvars,[]);
            
        case {'sims-zha','sz'}
            
            [abar,SIG,self.sampler]=vartools.sims_zha_prior(...
                self,y,warmup.Sigma,self.prior);
            
            B = reshape(abar,self.nvars,[]);
            
        case 'jeffrey'
            
            error('jeffrey priors not yet implemented')
            
        otherwise
            
            error(['unknown prior type ',self.prior.type])
            
    end
    
end

self.estim_param=[B(:);vech(SIG)];

end
