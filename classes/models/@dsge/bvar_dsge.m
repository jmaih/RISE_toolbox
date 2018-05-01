function [obj,retcode]=bvar_dsge(obj,varargin)
% bvar_dsge -- intermediary file for computing key elements for dsge-var
%
% ::
%
%
%   [obj,retcode]=bvar_dsge(obj,varargin)
%
% Args:
%
%    - **obj** [rise|dsge]: model object
%
%    - **varargin** [?]: ususal dsge options. The most important of which are:
%      - **dsgevar_lag** [integer|{4}]: number of lags in the VAR
%      - **dsgevar_constant** [false|{true}]: flag for having a constant in
%      the VAR
%      - **dsgevar_var_regime** [false|{true}]: use the VAR in simulations.
%      Otherwise use the DSGE
%      - **dsgevar_inner_param_uncertainty** [true|{false}]: make random draws
%      for the parameters around the mode for each simulation
%
% Returns:
%    :
%
%    - **obj** [rise|dsge]: model object
%
% Note:
%
%    - Because the BVAR-DSGE will have fewer variables than the DSGE, in
%    simulations, the missing variables will be stored as 0+1i in time series.
%    This is true for forecast, irf and simulate
%
% Example:
%
%    See also:

if isempty(obj)
    
    mydefaults=the_defaults();
    
    if nargout
        
        obj=mydefaults;
        
    else
        
        clear obj
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

if ~obj.is_dsge_var_model
    
    error('this function can only be used in the presence of a DSGE-VAR model')

end

if ~isempty(varargin)
    
    obj=set(obj,varargin{:});

end

if obj.markov_chains.regimes_number>1
    
    error([mfilename,':: dsge-var for markov switching not implemented yet'])

end

if max(obj.exogenous.shock_horizon(:))>0

    error([mfilename,':: dsge-var with anticipations not implemented yet'])

end

% initialize those and return if there is a problem
% note we take the negative of the penalty to maximize

dsge_var=create_dsge_var_tank(obj);

% load the elements computed in load_data, using the Schorfheide notation
%-------------------------------------------------------------------------
p=dsge_var.p;% the var order

T=dsge_var.T;% the sample size

const=dsge_var.constant; % flag for constant

n=dsge_var.n; % number of variables

k=const+n*p;

% theoretical autocovariances: impose that 
% 1- the model is not resolved otherwise we enter an infinite loop
%--------------------------------------------------------------------
[A,retcode]=theoretical_autocovariances(obj,'autocov_ar',p,...
    'autocov_model_resolve',false);

if ~retcode
    % VAR approximation of the DSGE model
    %-------------------------------------
    ids=obj.observables.state_id;
    
    steady_state=obj.solution.ss{1}(ids);
    
    [PHI_theta,SIG_theta,GXX,GYX,GXY,GYY,retcode]=var_approximation_to_the_dsge(steady_state,A(ids,ids,:),const);
    
    if ~retcode
        % the prior weight is given by the dsge model
        %---------------------------------------------
        lambda=obj.parameter_values(obj.dsge_prior_weight_id,1);
        
        if lambda*T>k+n
            % load the empirical moment matrices
            %-----------------------------------
            YY=dsge_var.YY;
            
            YX=dsge_var.YX;
            
            XX=dsge_var.XX;
            
            XY=dsge_var.XY;
            
            % the resulting Bayesian VAR combines prior(dsge) and the data
            % through the empirical moments
            [PHIb,SIGb,ltgxx,ltgxxi]=bvar_dsge_mode();
            
        else
            
            retcode=261;
            
        end
        
        store_dsge_var();
        
    end
    
end

    function store_dsge_var()

        if retcode
        
            dsge_var=[];
        
        else
            
            dsge_var.var_approx=struct('PHI',PHI_theta,...
                'SIG',SIG_theta,...
                'GXX',GXX);
            
            dsge_var.posterior.PHI=PHIb;
            
            dsge_var.posterior.SIG=SIGb;
            
            dsge_var.posterior.ZZi=ltgxxi;
            
            dsge_var.posterior.ZZ=ltgxx;
            
            dsge_var.posterior.inverse_wishart.df=fix((1+lambda)*T-k);
            
            dsge_var.lambda=lambda;
        
        end
        
        obj.dsge_var=dsge_var;
    
    end

    function [PHIb,SIGb,ltgxx,ltgxxi]=bvar_dsge_mode()
    
        if isinf(lambda)
        
            SIGb = SIG_theta;
            
            PHIb = PHI_theta;
            
            ltgxx=[];
            
            ltgxxi=[];
       
        else
            
            PHIb=(lambda/(1+lambda)*GXX+1/(1+lambda)*XX/T)\...
                (lambda/(1+lambda)*GXY+1/(1+lambda)*XY/T);
            
            ltgyx=lambda*T*GYX+YX;
        
            ltgxx=lambda*T*GXX+XX;
            
            ltgxxi=ltgxx\eye(size(ltgxx));
            
            SIGb=lambda*T*GYY+YY-ltgyx*ltgxxi*ltgyx'; % <---  SIGb=lambda*T*GYY+YY-ltgyx*(ltgxx\ltgyx');
            
            SIGb=SIGb/((1+lambda)*T);
        
        end
        
    end

    function [PHI,SIG,GXX,GYX,GXY,GYY,rcode]=var_approximation_to_the_dsge(varobs_steady,varobs_autocov,const)
        
        rcode=0;
        
        PHI=[];SIG=[];
        
        GYY=varobs_autocov(:,:,1);
        
        const=any(varobs_steady~=0)||const;
        
        GXX=nan(k);
        
        GYX=nan(n,k);
        
        if const
        
            GXX(:,1)=[1;repmat(varobs_steady,p,1)];
            
            GXX(1,2:end)=transpose(GXX(2:end,1));
            
            GYX(:,1)=varobs_steady;
        
        end
        
        for ii=1:p
        
            rr=const+((ii-1)*n+1:ii*n);
            
            for jj=ii:p
            
                cc=const+((jj-1)*n+1:jj*n);
                
                if ii==jj
                
                    GXX(rr,cc)=GYY;
                
                else
                    
                    GXX(rr,cc)=varobs_autocov(:,:,jj-ii+1);
                
                end
                
            end
            
            GYX(:,rr)=varobs_autocov(:,:,ii+1);
        
        end
        
        GXX=triu(GXX);
        
        GXX=GXX+triu(GXX,1)';
        
        GXY=transpose(GYX);
        
        Gxxi=GXX\eye(k);
        
        if any(any(isnan(Gxxi)))
        
            rcode=26;
        
        else
            
            PHI=Gxxi*GXY;
            
            SIG=GYY-GYX*Gxxi*GXY;
        
        end
        
    end

end

function dsge_var=create_dsge_var_tank(obj)

data=transpose(obj.data.y(:,obj.data.start:obj.data.finish));

const=obj.options.dsgevar_constant;

n=obj.observables.number(1); % endogenous observables

p=obj.options.dsgevar_lag;

Y=data(p+1:end,:);

smpl=size(Y,1);

X=nan(smpl,const+n*p);

if const
    
    X(:,1)=1;
    
end

for ii=1:p
    
    X(:,const+((ii-1)*n+1:ii*n))=data(p+1-ii:end-ii,:);
    
end

dsge_var=struct('YY',Y'*Y,'YX',Y'*X,...
    'XX',X'*X,'XY',X'*Y,'T',smpl,...
    'n',n,'p',p,'constant',const);

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;

d={'dsgevar_lag',4,@(x)num_fin_int(x) && x>=1,...
    'dsgevar_lag must be an integer >=1'% # lags
    
    'dsgevar_constant',true,@(x)islogical(x),...
    'dsgevar_constant must be true or false' % VAR admits constant
    
    'dsgevar_var_regime',true,@(x)islogical(x),...
    'dsgevar_var_regime must be true or false'% use the var for irf, forecasting and simulation
    
    'dsgevar_inner_param_uncertainty',false,@(x)islogical(x),...
    'dsgevar_inner_param_uncertainty must be true or false'
    }; %

end
