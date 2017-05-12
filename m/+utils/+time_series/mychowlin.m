function res=mychowlin(Ylow,Xhigh,s,aggreg_type,estim_method,ngrid)
% PURPOSE: Temporal disaggregation using the Chow-Lin method
% ------------------------------------------------------------
% SYNTAX: res=chowlin(Ylow,Xhigh,s,aggreg_type,estim_method);
% ------------------------------------------------------------
% INPUT: Ylow: Nx1 ---> vector of low frequency data
%        Xhigh: nxp ---> matrix of high frequency indicators (without intercept)
%        aggreg_type: type of disaggregation
%            aggreg_type=1 or 'flow' ---> sum (flow)
%            aggreg_type=2 or 'average' ---> average (index)
%            aggreg_type=3 or 'last'  ---> last element (stock) ---> interpolation
%            aggreg_type=4 or 'first' ---> first element (stock) ---> interpolation
%        s: number of high frequency data points for each low frequency data points 
%        estim_method: estimation method: 
%            estim_method=0 ---> weighted least squares 
%            estim_method=1 ---> maximum likelihood with grid
%            estim_method=2 ---> maximum likelihood with fmincon
% ------------------------------------------------------------
% OUTPUT: res: a structure
%           res.meth    ='Chow-Lin';
%           res.aggreg_type      = type of disaggregation
%           res.estim_method    = method of estimation
%           res.N       = nobs. of low frequency data
%           res.n       = nobs. of high-frequency data
%           res.pred    = number of extrapolations
%           res.s       = frequency conversion between low and high freq.
%           res.p       = number of regressors (including intercept)
%           res.Ylow       = low frequency data
%           res.Xhigh       = high frequency indicators
%           res.y       = high frequency estimate
%           res.y_dt    = high frequency estimate: standard deviation
%           res.y_lo    = high frequency estimate: sd - sigma
%           res.y_up    = high frequency estimate: sd + sigma
%           res.u       = high frequency residuals
%           res.U       = low frequency residuals
%           res.beta    = estimated model parameters
%           res.beta_sd = estimated model parameters: standard deviation
%           res.beta_t  = estimated model parameters: t ratios
%           res.rho     = innovational parameter
%           res.aic     = Information criterion: AIC
%           res.bic     = Information criterion: BIC
%           res.val     = Objective function used by the estimation method
%           res.r       = grid of innovational parameters used by the estimation method
% ------------------------------------------------------------
% SEE ALSO: 
% ------------------------------------------------------------

% Builds on LeSage's Econometrics toolbox

if nargin<6
    
    ngrid=[];
    
    if nargin<5
        
        estim_method=[];
        
        if nargin<4
            
            aggreg_type=[];
            
            if nargin<3
                
                error('s should be provided')
                
            end
            
        end
        
    end
    
end

if isempty(ngrid),ngrid=250; end       % Number of grid points

if isempty(aggreg_type),aggreg_type='average'; end

if isempty(estim_method),estim_method=0; end

% ------------------------------------------------------------
% Size of the problem

N = size(Ylow,1);    % Size of low-frequency input

[n,p] = size(Xhigh);    % Size of p high-frequency inputs (without intercept)

% ------------------------------------------------------------
% Preparing the X matrix: including an intercept

Xhigh=[ones(n,1),Xhigh];       % Expanding the regressor matrix

p=p+1;         % Number of p high-frequency inputs (plus intercept)

% ------------------------------------------------------------
% Generating the aggregation matrix

C = aggreg_matrix();

% -----------------------------------------------------------
% Expanding the aggregation matrix to perform
% extrapolation if needed

pred=max(0,n-s*N); % Number of required extrapolations

C=[C,zeros(N,pred)];

% -----------------------------------------------------------
% Temporal aggregation of the indicators

Xlow=C*Xhigh;

min_r=-0.99;

max_r=-min_r;

if estim_method==2

    rho0=0.5;
    
    [rho,minusFval]=fmincon(@(x)-objective(x),rho0,[],[],[],[],min_r,max_r,[],...
        struct('Display','iter'));
    
    val=-minusFval;
    
    r=rho;

else
    
    % Parameters of grid search
    val=zeros(1,ngrid);
    
    r=linspace(min_r,max_r,ngrid);
    
    % -----------------------------------------------------------
    % Evaluation of the objective function in the grid
    
    for h=1:ngrid
        
        val(h)=objective(r(h));
        
    end
    
    % -----------------------------------------------------------
    % Determination of optimal rho
    
    [~,hmax]=max(val);
    
    rho=r(hmax);

end

% -----------------------------------------------------------
% Final estimation with optimal rho

[Vh,~,iVl,beta,~,sigma_a,U]=one_rho(rho);

L=Vh*C'*iVl;                 % Filtering matrix

u=L*U;                     % High frequency residuals

% -----------------------------------------------------------
% Temporally disaggregated time series

y=Xhigh*beta+u;

% -----------------------------------------------------------
% Information criteria
% Note: p is expanded to include the innovational parameter

aic=log(sigma_a)+2*(p+1)/N;

bic=log(sigma_a)+log(N)*(p+1)/N;

% -----------------------------------------------------------
% VCV matrix of high frequency estimates

sigma_beta=sigma_a*(Xlow'*iVl*Xlow)\eye(p);

VCV_y=sigma_a*(eye(n)-L*C)*Vh+(Xhigh-L*Xlow)*sigma_beta*(Xhigh-L*Xlow)';

d_y=sqrt((diag(VCV_y)));   % Std. dev. of high frequency estimates
y_li=y-d_y;           % Lower lim. of high frequency estimates
y_ls=y+d_y;           % Upper lim. of high frequency estimates

% -----------------------------------------------------------
% -----------------------------------------------------------
% Loading the structure

res.meth='Chow-Lin';

% -----------------------------------------------------------
% Basic parameters 

res.aggreg_type       = aggreg_type;
res.estim_method      = estim_method;
res.N         = N;
res.n         = n;
res.pred      = pred;
res.s         = s;
res.p         = p;

% -----------------------------------------------------------
% Series

res.Ylow      = Ylow;
res.Xhigh     = Xhigh;
res.y         = y;
res.y_dt      = d_y;
res.y_lo      = y_li;
res.y_up      = y_ls;

% -----------------------------------------------------------
% Residuals

res.u         = u;
res.U         = U;

% -----------------------------------------------------------
% Parameters

res.beta      = beta;
res.beta_sd   = sqrt(diag(sigma_beta));
res.beta_t    = beta./sqrt(diag(sigma_beta));
res.rho       = rho;

% -----------------------------------------------------------
% Information criteria

res.aic       = aic;
res.bic       = bic;

% -----------------------------------------------------------
% Objective function

res.val       = val;
res.r         = r;

    function o=objective(r)
        
        [~,Vl,~,~,scp,sigma_ar]=one_rho(r);
        
        switch estim_method
            
            case 0
                
                o=-scp;   % Objective function = Weighted least squares
                
            case 1
                
                % Likelihood function
                loglik=-0.5*N*log(2*pi*sigma_ar)-0.5*log(det(Vl))-0.5*N;
                % loglik=(-N/2)*log(2*pi*sigma_a)-(1/2)*log(det(Vl))-(N/2);
                
                o=loglik;      % Objective function = Likelihood function
                
        end
        
    end

    function [Vh,Vl,iVl,beta,scp,sigma_a,U]=one_rho(rho)
        
        % High frequency VCV matrix (without sigma_a)
        Vh=1/(1-rho^2)*toeplitz(rho.^(0:n-1)); 
        
        % Low frequency VCV matrix (without sigma_a)
        Vl=C*Vh*C';                  
        
        iVl=Vl\eye(N);
        
        % beta estimator
        beta=(Xlow'*iVl*Xlow)\(Xlow'*iVl*Ylow);  
        
        % Low frequency residuals
        U=Ylow-Xlow*beta;                
        
        % Weighted least squares
        scp=U'*iVl*U;               
        
        % sigma_a estimator
        sigma_a=scp/(N-p);         
        
    end

    function [C]=aggreg_matrix()
        
        c = stud();
        
        C=kron(eye(N),c);
        
        function c=stud()
            % Generation of aggregation vector c
            
            switch aggreg_type
                
                case {1,'flow'} % sum
                    
                    c=ones(1,s);
                    
                case {2,'average','index'}
                    
                    c=ones(1,s)/s;
                    
                case {3,'last'}
                    
                    c=zeros(1,s);
                    
                    c(s)=1;
                    
                case {4,'first'}
                    
                    c=zeros(1,s);
                    
                    c(1)=1;
                    
                otherwise
                    
                    error ([' Unknown disaggregation type ',parser.any2str(aggreg_type)]);
                    
            end
            
        end
        
    end

end