function [Euler_errors,outsims]=accuracy(m,nxcuts,nsims,nburn,reset_params)
% Solves dsge models
%
% ::
%
%   Euler_errors=accuracy(m)
%   Euler_errors=accuracy(m,nxcuts)
%   Euler_errors=accuracy(m,nxcuts,nsims)
%   Euler_errors=accuracy(m,nxcuts,nsims,nburn)
%   Euler_errors=accuracy(m,nxcuts,nsims,nburn,reset_params)
%   [Euler_errors,outsims]=accuracy(...)
%
% Args:
%
%    m (rise | dsge): scalar or vector of model objects.
%
%    nxcuts ({10} | integer): Number of quantiles per discretized
%      continuous shock
%
%    nsims ({1000} | integer| 1x2 cell): Number of simulations to generate
%      the state of the endogenous variables. Alternatively, if the input
%      is a cell, no simulations will be performed. In this case nsims is
%      of the form nsims={Y,regs}, where:
%      - Y : simulations of endogenous variables
%      - regs : corresponding simulations for regimes
%
%    nburn ({100} | integer): Number of burned simulations
%
%    reset_params ({true} | false): reset the parameters to their correct
%      values in evaluating the model equations.
%
% Returns:
%    :
%
%    - **Euler_errors** [cell | struct]: if the number of models is greater
%      that one, the output is a cell array of structures. Each structure
%      has the following fields:
%      - EE [1 x nsols cell] : Euler equation errors for all the solutions
%      with each cell for one solution
%      - eqtn1, eqtn2,... : structures containing "mean" (mean of the euler
%      error across all simulations), "max" (maximum of the euler error
%      across all simulations) and "min" (minimum of the euler error across
%      all simulations)
%
%    - **outsims** [1 x 2 x nsols cell]: Simulations for the endogenous
%      variables (first element) and the regimes (second element) for each
%      solution (3rd dimension)
%
% Note:
%
%    - The errors are log10(abs(error))
%
%    - The smaller the standard deviation of the shocks the higher the
%      accuracy even if the solution is not very accurate in the first
%      place.
%
%    - The accuracy of the solution may critically depend on the tolerance
%      level used for solving the model in the first place. For instance,
%      if an iterative algorithm (e.g. mfi) is used, the solution might not
%      be 100% accurate even in a linear model
%
%    - The errors are provided in absolute value for all equations, rather
%      than being normalized as is customarily done in the literature
%
            
nmod=numel(m);

if nmod==0 % isempty(m)
    
    Euler_errors=cell(0,4);
    
    return
    
end

if nargin<5
    
    reset_params=[];
    
    if nargin<4
        
        nburn=[];
        
        if nargin<3
            
            nsims=[];
            
            if nargin<2
                
                nxcuts=[];
                
            end
            
        end
        
    end
    
end


if nmod>1
    
    Euler_errors=cell(1,nmod);
    
    nworkers=utils.parallel.get_number_of_workers();
    
    parfor(imod=1:nmod,nworkers)
        
        Euler_errors{imod}=accuracy(m(imod),nxcuts,nsims,nburn,reset_params);
        
    end
    
    return
    
end

[m,retcode]=solve(m);

if retcode
    
    error('unable to solve the model')
    
end

scurr=rng;

if isempty(reset_params),reset_params=true; end

sims=[];

if isempty(nsims)
    
    nsims=1000;
    
else
    
    if iscell(nsims)
        
        sims=nsims;
        
        nsims=size(sims{1},2);
        
    end
    
end

if isempty(nxcuts),nxcuts=10;end

if isempty(nburn),nburn=100;end

f01=m.routines.probs_times_dynamic;

ev=@utils.code.evaluate_functions;

order=m.options.solve_order;

[T,~,Qfunc,ss_logged_ordered,ov,isstate,log_vars]=load_solution(m,'ov');

iov=m.inv_order_var;

iov_log_vars=log_vars(iov);

[the_leads,the_lags]=...
    dsge_tools.create_endogenous_variables_indices(m.lead_lag_incidence.before_solve);

% get the parameters and definitions of the true model
%------------------------------------------------------
[p,d]=load_true_parameters_and_definitions();

% un-ordered sstate
s=m.solution.ss;

h=numel(s);

% initially all regimes equally likely
%-----------------------------------
pai0=1/h*ones(h,1);

% initial conditions
%---------------------
ylag=s;

for ireg=1:h
    
    ylag{ireg}(iov_log_vars)=log(ylag{ireg}(iov_log_vars));
    
end

y00=cell2mat(ylag)*pai0(:);

nv=m.endogenous.number;

nshks=sum(m.exogenous.number);

nsols=size(T,3);

Euler_errors=struct();

outsims=cell(1,2,nsols);

for isol=1:nsols
    
    [Y,regs,regs_distr]=simulate_data(y00);
    
    [w,shocks,regimes]=discretize_shocks(regs_distr);
    
    if isol==1
        
        nshocks=numel(regimes);
        
        K=nsims*nshocks;
        
        grid=cell(1,K);
        
        EE=zeros(nv,K,nsols);
        
    end
    
    regrid()
    
    for k=1:K
        
        r0=grid{k}.rt;
        
        xk=grid{k}.xt;
        
        ymk=grid{k}.yt;
        
        y0k=compute_forecast(ymk,xk,r0);
        
        def=d{r0};
        
        param=p(:,r0);
        
        ss=s{r0};
        
        for ishk=1:nshocks
            
            eki=0;
            
            xkp=shocks(:,ishk);
            
            for r1=1:h
                
                ypk=compute_forecast(y0k,xkp,r1);
                
                y=set_big_y(ymk,y0k,ypk);
                
                eki=eki+ev(f01,y,xk,ss,param,def,r0,r1);
                
            end
            
            EE(:,k,isol)=EE(:,k,isol)+w(ishk)*eki;
            
        end
        
    end
    
    outsims(:,:,isol)={Y,regs};
    
end

EE=log10(abs(EE));

cmax=zeros(1,nsols);

cmean=zeros(1,nsols);

cmin=zeros(1,nsols);

for iv=1:nv
    
    v=sprintf('eqtn%0.0f',iv);
    
    for iii=1:nsols
        
        cmax(iii)=max(EE(iv,:,iii));
        
        cmean(iii)=mean(EE(iv,:,iii));
        
        cmin(iii)=min(EE(iv,:,iii));
        
        if iv==1
            
            Euler_errors.EE{1,iii}=EE(:,:,iii);
            
        end
        
    end
    
    Euler_errors.(v).max=cmax;
    
    Euler_errors.(v).mean=cmean;
    
    Euler_errors.(v).min=cmin;
    
end

rng(scurr)

    function [w,shocks,regimes]=discretize_shocks(regs_distr)
        
        xw_discrete={1:h,regs_distr};
        
        n=nxcuts*ones(1,nshks);
        
        [x,w] = dsge_tools.quantization(n,xw_discrete);
        
        shocks=x(1:end-1,:);
        
        regimes=x(end,:);
        
    end

    function regrid()
        
        proto=struct('rt',[],'xt',[],'yt',[]);
        
        iter=0;
        
        for ii=1:nsims
            
            proto.yt=Y(:,ii);
            
            for jj=1:nshocks
                
                proto.xt=shocks(:,jj);
                
                iter=iter+1;
                
                proto.rt=regimes(jj);
                
                grid{iter}=proto;
                
            end
            
        end
        
    end

    function [Y,regs,regs_distr]=simulate_data(ylag)
        
        if isempty(sims)
            
            rng('default')
            
            shks=randn(nshks,nsims+nburn);
            
            % do all simulations in logs, exponentiate only when constructing
            % the large y vector
            %---> ylag(iov_log_vars)=log(ylag(iov_log_vars));
            
            [Y,regs]=utils.filtering.sim_engine(ylag(ov),shks,ss_logged_ordered,T,...
                Qfunc,isstate,pai0);
            
            Y=Y(iov,nburn+1:end);
            
            regs=regs(nburn+1:end);
            
            %---> Y(iov_log_vars,:)=exp(Y(iov_log_vars,:));
            
        else
            
            Y=sims{1};
            
            regs=sims{2};
            
        end
        
        % use an empirical distribution for the regimes as a workaround to
        % the fact that the transition probabilities can be time-varying
        %------------------------------------------------------------------
        regs_distr=pai0;
        
        for r=1:h
            
            regs_distr(r)=sum(regs==r)/nsims;
            
        end
        
    end

    function y=set_big_y(ylag,y0,ylead)
        
        [ylag_,y0_,ylead_]=log_varize(ylag,y0,ylead);
        
        y=[ylead_(the_leads);y0_;ylag_(the_lags)];
        
        function varargout=log_varize(varargin)
            
            varargout=varargin;
            
            if ~any(iov_log_vars)
                
                return
                
            end
            
            for ii=1:length(varargin)
                
                varargout{ii}(iov_log_vars)=exp(varargout{ii}(iov_log_vars));
                
            end
            
        end
        
    end

    function y1=compute_forecast(y0,shkt,rt)
        
        y0t=y0(ov)-ss_logged_ordered{rt};
        
        zt=[y0t(isstate);1;shkt];
        
        y1=ss_logged_ordered{rt};
        
        zkron=zt;
        
        for io=1:order
            
            y1=y1+T{io,rt,isol}*zkron;
            
            if io<order
                
                zkron=kron(zkron,zt);
                
            end
            
        end
        
        y1=y1(iov);
        
    end

    function [p,d]=load_true_parameters_and_definitions()
        
        m0=m;
        
        if reset_params
            
            m0=set(m0,'solve_perturbation_type','maih');
            
        end
        
        [defs,retcode_]=compute_definitions(m0);
        
        if retcode_
            
            error('unable to recompute definitions')
            
        end
        
        p=m0.parameter_values;
        
        d=defs; % m.solution.definitions;
        
    end

end