function [db,states,retcode] = simulate(obj,varargin)
% simulate - simulates a RISE model
%
% ::
%
%
%   [db,states,retcode] = simulate(obj,varargin)
%
% Args:
%
%    - **obj** [rfvar|dsge|rise|svar]: model object
%
%    - **varargin** : additional arguments including but not restricted to
%
%      - **simul_periods** [integer|{100}]: number of simulation periods
%
%      - **simul_burn** [integer|{100}]: number of burn-in periods. This
%      should not be confused with forecast_conditional_sampling_burnin, which
%      is used in the sampling from the truncated multivariate normal
%      distribution.
%
%      - **simul_historical_data** [ts|struct|{''}]: historical data from
%          which the simulations are based. If empty, the simulations start at
%          the steady state.
%
%      - **simul_history_end_date** [char|integer|serial date]: last date of
%          history
%
%      - **simul_regime** [integer|vector|function handle|{[]}]: regimes for
%          which the model is simulated. When it is a function handle, then it
%          should accept as inputs y (array over all regimes),
%          regimes_1_t_1(regimes from 1 to t-1), sims_1_t_1(the simulated
%          series up to t-1),varargin (possibly further arguments to the
%          function handle). The output is a logical vector that is true for
%          the columns that are acceptable/feasible and false otherwise.
%
%      - **simul_to_time_series** [{true}|false]: if true, the output is a
%          time series, else a cell array with a matrix and information on
%          elements that help reconstruct the time series.
%
%      - **simul_honor_constraints** [true|{false}]: honor restrictions during
%          simulations. If true, agents have to be able anticipate the future.
%
%      - **simul_frwrd_back_shoot** [true|{false}]: uses an algorithm that
%      checks the constraints are satisfied at all horizons instead of just
%      one at a time.
%
%      - **simul_shock_uncertainty** [{true}|false]: draw shocks over the
%      simulation horizon.
%
%      - **simul_honor_constraints_through_switch** [true|{false}]: if true,
%      constraints are honored through the switching mechanism. In that case
%      the number of regimes should be greater than 1. If false, constraints
%      are honored through an anticipatory behavior. In this case, there
%      should be shocks that are foreseen.
%
%      - **simul_anticipate_zero** [true|{false}]: When shocks are drawn, this
%      option allows to impose that agents continue to see only the
%      contemporaneous shocks.
%
%      - **simul_bgp_deviation** [true|{false}]: When the model is
%      nonstationary, a growth component appears in the solution. This option
%      enables or disables that component.
%
%      - **simul_restrictions_lb** [numeric|{-sqrt(eps)}]: effective lower
%      bound for restrictions in simulations.
%
% Returns:
%    :
%
%    - **db** [struct|cell array]: if **simul_to_time_series** is true, the
%      output is a time series, else a cell array with a matrix and
%      information on elements that help reconstruct the time series.
%
%    - **states** [vector]: history of the regimes over the forecast horizon
%
%    - **retcode** [integer]: if 0, the simulation went fine. Else something
%      got wrong. In that case one can understand the problem by running
%      decipher(retcode)
%
% Note:
%
%    - **simul_historical_data** contains the historical data as well as
%      conditional information over the forecast horizon. It may also include
%      as an alternative to **simul_regime**, a time series with name
%      **regime**, which indicates the regimes over the forecast horizon.
%
% Example:
%
%    See also:

if isempty(obj)
    
    mydefaults=the_defaults();
    
    if nargout
        
        db=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

nobj=numel(obj);

if nobj>1
    
    retcode=nan(1,nobj);
    
    states=cell(1,nobj);
    
    db=cell(1,nobj);
    
    for iobj=1:nobj
        
        [db{iobj},states{iobj},retcode(iobj)] = simulate(obj(iobj),varargin{:});
        
    end
    
    db=format_simulated_data_output(db);
    
    return
    
end

[obj,retcode]=solve(obj,varargin{:});
% note that the solving of the model may change the perturbation
% order. More explicitly optimal policy irfs will be computed for a
% perturbation of order 1 no matter what order the user chooses.
% This is because we have not solved the optimal policy problem
% beyond the first order.

db=[]; states=[];

if retcode
    
    return
    
end

% initial conditions
%-------------------
Initcond=set_simulation_initial_conditions(obj);

y0=Initcond.y;

if numel(y0)>1
    
    error('more than one initial conditions')
    
end

% use the steady state with possibly loglinear variables

[y,states,retcode,Qt,myshocks]=utils.forecast.multi_step(y0(1),...
    Initcond.log_var_steady_state,Initcond.T,...
    Initcond.state_vars_location,Initcond);

if retcode
    
    return %error(decipher(retcode))
    
end

% add initial conditions: (only the actual data on the first page)
% ONLY IF THERE WAS NO BURN-IN
%---------------------------------------------------------------
npy=size(y,3);

if Initcond.burn==0
    % we include history
    y=cat(2,y0(1).y(:,:,ones(npy,1)),y);
    
    y0cols=size(y0(1).y,2);
    
else
    % then we start after the end of history
    y0cols=0;
    
end

smply=size(y,2);

% recompose the output if we have a dsge-var
%-------------------------------------------

if Initcond.do_dsge_var
    
    endo_nbr=obj.endogenous.number;
    
    max_rows=obj.observables.number(1);
    
    relevant=obj.inv_order_var(obj.observables.state_id);
    
    yold=y;
    
    y=zeros(endo_nbr,smply)+1i;
    
    y(relevant,:)=yold(1:max_rows,:);
    
end

% exponentiate before doing anything: is_log_var is in the order_var order
%-------------------------------------------------------------------------
if isfield(Initcond,'is_log_var') && ~isempty(Initcond.is_log_var)
    
    y(Initcond.is_log_var,:,:)=exp(y(Initcond.is_log_var,:,:));
    
end

% put y in the correct order before storing
%------------------------------------------
y=re_order_output_rows(obj,y);

y=do_hp_filter(y,obj.options.simul_hpfilter_lambda);

start_date=date2serial(Initcond.simul_history_end_date)-y0cols+1;

nrs=size(states,1);

states=reshape(states,[nrs,1,npy]);

smpl_states=size(states,1);

[states_,markov_chains]=regimes2states(states);

nchains=numel(markov_chains);

[all_Q,all_names_i_j,nshifts]=set_transition_probabilities(Qt,obj.markov_chains);

[exo_nbr,smplx,npx]=size(myshocks);

myshocks=cat(2,zeros(exo_nbr,y0cols,npx),myshocks);

smplx=smplx+y0cols;

vnames=[obj.endogenous.name,obj.exogenous.name,'regime',markov_chains,all_names_i_j];

endo_nbr=obj.endogenous.number;

smpl=max(smplx,smply);

npxy=max(npx,npy);

% expand Q as necessary
% % % all_Q=[nan(nshifts,y0cols),all_Q];

all_Q=all_Q(:,:,ones(1,npxy));

yy=nan(smpl,endo_nbr+exo_nbr+1+nchains+nshifts,npxy);

yy(1:smply,1:endo_nbr,1:npy)=permute(y(1:endo_nbr,:,:),[2,1,3]);

yy(1:smplx,endo_nbr+(1:exo_nbr),1:npx)=permute(myshocks,[2,1,3]);

yy(y0cols+1:smply,endo_nbr+exo_nbr+(1:nchains+1),1:npy)=cat(2,states,states_);

yy(y0cols+1:smply,endo_nbr+exo_nbr+nchains+1+(1:nshifts),:)=permute(all_Q,[2,1,3]);

if obj.options.simul_to_time_series
    % store the simulations in a database: use the date for last observation in
    % history and not first date of forecast
    %--------------------------------------------------------------------------
    % store only the relevant rows in case we are dealing with a VAR with many
    % lags
    db=ts(start_date,yy,vnames);
    
    db=pages2struct(db);
    
else
    
    db={yy,...
        {
        '1=time, start_date=',start_date
        '2= endogenous names',vnames
        }
        };
end

    function [states,markov_chains]=regimes2states(regimes_history)
        
        regimes_tables=obj.markov_chains.regimes;
        
        markov_chains=regimes_tables(1,2:end);
        
        nchains=numel(markov_chains);
        
        reg_table=cell2mat(regimes_tables(2:end,2:end));
        
        states=nan(smpl_states,nchains,npy);
        
        for ip=1:npy
            % there is only one column and so the following should also
            % work max_reg=max(regimes_history(:,ip));
            max_reg=max(regimes_history(:,1,ip));
            
            for ireg=1:max_reg
                
                pos=regimes_history(:,ip)==ireg;
                
                n=sum(pos);
                
                states(pos,:,ip)=reg_table(ireg*ones(n,1),:);
                
            end
            
        end
        
    end

end

function [all_Q,all_names_i_j,nshifts]=set_transition_probabilities(Qt,MC)
nregs=MC.regimes_number;
spanq=size(Qt,3);
bigq=nan(nregs^2,spanq);
regimes_i_j=cell(1,nregs^2);
iter=0;
for i_reg=1:nregs
    for j_reg=1:nregs
        iter=iter+1;
        bigq(iter,:)=squeeze(Qt(i_reg,j_reg,:));
        regimes_i_j{iter}=sprintf('regime_%0.0f_%0.0f',i_reg,j_reg);
    end
end

regimes=MC.regimes;
Journal=MC.journal;
% Legacy
if ischar(Journal{1,1})
    journalize=@(x)eval(['[',x,']']);
    for jj=1:numel(Journal)
        Journal{jj}=journalize(Journal{jj});
    end
end

states_i_j=cell(1,MC.chains_number);
smallq=cell(MC.chains_number,1);
for ic=1:MC.chains_number
    cn=MC.chain_names{ic};
    loc=strcmp(cn,regimes(1,:));
    states_ic=cell2mat(regimes(2:end,loc));
    nstates=numel(unique(states_ic));
    states_i_j_k=cell(1,nstates^2);
    smallq_k=zeros(nstates^2,spanq);
    iter=0;
    for s0=1:nstates
        for r=1:nregs
            if Journal{r,1}(ic)==s0
                break
            end
        end
        Jr=Journal(r,:);
        
        for s1=1:nstates
            iter=iter+1;
            states_i_j_k{iter}=sprintf('%s_%0.0f_%0.0f',cn,s0,s1);
            for ii=1:nregs
                if Jr{ii}(2,ic)==s1
                    this_regime=sprintf('regime_%0.0f_%0.0f',r,ii);
                    bingo=strcmp(this_regime,regimes_i_j);
                    smallq_k(iter,:)=smallq_k(iter,:)+bigq(bingo,:);
                end
            end
        end
    end
    
    smallq{ic}=smallq_k;
    states_i_j{ic}=states_i_j_k;
end
states_i_j=[states_i_j{:}];
smallq=cell2mat(smallq);
% merge them all
%----------------
all_Q=[bigq;smallq];
all_names_i_j=[regimes_i_j,states_i_j];
nshifts=numel(all_names_i_j);
end

function y=do_hp_filter(y,lambda)

if isempty(lambda)
    
    return
    
end

y=permute(y,[2,1,3]);

for ii=1:size(y,3)
    
    y(:,:,ii)=utils.filtering.hpfilter(y(:,:,ii),lambda);
    
end

y=permute(y,[2,1,3]);

end

function sim_data=format_simulated_data_output(sim_data)

nobj=numel(sim_data);

if nobj==1
    
    sim_data=sim_data{1};
    
else
    
    if isempty(sim_data{1})||iscell(sim_data{1})
        
        return
        
    end
    
    sim_data=ts.concatenator(sim_data{:});
    
end

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;

num_vec_fin=@(x)isnumeric(x) && isvector(x) && all(isfinite(x));

num_vec_fin_int=@(x)num_vec_fin(x) && all(floor(x)==ceil(x)) && all(x>=0);

d={
    'simul_periods',100,@(x)num_fin_int(x),...
    'simul_periods must be a finite and positive integer'
    
    'simul_burn',100,@(x)num_fin_int(x),...
    'simul_burn must be a finite and positive integer'
    
    'simul_historical_data',ts.empty(0),@(x)isa(x,'ts')||isstruct(x),...
    'simul_historical_data must be a time series or a structure'
    
    'simul_history_end_date','',@(x)is_date(x)||is_serial(x),...
    'simul_history_end_date must be a valid date'
    
    'simul_regime',[],@(x)isa(x,'function_handle')||...
    (isa(x,'cell')&& isa(x{1},'function_handle'))||...
    (num_vec_fin_int(x) && all(x>=1)),...
    ['simul_regime must be a finite and positive vector of integers',...
    ' OR a function handle OR ',...
    'a cell array whose first element is a function handle']
    
    'simul_bgp_deviation',false,@(x)islogical(x),...
    'simul_bgp_deviation must be a logical'
    
    'simul_anticipate_zero',false,@(x)islogical(x),...
    'simul_anticipate_zero must be a logical'
    
    'simul_honor_constraints_through_switch',false,@(x)islogical(x),...
    'simul_honor_constraints_through_switch must be a logical'
    
    'simul_shock_uncertainty',true,@(x)islogical(x),...
    'simul_shock_uncertainty must be a logical'
    
    'simul_to_time_series',true,@(x)islogical(x),...
    'simul_to_time_series must be a logical'
    
    'simul_frwrd_back_shoot',false,@(x)islogical(x),...
    'simul_frwrd_back_shoot must be a logical'
    
    'simul_honor_constraints',false,@(x)islogical(x),...
    'simul_honor_constraints must be a logical'
    
    'simul_hpfilter_lambda',[],@(x)num_fin(x) && x>0,...
    ' simul_hpfilter_lambda must be >0'
    
    'simul_restrictions_lb',-sqrt(eps),@(x)num_fin(x),...
    ' simul_restrictions_lb must be a scalar and finite number'
    
    };

%     %         'simul_start_date','',... does not seem to be in use
end