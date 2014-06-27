function [db,State] = simulate(obj,varargin)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% the random numbers specifications and algorithms are specified outside
% this function. Need to check this again!!!
if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    db=struct('simul_periods',100,...
        'simul_burn',100,...
        'simul_algo','mt19937ar',...  % [{mt19937ar}| mcg16807| mlfg6331_64|mrg32k3a|shr3cong|swb2712]
        'simul_seed',0,...
        'simul_histval',[],...rise_time_series.empty(0)
        'simul_start_date','',...
        'simul_regime',[],...
        'simul_sig',1); % perturbation
    return
end

nobj=numel(obj);
if nobj>1
    State=cell(1,nobj);
    db=cell(1,nobj);
    for iobj=1:nobj
        [db{iobj},State{iobj}] = simulate(obj(iobj),varargin{:});
    end
    db=format_simulated_data_output(db);
    return
end
obj=set_options(obj,varargin{:});
simul_burn=obj.options.simul_burn;
simul_periods=obj.options.simul_periods;
simul_histval=obj.options.simul_histval;
simul_regime=obj.options.simul_regime;
simul_sig=obj.options.simul_sig;

% solution and system
%--------------------
[obj,retcode]=solve(obj);
if retcode
    error([mfilename,':: no solution: the model cannot be simulated'])
end
Q={obj.solution.Q,[],[]};
if  obj.is_endogenous_switching_model
    % take a handle to a private function. we can't access it otherwise
    Q{2}=obj.func_handles.transition_matrix;
    M=obj.parameter_values;
    % transition matrix should be invariant. And so, hopefully, the first
    % argument to Q{2} could be the updated data in any state.
    % order of the input arguments is y,x,param,ss
    Q{3}={[],mean(M,2),obj.solution.ss{1}}; % remaining arguments of Q{2} after the first one
end
exo_nbr=sum(obj.exogenous.number);
regimes_number=obj.markov_chains.regimes_number;
horizon=obj.options.solve_expect_order;

% measurement errors
%-------------------
measurement_errors_flag= ~isempty(obj.solution.H{1});
if measurement_errors_flag
    varobs_id= obj.observables.state_id;
    CS=nan(size(obj.solution.H{1}));
    for ii=1:regimes_number
        CS(:,:,ii)=diag(sqrt(diag(obj.solution.H{ii})));
    end
end

% initial conditions
%-------------------
endo_nbr=obj.endogenous.number(2);
PAI=1/regimes_number*ones(regimes_number,1); % ErgodicDistribution(Q)
if isempty(simul_regime)
    y0=0;
    % we start at the steady state
    for ireg=1:regimes_number
        y0=y0+obj.solution.ss{ireg};
    end
    y0=y0/regimes_number;
else
    y0=obj.solution.ss{simul_regime};
end

simul_start_date=0;
if ~isempty(simul_histval)
    if isstruct(simul_histval)
        simul_histval=rise_time_series.collect(simul_histval);
    end
    if ~isa(simul_histval,'rise_time_series')
        error('historical database must be a rise_time_series')
    end
    simul_start_date=obj.options.simul_start_date;
    if isempty(simul_start_date)
        simul_start_date=simul_histval.finish;
    end
    dbnames=simul_histval.varnames;
    dbdata=simul_histval(simul_start_date);
    db_locs=locate_variables(dbnames,obj.endogenous.name,true);
    for ivar=1:numel(dbnames)
        if ~isnan(db_locs(ivar))
            y0(db_locs(ivar))=dbdata(ivar);
        end
    end
    % for "serious" simulation, we do not burn anything. We consider how
    % the economy will evolve from the initial conditions.
    simul_burn=0;
end

State=nan(simul_periods+1,1);

for t=1:simul_burn+simul_periods
    % draw a state
    %-------------
    [regime,Q,PAI,retcode]=choose_state(simul_regime,Q,PAI,y0);
    if retcode
        return
    end
    % Simulate one step
    %------------------
    shocks=randn(exo_nbr,horizon);
    y1=simulation_engine(obj.solution,y0,shocks,simul_sig,regime,obj.options.solve_order,obj.options.irf_anticipate);
    % add the measurement errors
    %---------------------------
    if measurement_errors_flag
        y1(varobs_id)=y1(varobs_id)+CS(:,:,regime)*randn(obj.observables.number(1),1);
    end
    if t>simul_burn
        if t==simul_burn+1
            y=[y0,zeros(endo_nbr,simul_periods)];
        end
        % update probabilities
        y(:,t-simul_burn+1)=y1;
        State(t-simul_burn+1)=regime;
    end
    y0=y1;
end

% store the simulations in a database
%------------------------------------
db=rise_time_series(simul_start_date,y',obj.endogenous.name);
db=pages2struct(db);
end

function [st,Q,PAI,retcode]=choose_state(st,Q,PAI,y)
% st: state
% Q: transition matrix or ...
% PAI: current updated probabilities
% y: data for current period
retcode=0;
if isempty(st)
    endogenous_switching=~isempty(Q{2});
    Q0=Q{1};
    % update probabilities
    %---------------------
    if endogenous_switching
        % in the endogenous probability case the configuration of
        % the transition matrix will change
        shadow_transition_matrix=Q{2};
        Vargs=Q{3};
        [Q0,retcode]=online_function_evaluator(shadow_transition_matrix,y,Vargs{:},[],[],[]);
    end
    if retcode
        return
    end
    PAI=Q0'*PAI;
    csp=[0;cumsum(PAI)];
    st=find(csp>rand,1,'first')-1;
end
end


function sim_data=format_simulated_data_output(sim_data)
nobj=numel(sim_data);
if nobj==1
    sim_data=sim_data{1};
else
    if isempty(sim_data{1})
        return
    end
    sim_data=concatenate_series_from_different_models(sim_data);
end
end
