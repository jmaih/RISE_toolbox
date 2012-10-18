function [obj,State] = simulate(obj,varargin)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
% I could pre-allocate the simulation array and set initial conditions as
% an alternative to this simple setup. And if there are initial conditions
% there is no burn-in.
% Should also add a control for the parameter choice: mode, mean, calibration, etc.

Defaults=struct('simul_periods',100,...
    'simul_burn',100,...
	'simul_algo','mt19937ar',...  % [{mt19937ar}| mcg16807| mlfg6331_64|mrg32k3a|shr3cong|swb2712]
	'simul_seed',0);
% the random numbers specifications and algorithms are specified outside
% this function
if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    obj=Defaults;
    return
end

obj.options=mysetfield(obj.options,varargin{:});

simul_burn=obj.options.simul_burn;
simul_periods=obj.options.simul_periods;

[obj,retcode]=solve(obj);
if retcode
    error([mfilename,':: no solution: the model cannot be simulated'])
end
T=obj.T;
R=obj.R;
SS=obj.steady_state_and_balanced_growth_path(1:obj.NumberOfEndogenous(2),:);
Q={obj.Q,[],[]};
if  obj.is_endogenous_switching_model
    Q{2}=obj.shadow_transition_matrix;
    ss=obj.steady_state_and_balanced_growth_path;
    % collect the parameters for all regimes
    M=vertcat(obj.parameters.startval);
    x_ss=0;
    Q{3}={M(:,1),x_ss,ss(:,1),obj.evaluated_definitions(:,1)};
end
H=obj.H; % measurement errors
% locate the deterministic exogenous and remove them right here, right now
if obj.NumberOfObservables(2)
    badguys=[obj.varobs_exo.id];
    R(:,badguys,:,:)=[];
end
[endo_nbr,exo_nbr,horizon,NumberOfRegimes]=size(R);

measurement_errors_flag= ~isempty(H);
if measurement_errors_flag
    varobs_id=[obj.varobs.id];
    CS=nan(size(H));
    for ii=1:NumberOfRegimes
        CS(:,:,ii)=diag(sqrt(diag(H(:,:,ii))));
    end
end
y=zeros(endo_nbr,simul_periods);
y0=zeros(endo_nbr,1);
PAI=1/NumberOfRegimes*ones(NumberOfRegimes,1); % ErgodicDistribution(Q)
State=zeros(simul_periods,1);
good_cols=cell(1,NumberOfRegimes);
for st=1:NumberOfRegimes
    good_cols{st}=any(T(:,:,st));
end

endogenous_switching=~isempty(Q{2});
Q0=Q{1};
if endogenous_switching
    shadow_transition_matrix=Q{2};
    Vargs=Q{3};
end
clear Q

Qt=Q0;
for t=1:simul_burn+simul_periods
    % draw a state
    csp=[0;cumsum(PAI)];
    sj=find(csp>rand,1,'first')-1;
    statecols=good_cols{sj};
    % unconditional forecasts
    y0=T(:,statecols,sj)*y0(statecols);
    % add the shocks
    for hh=1:horizon
        y0=y0+R(:,:,hh,sj)*randn(exo_nbr,1);
    end
    % add the measurement errors
    if measurement_errors_flag
        y0(varobs_id)=y0(varobs_id)+CS(:,:,sj)*randn(obj.NumberOfObservables(1),1);
    end
    if t>simul_burn
        % update probabilities
        PAI=Qt'*PAI;
        y(:,t-simul_burn)=y0+SS(:,sj);
        State(t-simul_burn)=sj;
        if endogenous_switching
            [Qt,retcode]=transition_matrix_evaluation(...
                shadow_transition_matrix,y(:,t-simul_burn),Vargs{:});
            if retcode
                error([mfilename,':: not well behaved endogenous transition matrix'])
            end
        else
            Qt=Q0;
        end
    end
end

for ii=1:endo_nbr
    obj.varendo(ii)=obj.varendo(ii).set_properties('value',y(ii,:));
end

end

