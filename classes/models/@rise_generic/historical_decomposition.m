function [Histdec,obj]=historical_decomposition(obj,varargin)
% historical_decomposition Computes historical decompositions of a DSGE model
%
% Syntax
% -------
% ::
%
%   [Histdec,obj]=history_dec(obj)
%   [Histdec,obj]=history_dec(obj,varargin)
%
% Inputs
% -------
%
% - obj : [rise|dsge|rfvar|svar] model(s) for which to compute the
%   decomposition. obj could be a vector of models
%
% - varargin : standard optional inputs **coming in pairs**. Among which:
%   - **histdec_start_date** : [char|numeric|{''}] : date at which the
%     decomposition starts. If empty, the decomposition starts at he
%     beginning of the history of the dataset
%
% Outputs
% --------
%
% - Histdec : [struct|cell array] structure or cell array of structures
%   with the decompositions in each model. The decompositions are given in
%   terms of:
%   - the exogenous variables
%   - **InitialConditions** : the effect of initial conditions
%   - **risk** : measure of the effect of non-certainty equivalence
%   - **switch** : the effect of switching (which is also a shock!!!)
%   - **steady_state** : the contribution of the steady state
%
% Remarks
% --------
%
% - the elements that do not contribute to any of the variables are
%   automatically discarded.
%
% - **N.B** : a switching model is inherently nonlinear and so, strictly
%   speaking, the type of decomposition we do for linear/linearized
%   constant-parameter models is not feasible. RISE takes an approximation
%   in which the variables, shocks and states matrices across states are
%   averaged. The averaging weights are the smoothed probabilities.
%
% Examples
% ---------
%
% See also: 

if isempty(obj)
    Histdec=struct('histdec_start_date','');
    return
end

nobj=numel(obj);
if nobj>1
    Histdec=cell(1,nobj);
    for iobj=1:nobj
        [Histdec{iobj},obj(iobj)]=historical_decomposition(obj(iobj),varargin{:});
    end
    return
end

obj=set(obj,varargin{:});

[obj,~,~,retcode]=obj.filter();	% by default, this will also smooth....

% extract the state and transition matrices
if retcode
    error([mfilename,':: model could not be solved'])
end

histdec_start_date=obj.options.histdec_start_date;

endo_names=obj.endogenous.name;
exo_names=obj.exogenous.name;
endo_nbr=obj.endogenous.number;
% number of shocks: this is potentially a problem if there are exogenous
% deterministic variables
exo_nbr = numel(exo_names);
reg_nbr=obj.markov_chains.regimes_number;

% collect the state matrices
iov=obj.inv_order_var;
T=cell(1,reg_nbr);
R=cell(1,reg_nbr);
Tsig=cell(1,reg_nbr);
T0=zeros(endo_nbr);
for ireg=1:reg_nbr
    tmp=obj.solution.Tz{ireg};
    Ty=T0;
    Ty(:,obj.locations.after_solve.t.pb)=tmp(:,obj.locations.after_solve.z.pb);
    T{ireg}=Ty(:,iov);
    R{ireg}=tmp(:,obj.locations.after_solve.z.e_0(1):end);
    Tsig{ireg}=tmp(:,obj.locations.after_solve.z.sig);
end

smoothed_variables=ts.collect(obj.filtering.Expected_smoothed_variables);
smoothed_shocks=ts.collect(obj.filtering.Expected_smoothed_shocks);
smoothed_probabilities=ts.collect(obj.filtering.smoothed_regime_probabilities);
varList=smoothed_variables.varnames;
shockList=smoothed_shocks.varnames;

hist_start_date=smoothed_variables.date_numbers(1);
hist_end_date=smoothed_variables.date_numbers(end);
if isempty(histdec_start_date)
    histdec_start_date=hist_start_date;
else
    histdec_start_date=date2serial(histdec_start_date);
end
if histdec_start_date<hist_start_date || histdec_start_date>hist_end_date
    error([mfilename,':: the decomposition start date must lie between ',...
        serial2date(hist_start_date),' and ',serial2date(hist_end_date)])
end
smoothed_variables=smoothed_variables(hist_start_date:hist_end_date,:);
smoothed_variables=double(smoothed_variables);
smoothed_variables=permute(smoothed_variables,[2,1]);
NumberOfObservations=size(smoothed_variables,2);
smoothed_probabilities=smoothed_probabilities(hist_start_date:hist_end_date,:);
smoothed_probabilities=permute(double(smoothed_probabilities),[2,1]);
% apply the smoothed probabilities to SS
Probs=smoothed_probabilities;

% smoothed shocks
%-----------------
smoothed_shocks=double(smoothed_shocks(hist_start_date:hist_end_date,:,:));
smoothed_shocks=smoothed_shocks(:,:);% expand the pages immediately
% do not transpose shocks!!!

% re-order the smoothed variables
%---------------------------------
varlocs=locate_variables(varList,endo_names);
smoothed_variables=smoothed_variables(varlocs,:);
shocklocs=locate_variables(shockList,exo_names);
nshocks=numel(shockList);
horizon=max(obj.exogenous.shock_horizon(:))+1;
shocklocs=repmat({shocklocs(:).'},1,horizon);
offset=0;
for ih=2:horizon
    offset=offset+exo_nbr;
    shocklocs{ih}=shocklocs{ih}+offset;
end
shocklocs=cell2mat(shocklocs);
smoothed_shocks=smoothed_shocks(:,shocklocs);

% initialize shocks partitions
%------------------------------
proto=(1:nshocks:size(smoothed_shocks,2))-1;

% exogenous--initial condition--perturbation--switch--steady state
%------------------------------------------------------------------
z = zeros(endo_nbr,exo_nbr+4,NumberOfObservations);

init_id=exo_nbr+1;
risk_id=exo_nbr+2;
switch_id=exo_nbr+3;
ss_id=exo_nbr+4;
exo_plus_risk_plus_steady=[(1:exo_nbr),risk_id,ss_id];
exo_plus_init_plus_risk_plus_steady=[(1:exo_nbr),init_id,risk_id,ss_id];
ContributingNames=[exo_names,'InitialConditions','risk','switch','steady_state'];

z=HistoricalDecompositionEngine(z);

if max(max(abs(squeeze(sum(z,2))-smoothed_variables)))>1e-9
    error([mfilename,':: Decomposition failed'])
end

% discard non-contributing names : could be set as an option
%------------------------------------------------------------
good = any(squeeze(any(z,1)),2);

z=z(:,good,:);
discarded=ContributingNames(~good);
ContributingNames=ContributingNames(good);

if ~isempty(discarded)
    disp(discarded)
    disp('The above variables do not contribute to any of the endogenous and so have been discarded')
end

ncontribs=numel(ContributingNames);

Histdec=struct();
for ii=1:endo_nbr
    theData=transpose(squeeze(z(ii,1:end,:)));
    Histdec.(endo_names{ii})=ts(serial2date(histdec_start_date),...
        theData(:,1:ncontribs),ContributingNames(1:ncontribs));
end

    function [z,SS]=HistoricalDecompositionEngine(z)
        [SS,Tsig_t]=smooth_vector(obj.solution.ss,Tsig);
        z(:,ss_id,:)=SS;
        
        for t=1:NumberOfObservations
            [A,B]=expected_state_matrices(t);
            if t > 1
                % evolve shocks
                %---------------
                z(:,1:exo_nbr,t) = A*z(:,1:exo_nbr,t-1);
                % evolve risk
                %-------------
                z(:,risk_id,t) = A*z(:,risk_id,t-1);
            end
            
            % add new shocks
            %---------------
            z(:,1:exo_nbr,t) = z(:,1:exo_nbr,t) + bsxfun_shocks();
            
            % add new risk
            %---------------
            z(:,risk_id,t) = z(:,risk_id,t)+Tsig_t(:,t);
            
            % initial conditions
            %--------------------
            if t==1
                z(:,init_id,t) = smoothed_variables(:,t) - sum(z(:,exo_plus_risk_plus_steady,t),2);
            else
                % initial condition are expected to fade if A is
                % well-behaved
                %----------------------------------------------------------
                z(:,init_id,t)=A*z(:,init_id,t-1);
                % switch
                %-------
                z(:,switch_id,t)=smoothed_variables(:,t) - sum(z(:,exo_plus_init_plus_risk_plus_steady,t),2);
            end
            if obj.options.debug
                discrep=max(abs(sum(z(:,1:end,t),2)-smoothed_variables(:,t)));
                fprintf(1,'t=%0.0f discrep=%0.15f \n',t,discrep);
            end
        end
        if obj.options.debug
            keyboard
        end
        
        function B_S=bsxfun_shocks()
            tmp_=bsxfun(@times,B,smoothed_shocks(t,:));
            B_S=cell(1,nshocks);
            % partition shocks: each shock will have many columns in case
            % of anticipation
            for ishock=1:nshocks
                B_S{ishock}=sum(tmp_(:,proto+ishock),2);
            end
            B_S=cell2mat(B_S);
        end
        
        function [A,B]=expected_state_matrices(t)
            probs_t=Probs(:,t);
            A=probs_t(1)*T{1};
            B=probs_t(1)*R{1};
            for st=2:reg_nbr
                A=A+probs_t(st)*T{st};
                B=B+probs_t(st)*R{st};
            end
        end
    end

    function varargout=smooth_vector(varargin)
        n=nargin;
        varargout=repmat({zeros(endo_nbr,NumberOfObservations)},1,n);
        for t=1:NumberOfObservations
            probs_t=Probs(:,t);
            for st=1:reg_nbr
                for item=1:n
                    varargout{item}(:,t)=varargout{item}(:,t)+probs_t(st)*varargin{item}{st};
                end
            end
        end
    end
end


%{
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 


% PURPOSE: Computes historical decompositions of a DSGE model
%-------------------------------------------------------------
% USAGE:
% where:
%
%
%
%
%-------------------------------------------------------------
% EXAMPLE:
%-------------------------------------------------------------
% LOG:
% 1.
%
%
%-------------------------------------------------------------
% SEE ALSO:
%-------------------------------------------------------------

% written [December 28, 2010] by:
% Junior Maih, Dept of Economics
% The Central Bank of Norway
% junior.maih@norges-bank.no
% this update [July 18, 2014]

if isempty(obj)
    Histdec=struct('histdec_start_date','',...
        'histdec_add_constant',false,...
        'histdec_method',0);
    return
end


nobj=numel(obj);
if nobj>1
    Histdec=cell(1,nobj);
    for iobj=1:nobj
        [Histdec{iobj},obj(iobj)]=historical_decomposition(obj(iobj),varargin{:});
    end
    return
end

if ~isfield(obj.options,'histdec_add_constant')
    obj.options.histdec_add_constant=true;
end

obj=set(obj,varargin{:});

[obj,~,~,retcode]=obj.filter();	% by default, this will also smooth....

% extract the state and transition matrices
if retcode
    error([mfilename,':: model could not be solved'])
end

histdec_start_date=obj.options.histdec_start_date;
histdec_method  =obj.options.histdec_method  ;

endo_names=obj.endogenous.name;
exo_names=obj.exogenous.name;
endo_nbr=obj.endogenous.number;
% number of shocks: this is potentially a problem if there are exogenous
% deterministic variables
exo_nbr = numel(exo_names);
reg_nbr=obj.markov_chains.regimes_number;

% collect the state matrices
[T,R]=set_solution_to_companion(obj);
SS=obj.solution.ss;

smoothed_variables=ts.collect(obj.filtering.Expected_smoothed_variables);
smoothed_shocks=ts.collect(obj.filtering.Expected_smoothed_shocks);
smoothed_probabilities=ts.collect(obj.filtering.smoothed_regime_probabilities);
varList=smoothed_variables.varnames;
shockList=smoothed_shocks.varnames;

hist_start_date=smoothed_variables.date_numbers(1);
hist_end_date=smoothed_variables.date_numbers(end);
if isempty(histdec_start_date)
    histdec_start_date=hist_start_date;
else
    histdec_start_date=date2serial(histdec_start_date);
end
if histdec_start_date<hist_start_date || histdec_start_date>hist_end_date
    error([mfilename,':: the decomposition start date must lie between ',...
        serial2date(hist_start_date),' and ',serial2date(hist_end_date)])
end
smoothed_variables=smoothed_variables(hist_start_date:hist_end_date,:);
smoothed_variables=double(smoothed_variables);
smoothed_variables=permute(smoothed_variables,[2,1]);
NumberOfObservations=size(smoothed_variables,2);
smoothed_probabilities=smoothed_probabilities(hist_start_date:hist_end_date,:);
smoothed_probabilities=permute(double(smoothed_probabilities),[2,1]);
% apply the smoothed probabilities to SS
Probs=smoothed_probabilities;

smoothed_shocks=smoothed_shocks(hist_start_date:hist_end_date,:,:);
smoothed_shocks=permute(double(smoothed_shocks),[2,1,3]);

% re-order the smoothed variables
varlocs=locate_variables(varList,endo_names);
smoothed_variables=smoothed_variables(varlocs,:);
varlocs=locate_variables(shockList,exo_names);
smoothed_shocks=smoothed_shocks(varlocs,:);

% exogenous--initial condition--steady state
%-------------------------------------------
z = zeros(endo_nbr,exo_nbr+2,NumberOfObservations);

z=HistoricalDecompositionEngine(z,smoothed_shocks);

if max(max(abs(squeeze(sum(z,2))-smoothed_variables)))>1e-9
    error([mfilename,':: Decomposition failed'])
end

ContributingNames=[exo_names,'InitialConditions','steady_state'];

ncontribs=numel(ContributingNames)-(1-obj.options.histdec_add_constant);

Histdec=struct();
for ii=1:endo_nbr
    theData=transpose(squeeze(z(ii,1:end,:)));
    Histdec.(endo_names{ii})=ts(serial2date(histdec_start_date),...
        theData(:,1:ncontribs),ContributingNames(1:ncontribs));
end

    function z=HistoricalDecompositionEngine(z,epsilon)
        exo_plus_steady=[(1:exo_nbr),exo_nbr+2];
        NumberOfAnticipatedSteps=size(R{1},3);
        
        for t=1:NumberOfObservations
            [A,B,SS_t]=expected_state_matrices(t);
            z(:,end,t)=SS_t;
            if histdec_method
                if t==1
                    % Collect the shocks
                    for j=1:NumberOfAnticipatedSteps
                        z(:,1:exo_nbr,t) = z(:,1:exo_nbr,t) + ...
                            B(:,:,j).*repmat(epsilon(:,t,j)',endo_nbr,1);
                    end
                    % Find initial condition contribution
                    z(:,exo_nbr+1,t) = smoothed_variables(:,t)-SS_t - ...
                        sum(z(:,exo_plus_steady,t),2);
                else
                    % evolve initial condition
                    z(:,exo_nbr+1,t)=A*z(:,exo_nbr+1,t-1);
                    % collect the shocks
                    z(:,1:exo_nbr,t)=A*z(:,1:exo_nbr,t-1);
                    for j=1:NumberOfAnticipatedSteps
                        z(:,1:exo_nbr,t) = z(:,1:exo_nbr,t) + ...
                            B(:,:,j).*repmat(epsilon(:,t,j)',endo_nbr,1);
                    end
                end
            else
                if t > 1
                    z(:,1:exo_nbr,t) = A*z(:,1:exo_nbr,t-1);
                end
                
                % collect the shocks
                for j=1:NumberOfAnticipatedSteps
                    z(:,1:exo_nbr,t) = z(:,1:exo_nbr,t) + ...
                        B(:,:,j).*repmat(epsilon(:,t,j)',endo_nbr,1);
                end
                % initial conditions
                z(:,exo_nbr+1,t) = smoothed_variables(:,t) - sum(z(:,exo_plus_steady,t),2);
            end
            if obj.options.debug
                discrep=max(abs(sum(z(:,1:end,t),2)-smoothed_variables(:,t)));
                fprintf(1,'t=%0.0f discrep=%0.15f \n',t,discrep);
            end
        end
        if obj.options.debug
            keyboard
        end
        
        function [A,B,SS_t]=expected_state_matrices(t)
            probs_t=Probs(:,t);
            A=probs_t(1)*T{1};
            B=probs_t(1)*R{1};
            SS_t=probs_t(1)*SS{1};
            for st=2:reg_nbr
                A=A+probs_t(st)*T{st};
                B=B+probs_t(st)*R{st};
                SS_t=SS_t+probs_t(st)*SS{st};
            end
        end
    end

end
%}
