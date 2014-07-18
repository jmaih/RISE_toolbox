function [Histdec,obj]=historical_decomposition(obj,varargin)

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
endo_nbr=obj.endogenous.number(end);
% number of shocks: this is potentially a problem if there are exogenous
% deterministic variables
exo_nbr = numel(exo_names);
reg_nbr=obj.markov_chains.regimes_number;

% collect the state matrices
[T,R]=inv_order_var_solution();
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

z = zeros(endo_nbr,exo_nbr+1,NumberOfObservations);

z=HistoricalDecompositionEngine(z,smoothed_shocks);

if max(max(abs(squeeze(sum(z,2))-smoothed_variables)))>1e-9
    error([mfilename,':: Decomposition failed'])
end

ContributingNames=[exo_names,'InitialConditions'];

Histdec=struct();
for ii=1:endo_nbr
    theData=transpose(squeeze(z(ii,1:end,:)));
    Histdec.(endo_names{ii})=ts(serial2date(histdec_start_date),...
        theData,ContributingNames);
end

    function [Tz,Re]=inv_order_var_solution()
        %         ov=obj.order_var.after_solve;
        iov=obj.inv_order_var.after_solve;
        z_pb=obj.locations.after_solve.z.pb;
        t_pb=obj.locations.after_solve.t.pb;
        e_0=obj.locations.after_solve.z.e_0;
        Re=cell(1,reg_nbr);
        Tz=Re;
        tmp=zeros(endo_nbr);
        for isol=1:reg_nbr
            tmp(:,t_pb)=obj.solution.Tz{isol}(:,z_pb);
            % separate autoregressive part from shocks
            %----------------------------------------
            Tz{isol}=tmp(iov,iov);
            Re{isol}=obj.solution.Tz{isol}(:,e_0(1):end);
            if isol==1
                npges=size(Re{isol},2)/exo_nbr;
            end
            Re{isol}=reshape(Re{isol},[endo_nbr,exo_nbr,npges]);
        end
    end

    function z=HistoricalDecompositionEngine(z,epsilon)
        
        NumberOfAnticipatedSteps=size(R,3);
        
        for t=1:NumberOfObservations
            [A,B,SS_t]=expected_state_matrices(t);
            smoothed_variables(:,t)=smoothed_variables(:,t)-SS_t;
            if histdec_method
                if t==1
                    % Collect the shocks
                    for j=1:NumberOfAnticipatedSteps
                        z(:,1:exo_nbr,t) = z(:,1:exo_nbr,t) + ...
                            B(:,:,j).*repmat(epsilon(:,t,j)',endo_nbr,1);
                    end
                    % Find initial condition contribution
                    z(:,exo_nbr+1,t) = smoothed_variables(:,t) - ...
                        sum(z(:,1:exo_nbr,t),2);
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
                z(:,exo_nbr+1,t) = smoothed_variables(:,t) - sum(z(:,1:exo_nbr,t),2);
            end
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

