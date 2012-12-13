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
% this update [December 09, 2012]

Defaults=struct('histdec_start_date','',...
    'histdec_method',0);

if isempty(obj)
    Histdec=Defaults;
    return
end

if numel(obj)>1
    error([mfilename,':: Historical decomposition for multiple models not permitted'])
end

obj=set_options(obj,varargin{:});

[obj,~,~,retcode]=obj.filter();	% by default, this will also smooth....

% extract the state and transition matrices
if retcode
    error([mfilename,':: model could not be solved'])
end

oldoptions=obj.options;
thisDefault=Defaults;
Fields=fieldnames(thisDefault);
for ifield=1:numel(Fields)
    vi=Fields{ifield};
    if isfield(oldoptions,vi)
        thisDefault.(vi)=oldoptions.(vi);
    end
end
histdec_start_date=thisDefault.histdec_start_date;
histdec_method  =thisDefault.histdec_method  ;

% collect the state matrices
T=obj.T;
R=obj.R;
SS=obj.steady_state_and_balanced_growth_path(1:obj.NumberOfEndogenous(2),:);

smoothed_variables=rise_time_series.collect(obj.Filters.Expected_smoothed_variables);
smoothed_shocks=rise_time_series.collect(obj.Filters.Expected_smoothed_shocks);
smoothed_probabilities=rise_time_series.collect(obj.Filters.smoothed_probabilities);
varList=smoothed_variables.varnames;
shockList=smoothed_shocks.varnames;
probList=smoothed_probabilities.varnames;
hist_start_date=smoothed_variables.TimeInfo(1);
hist_end_date=smoothed_variables.TimeInfo(end);
if isempty(histdec_start_date)
    histdec_start_date=hist_start_date;
else
    histdec_start_date=rise_date(histdec_start_date);
end
if histdec_start_date<hist_start_date || histdec_start_date>hist_end_date
    error([mfilename,':: the decomposition start date must lie between ',...
        hist_start_date.date,' and ',hist_end_date.date])
end
smoothed_variables=window(smoothed_variables,hist_start_date);
smoothed_variables=double(smoothed_variables);
NumberOfObservations=size(smoothed_variables,1);
smoothed_probabilities=window(smoothed_probabilities,hist_start_date);
smoothed_probabilities=permute(double(smoothed_probabilities),[2,1]);
% apply the smoothed probabilities to SS
Regimes=obj.Regimes;
[nreg,nchains]=size(Regimes);
Probs=ones(nreg,NumberOfObservations);
chainNames={obj.markov_chains.name};
for ireg=1:nreg
    reg_i=Regimes(ireg,:);
    for ichain=1:nchains
        ch_reg=[chainNames{ichain},'_',int2str(reg_i(ichain))];
        loc=strcmp(ch_reg,probList);
        Probs(ireg,:)=Probs(ireg,:).*smoothed_probabilities(loc,:);
    end
end

% SS=SS(:,1); % to be adjusted for markov chain later
% smoothed_variables=bsxfun(@minus,smoothed_variables,transpose(SS));
smoothed_variables=permute(smoothed_variables,[2,1]);
smoothed_shocks=window(smoothed_shocks,hist_start_date);
smoothed_shocks=permute(double(smoothed_shocks),[2,1,3]);

endo_names={obj.varendo.name};
exo_names={obj.varexo.name};

% re-order the smoothed variables
varlocs=locate_variables(varList,endo_names);
smoothed_variables=smoothed_variables(varlocs,:);
varlocs=locate_variables(shockList,exo_names);
smoothed_shocks=smoothed_shocks(varlocs,:);

% number of variables
endo_nbr = numel(endo_names);

% number of shocks: this is potentially a problem if there are exogenous
% deterministic variables 
exo_nbr = numel(exo_names);

% orig_smoothed_variables=smoothed_variables;
z = zeros(endo_nbr,exo_nbr+1,NumberOfObservations);

z=HistoricalDecompositionEngine(z,smoothed_shocks);

if max(max(abs(squeeze(sum(z,2))-smoothed_variables)))>1e-9
    error([mfilename,':: Decomposition failed'])
end

ContributingNames=[exo_names,'InitialConditions'];

Histdec=struct();
for ii=1:endo_nbr
    theData=transpose(squeeze(z(ii,1:end,:)));
    if ii==1
        Histdec.(endo_names{ii})=rise_time_series(histdec_start_date,...
            theData,ContributingNames);
        TimeInfo=Histdec.(endo_names{ii}).TimeInfo;
    else
        Histdec.(endo_names{ii})=rise_time_series(TimeInfo,theData,...
            ContributingNames);
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
            A=probs_t(1)*T(:,:,1);
            B=probs_t(1)*R(:,:,:,1);
            SS_t=probs_t(1)*SS(:,1);
            for st=2:nreg
                A=A+probs_t(st)*T(:,:,st);
                B=B+probs_t(st)*R(:,:,:,st);
                SS_t=SS_t+probs_t(st)*SS(:,st);
            end
        end
    end

end
%%%%function Packeddec_db=SelectionAndRepackagingEngine(Histdec,Packages,var_list)
%%%%
%%%%endo_names=fieldnames(Histdec);
%%%%% indices of endogenous variables
%%%%i_var=locate_variables(var_list,char(endo_names));
%%%%out=setdiff(1:numel(endo_names),i_var);
%%%%Packeddec_db=rmfield(Histdec,endo_names(out));
%%%%
%%%%NumberOfPackages=size(Packages,2);
%%%%if NumberOfPackages==0
%%%%    return
%%%%end
%%%%%
%%%%ContributingNames=Packeddec_db.(deblank(var_list(1,:))).varnames;
%%%%startdate=Packeddec_db.(deblank(var_list(1,:))).start;
%%%%
%%%%% Now Repackage things up
%%%%DejaVu=[];
%%%%PackNames='';
%%%%Package_id=cell(1,NumberOfPackages);
%%%%for j=1:NumberOfPackages
%%%%    PackNames=strvcat(PackNames,Packages{1,j}); %#ok<VCAT>
%%%%    Package_id{j} = locate_variables(Packages{2,j},ContributingNames);
%%%%    tmp=intersect(DejaVu,Package_id{j});
%%%%    if ~isempty(tmp)
%%%%        tmp=mat2str(transpose(tmp(:)));
%%%%        error([mfilename,':: shocks ',tmp,' have been declared at least twice'])
%%%%    end
%%%%    DejaVu=union(DejaVu,Package_id{j});
%%%%end
%%%%Rest=setdiff(1:Packeddec_db.(deblank(var_list(1,:))).NumberOfVariables,DejaVu);
%%%%if ~isempty(Rest)
%%%%    PackNames=strvcat(PackNames,'Rest'); %#ok<VCAT>
%%%%    Package_id{NumberOfPackages+1}=Rest;
%%%%    NumberOfPackages=NumberOfPackages+1;
%%%%end
%%%%for j=1:numel(i_var)
%%%%    db=Packeddec_db.(deblank(var_list(j,:)));
%%%%    data=cell2mat(db.data(2:end,2:end));
%%%%    Newdata=nan(db.NumberOfObservations,NumberOfPackages);
%%%%    for p=1:NumberOfPackages
%%%%        Newdata(:,p)=sum(data(:,Package_id{p}),2);
%%%%    end
%%%%    Packeddec_db.(deblank(var_list(j,:)))=rise_time_series(startdate,Newdata,PackNames);
%%%%end
