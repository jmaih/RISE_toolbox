function obj=historical_decomposition(obj,param_type,dec_start_date,method)

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
% this update [Mai 03, 2011]
if isempty(obj)
    obj=struct();
    return
end

if nargin<4
    method=0;
    if nargin<3
        dec_start_date=[];
        if nargin<2
            param_type=[];
            if nargin <1
                error([mfilename,':: At least a model object should be provided'])
            end
        end
    end
end

mode_=vertcat(obj.estimated_parameters.mode);
mean_=vertcat(obj.estimated_parameters.mean);
startval=vertcat(obj.estimated_parameters.startval);
params=[];
if ~isempty(param_type)
    if strcmp(param_type,'mode')
        if isempty(mode_)
            warning([mfilename,':: mode has not been estimated']) %#ok<WNTAG>
        end
        params=mode_;
    elseif strcmp(param_type,'mean')
        if isempty(mean_)
            if isempty(mode_)
                warning([mfilename,':: model has not been estimated, using calibration']) %#ok<WNTAG>
            else
                warning([mfilename,':: mean has not been estimated, using the mode']) %#ok<WNTAG>
            end
            params=mode_;
        else
            params=mean_;
        end
    elseif strcmp(param_type,'calibration')
        params=startval;
    else
        error([mfilename,':: unrecognized parameter type ',param_type])
    end
end

if obj.NumberOfRegimes>1
    error([mfilename,':: Historical decomposition for multiple regime models not yet implemented'])
end
[obj,junk,junk,retcode]=obj.filter('evaluate_params',params);	% by default, this will also smooth....

% extract the state and transition matrices
if retcode
    error([mfilename,':: model could not be solved'])
end
% collect the state matrices
T=obj.T;
R=obj.R;
SS=obj.steady_state_and_balanced_growth_path(1:obj.NumberOfEndogenous(2),:);

hist_start_date=obj.Filters.smoothed_variables.regime_1.TimeInfo(1);
hist_end_date=obj.Filters.smoothed_variables.regime_1.TimeInfo(end);
if isempty(dec_start_date)
    dec_start_date=hist_start_date;
elseif ischar(dec_start_date)
    dec_start_date=rise_date(dec_start_date);
end
if dec_start_date<hist_start_date || dec_start_date>hist_end_date
    error([mfilename,':: the decomposition start date must lie between ',hist_start_date.date,' and ',hist_end_date.date])
end
% hist_start=1;
% hist_finish=obj.options.data.NumberOfObservations;
% dec_start=Date2Observation(hist_start_date,dec_start_date);

% number of variables
endo_nbr = obj.NumberOfEndogenous(2);

% number of shocks
exo_nbr = obj.NumberOfExogenous;

smoothed_vars=bsxfun(@minus,...
    cell2mat(obj.Filters.smoothed_variables.regime_1.data(2:end,2:end)),...
    transpose(SS(:,1)));
smoothed_shocks=cell2mat(obj.Filters.smoothed_shocks.regime_1.data(2:end,2:end,:));

% keep only the relevant part of the sample
dec_start=hist_start_date.date_2_observation(dec_start_date);
smoothed_vars=smoothed_vars(:,dec_start:end);
smoothed_shocks=smoothed_shocks(:,dec_start:end,:);
% reshape both into var x time x ...
smoothed_vars=transpose(smoothed_vars);
smoothed_shocks=permute(smoothed_shocks,[2,1,3]);
NumberOfObservations=size(smoothed_shocks,2);

z = zeros(endo_nbr,exo_nbr+2,NumberOfObservations);
z(:,end,:)=smoothed_vars;

z=HistoricalDecompositionEngine(T,R,z,smoothed_shocks,method);

if max(max(abs(squeeze(sum(z(:,1:end-1,:),2))-squeeze(z(:,end,:)))))>1e-9
    error([mfilename,':: Decomposition failed'])
end

ContributingNames=char(char(obj.varexo.name),'InitialConditions');
Histdec=rise_time_series(dec_start_date,permute(z(:,1:end-1,:),[3,1,2]),...
    char(obj.varendo.name),ContributingNames);
obj=obj.set_properties('hist_dec',Histdec);

% for j=1:endo_nbr
%     vj=obj.varendo(j).name;
%     Histdec.(vj)=rise_time_series(dec_start_date,squeeze(zz(:,1:end-1,j)),ContributingNames);
% end

function z=HistoricalDecompositionEngine(A,B,z,epsilon,method)

[endo_nbr,exo_nbrPlus2,NumberOfObservations]=size(z);
exo_nbr=exo_nbrPlus2-2;

NumberOfAnticipatedSteps=size(B,3);

for t=1:NumberOfObservations
    if method
        if t==1
            % Collect the shocks
            for j=1:NumberOfAnticipatedSteps
                z(:,1:exo_nbr,t) = z(:,1:exo_nbr,t) + B(:,:,j).*repmat(epsilon(:,t,j)',endo_nbr,1);
            end
            % Find initial condition contribution
            z(:,exo_nbr+1,t) = z(:,exo_nbr+2,t) - sum(z(:,1:exo_nbr,t),2);
        else
            % evolve initial condition
            z(:,exo_nbr+1,t)=A*z(:,exo_nbr+1,t-1);
            % collect the shocks
            z(:,1:exo_nbr,t)=A*z(:,1:exo_nbr,t-1);
            for j=1:NumberOfAnticipatedSteps
                z(:,1:exo_nbr,t) = z(:,1:exo_nbr,t) + B(:,:,j).*repmat(epsilon(:,t,j)',endo_nbr,1);
            end
        end
    else
        if t > 1
            z(:,1:exo_nbr,t) = A*z(:,1:exo_nbr,t-1);
        end
        
        % collect the shocks
        for j=1:NumberOfAnticipatedSteps
            z(:,1:exo_nbr,t) = z(:,1:exo_nbr,t) + B(:,:,j).*repmat(epsilon(:,t,j)',endo_nbr,1);
        end
        % initial conditions
        z(:,exo_nbr+1,t) = z(:,exo_nbr+2,t) - sum(z(:,1:exo_nbr,t),2);
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
