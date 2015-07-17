function [PbParams,M1,M2,RM2i,OMG,DYbar,UU,SS,Yf,MU,LB,UB,R,retcode]=projection_sub_engine(Hypothesis,Nsteps,H,G,Y0,...
    EndogenousConditions,ShocksConditions,verbose) 
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

if nargin<8
    verbose=false;
    if nargin<7
        ShocksConditions=[];
        if nargin<6
            EndogenousConditions=[];
        end
    end
end

[endo_nbr,exo_nbr,NumberOfAnticipatedPeriods]=size(G);

[MUy,OMGy,LBy,UBy,yrest_id,ncv,ncpe]=ExtractConditions(EndogenousConditions);

[MUx,OMGx,LBx,UBx,xrest_id,ncvx,ncpx]=ExtractConditions(ShocksConditions);

ncp=max(ncpe,ncpx);

if NumberOfAnticipatedPeriods<1
    error([mfilename,':: Number of anticipated periods cannot be less than 1'])
end
if NumberOfAnticipatedPeriods==1
    Hypothesis=0;
end

if NumberOfAnticipatedPeriods>ncp+1
    error([mfilename,':: Number of anticipated cannot exceed the number of conditioning periods'])
end

% expected number of iterations
maxiter=max(Nsteps,ncp);

% build restrictions on shocks i.e. genuine endogenous + pure shock restrictions
[R,HHRestr,NumberOfCondShocksPeriods]=utils.forecast.conditional.build_shock_restrictions(...
    H,G,yrest_id,xrest_id,ncp,NumberOfAnticipatedPeriods,Hypothesis); %#ok<ASGLU>

% Unconditional forecasts
Yf=utils.forecast.conditional.initialize_array(@zeros,endo_nbr,maxiter,1,Y0);
Yf=utils.forecast.conditional.compute_forecasts(Yf,H,G,[],maxiter,0);

% Uconditional distribution of endogenous and possibly shocks
% over the forecast horizon:
[DYbar,OMG]=utils.forecast.conditional.unconditional_distribution_of_endogenous_and_shocks(R,Yf(yrest_id,1+(1:ncp)));
EndoCond={MUy,OMGy,LBy,UBy,ncv,ncp};
ShocksCond={MUx,OMGx,LBx,UBx,ncvx,ncpx};
[MU,OMG,LB,UB]=utils.forecast.conditional.conditional_distribution_of_endogenous_and_shocks(DYbar,OMG,EndoCond,ShocksCond);

% Trim everything by removing unconstrained locations
% put covariance matrices in cell to indicate that it is a covariance
% matrix in case OMG is a scalar
[MU,OMG,LB,UB,R,DYbar]=utils.forecast.conditional.remove_holes(MU,{OMG},LB,UB,R,DYbar);

[M1,M2,RM2i,qq,UU,SS,retcode]=utils.forecast.conditional.null_and_column_spaces(R,verbose);

PbParams.maxiter=maxiter;
PbParams.endo_nbr=endo_nbr;
PbParams.exo_nbr=exo_nbr;
PbParams.NumberOfAnticipatedPeriods=NumberOfAnticipatedPeriods;
PbParams.NumberOfCondShocksPeriods=NumberOfCondShocksPeriods;
PbParams.qq=qq;
PbParams.kk=size(R,2);

end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function [MU,OMG,LB,UB,rest_id,ncv,ncp]=ExtractConditions(Conditions)
if isempty(Conditions)
    MU=zeros(0,1);
    OMG=zeros(0);
    LB=zeros(0,1);
    UB=zeros(0,1);
    rest_id=zeros(0,1);
    ncv=0;
    ncp=0;
    return
end

MU		=Conditions{1};
OMG		=Conditions{2};
LB    	=Conditions{3}.LB;
UB    	=Conditions{3}.UB;
rest_id	=Conditions{4};

% % % % % % % trim from the end: This is seems to create more problems
% % % % % % %------------------- should be revisited.
% % % % % % discard=false(1,size(MU,2));
% % % % % % not_nan_mu=any(~isnan(MU),1);
% % % % % % good=find(not_nan_mu,1,'last');
% % % % % % if ~isempty(good)
% % % % % %     discard(good+1:end)=true;
% % % % % % end
% % % % % % MU=MU(:,~discard);
% % % % % % LB=LB(:,~discard);
% % % % % % UB=UB(:,~discard);
if isempty(MU)
    [MU,OMG,LB,UB,rest_id,ncv,ncp]=ExtractConditions([]);
else
    [ncv,ncp]=size(MU);
    if ~isequal(size(LB),[ncv,ncp])||~isequal(size(UB),[ncv,ncp])
        error([mfilename,':: LB and UB must have the same size as the Central Tendency'])
    end
    if ~isequal(ncv,numel(rest_id))
        error([mfilename,':: Number of conditioning variables inconsistent with the bounds'])
    end
    MU=MU(:);
    LB=LB(:);
    UB=UB(:);
end

end
