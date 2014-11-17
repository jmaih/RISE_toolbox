function [Yf,CYfMean,Emean,CYf,E]=forecast_engine(Y0,H,G,EndogenousConditions,...
    ShocksConditions,Nsteps,NumberOfSimulations,Hypothesis,method)
% Description:
% Usage:
% Arguments:
% Value:
% Author(s): Junior Maih (junior.maih@norges-bank.no)
% See Also:
% Examples:

if nargin<9
    method='ghk';
    if nargin <8
        Hypothesis=[];
    end
end
if isempty(Hypothesis)
    Hypothesis=0;
end

[PbParams,M1,M2,RM2i,OMG,DYbar,UU,SS,Yf,MU,LB,UB]=...
    utils.forecast.conditional.projection_sub_engine(Hypothesis,Nsteps,H,G,Y0,...
    EndogenousConditions,ShocksConditions);
Yf=format_forecast_output(Yf);
if nargout>1
    [rstar,rlow,rhigh,OMG]=ConditionalDistributionOfShocks(MU,OMG,LB,UB,DYbar);
    
    % check whether there are more holes
    if any(isnan(rlow))||any(isnan(rhigh))
        error([mfilename,':: LB and UB cannot have holes at other places than the holes in the Central tendency'])
    end
    if any(rlow>rhigh)
        error([mfilename,':: some elements in LB greater than their counterparts in UB'])
    end
    
    CSOMG=transpose(chol(OMG));
    % Compute the theoretical mean (doesn't change with the tightening in a symmetric
    % distribution)
    rmean=TruncatedMultivariateNormalMean(rstar,rlow,rhigh,CSOMG);
    gam2mean=RM2i*rmean;
    epsilon=M2*gam2mean;
    
    ShockHorizon=PbParams.NumberOfAnticipatedPeriods-1+PbParams.maxiter;
    Emean=utils.forecast.conditional.initialize_array(@zeros,PbParams.exo_nbr,ShockHorizon);
    Emean(:,1:PbParams.NumberOfCondShocksPeriods)=reshape(epsilon,PbParams.exo_nbr,[]);
    % Compute conditional forecasts
    CYfMean=utils.forecast.conditional.initialize_array(@zeros,PbParams.endo_nbr,PbParams.maxiter,1,Y0);
    CYfMean=utils.forecast.conditional.compute_forecasts(CYfMean,H,G,Emean,Nsteps,PbParams.NumberOfAnticipatedPeriods);
    
    % keep only the relevant part
    %-----------------------------
    CYfMean=format_forecast_output(CYfMean);
    if nargout>2
        ShockHorizonFinal=PbParams.NumberOfAnticipatedPeriods-1+Nsteps;
        Emean=format_shock_output(Emean);
        if nargout>3
            seed=[];
            r=utils.forecast.conditional.truncated_mv_normal_rnd(rstar,OMG,rlow,rhigh,NumberOfSimulations,method,seed);
            % draw gam2
            gam2=RM2i*r;
            clear r
            % draw gam1
            gam1=randn(PbParams.kk-PbParams.qq,NumberOfSimulations);
            % form epsilon
            epsilon=M1*gam1+M2*gam2;
            clear gam1 gam2
            
            % Build E
            E=utils.forecast.conditional.initialize_array(@randn,PbParams.exo_nbr,ShockHorizon,NumberOfSimulations);
            if ~isequal(Hypothesis,0)
                E=0*E;
            end
            E(:,1:PbParams.NumberOfCondShocksPeriods,:)=reshape(epsilon,[PbParams.exo_nbr,PbParams.NumberOfCondShocksPeriods,NumberOfSimulations]);
            clear epsilon
            
            CYf=utils.forecast.conditional.initialize_array(@zeros,PbParams.endo_nbr,PbParams.maxiter,NumberOfSimulations,Y0);
            for isim=1:NumberOfSimulations
                % Compute conditional forecasts
                CYf(:,:,isim)=utils.forecast.conditional.compute_forecasts(CYf(:,:,isim),H,G,E(:,:,isim),Nsteps,PbParams.NumberOfAnticipatedPeriods);
            end
            % Reset the forecasts to their normal length
            CYf=format_forecast_output(CYf);
            if nargout>4
                % and also the shocks (retain only the ones that matter)
                E=format_shock_output(E);
                
                % Decide whether to perform a compatibility test
                if isequal(rlow,rhigh)&& ~isempty(rlow)
                    GuerreroPena(SS,rstar,UU,RM2i,Hypothesis)
                end
            end
        end
    end
end

    function E=format_shock_output(E)
        E=E(:,1:ShockHorizonFinal,:);
    end

    function FK=format_forecast_output(FK)
        FK=FK(:,1:Nsteps+1,:);
    end

end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function GuerreroPena(SS,rstar,UU,RM2i,Hypothesis)
N=200;
qq=size(RM2i,1);
Di=diag(1./diag(SS)); % <--- D=SS(:,1:qq); Di=diag(1./diag(D));
GPstat=rstar'*(UU*Di*RM2i)*rstar;
xx=linspace(0,1,N)';
CritValDistr=chi2inv(1-xx,qq);
MaxDist=find(CritValDistr==inf,1,'last')+1;
GP=linspace(0,CritValDistr(MaxDist),N)';
GP=GP(2:end);
Dist=chi2pdf(GP,qq);
CumDist=chi2cdf(GP,qq);
CritVals=chi2inv(1-[.05,.01],qq);

figure('name',['Chi square Compatibility test for hard conditions (',upper(Hypothesis),' hypothesis)'])
subplot(2,2,1)
plot(CritValDistr,xx,[GPstat,GPstat],[0,1],[CritVals(1),CritVals(1)],[0,1],[CritVals(2),CritVals(2)],[0,1],'linewidth',2)
title(['\chi_{(',int2str(qq),')}^{2} inverse cdf'])
legend('distribution','Statistic','Critical 5%','Critical 1%')%,'Location','SouthOutside'
subplot(2,2,2)
plot(GP,Dist,[GPstat,GPstat],[0,max(Dist)],[CritVals(1),CritVals(1)],[0,max(Dist)],[CritVals(2),CritVals(2)],[0,max(Dist)],'linewidth',2)
title(['\chi_{(',int2str(qq),')}^{2} pdf'])
subplot(2,2,3),plot(GP,CumDist,[GPstat,GPstat],[0,1],[CritVals(1),CritVals(1)],[0,1],[CritVals(2),CritVals(2)],[0,1],'linewidth',2)
title(['\chi_{(',int2str(qq),')}^{2} cdf'])
Pvalue=1-chi2cdf(GPstat,qq);
disp(['p-value for compatibility test under hard conditions:: ',num2str(Pvalue)])
end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function [rstar,rlow,rhigh,OMG]=ConditionalDistributionOfShocks(MU,OMG,LB,UB,DTy)
rstar=MU-DTy;
rlow=LB-DTy;
rhigh=UB-DTy;

end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function xx=TruncatedMultivariateNormalMean(m,lb,ub,CholCov)
npar=size(CholCov,1);
xx=nan(npar,1);
BOUNDS=[lb-m,ub-m];
for i=1:npar
    tmp=(BOUNDS(i,:)-CholCov(i,1:i-1)*xx(1:i-1,1))/CholCov(i,i);
    xx(i)=TruncatedNormalMean(tmp(1),tmp(2));
end
xx=m+CholCov*xx;
end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function y=TruncatedNormalMean(lb,ub)
PHIl = NormalCumulativeDistribution(lb);
PHIr = NormalCumulativeDistribution(ub);
crit=1e-8;

finite_lb_and_ub=isfinite(lb) & isfinite(ub);
same=finite_lb_and_ub &  abs(ub-lb)<crit;
tails=abs(PHIr-PHIl)<crit & ~ same;
good_tails=tails & finite_lb_and_ub;
bad_tails=tails & ~finite_lb_and_ub;
others=~same & ~tails; clear tails

if same % same, no problems
    y=ub;
elseif good_tails% assume a uniform distribution in the tails
    y=.5*(ub+lb);
elseif others% normal distribution for nice ones
    y=-(NormalDensity(ub)-NormalDensity(lb))/(PHIr-PHIl);
elseif bad_tails% Nasty ones
    y=0;
end
% if abs(lb-ub)<1e-9
%     y=ub;
% else
%     y=-(NormalDensity(ub)-NormalDensity(lb))/...
%         (NormalCumulativeDistribution(ub)-NormalCumulativeDistribution(lb));
% end
end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function cdf=NormalCumulativeDistribution(x)
cdf=.5*(1+erf(x/sqrt(2)));
end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function pdf=NormalDensity(x)
pdf=1/sqrt(2*pi)*exp(-.5*x^2);
end