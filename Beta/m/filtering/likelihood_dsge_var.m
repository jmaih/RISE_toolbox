function [LogLik,Incr,retcode,obj]=likelihood_dsge_var(params,obj)
% N.B: This function assumes the likelihood is to be maximized

if nargin~=2
    error([mfilename,':: Number of arguments must be 2'])
end

if obj.NumberOfRegimes>1
    error([mfilename,':: dsge-var for markov switching not implemented yet'])
end
if obj.options.order>1
    error([mfilename,':: dsge-var with anticipations not implemented yet'])
end
%% this important output is not created yet
Incr=[];

%% solve the dsge model
[obj,retcode]=solve(obj,'evaluate_params',params);

%% initialize those and return if there is a problem
% note we take the negative of the penalty in order to maximize
LogLik=-obj.options.Penalty;
if ~retcode
    %% load the elements computed in load_data, using the Schorfheide notation
    %  the fields are: 'YY','YX','XX','XY','T','n','p','constant'
    dsge_var=obj.dsge_var;
    p=dsge_var.p;% the var order
    %% theoretical autocovariances
    [A,info]=theoretical_autocovariances(obj,'ar',p);
    if ~info
        %% load further elements computed in load_data
        T=dsge_var.T;% the sample size
        const=dsge_var.constant; % flag for constant
        n=dsge_var.n; % number of variables
        
        %% VAR approximation of the DSGE model
        ids=[obj.varobs.id];
        steady_state=vertcat(obj.varendo(ids).det_steady_state);
        [PHIstar,SIGstar,GXX,GYX,GXY,GYY,retcode]=var_approximation_to_the_dsge(steady_state,A(ids,ids,:),const);
        if ~retcode
            k=const+n*p;
            %% the prior weight is given by the dsge model
%             lambda=obj.parameters(obj.dsge_prior_weight_id).startval;
            lambda=obj.parameters_image{2,end}(obj.dsge_prior_weight_id);
            if lambda*T>k+n
                %% load the empirical moment matrices
                YY=dsge_var.YY;
                YX=dsge_var.YX;
                XX=dsge_var.XX;
                XY=dsge_var.XY;
                
                %% the resulting Bayesian VAR combines prior(dsge) and the data
                % through the empirical moments
                [PHIb,SIGb,ltgxx,ltgxxi]=bvar_dsge_mode(PHIstar,SIGstar,GXX,XX,GXY,XY,GYX,YX,GYY,YY,T,lambda);
                
                % compute likelihood
                LogLik=bvar_dsge_likelihood(PHIb,SIGb,SIGstar,YY,YX,XY,XX,ltgxx,GXX,lambda,n,T,k);
                %	            dsge_var.GXX=GXX;
                %             dsge_var.GYY=GYY;
                %             dsge_var.GYX=GYX;
                dsge_var.var_approx.PHI=PHIstar;
                dsge_var.var_approx.SIG=SIGstar;
                dsge_var.posterior.PHI=PHIb;
                dsge_var.posterior.SIG=SIGb;
                dsge_var.posterior.ZZi=ltgxxi;
                dsge_var.posterior.inverse_wishart.df=[fix((1+lambda)*T-k),n];
                obj=obj.set_properties('dsge_var',dsge_var);
            end
            %% evaluate the posterior
            if obj.options.debug
                disp(LogLik)
            end
            if obj.options.kf_filtering_level
                % now we filter the data. Hopefully, the parameters
                % estimated using the dsge-var do not have a low density as
                % seen from the point of view of the pure dsge.
                obj=filter(obj);
            end % kf_filtering_level
        end
    end % if info
end % if retcode

% lT=lambda*smpl;
% V=lT*GXX+Mxx;
% Vi=V\eye(k);
% tmp=lT*GYX+Myx;
% PHIpost=Vi*tmp';
% SIGpost=1/((1+lambda)*smpl)*(lT*GYY+Myy-tmp*Vi*tmp');
% function []=var_prior()
%
%
% function []=var_posterior()

function [PHIb,SIGb,ltgxx,ltgxxi]=bvar_dsge_mode(PHI_theta,SIG_theta,GXX,XX,GXY,XY,GYX,YX,GYY,YY,T,lambda)
if isinf(lambda)
    SIGb = SIG_theta;
    PHIb = PHI_theta;
    ltgxx=[];
    ltgxxi=[];
else
    PHIb=(lambda/(1+lambda)*GXX+1/(1+lambda)*XX/T)\...
        (lambda/(1+lambda)*GXY+1/(1+lambda)*XY/T);
    
    ltgyx=lambda*T*GYX+YX;
    ltgxx=lambda*T*GXX+XX;
    ltgxxi=ltgxx\eye(size(ltgxx));
    SIGb=lambda*T*GYY+YY-ltgyx*ltgxxi*ltgyx'; % <---  SIGb=lambda*T*GYY+YY-ltgyx*(ltgxx\ltgyx');
    SIGb=SIGb/((1+lambda)*T);
end


function [PHI,SIG,GXX,GYX,GXY,GYY,rcode]=var_approximation_to_the_dsge(varobs_steady,varobs_autocov,const)
rcode=0;
PHI=[];SIG=[];
GYY=varobs_autocov(:,:,1);
const=any(varobs_steady~=0)||const;
[~,n,p_plus_one]=size(varobs_autocov);
p=p_plus_one-1;
k=const+n*p;
GXX=nan(k);
GYX=nan(n,k);
if const
    GXX(:,1)=[1;repmat(varobs_steady,p,1)];
    GXX(1,2:end)=transpose(GXX(2:end,1));
    GYX(:,1)=varobs_steady;
end

for ii=1:p
    rr=const+((ii-1)*n+1:ii*n);
    for jj=ii:p
        cc=const+((jj-1)*n+1:jj*n);
        if ii==jj
            GXX(rr,cc)=GYY;
        else
            GXX(rr,cc)=varobs_autocov(:,:,jj-ii+1);
        end
    end
    GYX(:,rr)=varobs_autocov(:,:,ii+1);
end
GXX=triu(GXX);
GXX=GXX+triu(GXX,1)';
GXY=transpose(GYX);
Gxxi=GXX\eye(k);
if any(any(isnan(Gxxi)))
    rcode=26;
else
    PHI=Gxxi*GXY;
    SIG=GYY-GYX*Gxxi*GXY;
end



function LogLik=bvar_dsge_likelihood(PHIb,SIGb,SIGstar,YY,YX,XY,XX,ltgxx,GXX,lambda,n,T,k)
if isinf(lambda)
    LogLik = -.5*T*log(det(SIGb))...
        -.5*n*T*log(2*pi)...
        -.5*trace(SIGb\(YY-YX*PHIb-PHIb'*XY+PHIb'*XX*PHIb));
else
    LogLik=-.5*n*log(det(ltgxx))...
        -.5*((1+lambda)*T-k)*log(det((1+lambda)*T*SIGb))...
        +.5*n*log(det(lambda*T*GXX))...
        +.5*(lambda*T-k)*log(det(lambda*T*SIGstar))...
        -.5*n*T*log(2*pi)...
        +.5*n*T*log(2)...
        +sum(gammaln(.5*((1+lambda)*T-k+1-(1:n))))...
        -sum(gammaln(.5*(lambda*T-k+1-(1:n))));
end

