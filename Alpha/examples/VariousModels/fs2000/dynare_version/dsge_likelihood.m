function [fval,DLIK,Hess,exit_flag,ys,trend_coeff,info,Model,DynareOptions,...
    BayesInfo,DynareResults] = dsge_likelihood(xparam1,DynareDataset,...
    DynareOptions,Model,EstimatedParameters,BayesInfo,DynareResults,...
    derivatives_info)
DEBUG=false;

% Evaluates the posterior kernel of a dsge model.

%@info:
%! @deftypefn {Function File} {[@var{fval},@var{exit_flag},@var{ys},@var{trend_coeff},@var{info},@var{Model},@var{DynareOptions},@var{BayesInfo},@var{DynareResults},@var{DLIK},@var{AHess}] =} dsge_likelihood (@var{xparam1},@var{DynareDataset},@var{DynareOptions},@var{Model},@var{EstimatedParameters},@var{BayesInfo},@var{DynareResults},@var{derivatives_flag})
%! @anchor{dsge_likelihood}
%! @sp 1
%! Evaluates the posterior kernel of a dsge model.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item xparam1
%! Vector of doubles, current values for the estimated parameters.
%! @item DynareDataset
%! Matlab's structure describing the dataset (initialized by dynare, see @ref{dataset_}).
%! @item DynareOptions
%! Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%! @item Model
%! Matlab's structure describing the Model (initialized by dynare, see @ref{M_}).
%! @item EstimatedParamemeters
%! Matlab's structure describing the estimated_parameters (initialized by dynare, see @ref{estim_params_}).
%! @item BayesInfo
%! Matlab's structure describing the priors (initialized by dynare, see @ref{bayesopt_}).
%! @item DynareResults
%! Matlab's structure gathering the results (initialized by dynare, see @ref{oo_}).
%! @item derivates_flag
%! Integer scalar, flag for analytical derivatives of the likelihood.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item fval
%! Double scalar, value of (minus) the likelihood.
%! @item exit_flag
%! Integer scalar, equal to zero if the routine return with a penalty (one otherwise).
%! @item ys
%! Vector of doubles, steady state level for the endogenous variables.
%! @item trend_coeffs
%! Matrix of doubles, coefficients of the deterministic trend in the measurement equation.
%! @item info
%! Integer scalar, error code.
%! @table @ @code
%! @item info==0
%! No error.
%! @item info==1
%! The model doesn't determine the current variables uniquely.
%! @item info==2
%! MJDGGES returned an error code.
%! @item info==3
%! Blanchard & Kahn conditions are not satisfied: no stable equilibrium.
%! @item info==4
%! Blanchard & Kahn conditions are not satisfied: indeterminacy.
%! @item info==5
%! Blanchard & Kahn conditions are not satisfied: indeterminacy due to rank failure.
%! @item info==6
%! The jacobian evaluated at the deterministic steady state is complex.
%! @item info==19
%! The steadystate routine thrown an exception (inconsistent deep parameters).
%! @item info==20
%! Cannot find the steady state, info(2) contains the sum of square residuals (of the static equations).
%! @item info==21
%! The steady state is complex, info(2) contains the sum of square of imaginary parts of the steady state.
%! @item info==22
%! The steady has NaNs.
%! @item info==23
%! M_.params has been updated in the steadystate routine and has complex valued scalars.
%! @item info==24
%! M_.params has been updated in the steadystate routine and has some NaNs.
%! @item info==30
%! Ergodic variance can't be computed.
%! @item info==41
%! At least one parameter is violating a lower bound condition.
%! @item info==42
%! At least one parameter is violating an upper bound condition.
%! @item info==43
%! The covariance matrix of the structural innovations is not positive definite.
%! @item info==44
%! The covariance matrix of the measurement errors is not positive definite.
%! @item info==45
%! Likelihood is not a number (NaN).
%! @item info==45
%! Likelihood is a complex valued number.
%! @end table
%! @item Model
%! Matlab's structure describing the model (initialized by dynare, see @ref{M_}).
%! @item DynareOptions
%! Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%! @item BayesInfo
%! Matlab's structure describing the priors (initialized by dynare, see @ref{bayesopt_}).
%! @item DynareResults
%! Matlab's structure gathering the results (initialized by dynare, see @ref{oo_}).
%! @item DLIK
%! Vector of doubles, score of the likelihood.
%! @item AHess
%! Matrix of doubles, asymptotic hessian matrix.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{dynare_estimation_1}, @ref{mode_check}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{dynare_resolve}, @ref{lyapunov_symm}, @ref{schur_statespace_transformation}, @ref{kalman_filter_d}, @ref{missing_observations_kalman_filter_d}, @ref{univariate_kalman_filter_d}, @ref{kalman_steady_state}, @ref{getH}, @ref{kalman_filter}, @ref{score}, @ref{AHessian}, @ref{missing_observations_kalman_filter}, @ref{univariate_kalman_filter}, @ref{priordens}
%! @end deftypefn
%@eod:

% Copyright (C) 2004-2011 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT FR

% Declaration of the penalty as a persistent variable.

% Persistent variable 'penalty' is used to compute an endogenous penalty to
% the value 'fval' when various conditions are encountered. These conditions
% set also 'exit_flag' equal to 0 instead of 1.  It is only when
% dsge_likelihood() is called by an optimizer called by
% dynare_estimation_1() that 'exit_flag' is ignored and penalized 'fval' is
% actually used.
% In that case, 'penalty' is properly initialized, at the very end of the
% present function, by a call to dsge_likelihood() made in
% initial_estimation_checks(). If a condition triggers exit_flag ==
% 0, initial_estimation_checks() triggers an error.
% In summary, an initial call to the present function, without triggering
% any condition, guarantees that 'penalty' is properly initialized when needed.

persistent penalty

% Initialization of the returned variables and others...
fval        = [];
ys          = [];
trend_coeff = [];
exit_flag   = 1;
info        = 0;
singularity_flag = 0;
DLIK        = [];
Hess       = [];

if DynareOptions.estimation_dll
    [fval,exit_flag,ys,trend_coeff,info,params,H,Q] ...
        = logposterior(xparam1,DynareDataset, DynareOptions,Model, ...
        EstimatedParameters,BayesInfo,DynareResults);
    Model.params = params;
    if ~isequal(Model.H,0)
        Model.H = H;
    end
    Model.Sigma_e = Q;
    DynareResults.dr.ys = ys;
    return
end

% Set flag related to analytical derivatives.
analytic_derivation = DynareOptions.analytic_derivation;
if nargout==1,
    analytic_derivation=0;
end

%------------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%------------------------------------------------------------------------------

% Return, with endogenous penalty, if some parameters are smaller than the lower bound of the prior domain.
if ~isequal(DynareOptions.mode_compute,1) && any(xparam1<BayesInfo.lb)
    k = find(xparam1<BayesInfo.lb);
    fval = penalty+sum((BayesInfo.lb(k)-xparam1(k)).^2);
    exit_flag = 0;
    info = 41;
    if analytic_derivation,
        DLIK=ones(length(xparam1),1);
    end
    return
end

% Return, with endogenous penalty, if some parameters are greater than the upper bound of the prior domain.
if ~isequal(DynareOptions.mode_compute,1) && any(xparam1>BayesInfo.ub)
    k = find(xparam1>BayesInfo.ub);
    fval = penalty+sum((xparam1(k)-BayesInfo.ub(k)).^2);
    exit_flag = 0;
    info = 42;
    if analytic_derivation,
        DLIK=ones(length(xparam1),1);
    end
    return
end

% Get the diagonal elements of the covariance matrices for the structural innovations (Q) and the measurement error (H).
Q = Model.Sigma_e;
H = Model.H;
for i=1:EstimatedParameters.nvx
    k =EstimatedParameters.var_exo(i,1);
    Q(k,k) = xparam1(i)*xparam1(i);
end
offset = EstimatedParameters.nvx;
if EstimatedParameters.nvn
    for i=1:EstimatedParameters.nvn
        k = EstimatedParameters.var_endo(i,1);
        H(k,k) = xparam1(i+offset)*xparam1(i+offset);
    end
    offset = offset+EstimatedParameters.nvn;
else
    H = zeros(DynareDataset.info.nvobs);
end

% Get the off-diagonal elements of the covariance matrix for the structural innovations. Test if Q is positive definite.
if EstimatedParameters.ncx
    for i=1:EstimatedParameters.ncx
        k1 =EstimatedParameters.corrx(i,1);
        k2 =EstimatedParameters.corrx(i,2);
        Q(k1,k2) = xparam1(i+offset)*sqrt(Q(k1,k1)*Q(k2,k2));
        Q(k2,k1) = Q(k1,k2);
    end
    % Try to compute the cholesky decomposition of Q (possible iff Q is positive definite)
    [CholQ,testQ] = chol(Q);
    if testQ
        % The variance-covariance matrix of the structural innovations is not definite positive. We have to compute the eigenvalues of this matrix in order to build the endogenous penalty.
        a = diag(eig(Q));
        k = find(a < 0);
        if k > 0
            fval = penalty+sum(-a(k));
            exit_flag = 0;
            info = 43;
            return
        end
    end
    offset = offset+EstimatedParameters.ncx;
end

% Get the off-diagonal elements of the covariance matrix for the measurement errors. Test if H is positive definite.
if EstimatedParameters.ncn
    for i=1:EstimatedParameters.ncn
        k1 = DynareOptions.lgyidx2varobs(EstimatedParameters.corrn(i,1));
        k2 = DynareOptions.lgyidx2varobs(EstimatedParameters.corrn(i,2));
        H(k1,k2) = xparam1(i+offset)*sqrt(H(k1,k1)*H(k2,k2));
        H(k2,k1) = H(k1,k2);
    end
    % Try to compute the cholesky decomposition of H (possible iff H is positive definite)
    [CholH,testH] = chol(H);
    if testH
        % The variance-covariance matrix of the measurement errors is not definite positive. We have to compute the eigenvalues of this matrix in order to build the endogenous penalty.
        a = diag(eig(H));
        k = find(a < 0);
        if k > 0
            fval = penalty+sum(-a(k));
            exit_flag = 0;
            info = 44;
            return
        end
    end
    offset = offset+EstimatedParameters.ncn;
end

% Update estimated structural parameters in Mode.params.
if EstimatedParameters.np > 0
    Model.params(EstimatedParameters.param_vals(:,1)) = xparam1(offset+1:end);
end

% Update Model.Sigma_e and Model.H.
Model.Sigma_e = Q;
Model.H = H;

%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------

% Linearize the model around the deterministic sdteadystate and extract the matrices of the state equation (T and R).
[T,R,SteadyState,info,Model,DynareOptions,DynareResults] = dynare_resolve(Model,DynareOptions,DynareResults,'restrict');

% Return, with endogenous penalty when possible, if dynare_resolve issues an error code (defined in resol).
if info(1) == 1 || info(1) == 2 || info(1) == 5 || info(1) == 7 || info(1) == 22 || info(1) == 24
    fval = penalty+1;
    info = info(1);
    exit_flag = 0;
    if analytic_derivation,
        DLIK=ones(length(xparam1),1);
    end
    return
elseif info(1) == 3 || info(1) == 4 || info(1)==6 ||info(1) == 19 || info(1) == 20 || info(1) == 21  || info(1) == 23
    fval = penalty+info(2);
    info = info(1);
    exit_flag = 0;
    if analytic_derivation,
        DLIK=ones(length(xparam1),1);
    end
    return
end

% Define a vector of indices for the observed variables. Is this really usefull?...
BayesInfo.mf = BayesInfo.mf1;

% Define the constant vector of the measurement equation.
if DynareOptions.noconstant
    constant = zeros(DynareDataset.info.nvobs,1);
else
    if DynareOptions.loglinear
        constant = log(SteadyState(BayesInfo.mfys));
    else
        constant = SteadyState(BayesInfo.mfys);
    end
end

% Define the deterministic linear trend of the measurement equation.
if BayesInfo.with_trend
    trend_coeff = zeros(DynareDataset.info.nvobs,1);
    t = DynareOptions.trend_coeffs;
    for i=1:length(t)
        if ~isempty(t{i})
            trend_coeff(i) = evalin('base',t{i});
        end
    end
    trend = repmat(constant,1,DynareDataset.info.ntobs)+trend_coeff*[1:DynareDataset.info.ntobs];
else
    trend = repmat(constant,1,DynareDataset.info.ntobs);
end

% Get needed informations for kalman filter routines.
start = DynareOptions.presample+1;
Z = BayesInfo.mf; % old mf
no_missing_data_flag = ~DynareDataset.missing.state;
mm = length(T); % old np
pp = DynareDataset.info.nvobs;
rr = length(Q);
kalman_tol = DynareOptions.kalman_tol;
riccati_tol = DynareOptions.riccati_tol;
Y   = DynareDataset.data-trend;

%------------------------------------------------------------------------------
% 3. Initial condition of the Kalman filter
%------------------------------------------------------------------------------
kalman_algo = DynareOptions.kalman_algo;

% resetting measurement errors covariance matrix for univariate filters
if (kalman_algo == 2) || (kalman_algo == 4)
    if isequal(H,0)
        H = zeros(pp,1);
        mmm = mm;
    else
        if all(all(abs(H-diag(diag(H)))<1e-14))% ie, the covariance matrix is diagonal...
            H = diag(H);
            mmm = mm;
        else
            Z = [Z, eye(pp)];
            T = blkdiag(T,zeros(pp));
            Q = blkdiag(Q,H);
            R = blkdiag(R,eye(pp));
            Pstar = blkdiag(Pstar,H);
            Pinf  = blckdiag(Pinf,zeros(pp));
            H = zeros(pp,1);
            mmm   = mm+pp;
        end
    end
end


diffuse_periods = 0;
correlated_errors_have_been_checked = 0;
singular_diffuse_filter = 0;
switch DynareOptions.lik_init
    case 1% Standard initialization with the steady state of the state equation.
        if kalman_algo~=2
            % Use standard kalman filter except if the univariate filter is explicitely choosen.
            kalman_algo = 1;
        end
        if DynareOptions.lyapunov_fp == 1
            Pstar = lyapunov_symm(T,Q,DynareOptions.lyapunov_fixed_point_tol,DynareOptions.lyapunov_complex_threshold, 3, R);
        elseif DynareOptions.lyapunov_db == 1
            Pstar = disclyap_fast(T,R*Q*R',DynareOptions.lyapunov_doubling_tol);
        elseif DynareOptions.lyapunov_srs == 1
            Pstar = lyapunov_symm(T,Q,DynareOptions.lyapunov_fixed_point_tol,DynareOptions.lyapunov_complex_threshold, 4, R);
        else
            Pstar = lyapunov_symm(T,R*Q*R',DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold);
        end;
        Pinf  = [];
        a     = zeros(mm,1);
        Zflag = 0;
    case 2% Initialization with large numbers on the diagonal of the covariance matrix if the states (for non stationary models).
        if kalman_algo ~= 2
            % Use standard kalman filter except if the univariate filter is explicitely choosen.
            kalman_algo = 1;
        end
        Pstar = DynareOptions.Harvey_scale_factor*eye(mm);
        Pinf  = [];
        a     = zeros(mm,1);
        Zflag = 0;
    case 3% Diffuse Kalman filter (Durbin and Koopman)
        % Use standard kalman filter except if the univariate filter is explicitely choosen.
        if kalman_algo == 0
            kalman_algo = 3;
        elseif ~((kalman_algo == 3) || (kalman_algo == 4))
            error(['diffuse filter: options_.kalman_algo can only be equal ' ...
                'to 0 (default), 3 or 4'])
        end
        
        [Z,T,R,QT,Pstar,Pinf] = schur_statespace_transformation(Z,T,R,Q,DynareOptions.qz_criterium);
        Zflag = 1;
        % Run diffuse kalman filter on first periods.
        if (kalman_algo==3)
            % Multivariate Diffuse Kalman Filter
            if no_missing_data_flag
                [dLIK,dlik,a,Pstar] = kalman_filter_d(Y, 1, size(Y,2), ...
                    zeros(mm,1), Pinf, Pstar, ...
                    kalman_tol, riccati_tol, DynareOptions.presample, ...
                    T,R,Q,H,Z,mm,pp,rr);
            else
                [dLIK,dlik,a,Pstar] = missing_observations_kalman_filter_d(DynareDataset.missing.aindex,DynareDataset.missing.number_of_observations,DynareDataset.missing.no_more_missing_observations, ...
                    Y, 1, size(Y,2), ...
                    zeros(mm,1), Pinf, Pstar, ...
                    kalman_tol, riccati_tol, DynareOptions.presample, ...
                    T,R,Q,H,Z,mm,pp,rr);
            end
            diffuse_periods = length(dlik);
            if isinf(dLIK)
                % Go to univariate diffuse filter if singularity problem.
                singular_diffuse_filter = 1;
            end
        end
        if singular_diffuse_filter || (kalman_algo==4)
            % Univariate Diffuse Kalman Filter
            if isequal(H,0)
                H1 = zeros(pp,1);
                mmm = mm;
            else
                if all(all(abs(H-diag(diag(H)))<1e-14))% ie, the covariance matrix is diagonal...
                    H1 = diag(H);
                    mmm = mm;
                else
                    Z = [Z, eye(pp)];
                    T = blkdiag(T,zeros(pp));
                    Q = blkdiag(Q,H);
                    R = blkdiag(R,eye(pp));
                    Pstar = blkdiag(Pstar,H);
                    Pinf  = blckdiag(Pinf,zeros(pp));
                    H1 = zeros(pp,1);
                    mmm   = mm+pp;
                end
            end
            % no need to test again for correlation elements
            correlated_errors_have_been_checked = 1;
            
            [dLIK,dlik,a,Pstar] = univariate_kalman_filter_d(DynareDataset.missing.aindex,...
                DynareDataset.missing.number_of_observations,...
                DynareDataset.missing.no_more_missing_observations, ...
                Y, 1, size(Y,2), ...
                zeros(mmm,1), Pinf, Pstar, ...
                kalman_tol, riccati_tol, DynareOptions.presample, ...
                T,R,Q,H1,Z,mmm,pp,rr);
            diffuse_periods = length(dlik);
        end
    case 4% Start from the solution of the Riccati equation.
        if kalman_algo ~= 2
            kalman_algo = 1;
        end
        if isequal(H,0)
            [err,Pstar] = kalman_steady_state(transpose(T),R*Q*transpose(R),transpose(build_selection_matrix(Z,np,length(Z))));
        else
            [err,Pstar] = kalman_steady_state(transpose(T),R*Q*transpose(R),transpose(build_selection_matrix(Z,np,length(Z))),H);
        end
        if err
            disp(['dsge_likelihood:: I am not able to solve the Riccati equation, so I switch to lik_init=1!']);
            DynareOptions.lik_init = 1;
            Pstar = lyapunov_symm(T,R*Q*R',DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold);
        end
        Pinf  = [];
    otherwise
        error('dsge_likelihood:: Unknown initialization approach for the Kalman filter!')
end

if analytic_derivation
    no_DLIK = 0;
    full_Hess = analytic_derivation==2;
    asy_Hess = analytic_derivation==-2;
    if asy_Hess,
        analytic_derivation=1;
    end
    DLIK = [];
    AHess = [];
    if nargin<8 || isempty(derivatives_info)
        [A,B,nou,nou,Model,DynareOptions,DynareResults] = dynare_resolve(Model,DynareOptions,DynareResults);
        if ~isempty(EstimatedParameters.var_exo)
            indexo=EstimatedParameters.var_exo(:,1);
        else
            indexo=[];
        end
        if ~isempty(EstimatedParameters.param_vals)
            indparam=EstimatedParameters.param_vals(:,1);
        else
            indparam=[];
        end
        
        if full_Hess,
            [dum, DT, DOm, DYss, dum2, D2T, D2Om, D2Yss] = getH(A, B, Model,DynareResults,0,indparam,indexo);
        else
            [dum, DT, DOm, DYss] = getH(A, B, Model,DynareResults,0,indparam,indexo);
        end
    else
        DT = derivatives_info.DT;
        DOm = derivatives_info.DOm;
        DYss = derivatives_info.DYss;
        if isfield(derivatives_info,'full_Hess'),
            full_Hess = derivatives_info.full_Hess;
        end
        if full_Hess,
            D2T = derivatives_info.D2T;
            D2Om = derivatives_info.D2Om;
            D2Yss = derivatives_info.D2Yss;
        end
        if isfield(derivatives_info,'no_DLIK'),
            no_DLIK = derivatives_info.no_DLIK;
        end
        clear('derivatives_info');
    end
    iv = DynareResults.dr.restrict_var_list;
    DYss = [zeros(size(DYss,1),offset) DYss];
    DT = DT(iv,iv,:);
    DOm = DOm(iv,iv,:);
    DYss = DYss(iv,:);
    DH=zeros([size(H),length(xparam1)]);
    DQ=zeros([size(Q),length(xparam1)]);
    DP=zeros([size(T),length(xparam1)]);
    if full_Hess,
        for j=1:size(D2Yss,1),
            tmp(j,:,:) = blkdiag(zeros(offset,offset), squeeze(D2Yss(j,:,:)));
        end
        D2Yss = tmp;
        D2T = D2T(iv,iv,:,:);
        D2Om = D2Om(iv,iv,:,:);
        D2Yss = D2Yss(iv,:,:);
        D2H=zeros([size(H),length(xparam1),length(xparam1)]);
        D2P=zeros([size(T),length(xparam1),length(xparam1)]);
    end
    for i=1:EstimatedParameters.nvx
        k =EstimatedParameters.var_exo(i,1);
        DQ(k,k,i) = 2*sqrt(Q(k,k));
        dum =  lyapunov_symm(T,DOm(:,:,i),DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold);
        kk = find(abs(dum) < 1e-12);
        dum(kk) = 0;
        DP(:,:,i)=dum;
        if full_Hess
            for j=1:i,
                dum =  lyapunov_symm(T,D2Om(:,:,i,j),DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold);
                kk = (abs(dum) < 1e-12);
                dum(kk) = 0;
                D2P(:,:,i,j)=dum;
                D2P(:,:,j,i)=dum;
            end
        end
    end
    offset = EstimatedParameters.nvx;
    for i=1:EstimatedParameters.nvn
        k = EstimatedParameters.var_endo(i,1);
        DH(k,k,i+offset) = 2*sqrt(H(k,k));
        if full_Hess
            D2H(k,k,i+offset,i+offset) = 2;
        end
    end
    offset = offset + EstimatedParameters.nvn;
    for j=1:EstimatedParameters.np
        dum =  lyapunov_symm(T,DT(:,:,j+offset)*Pstar*T'+T*Pstar*DT(:,:,j+offset)'+DOm(:,:,j+offset),DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold);
        kk = find(abs(dum) < 1e-12);
        dum(kk) = 0;
        DP(:,:,j+offset)=dum;
        if full_Hess
            DTj = DT(:,:,j+offset);
            DPj = dum;
            for i=1:j,
                DTi = DT(:,:,i+offset);
                DPi = DP(:,:,i+offset);
                D2Tij = D2T(:,:,i,j);
                D2Omij = D2Om(:,:,i,j);
                tmp = D2Tij*Pstar*T' + T*Pstar*D2Tij' + DTi*DPj*T' + DTj*DPi*T' + T*DPj*DTi' + T*DPi*DTj' + DTi*Pstar*DTj' + DTj*Pstar*DTi' + D2Omij;
                dum = lyapunov_symm(T,tmp,DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold);
                dum(abs(dum)<1.e-12) = 0;
                D2P(:,:,i+offset,j+offset) = dum;
                D2P(:,:,j+offset,i+offset) = dum;
            end
        end
    end
    if analytic_derivation==1,
        analytic_deriv_info={analytic_derivation,DT,DYss,DOm,DH,DP};
    else
        analytic_deriv_info={analytic_derivation,DT,DYss,DOm,DH,DP,D2T,D2Yss,D2Om,D2H,D2P};
    end
else
    analytic_deriv_info={0};
end

%------------------------------------------------------------------------------
% 4. Likelihood evaluation
%------------------------------------------------------------------------------

if ((kalman_algo==1) || (kalman_algo==3))% Multivariate Kalman Filter
    if no_missing_data_flag
        if DynareOptions.block == 1
            [err, LIK] = block_kalman_filter(T,R,Q,H,Pstar,Y,start,Z,kalman_tol,riccati_tol, Model.nz_state_var, Model.n_diag, Model.nobs_non_statevar);
            mexErrCheck('block_kalman_filter', err);
        else
            [LIK,lik] = kalman_filter(Y,diffuse_periods+1,size(Y,2), ...
                a,Pstar, ...
                kalman_tol, riccati_tol, ...
                DynareOptions.presample, ...
                T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods, ...
                analytic_deriv_info{:});
            
        end
    else
        [LIK,lik] = missing_observations_kalman_filter(DynareDataset.missing.aindex,DynareDataset.missing.number_of_observations,DynareDataset.missing.no_more_missing_observations,Y,diffuse_periods+1,size(Y,2), ...
            a, Pstar, ...
            kalman_tol, DynareOptions.riccati_tol, ...
            DynareOptions.presample, ...
            T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods);
    end
    if analytic_derivation,
        LIK1=LIK;
        LIK=LIK1{1};
    end
    if isinf(LIK)
        if kalman_algo == 1
            kalman_algo = 2;
        else
            kalman_algo = 4;
        end
    else
        if DynareOptions.lik_init==3
            LIK = LIK + dLIK;
            if analytic_derivation==0 && nargout==2,
                lik = [dlik; lik];
            end
        end
    end
end

if (kalman_algo==2) || (kalman_algo==4)
    % Univariate Kalman Filter
    % resetting measurement error covariance matrix when necessary                                                           %
    if ~correlated_errors_have_been_checked
        if isequal(H,0)
            H = zeros(pp,1);
            mmm = mm;
            if analytic_derivation,
                DH = zeros(pp,length(xparam1));
            end
        else
            if all(all(abs(H-diag(diag(H)))<1e-14))% ie, the covariance matrix is diagonal...
                H = diag(H);
                mmm = mm;
                if analytic_derivation,
                    for j=1:pp,
                        tmp(j,:)=DH(j,j,:);
                    end
                    DH=tmp;
                end
            else
                Z = [Z, eye(pp)];
                T = blkdiag(T,zeros(pp));
                Q = blkdiag(Q,H);
                R = blkdiag(R,eye(pp));
                Pstar = blkdiag(Pstar,H);
                Pinf  = blckdiag(Pinf,zeros(pp));
                H = zeros(pp,1);
                mmm   = mm+pp;
            end
        end
        if analytic_derivation,
            analytic_deriv_info{5}=DH;
        end
    end
    
    [LIK, lik] = univariate_kalman_filter(DynareDataset.missing.aindex,DynareDataset.missing.number_of_observations,DynareDataset.missing.no_more_missing_observations,Y,diffuse_periods+1,size(Y,2), ...
        a,Pstar, ...
        DynareOptions.kalman_tol, ...
        DynareOptions.riccati_tol, ...
        DynareOptions.presample, ...
        T,Q,R,H,Z,mmm,pp,rr,Zflag,diffuse_periods,analytic_deriv_info{:});
    if analytic_derivation,
        LIK1=LIK;
        LIK=LIK1{1};
    end
    if DynareOptions.lik_init==3
        LIK = LIK+dLIK;
        if analytic_derivation==0 && nargout==2,
            lik = [dlik; lik];
        end
    end
end

if analytic_derivation
    if no_DLIK==0
        DLIK = LIK1{2};
        %                 [DLIK] = score(T,R,Q,H,Pstar,Y,DT,DYss,DOm,DH,DP,start,Z,kalman_tol,riccati_tol);
    end
    if full_Hess,
        Hess = -LIK1{3};
        %                     [Hess, DLL] = get_Hessian(T,R,Q,H,Pstar,Y,DT,DYss,DOm,DH,DP,D2T,D2Yss,D2Om,D2H,D2P,start,Z,kalman_tol,riccati_tol);
        %                     Hess0 = getHessian(Y,T,DT,D2T, R*Q*transpose(R),DOm,D2Om,Z,DYss,D2Yss);
    end
    if asy_Hess,
        [Hess] = AHessian(T,R,Q,H,Pstar,Y,DT,DYss,DOm,DH,DP,start,Z,kalman_tol,riccati_tol);
    end
end


if isnan(LIK)
    info = 45;
    exit_flag = 0;
    return
end
if imag(LIK)~=0
    likelihood = penalty;
else
    likelihood = LIK;
end

% ------------------------------------------------------------------------------
% 5. Adds prior if necessary
% ------------------------------------------------------------------------------
if analytic_derivation
    if full_Hess,
        [lnprior, dlnprior, d2lnprior] = priordens(xparam1,BayesInfo.pshape,BayesInfo.p6,BayesInfo.p7,BayesInfo.p3,BayesInfo.p4,false,DEBUG);
        Hess = Hess - d2lnprior;
    else
        [lnprior, dlnprior] = priordens(xparam1,BayesInfo.pshape,BayesInfo.p6,BayesInfo.p7,BayesInfo.p3,BayesInfo.p4,false,DEBUG);
    end
    if no_DLIK==0
        DLIK = DLIK - dlnprior';
    end
else
    lnprior = priordens(xparam1,BayesInfo.pshape,BayesInfo.p6,BayesInfo.p7,BayesInfo.p3,BayesInfo.p4,false,DEBUG);
end
fval    = (likelihood-lnprior);

if DEBUG
    disp([mfilename,':: likelihood ',num2str(likelihood)])
    disp([mfilename,':: lnprior ',num2str(lnprior)])
    disp([mfilename,':: logpost ',num2str(fval)])
    keyboard
end
% Update DynareOptions.kalman_algo.
DynareOptions.kalman_algo = kalman_algo;

% Update the penalty.
penalty = fval;

if analytic_derivation==0 && nargout==2,
    lik=lik(start:end,:);
    DLIK=[-lnprior; lik(:)];
end