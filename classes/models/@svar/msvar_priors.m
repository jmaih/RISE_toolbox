function obj=msvar_priors(obj,estim_names)
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


if isempty(obj)
    obj=struct('vp_mnst_overall_tightness',3,...1
        'vp_mnst_relative_tightness_lags',1,...
        'vp_mnst_relative_tightness_constant',0.1,...
        'vp_mnst_tightness_on_lag_decay',0.5,...1.2
        'vp_mnst_unit_root_vars','all',... % {'all'},'none',cellarray={'v1','v5','v10'}
        'vp_mnst_stationary_var_mean',0.5,... % unit root = 1
        'vp_natconj_normwish_variance',10,...
        'vp_analytical_post_mode',true,... % compute the posterior mode analytically if possible
        'vp_gls_ar1_processes',true,... % use the covariance formed by the ar1 processes when forming the posterior mode( else use the covariance of the OLS)
        'vp_prior_type','minnesota'); % {minnesota},none,normal_wishart,indep_normal_wishart,jeffrey(diffuse) 
    return
end
s=quick_ar1_processes();
obj.miscellaneous.constant_var.sigma_=diag(s.^2);

prior_type=obj.options.vp_prior_type;
use_priors=~strcmp(prior_type,'none');
if ~use_priors
    % use the minnesota to build the uniform... for now
    prior_type='minnesota';
end

nvar=numel(s);
if strcmp(prior_type,'minnesota')
    W_mnst=eye(nvar);
    W_mnst(W_mnst==0)=obj.options.vp_mnst_relative_tightness_lags;
        theta=obj.options.vp_mnst_overall_tightness;
        phi=obj.options.vp_mnst_tightness_on_lag_decay;
        if ~(0<=phi && phi<=1)
            error('thightness on lag decay expected to be between 0 and 1')
        end
        lam_3=1/obj.options.vp_mnst_relative_tightness_constant;
        is_unit_root=true(1,nvar);
        unit_root_vars=obj.options.vp_mnst_unit_root_vars;
        if ischar(unit_root_vars)
            unit_root_vars=cellstr(unit_root_vars);
        end
        if numel(unit_root_vars)==1
            if strcmp(unit_root_vars{1},'none')
                is_unit_root=~is_unit_root;
            end
        else
            is_unit_root(~ismember(obj.endogenous.name,unit_root_vars))=false;
        end
end

MyPriors=set_priors_structure();

obj=setup_priors(obj,MyPriors);

    function sd=quick_ar1_processes()
        nvar=numel(obj.data.varobs_id);
        sd=nan(nvar,1);
        for ivar=1:nvar
            y=obj.data.y(ivar,:)';
            y(isnan(y))=[];
            x=y(1:end-1);
            y=y(2:end);
            coef=x\y;
            resid=y-coef*x;
            sd(ivar)=std(resid);
        end
    end

    function p=set_priors_structure()
        % lag matrices a0, a1,...,ap
        %---------------------------
        [lag_names]=vartools.select_parameter_type(estim_names,'lag_coef');
        decomp=regexp(lag_names,'[a-zA-Z]+(?<lag>\d+)_(?<eqtn>\d+)_(?<vrbl>\d+)','names');
        decomp=[decomp{:}];
        p=struct();
        for ip=1:numel(lag_names)
            pname=lag_names{ip};
            lag=str2double(decomp(ip).lag);
            eqtn=str2double(decomp(ip).eqtn);
            vn=str2double(decomp(ip).vrbl);
            
            [m,sd]=set_var_prior(eqtn,vn,lag);
            
            % set the hyperparameters directly
            %---------------------------------
            if use_priors
                p.(pname)={m,m,sd,'normal_pdf'};
            else
                p.(pname)={m,m-3*sd,m+3*sd};
            end
        end
        
        % deterministic terms
        %--------------------
        [determ_names]=vartools.select_parameter_type(estim_names,'det_coef'); 
        for ip=1:numel(determ_names)
            pname=determ_names{ip};
            eqtn=str2double(regexprep(pname,'[a-z0-9A-Z]+_([a-z0-9A-Z]+)_.+','$1')); % eqtn=str2double(regexprep(pname,'\w+_(\w+)_.+','$1')); % eqtn=str2double(pname(3));
            [m,sd]=set_var_prior(eqtn,nan,nan);
            if use_priors
                p.(pname)={m,m,sd,'normal_pdf'};
            else
                p.(pname)={m,m-3*sd,m+3*sd};
            end
        end
        
        % standard deviations and correlations
        %-------------------------------------
        [std_names]=vartools.select_parameter_type(estim_names,'stdev_corr'); 
        decomp=regexp(std_names,'[a-zA-Z]+_(?<eqtn>\d+)_(?<vrbl>\d+)','names');
        decomp=[decomp{:}];
        for ip=1:numel(std_names)
            pname=std_names{ip};
            eqtn=str2double(decomp(ip).eqtn);
            v=str2double(decomp(ip).vrbl);
            if v==eqtn
                m=s(eqtn);
                sd=5;
                if use_priors
                    p.(pname)={m,m,sd,'inv_gamma_pdf'};
                else
                    p.(pname)={m,0,m+3*sd};
                end
            else
                % correlation terms truncated to lie inside [-1 1]
                m=0;
                sd=1;
                if use_priors
                    p.(pname)={m,m,sd,'normal_pdf',-1,1};
                else
                    p.(pname)={m,-1,1};
                end
            end
        end
        
        % standard deviations of various processes
        %-----------------------------------------
        [std_names_theta]=vartools.select_parameter_type(estim_names,'theta_coef'); 
        for ip=1:numel(std_names_theta)
            pname=std_names_theta{ip};
            m=0.01;
            sd=100;
            if prior_type
                p.(pname)={m,m,sd,'inv_gamma_pdf'};
            else
                p.(pname)={m,0,m+3*sd};
            end
        end
        
        % AR coefficients on processes
        %-----------------------------
        [rho_names]=vartools.select_parameter_type(estim_names,'ar_coef'); 
        for ip=1:numel(rho_names)
            pname=rho_names{ip};
            m=1;
            sd=1;
            if use_priors
                p.(pname)={m,m,sd,'normal_pdf'};
            else
                p.(pname)={m,0,m+3*sd};
            end
        end
        % transition probabilities
        %-------------------------
        [transprob_names]=vartools.select_parameter_type(estim_names,'trans_probs'); 
        markov_chains=obj.construction_data.markov_chains;
        if ~isempty(markov_chains)
            chain_names=(markov_chains.name);
            for ip=1:numel(transprob_names)
                pname=transprob_names{ip};
                underscores=find(pname=='_');
                cname=pname(1:underscores(1)-1);
                cloc= strcmp(cname,chain_names);
                state=str2double(pname(underscores(2)+1:underscores(3)-1));
                duration=markov_chains(cloc).states_expected_duration;
                h=numel(duration);
                d=real(duration(state));
                abar=imag(duration(state));
                qii=1-1/d;
                a=qii*(h-1)*abar/(1-qii);
                a0=a+(h-1)*abar;
                if h==2
                    m=1-qii;
                    sd=sqrt(abar*(a0-abar)/(a0^2*(a0+1)));
                    if use_priors
                        p.(pname)={m,m,sd,'beta_pdf'};
                    else
                        p.(pname)={m,0,1};
                    end
                else
                    if use_priors
                        error('please remind junior.maih@gmail.com to implement the dirichlet')
                        hyperparams=abar(ones(1,h));
                        hyperparams(state)=a;
                        % find all the parameters belonging to this dirichlet and
                        % set them as a package. The diagonal elements will be
                        % computed as 1-sum(qij)
                        %--------------------------------------------------------
                        p.something=dirichlet(hyperparams);
                    else
                        m=abar/a0;
                        p.(pname)={m,0,1};
                    end
                end
            end
        end
        % put back in the original order
        %-------------------------------
        pp=struct();
        for iname=1:numel(estim_names)
            pp.(estim_names{iname})=p.(estim_names{iname});
        end
        p=pp; clear pp
    end

    function [m,sd]=set_var_prior(eqtn,vn,lag)
        m=0;
        switch prior_type
            case {'diffuse','jeffrey', 'none'}
                sd = sqrt(realmax);
            case 'minnesota'
                is_variable=~isnan(vn);
                if is_variable
                    if lag==1 && vn==eqtn
                        if is_unit_root(vn)
                            m=1;
                        else
                            m=obj.options.vp_mnst_stationary_var_mean;
                        end
                    end
                    % standard deviation
                    %-------------------
                    sd=theta*W_mnst(eqtn,vn)*max(1,lag)^(-phi)*s(vn)/s(eqtn);
                else
                    sd=lam_3*s(eqtn);
                end
            case {'normal_wishart','indep_normal_wishart'}
                sd = sqrt(obj.options.vp_natconj_normwish_variance);
                % Maybe allow the user to specify their own matrix here
                % p.v_prior = nvar + 1; p.S_prior = eye(nvar);
            otherwise
                error(['unknown prior type "',prior_type,'"'])
        end
    end
end