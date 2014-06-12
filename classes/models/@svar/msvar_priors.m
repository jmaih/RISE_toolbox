function obj=msvar_priors(obj)

if isempty(obj)
    obj=struct('minnesota_overall_tightness',3,...1
        'minnesota_relative_tightness_lags',1,...
        'minnesota_relative_tightness_constant',0.1,...
        'minnesota_tightness_on_lag_decay',0.5,...1.2
        'minnesota_weight_on_nvar_sum_coef',1,...
        'minnesota_weight_on_single_dummy_initial',1,...
        'minnesota_co_persistence',5,...
        'minnesota_own_persistence',2,...
        'minnesota_weight_on_variance_covariance',1,...
        'minnesota_flat',0,...
        'minnesota_use_priors',true);
    return
end
estim_names=obj.estimated_parameters_list;
s=quick_ar1_processes();

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
        use_priors=obj.options.minnesota_use_priors;
        % lag matrices a0, a1,...,ap
        %---------------------------
        lag_names=regexp(estim_names,'(?<!w+)a\d+_\d+_\d+(_\w+_\d+)?(?!\w+)','match');
        lag_names=[lag_names{:}];
        % lag_locs=locate_variables(lag_names,estim_names);
        p=struct();
        nvar=numel(s);
        W=eye(nvar);
        W(W==0)=obj.options.minnesota_relative_tightness_lags;
        theta=obj.options.minnesota_overall_tightness;
        phi=obj.options.minnesota_tightness_on_lag_decay;
        if ~(0<=phi && phi<=1)
            error('thightness on lag decay expected to be between 0 and 1')
        end
        for ip=1:numel(lag_names)
            pname=lag_names{ip};
            lag=str2double(pname(2));
            eqtn=str2double(pname(4));
            vn=str2double(pname(6));
            % mean
            %-----
            m=0;
            if lag==1 && vn==eqtn % && is_unit_root(vn)
                m=1;
            end
            % standard deviation
            %-------------------
            sd=theta*W(eqtn,vn)*max(1,lag)^(-phi)*s(vn)/s(eqtn);
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
        determ_names=regexp(estim_names,'c_\d+_\d+(_\w+_\d+)?(?!\w+)','match');
        determ_names=[determ_names{:}];
        % determ_locs=locate_variables(determ_names,estim_names);
        lam_3=1/obj.options.minnesota_relative_tightness_constant;
        for ip=1:numel(determ_names)
            pname=determ_names{ip};
            eqtn=str2double(pname(3));
            m=0;
            sd=lam_3*s(eqtn);
            if use_priors
                p.(pname)={m,m,sd,'normal_pdf'};
            else
                p.(pname)={m,m-3*sd,m+3*sd};
            end
        end
        
        % standard deviations and correlations
        %-------------------------------------
        std_names=regexp(estim_names,'(sig|omg)_\d+_\d+(_\w+_\d+)?(?!\w+)','match');
        std_names=[std_names{:}];
        for ip=1:numel(std_names)
            pname=std_names{ip};
            eqtn=str2double(pname(5));
            v=str2double(pname(7));
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
        std_names_theta=regexp(estim_names,'theta_(a\d+|c|sig|omg)_\d+_\d+(?!\w+)','match');
        std_names_theta=[std_names_theta{:}];
        for ip=1:numel(std_names_theta)
            pname=std_names_theta{ip};
            m=0.01;
            sd=100;
            if use_priors
                p.(pname)={m,m,sd,'inv_gamma_pdf'};
            else
                p.(pname)={m,0,m+3*sd};
            end
        end
        
        % AR coefficients on processes
        %-----------------------------
        rho_names=regexp(estim_names,'rho_(a\d+|c|sig|omg)_\d+_\d+(?!\w+)','match');
        rho_names=[rho_names{:}];
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
        transprob_names=regexp(estim_names,'\w+_tp_\d+_\d+','match');
        transprob_names=[transprob_names{:}];
        % transprob_locs=locate_variables(transprob_names,estim_names);
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
end