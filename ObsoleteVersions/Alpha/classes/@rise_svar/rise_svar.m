classdef rise_svar
    %Based on
    % 1- Juan F. Rubio-Ramirez & Daniel F. Waggoner & Tao Zha, 2010.
    % "Structural Vector Autoregressions: Theory of Identification and
    % Algorithms for Inference," Review of Economic Studies, Oxford
    % University Press, vol. 77(2), pages 665-696.
    % 2- Andrew Binning, 2013: "Underidentified SVAR models: A framework for
    % combining short and long-run restrictions with sign-restrictions,"
    % Working Paper 2013/14, Norges Bank.
    properties
        constant
        nlags
        lag_structure_restrictions
        short_run_restrictions
        long_run_restrictions
        sign_restrictions
        endo_names
        exo_names
        data
        estim_start_date
        estim_end_date
        irf_periods
        draws
        prior
        use_priors
        nsteps
    end
    properties(SetAccess=protected)
        reduced_form
        number_of_observations
        sf_A0
        sf_Alags
        endo_nbr
        exo_nbr
        nstates=1;
        R
    end
    properties(Dependent=true)
        f
        sr
        Q
        index
        exact_identification_status
        exact_identification_flag
        is_globally_identified
        Q_index_ident
        is_sign_restriction
        restrictions_structure
    end
    
    methods
        function obj=rise_svar(r,varargin)
            obj.endo_names=sort(r.endogenous(:,1))';
            obj.exo_names=sort(r.exogenous(:,1))';
            obj.endo_nbr=numel(obj.endo_names);
            obj.exo_nbr=numel(obj.exo_names);
            if obj.exo_nbr~=obj.endo_nbr
                error('number of exogenous must be same as number of endogenous')
            end
            quick_fill={'nlags','constant','short_run_restrictions',...
                'long_run_restrictions','lag_structure_restrictions',...
                'sign_restrictions','irf_periods','prior','use_priors'};
            for ifill=1:numel(quick_fill)
                ff=quick_fill{ifill};
                obj.(ff)=r.(ff);
            end
            obj=set(obj,varargin{:});
        end
        function obj=set(obj,varargin)
            for iarg=1:2:length(varargin)
                obj.(varargin{iarg})=varargin{iarg+1};
            end
        end
        function exact_identification_flag=get.exact_identification_flag(obj)
            exact_identification_flag=obj.Q_index_ident{4};
        end
        function restrictions_structure=get.restrictions_structure(obj)
            restrictions_structure=struct('short',[],'lag_structure',[],'long',[]);
            if ~isempty(obj.short_run_restrictions)
                lagrest=obj.short_run_restrictions;
                tmp=regexp(lagrest(:,1),'@','split');
                tmp=vertcat(tmp{:});
                nleadlags=get_max_lead_lag(tmp(:,1));
                restrictions_structure.short=0:nleadlags;
            end
            if ~isempty(obj.lag_structure_restrictions)
                lagrest=obj.lag_structure_restrictions;
                tmp=regexp(lagrest(:,1),'@','split');
                tmp=vertcat(tmp{:});
                nleadlags=get_max_lead_lag(tmp(:,1));
                restrictions_structure.lag_structure=0:nleadlags;
            end
            if ~isempty(obj.long_run_restrictions)
                restrictions_structure.long=inf;
            end
        end
        function sf_residuals=structural_form_residuals(obj)
            sf_residuals=[];
            if ~isempty(obj.rf_residuals) && ~isempty(obj.R)
                sf_residuals_=nan(obj.exo_nbr,obj.number_of_observations,obj.nstates);
                init_date=obs2date(obj.estim_start_date,obj.nlags+1);
                sf_residuals=struct();
                for istate=1:obj.nstates
                    sf_residuals_(:,:,istate)=obj.R{istate}*obj.rf_residuals;
                    if obj.nstates==1
                    sf_residuals=...
                        pages2struct(rise_time_series(init_date,sf_residuals_(:,:,istate)',obj.exo_names));
                    else
                        regime_name=sprintf('regime_%0.0f',istate);
                        sf_residuals.(regime_name)=...
                            page2struct(rise_time_series(init_date,sf_residuals_(:,:,istate)',obj.exo_names));
                    end
                end
            end
        end
        function f=get.f(obj)
            % preallocate f for long-run
            %---------------------------
            restrict=obj.restrictions_structure;
            fshort=ones(obj.endo_nbr*numel(restrict.short),obj.exo_nbr);
            % process short-run restrictions: now we are dealing with
            % future impacts. Therefore signs should be positive and there
            % is no restriction on the leads length
            %----------------------------------------------------------
            
            % process short-run restrictions at all leads
            %--------------------------------------------
            if ~isempty(obj.short_run_restrictions)
                rhs_type='exo';
                lagrest=obj.short_run_restrictions;
                tmp=regexp(lagrest(:,1),'@','split');
                tmp=vertcat(tmp{:});
                for irow=1:size(lagrest,1)
                    v=tmp{irow,1};
                    shockname=tmp{irow,2};
                    fshort=process_lag_structure(fshort,v,0,shockname,rhs_type);
                end
            end

            % process lag structure restrictions at all lags
            %------------------------------------------------
            flag=ones(obj.endo_nbr*numel(restrict.lag_structure),obj.endo_nbr);
            if ~isempty(obj.lag_structure_restrictions)
                rhs_type='endo';
                lagrest=obj.lag_structure_restrictions;
                tmp=regexp(lagrest(:,1),'@','split');
                tmp=vertcat(tmp{:});
                vals=cell2mat(lagrest(:,2));
                for irow=1:size(lagrest,1)
                    v=tmp{irow,1};
                    rhsvarname=tmp{irow,2};
                    flag=process_lag_structure(flag,v,vals(irow),rhsvarname,rhs_type);
                end
            end

            % process long-run restrictions
            %------------------------------
            flong=ones(obj.endo_nbr*numel(restrict.long),obj.exo_nbr);
            if ~isempty(obj.long_run_restrictions)
                SIZ=size(flong);
                tmp=regexp(obj.long_run_restrictions(:,1),'@','split');
                tmp=vertcat(tmp{:});
                % check there is no time subscript in the names
                for irow=1:size(tmp,1)
                    if any(tmp{irow,1}=='(')||any(tmp{irow,1}=='{')
                        error(['time subscripts not allowed in long-run restrictions:: ',tmp{irow,1}])
                    end
                end
                endo_ids=locate_variables(tmp(:,1),obj.endo_names);
                exo_ids=locate_variables(tmp(:,2),obj.exo_names);
                IND = sub2ind(SIZ,endo_ids,exo_ids);
                flong(IND)=0;
            end
            f=[fshort;flag;flong];
            
            function f=process_lag_structure(f,v,vval,shockname,rhs_type)
                if nargin<5
                    rhs_type='exo';
                end
                lag=0;
                if any(v=='(')||any(v=='{')
                    vorig=v;
                    if ~(any(v==')')||any(v=='}'))
                        error(['wrong variable name ',vorig])
                    end
                    lag=regexprep(v,'(\w+)({|\()(\d+|\-\d+|\+\d+)(}|\))','$3');
                    loc=find(v=='(');if isempty(loc),loc=find(v=='{');end
                    v=v(1:loc-1);
                    lag=str2double(lag);
                end
                if strcmp(rhs_type,'endo')&& abs(lag)>obj.nlags
                    error(['for restrictions on the lag structure, the lag cannot exceed the maximum number of lags:: ',vorig])
                end
                pos=abs(lag)+1;
                endo_id=locate_variables(v,obj.endo_names);
                rhs_id=locate_variables(shockname,obj.([rhs_type,'_names']));
                f((pos-1)*obj.endo_nbr+endo_id,rhs_id)=vval;
            end
        end
        function sr=get.sr(obj)
            sr=nan(obj.endo_nbr,obj.exo_nbr);
            if ~isempty(obj.sign_restrictions)
                tmp=regexp(obj.sign_restrictions(:,1),'@','split');
                tmp=vertcat(tmp{:});
                vals=obj.sign_restrictions(:,2);
                endo_ids=locate_variables(tmp(:,1),obj.endo_names);
                exo_ids=locate_variables(tmp(:,2),obj.exo_names);
                for isign=1:numel(endo_ids)
                    vi=vals{isign};
                    if isa(vi,'double')
                        ss=logical(vi);
                    else
                        ss=strcmp(vi,'+')-strcmp(vi,'-');
                    end
                    sr(endo_ids(isign),exo_ids(isign))=ss;
                end
            end
        end
        function is_sign_restriction=get.is_sign_restriction(obj)
            is_sign_restriction=~all(isnan(obj.sr(:)));
        end
        function Q_index_ident=get.Q_index_ident(obj)
            %==============================================================
            % finds the Q matrices that describe the linear restrictions on
            % the shock impact matrix.
            %
            % inputs:
            %       endo_nbr = number of dependent variables
            %       f = matrix of short run, lag-structure and long run
            %       restrictions 
            %
            % outputs:
            %       Q = a cell that contains the linear restrictions for
            %           each equation
            %       index = the original column order of the matrix of
            %           restrictions
            %       flag = indicates whether the model is over, under or
            %           exactly identified
            %==============================================================
            E = eye(obj.endo_nbr);
            Q_init = cell(obj.endo_nbr,2);
            for ii = 1:obj.endo_nbr
                Q_init{ii,1} = double(diag(obj.f*E(:,ii)==0));
                Q_init{ii,2} = rank(Q_init{ii,1});
            end
            for ii = 1:obj.endo_nbr
                temp = Q_init{ii,1};
                Q_init{ii,1} = temp(logical(sum(temp,2)),:);
            end
            [new,ord] = sort([Q_init{:,2}],2,'descend');
            
            % Check exact identification
            %---------------------------
            
            if any(new - (obj.endo_nbr - (1:obj.endo_nbr)) > 0)  % over-identified
                flag = 1;
                ident='over identified';
            elseif all(new - (obj.endo_nbr - (1:obj.endo_nbr)) == 0) % exactly identified
                flag = 0;
                ident='exactly identified';
            elseif any(new - (obj.endo_nbr - (1:obj.endo_nbr)) < 0) % under-identified
                flag = -1;
                ident='under identified';
            end
            % reorder the Qs
            %---------------
            index = nan(obj.endo_nbr,1);
            index(ord) = 1:obj.endo_nbr;
            n=obj.endo_nbr;
            Q = cell(n,1);
            is_globally_identified=true;
            for jj = 1:obj.endo_nbr
                Q{jj} = Q_init{ord(jj),1};
                if is_globally_identified
                    M=[Q{jj}*obj.f
                        eye(jj),zeros(jj,n-jj)];
                    is_globally_identified=rank(M)==n;
                end
            end
            Q_index_ident={Q,index,ident,flag,is_globally_identified};
        end
        function Q=get.Q(obj)
            Q=obj.Q_index_ident{1};
        end
        function index=get.index(obj)
            index=obj.Q_index_ident{2};
        end
        function exact_identification_status=get.exact_identification_status(obj)
            exact_identification_status=obj.Q_index_ident{3};
        end
        function is_globally_identified=get.is_globally_identified(obj)
            is_globally_identified=obj.Q_index_ident{5};
        end
        function obj=estimate_reduced_form(obj)
            if isempty(obj.data)
                error('data not provided')
            end
            [data_values,obj.estim_start_date,obj.estim_end_date]=...
                data_request(obj.data,obj.endo_names,obj.estim_start_date,obj.estim_end_date);
            if any(any(isnan(data_values)))
                error('missing observations in the VAR, estimation cannot take place')
            end
            do_reduced_form();
            function do_reduced_form()
                obj = bvar_reduced_form(obj);
                return
                data_values=transpose(data_values);
                [nSample,nvar] = size(data_values);
                nexo=nvar-obj.endo_nbr;
                % number of coefficients in *each* equation, RHS coefficients only.
                ncoef = obj.endo_nbr*obj.nlags+nexo+obj.constant;
                Ydum=[];
                Xdum=[];
                ndobs=0;
                % construct X for Y = X*B + U where X = X: (T-obj.nlags)*k, Y: (T-obj.nlags)*obj.endo_nbr
                %    columns: k = # of [obj.endo_nbr for 1st lag, ..., obj.endo_nbr for last lag, exo var, const]
                if ~isempty(obj.prior)
                    ndobs=obj.endo_nbr+1;         % number of dummy observations
                    % Dummies on right-hand side
                    %---------------------------
                    Xdum = zeros(ndobs,ncoef);
                    % constant term to be put at the end
                    %-----------------------------------
                    if obj.constant
                        Xdum(end,end) = 1;  % the first obj.endo_nbr periods: no or zero constant!
                    end
                    
                    xdgelint = mean(data_values(1:obj.nlags,1:obj.endo_nbr),1); % mean of the first obj.nlags initial conditions
                    % rows 1:obj.endo_nbr
                    Xdum(1:obj.endo_nbr,1:obj.endo_nbr*obj.nlags)=repmat(diag(xdgelint),1,obj.nlags);
                    % last row
                    tmp=repmat(xdgelint',1,obj.nlags);
                    Xdum(ndobs,1:obj.endo_nbr*obj.nlags) = tmp(:)';
                    
                    % Dummies on left hand side
                    %--------------------------
                    Ydum = zeros(ndobs,obj.endo_nbr);
                    for k=1:obj.endo_nbr
                        Ydum(ndobs,k) = xdgelint(k);
                        Ydum(k,k) = xdgelint(k);
                    end
                    % multiply dummies with hyperparameter
                    %-------------------------------------
                    Xdum(1:obj.endo_nbr,:) = 1*obj.prior.weight_on_nvar_sum_coef*Xdum(1:obj.endo_nbr,:);
                    Ydum(1:obj.endo_nbr,:) = obj.prior.weight_on_nvar_sum_coef*Ydum(1:obj.endo_nbr,:);
                    Xdum(obj.endo_nbr+1,:) = obj.prior.weight_on_single_dummy_initial*Xdum(obj.endo_nbr+1,:);
                    Ydum(obj.endo_nbr+1,:) = obj.prior.weight_on_single_dummy_initial*Ydum(obj.endo_nbr+1,:);
                end
                nobs=nSample-obj.nlags;
                
                Xreal = zeros(nobs,ncoef);
                if obj.constant
                    Xreal(:,end) = 1;
                end
                % add data on other exogenous variables (OTHER THAN THE CONSTANT TERM)
                %---------------------------------------------------------------------
                Xreal(:,obj.endo_nbr*obj.nlags+(1:nexo)) = data_values(obj.nlags+1:end,obj.endo_nbr+(1:nexo));
                % add data on endogenous variables
                %---------------------------------
                for k=1:obj.nlags
                    Xreal(:,obj.endo_nbr*(k-1)+1:obj.endo_nbr*k) = data_values(obj.nlags+1-k:nSample-k,1:obj.endo_nbr);
                end
                Yreal = data_values(obj.nlags+1:nSample,1:obj.endo_nbr);
                
                % concatenate dummies and real data
                %----------------------------------
                Y=[Ydum;Yreal];
                X=[Xdum;Xreal];
                
                % reduced form
                %-------------
                [~,xr]=qr(X,0);
%                 xtx=xr'*xr;
                xty=X'*Y;
%                 [~,yr]=qr(Y,0);
%                 yty=yr'*yr;
                Bh = xr\(xr'\xty);
                e=Y-X*Bh;
                % transpose everything
                %---------------------
%                 xty=xty';
%                 xtx=xtx';
                obj.number_of_observations = nobs+ndobs;
                obj.rf_residuals=e';
                obj.rf_B={Bh'};
                obj.rf_Sigma ={obj.rf_residuals*obj.rf_residuals'/(nobs-ncoefs)};
                obj.rf_C = {chol(obj.rf_Sigma{1},'lower')};
                data_values=transpose(data_values);
            end
        end
        function myirfs=irf(obj)
            redraw_BC=obj.use_priors;
            if redraw_BC||(obj.is_sign_restriction && obj.exact_identification_flag==-1)
                ndraws=obj.draws;
            else
                ndraws=1;
            end
            len=obj.irf_periods;
            n=obj.endo_nbr;
            myirfs=struct();
            B=[];
            C=[];
            %             shocks = eye(obj.exo_nbr);
            idraw=0;
            while idraw<ndraws
                [obj,B,C] = find_structural_form(obj,B,C,redraw_BC);
                is_break=false;
                attempt=idraw+1;
                for istate=1:obj.nstates
                    T_st=companion_form(B{istate},obj.reduced_form.posterior.nx,obj.nlags);
                    R_st=obj.R{istate};
                    for ishock=1:obj.exo_nbr
                        if obj.is_sign_restriction
                          sr_index = ~isnan(obj.sr(:,ishock)); 
                          if any(sr_index)
                              chk=sign(R_st(sr_index,ishock))-obj.sr(sr_index,ishock);
                              if any(chk~=0)
                                  % sign violation
                                  is_break=true;
                                  break
                              end
                          end
                        end
                        tmp=zeros(n,len);
                        y0=[R_st(:,ishock)
                            zeros((obj.nlags-1)*n,1)];
                        tmp(:,1)=y0(1:n,1);
                        for t=2:len
                            y0=T_st*y0;
                            tmp(:,t)=y0(1:n,1);
                        end
                        myirfs.(obj.exo_names{ishock})(:,:,attempt,istate)=tmp;
                    end
                end
                if ~is_break
                    idraw=idraw+1;
                    fprintf(1,'%0.0f in %0.0f\n',idraw,ndraws);
                end
            end
            for ishock=1:obj.exo_nbr
                data_=myirfs.(obj.exo_names{ishock});
                myirfs.(obj.exo_names{ishock})=rise_time_series(1,permute(data_,[2,1,3]),obj.endo_names);
            end
        end
        function [obj,B,C] = find_structural_form(obj,B,C,redraw_BC)
            if nargin<4
                redraw_BC=[];
                if nargin<3
                    C=[];
                    if nargin<2
                        B=[];
                    end
                end
            end
            if isempty(redraw_BC)
                redraw_BC=obj.use_priors;
            end
            if isempty(B)||isempty(C)
                obj=estimate_reduced_form(obj);
            end
            if redraw_BC
                Sigma = rand_inverse_wishart(obj.endo_nbr,...
                    obj.reduced_form.posterior.df,...
                    obj.reduced_form.posterior.S_inv_upper_chol);
                Sigma_upper_chol = chol(Sigma);
                C = transpose(Sigma_upper_chol);
                % Get the Autoregressive matrices from a matrix variate distribution.
                Phi = rand_matrix_normal(obj.reduced_form.posterior.ncoef_per_eqtn,...
                    obj.endo_nbr, obj.reduced_form.posterior.PhiHat,...
                    C, obj.reduced_form.posterior.XXi_lower_chol);
                B={transpose(Phi)};
                C={C};
            else
                if isempty(B)
                    B={transpose(obj.reduced_form.posterior.PhiHat)};
                end
                if isempty(C)
                    C={chol(obj.reduced_form.posterior.S/obj.reduced_form.posterior.df,'lower')};
                end
            end
            for st=1:obj.nstates
                % generate a draw for the impact matrix
                %--------------------------------------
                if redraw_BC
                    CC_st=C{st};
                else
                    CC_st=generate_draw(C{st});
                end
                % find a rotation matrix of C consistent with the
                % restrictions
                %-------------------------------------------------
                P=compute_rotation(CC_st);
                % Compute corresponding short-run impact
                %---------------------------------------
                obj.R{st} = transpose(CC_st*P);
                % implied state matrices
                %-----------------------
                obj.sf_A0{st}=obj.R{st}\eye(obj.endo_nbr);
                obj.sf_Alags{st}=obj.sf_A0{st}*B{st};
            end
            function C=generate_draw(C)
                k=obj.endo_nbr;
                newmatrix = randn(k);
                [Q,RR] = qr(newmatrix);
                for ii = 1:k
                    if RR(ii,ii)<0
                        Q(:,ii) = -Q(:,ii);
                    end
                end
                C = C*Q;
            end
            function P=compute_rotation(C)
                L0 = C;                
                if obj.constant
                    beta_temp = B{st}(:,1:end-1);
                else
                    beta_temp = B{st};
                end
                restrictions=obj.restrictions_structure;
                % process short-run restrictions
                %-------------------------------
                nrest=numel(restrictions.short);
                nrows=obj.endo_nbr*nrest;
                FS=zeros(nrows,obj.exo_nbr);
                if nrest
                    max_order=restrictions.short(end);
                    FS(1:obj.endo_nbr,:)=L0;
                    T=companion_form(B{st},obj.reduced_form.posterior.nx,obj.nlags);
                    G=[L0;zeros(obj.endo_nbr*(obj.nlags-1),obj.exo_nbr)];
                    y0=G;
                    for iorder=1:max_order
                        y1=T*y0;
                        FS((1:obj.endo_nbr)+obj.endo_nbr*iorder,:)=y1(1:obj.endo_nbr,:);
                        y0=y1;
                    end
                end 
                % process lag structure restrictions
                %-----------------------------------
                nrest=numel(restrictions.lag_structure);
                nrows=obj.endo_nbr*nrest;
                FLAG=zeros(nrows,obj.exo_nbr);  
                if nrest
                    max_order=restrictions.lag_structure(end);
                    iC=C\eye(obj.endo_nbr);
                    FLAG(1:obj.endo_nbr,:)=iC';
                    Ai=iC*beta_temp(:,1:max_order*obj.endo_nbr);
                    for iorder=1:max_order
                        FLAG((1:obj.endo_nbr)+obj.endo_nbr*iorder,:)=Ai(:,1:obj.endo_nbr)';
                        Ai=Ai(:,obj.endo_nbr+1:end);
                    end
                end
                
                % process long run restrictions
                %------------------------------
                nrest=numel(restrictions.long);
                nrows=obj.endo_nbr*nrest;
                FLONG=zeros(nrows,obj.exo_nbr);  
                if nrest
                    beta = zeros(obj.endo_nbr,obj.endo_nbr);
                    for ii = 1:obj.nlags
                        beta_i=beta_temp(:,(ii-1)*obj.endo_nbr+1:ii*obj.endo_nbr);
                        beta = beta + beta_i;
                    end
                    Linf = (eye(obj.endo_nbr)-beta)\L0;
                    FLONG = Linf;
                end
                F=[FS;FLAG;FLONG];
                
                P = zeros(obj.endo_nbr,obj.endo_nbr);
                for ii = 1:obj.endo_nbr
                    if ii == 1
                        Qtilde = obj.Q{ii}*F;
                    else
                        Qtilde = [obj.Q{ii}*F;P'];
                    end
                    [QQ,~] = qr(Qtilde');
                    P_temp = QQ(:,end);
                    P(:,ii) = P_temp;
                end
                P = P(:,obj.index);
            end
        end
        function [sims_no_shock,explosive_vars]=forecast(obj,varargin)
            obj=set(obj,varargin{:});
            obj=estimate_reduced_form(obj);
            [data_values,obj.estim_start_date,obj.estim_end_date]=...
                data_request(obj.data,obj.endo_names,obj.estim_start_date,obj.estim_end_date);
            
            sims_no_shock = NaN(obj.endo_nbr, obj.nsteps, obj.draws);
            % sims_with_shocks = NaN(obj.endo_nbr, nsteps, obj.draws);
            
            % Declaration of the companion matrix
            Companion_matrix = diag(ones(obj.endo_nbr*(obj.nlags-1),1),-obj.endo_nbr);
            
            % Number of explosive VAR models sampled
            explosive_vars = 0;
            % Loop counter initialization
            d = 0;
            while d < obj.draws
                test=0;
                if obj.use_priors
                    Sigma = rand_inverse_wishart(obj.endo_nbr,...
                        obj.reduced_form.posterior.df,...
                        obj.reduced_form.posterior.S_inv_upper_chol);
                    Sigma_upper_chol = chol(Sigma);
                    C = transpose(Sigma_upper_chol);
                    % Get the Autoregressive matrices from a matrix variate distribution.
                    Phi = rand_matrix_normal(obj.reduced_form.posterior.ncoef_per_eqtn,...
                        obj.endo_nbr, obj.reduced_form.posterior.PhiHat,...
                        C, obj.reduced_form.posterior.XXi_lower_chol);
                    B={transpose(Phi)};
                    C={C};
                    % All the eigenvalues of the companion matrix have to be on or inside the unit circle
                    Companion_matrix(1:obj.endo_nbr,:) = Phi(1:obj.endo_nbr*obj.nlags,:)';
                    test = (abs(eig(Companion_matrix)));
                else
                    if isempty(B)
                        B={transpose(obj.reduced_form.posterior.PhiHat)};
                    end
                    if isempty(C)
                        C={chol(obj.reduced_form.posterior.S/obj.reduced_form.posterior.df,'lower')};
                    end
                end
                
                if any(test>1.0000000000001)
                    explosive_vars = explosive_vars+1;
                end
                d = d+1;
                
                % Without shocks
                lags_data = fliplr(data_values(1:obj.endo_nbr,end-obj.nlags+1:end));
                for t = 1:obj.nsteps
                    X = lags_data(:);
                    if obj.constant
                        X=[X;1];
                    end
                    y = B{1}*X;
                    lags_data=[y,lags_data(:,1:end-1)];
                    sims_no_shock(:, t, d) = y;
                end
                
                %    % With shocks
                %    lags_data = fliplr(data_values(1:obj.endo_nbr,end-obj.nlags+1:end));
                %    for t = 1:nsteps
                %        X = [ reshape(flipdim(lags_data, 1)', 1, obj.endo_nbr*obj.nlags) forecast_data.xdata(t, :) ];
                %        shock = (Sigma_lower_chol * randn(obj.endo_nbr, 1))';
                %        y = X * Phi + shock;
                %        lags_data(1:end-1,:) = lags_data(2:end, :);
                %        lags_data(end,:) = y;
                %        sims_with_shocks(:, t, d) = y;
                %    end
            end
            % add initial conditions
            sims_no_shock=cat(2,repmat(data_values(1:obj.endo_nbr,end-obj.nlags+1:end),1,1,obj.draws),sims_no_shock);
            % put to time series format
            sims_no_shock=rise_time_series(obs2date(obj.estim_end_date,-obj.nlags),permute(sims_no_shock,[2,1,3]),obj.endo_names);
% % %             sims_no_shock=pages2struct(sims_no_shock);
        end
    end
    methods(Static)
        function r=template()
            r=struct();
            r.model='svar';
            r.constant=true;
            r.nlags=4;
            r.endogenous={};
            r.exogenous={};
            r.long_run_restrictions={};
            r.short_run_restrictions={};
            r.lag_structure_restrictions={};
            r.sign_restrictions={};
            r.draws = 1000; % Number of draws
            r.irf_periods=40;
            r.use_priors=true;

            r.prior=struct('overall_tightness',3,...1
                'relative_tightness_lags',1,...
                'relative_tightness_constant',0.1,...
                'tightness_on_lag_decay',0.5,...1.2
                'weight_on_nvar_sum_coef',1,...
                'weight_on_single_dummy_initial',1,...
                'co_persistence',5,...
                'own_persistence',2,...
                'weight_on_variance_covariance',1,...
                'flat',0);
            r.markov_chains=struct('name',{},'nstates',{},'controled_parameters',{},'transition_matrix',{});
        end
        function [ydum,xdum,breaks]=varprior(nv,nx,lags,mnprior,vprior)
            %function [ydum,xdum,breaks]=varprior(nv,nx,lags,mnprior,vprior)
            % ydum, xdum:   dummy observation data that implement the prior
            % breaks:       vector of points in the dummy data after which new dummy obs's start
            %                   Set breaks=T+[0;breaks], ydata=[ydata;ydum], xdum=[xdata;xdum], where
            %                   actual data matrix has T rows, in preparing input for rfvar3
            % nv,nx,lags: VAR dimensions
            % mnprior.tight:Overall tightness of Minnesota prior
            % mnprior.decay:Standard deviations of lags shrink as lag^(-decay)
            % vprior.sig:   Vector of prior modes for diagonal elements of r.f. covariance matrix
            % vprior.w:     Weight on prior on vcv.  1 corresponds to "one dummy observation" weight
            %                   Should be an integer, and will be rounded if not.  vprior.sig is needed
            %                   to scale the Minnesota prior, even if the prior on sigma is not used itself.
            %                   Set vprior.w=0 to achieve this.
            % Note:         The original Minnesota prior treats own lags asymmetrically, and therefore
            %                   cannot be implemented entirely with dummy observations.  It is also usually
            %                   taken to include the sum-of-coefficients and co-persistence components
            %                   that are implemented directly in rfvar3.m.  The diagonal prior on v, combined
            %                   with sum-of-coefficients and co-persistence components and with the unit own-first-lag
            %                   prior mean generates larger prior variances for own than for cross-effects even in
            %                   this formulation, but here there is no way to shrink toward a set of unconstrained
            %                   univariate AR's.
            
            % Original file downloaded from:
            % http://sims.princeton.edu/yftp/VARtools/matlab/varprior.m
            
            if ~isempty(mnprior)
                xdum = zeros(lags+1,nx,lags,nv);
                ydum = zeros(lags+1,nv,lags,nv);
                for il = 1:lags
                    ydum(il+1,:,il,:) = il^mnprior.decay*diag(vprior.sig);
                end
                ydum(1,:,1,:) = diag(vprior.sig);
                ydum = mnprior.tight*reshape(ydum,[lags+1,nv,lags*nv]);
                ydum = flipdim(ydum,1);
                xdum = mnprior.tight*reshape(xdum,[lags+1,nx,lags*nv]);
                xdum = flipdim(xdum,1);
                breaks = (lags+1)*[1:(nv*lags)]';
                lbreak = breaks(end);
            else
                ydum = [];
                xdum = [];
                breaks = [];
                lbreak = 0;
            end
            if ~isempty(vprior) && vprior.w>0
                ydum2 = zeros(lags+1,nv,nv);
                xdum2 = zeros(lags+1,nx,nv);
                ydum2(end,:,:) = diag(vprior.sig);
                for i = 1:vprior.w
                    ydum = cat(3,ydum,ydum2);
                    xdum = cat(3,xdum,xdum2);
                    breaks = [breaks;(lags+1)*[1:nv]'+lbreak];
                    lbreak = breaks(end);
                end
            end
            dimy = size(ydum);
            ydum = reshape(permute(ydum,[1 3 2]),dimy(1)*dimy(3),nv);
            xdum = reshape(permute(xdum,[1 3 2]),dimy(1)*dimy(3),nx);
            breaks = breaks(1:(end-1));
        end
    end
end
function maxleadlag=get_max_lead_lag(names)
maxleadlag=0;
for ii=1:numel(names)
    if any(names{ii}=='(')
        opening=strfind(names{ii},'(');
        closing=strfind(names{ii},')');
    elseif any(names{ii}=='}')
        opening=strfind(names{ii},'{');
        closing=strfind(names{ii},'}');
    else
        continue
    end
    maxleadlag=max(maxleadlag,abs(str2double(names{ii}(opening+1:closing-1))));
end
end

function T=companion_form(B,nx,nlags)
B = B(:,1:end-nx);
n=size(B,1);
dd=n*(nlags-1);
T=[B;
    eye(dd),zeros(dd,n)];
end

function obj = bvar_reduced_form(obj)
%function obj = bvar_reduced_form(obj)
% bvar_reduced_form  Routines shared between BVAR methods
% Computes several things for the estimations of a BVAR(obj)
%
% INPUTS:
%     obj: rise_var object
%
% OUTPUTS:obj.reduced_form
%    posterior:     a structure describing the posterior distribution (which is
%                   normal-Inverse-Wishart)
%                   Its fields are:
%                   - df: degrees of freedom of the inverse-Wishart distribution
%                   - S: matrix parameter for the inverse-Wishart distribution
%                   - XXi: first component of the VCV of the matrix-normal
%                     distribution (the other one being drawn from the
%                     inverse-Wishart)
%                   - PhiHat: mean of the matrix-normal distribution
%    prior:         a structure describing the prior distribution
%                   Its fields are the same than for the posterior

% Bring in data
%--------------
train=0;
if isempty(obj.data)
    error('data not provided')
end
[data_values,obj.estim_start_date,obj.estim_end_date]=...
    data_request(obj.data,obj.endo_names,obj.estim_start_date,obj.estim_end_date);
if any(any(isnan(data_values)))
    error('missing observations in the VAR, estimation cannot take place')
end
data_values=transpose(data_values);
nvar = size(data_values,2);
nexo=nvar-obj.endo_nbr;
nx=nexo+obj.constant;

% Set priors
%--------------
mnprior=[];
vprior=[];
flat=0;
if obj.use_priors && ~isempty(obj.prior)
    mnprior.tight = obj.prior.overall_tightness;
    mnprior.decay = obj.prior.tightness_on_lag_decay;
    
    lambda = obj.prior.co_persistence;
    mu = obj.prior.own_persistence;
    flat = obj.prior.flat;
    
    % Use only initializations obj.nlags for the variance prior
    vprior.sig = std(data_values(1:obj.nlags,:))';
    vprior.w = obj.prior.weight_on_variance_covariance;
end
is_minnesota_prior=~isempty(mnprior);
is_varcov_prior=~isempty(vprior) && vprior.w>0;
% Compute reduced form
%---------------------
[ydum, xdum, pbreaks] = varprior(nx);

ydata = data_values(:,1:obj.endo_nbr);
T = size(ydata, 1);
xdata = ones(T,nx);

% Posterior density
%------------------
var = rfvar3([ydata; ydum], [xdata; xdum], [T; T+pbreaks]);
Tu = size(var.u, 1);

posterior.ncoef_per_eqtn=obj.endo_nbr*obj.nlags+nx;
posterior.df = Tu -posterior.ncoef_per_eqtn- flat*(obj.endo_nbr+1);
posterior.S = var.u' * var.u;
posterior.XXi = var.xxi;
posterior.PhiHat = var.B;
posterior.S_inv_upper_chol=chol(inv(posterior.S));
posterior.XXi_lower_chol=chol(posterior.XXi,'lower');
posterior.nx=nx;

% Prior density
Tp = train + obj.nlags;
if nx
    xdata = xdata(1:Tp, :);
else
    xdata = [];
end
prior=[];
if is_minnesota_prior||is_varcov_prior
    varp = rfvar3([ydata(1:Tp,:);ydum],[xdata; xdum],[Tp;Tp+pbreaks]);
    Tup = size(varp.u,1);
    
    prior.df = Tup - obj.endo_nbr*obj.nlags - nx - flat*(obj.endo_nbr+1);
    prior.S = varp.u' * varp.u;
    prior.XXi = varp.xxi;
    prior.PhiHat = varp.B;
    
    if prior.df < obj.endo_nbr
        error('Too few degrees of freedom in the inverse-Wishart part of prior distribution. You should increase training sample size.')
    end
end
obj.reduced_form=struct('prior',prior,'posterior',posterior);

    function [ydum,xdum,breaks]=varprior(nx)
        %function [ydum,xdum,breaks]=varprior(obj.endo_nbr,nx,obj.nlags,mnprior,vprior)
        % ydum, xdum:   dummy observation data that implement the prior
        % breaks:       vector of points in the dummy data after which new dummy obs's start
        %                   Set breaks=T+[0;breaks], ydata=[ydata;ydum], xdum=[xdata;xdum], where
        %                   actual data matrix has T rows, in preparing input for rfvar3
        % obj.endo_nbr,nx,obj.nlags: VAR dimensions
        % mnprior.tight:Overall tightness of Minnesota prior
        % mnprior.decay:Standard deviations of obj.nlags shrink as lag^(-decay)
        % vprior.sig:   Vector of prior modes for diagonal elements of r.f. covariance matrix
        % vprior.w:     Weight on prior on vcv.  1 corresponds to "one dummy observation" weight
        %                   Should be an integer, and will be rounded if not.  vprior.sig is needed
        %                   to scale the Minnesota prior, even if the prior on sigma is not used itself.
        %                   Set vprior.w=0 to achieve this.
        % Note:         The original Minnesota prior treats own obj.nlags asymmetrically, and therefore
        %                   cannot be implemented entirely with dummy observations.  It is also usually
        %                   taken to include the sum-of-coefficients and co-persistence components
        %                   that are implemented directly in rfvar3.m.  The diagonal prior on v, combined
        %                   with sum-of-coefficients and co-persistence components and with the unit own-first-lag
        %                   prior mean generates larger prior variances for own than for cross-effects even in
        %                   this formulation, but here there is no way to shrink toward a set of unconstrained
        %                   univariate AR's.
        
        % Original file downloaded from:
        % http://sims.princeton.edu/yftp/VARtools/matlab/varprior.m
        if is_minnesota_prior
            xdum = zeros(obj.nlags+1,nx,obj.nlags,obj.endo_nbr);
            ydum = zeros(obj.nlags+1,obj.endo_nbr,obj.nlags,obj.endo_nbr);
            for il = 1:obj.nlags
                ydum(il+1,:,il,:) = il^mnprior.decay*diag(vprior.sig);
            end
            ydum(1,:,1,:) = diag(vprior.sig);
            ydum = mnprior.tight*reshape(ydum,[obj.nlags+1,obj.endo_nbr,obj.nlags*obj.endo_nbr]);
            ydum = flipdim(ydum,1);
            xdum = mnprior.tight*reshape(xdum,[obj.nlags+1,nx,obj.nlags*obj.endo_nbr]);
            xdum = flipdim(xdum,1);
            breaks = (obj.nlags+1)*(1:(obj.endo_nbr*obj.nlags))';
            lbreak = breaks(end);
        else
            ydum = [];
            xdum = [];
            breaks = [];
            lbreak = 0;
        end
        if is_varcov_prior
            ydum2 = zeros(obj.nlags+1,obj.endo_nbr,obj.endo_nbr);
            xdum2 = zeros(obj.nlags+1,nx,obj.endo_nbr);
            ydum2(end,:,:) = diag(vprior.sig);
            for i = 1:vprior.w
                ydum = cat(3,ydum,ydum2);
                xdum = cat(3,xdum,xdum2);
                breaks = [breaks;(obj.nlags+1)*(1:obj.endo_nbr)'+lbreak]; 
                lbreak = breaks(end);
            end
        end
        if is_varcov_prior||is_minnesota_prior
            dimy = size(ydum);
            ydum = reshape(permute(ydum,[1 3 2]),dimy(1)*dimy(3),obj.endo_nbr);
            xdum = reshape(permute(xdum,[1 3 2]),dimy(1)*dimy(3),nx);
            breaks = breaks(1:(end-1));
        end
    end

    function var=rfvar3(ydata,xdata,breaks)
        %function var=rfvar3(ydata,obj.nlags,xdata,breaks,lambda,mu)
        % This algorithm goes for accuracy without worrying about memory requirements.
        % ydata:   dependent variable data matrix
        % xdata:   exogenous variable data matrix
        % obj.nlags:    number of obj.nlags
        % breaks:  rows in ydata and xdata after which there is a break.  This allows for
        %          discontinuities in the data (e.g. war years) and for the possibility of
        %          adding dummy observations to implement a prior.  This must be a column vector.
        %          Note that a single dummy observation becomes obj.nlags+1 rows of the data matrix,
        %          with a break separating it from the rest of the data.  The function treats the
        %          first obj.nlags observations at the top and after each "break" in ydata and xdata as
        %          initial conditions.
        % lambda:  weight on "co-persistence" prior dummy observations.  This expresses
        %          belief that when data on *all* y's are stable at their initial levels, they will
        %          tend to persist at that level.  lambda=5 is a reasonable first try.  With lambda<0,
        %          constant term is not included in the dummy observation, so that stationary models
        %          with means equal to initial ybar do not fit the prior mean.  With lambda>0, the prior
        %          implies that large constants are unlikely if unit roots are present.
        % mu:      weight on "own persistence" prior dummy observation.  Expresses belief
        %          that when y_i has been stable at its initial level, it will tend to persist
        %          at that level, regardless of the values of other variables.  There is
        %          one of these for each variable.  A reasonable first guess is mu=2.
        %      The program assumes that the first obj.nlags rows of ydata and xdata are real data, not dummies.
        %      Dummy observations should go at the end, if any.  If pre-sample x's are not available,
        %      repeating the initial xdata(obj.nlags+1,:) row or copying xdata(obj.nlags+1:2*obj.nlags,:) into
        %      xdata(1:obj.nlags,:) are reasonable subsititutes.  These values are used in forming the
        %      persistence priors.
        
        % Original file downloaded from:
        % http://sims.princeton.edu/yftp/VARtools/matlab/rfvar3.m
        
        [T,nvar_] = size(ydata);
        nox = isempty(xdata);
        if ~nox
            [T2,nx] = size(xdata);
        else
            T2 = T;
            nx = 0;
            xdata = zeros(T2,0);
        end
        % note that x must be same length as y, even though first part of x will not be used.
        % This is so that the obj.nlags parameter can be changed without reshaping the xdata matrix.
        if T2 ~= T, error('Mismatch of x and y data lengths'),end
        if nargin < 4
            nbreaks = 0;
            breaks = [];
        else
            nbreaks = length(breaks);
        end
        breaks = [0;breaks;T];
        smpl = [];
        for nb = 1:nbreaks+1
            smpl = [smpl;(breaks(nb)+obj.nlags+1:breaks(nb+1))']; %#ok<*AGROW>
        end
        Tsmpl = size(smpl,1);
        X = zeros(Tsmpl,nvar_,obj.nlags);
        for is = 1:length(smpl)
            X(is,:,:) = ydata(smpl(is)-(1:obj.nlags),:)';
        end
        X = [X(:,:) xdata(smpl,:)];
        y = ydata(smpl,:);
        % Everything now set up with input data for y=Xb+e
        
        % Add persistence dummies
        if is_minnesota_prior && (lambda ~= 0 || mu > 0)
            ybar = mean(ydata(1:obj.nlags,:),1);
            if ~nox
                xbar = mean(xdata(1:obj.nlags,:),1);
            else
                xbar = [];
            end
            if lambda ~= 0
                if lambda>0
                    xdum_ = lambda*[repmat(ybar,1,obj.nlags) xbar];
                else
                    lambda = -lambda;
                    xdum_ = lambda*[repmat(ybar,1,obj.nlags) zeros(size(xbar))];
                end
                ydum_ = zeros(1,nvar_);
                ydum_(1,:) = lambda*ybar;
                y = [y;ydum_];
                X = [X;xdum_];
            end
            if is_minnesota_prior && mu>0
                xdum_ = [repmat(diag(ybar),1,obj.nlags) zeros(nvar_,nx)]*mu;
                ydum_ = mu*diag(ybar);
                X = [X;xdum_];
                y = [y;ydum_];
            end
        end
        
        % Compute OLS regression and residuals
        [vl,d,vr] = svd(X,0);
        di = 1./diag(d);
        B = (vr.*repmat(di',nvar_*obj.nlags+nx,1))*vl'*y;
        u = y-X*B;
        xxi = vr.*repmat(di',nvar_*obj.nlags+nx,1);
        xxi = xxi*xxi';
        
        var.B = B;
        var.u = u;
        var.xxi = xxi;
    end
end

function G = rand_inverse_wishart(m, v, H_inv_upper_chol)
% function G = rand_inverse_wishart(m, v, H_inv_upper_chol)
% rand_inverse_wishart  Pseudo random matrices drawn from an
% inverse Wishart distribution
% G = rand_inverse_wishart(m, v, H_inv_upper_chol)
% Returns an m-by-m matrix drawn from an inverse-Wishart distribution.
%
% INPUTS:
%     m:          dimension of G and H_inv_upper_chol.
%     v:          degrees of freedom, greater or equal than m.
%     H_inv_chol: upper cholesky decomposition of the inverse of the
%                 matrix parameter.
%                 The upper cholesky of the inverse is requested here
%                 in order to avoid to recompute it at every random draw.
%                 H_inv_upper_chol = chol(inv(H))
% OUTPUTS:
%     G:          G ~ IW(m, v, H) where H = inv(H_inv_upper_chol'*H_inv_upper_chol)
%                 or, equivalently, using the correspondence between Wishart and
%                 inverse-Wishart: inv(G) ~ W(m, v, S) where 
%                 S = H_inv_upper_chol'*H_inv_upper_chol = inv(H)
%  
% SPECIAL REQUIREMENT
%     none
%    

X = randn(v, m) * H_inv_upper_chol; 


% At this point, X'*X is Wishart distributed
% G = inv(X'*X);

% Rather compute inv(X'*X) using the SVD
[~,S,V] = svd(X,0);
SSi = 1 ./ (diag(S) .^ 2);
G = (V .* repmat(SSi', m, 1)) * V';
end

function B = rand_matrix_normal(n, p, M, Omega_lower_chol, Sigma_lower_chol)

% function B = rand_matrix_normal(n, p, M, Omega_lower_chol, Sigma_lower_chol)
% Pseudo random matrices drawn from a matrix-normal distribution
% B ~ MN_n*p(M, Omega, Sigma) 
% Equivalent to vec(B) ~ N(vec(Mu), kron(Omega, Sigma))
%
% INPUTS
%    n:                 row
%    p:                 column
%    M:                 (n*p) matrix, mean
%    Omega_lower_chol:  (p*p), lower Cholesky decomposition of Omega,
%                       (Omega_lower_chol = chol(Omega, 'lower'))
%    Sigma_lower_chol:  (n*n), lower Cholesky decomposition of Sigma,
%                       (Sigma_lower_chol = chol(Sigma, 'lower'))
%    
% OUTPUTS
%    B:                 (n*p) matrix drawn from a Matrix-normal distribution
%        
% SPECIAL REQUIREMENTS
%    Same notations than: http://en.wikipedia.org/wiki/Matrix_normal_distribution

B1 = randn(n * p, 1);
B2 = kron(Omega_lower_chol, Sigma_lower_chol) * B1;
B3 = reshape(B2, n, p);
B = B3 + M;
end
%             % process the lag structure
%             %--------------------------
%             the_sign=-1;
%             if isempty(obj.lag_structure_restrictions)
%                 f_lagrest=ones(0,obj.exo_nbr);
%                 lags=[]; % from 0 to nlags
%             else
%                 lags=false(1,obj.nlags+1); % from 0 to nlags
%                 lagrest=obj.lag_structure_restrictions;
%                 tmp=regexp(lagrest(:,1),'@','split');
%                 tmp=vertcat(tmp{:});
%                 f_lagrest=ones(obj.endo_nbr*(obj.nlags+1),obj.exo_nbr);
%                 vals=cell2mat(lagrest(:,2));
%                 for irow=1:size(lagrest,1)
%                     v=tmp{irow,1};
%                     shockname=tmp{irow,2};
%                     [f_lagrest,lags]=process_lag_structure(f_lagrest,lags,v,vals(irow),shockname);
%                 end
%                 lags=find(lags);
%                 % discard the lags that do not contribute
%                 %----------------------------------------
%                 discard=false(obj.endo_nbr*(obj.nlags+1),1);
%                 for ilag=1:numel(lags)
%                     thislag=lags(ilag);
%                     discard((thislag-1)*obj.endo_nbr+1:thislag*obj.endo_nbr)=true;
%                 end
%                 f_lagrest(discard,:)=[];
%             end
%             f=[f;f_lagrest];
