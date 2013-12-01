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
%         lag_structure_restrictions
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
    end
    properties(SetAccess=protected)
        rf_B
        rf_Sigma
        rf_C
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
            obj.nlags=r.nlags;
            obj.constant=r.constant;
            obj.short_run_restrictions=r.short_run_restrictions;
            obj.long_run_restrictions=r.long_run_restrictions;
            obj.sign_restrictions=r.sign_restrictions;
            obj.irf_periods=r.irf_periods;
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
            restrictions_structure=[];
            if ~isempty(obj.short_run_restrictions)
                lagrest=obj.short_run_restrictions;
                tmp=regexp(lagrest(:,1),'@','split');
                tmp=vertcat(tmp{:});
                nleadlags=get_max_lead_lag(tmp(:,1));
                restrictions_structure=0:nleadlags;
            end
            if ~isempty(obj.long_run_restrictions)
                restrictions_structure=[restrictions_structure,inf];
            end
        end
        function f=get.f(obj)
            % preallocate f for long-run
            %---------------------------
            leads=obj.restrictions_structure;
            f=ones(obj.endo_nbr*numel(leads),obj.exo_nbr);
            % process short-run restrictions: now we are dealing with
            % future impacts. Therefore signs should be positive and there
            % is no restriction on the leads length
            %----------------------------------------------------------
            
            % process short-run restrictions at all leads
            %--------------------------------------------
            if ~isempty(obj.short_run_restrictions)
                lagrest=obj.short_run_restrictions;
                tmp=regexp(lagrest(:,1),'@','split');
                tmp=vertcat(tmp{:});
                vals=cell2mat(lagrest(:,2));
                for irow=1:size(lagrest,1)
                    v=tmp{irow,1};
                    shockname=tmp{irow,2};
                    f=process_lag_structure(f,v,vals(irow),shockname);
                end
            end

            % process long-run restrictions
            %------------------------------
            if ~isempty(obj.long_run_restrictions)
                nprev=numel(leads);
                SIZ=size(f);
                tmp=regexp(obj.long_run_restrictions(:,1),'@','split');
                tmp=vertcat(tmp{:});
                % check there is no time subscript in the names
                for irow=1:size(tmp,1)
                    if any(tmp{irow,1}=='(')||any(tmp{irow,1}=='{')
                        error(['time subscripts not allowed in long-run restrictions:: ',tmp{irow,1}])
                    end
                end
                vals=cell2mat(obj.long_run_restrictions(:,2));
                endo_ids=locate_variables(tmp(:,1),obj.endo_names);
                exo_ids=locate_variables(tmp(:,2),obj.exo_names);
                IND = sub2ind(SIZ,(nprev-1)*obj.endo_nbr+endo_ids,exo_ids);
                f(IND)=vals;
            end
            function f=process_lag_structure(f,v,vval,shockname)
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
                pos=abs(lag)+1;
                endo_id=locate_variables(v,obj.endo_names);
                exo_id=locate_variables(shockname,obj.exo_names);
                f((pos-1)*obj.endo_nbr+endo_id,exo_id)=vval;
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
            %       f = matrix of short and long run restrictions
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
            for ii = 1:obj.endo_nbr
                index(ord(ii)) = ii;
            end
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
                Y=[];
            else
                [data_values,obj.estim_start_date,obj.estim_end_date]=...
                    data_request(obj.data,obj.endo_names,obj.estim_start_date,obj.estim_end_date);
                if any(any(isnan(data_values)))
                    error('missing observations in the VAR, estimation cannot take place')
                end
                Y=data_values(:,obj.nlags+1:end);
                [nvobs,nobs]=size(Y);
                X=zeros(nvobs*obj.nlags,nobs);
                for ilag=1:obj.nlags
                    X((ilag-1)*nvobs+1:ilag*nvobs,:)=data_values(:,(obj.nlags+1:end)-ilag);
                end
                if obj.constant
                    X=[ones(1,nobs);X];
                end
            end
            if isempty(Y)
                error('data not provided')
            end
            obj.rf_B = {Y/X};
            u = Y - obj.rf_B{1}*X;
            nrows=size(obj.rf_B{1},1);
            obj.rf_Sigma ={u*u'/(nobs-nrows)};
            obj.rf_C = {chol(obj.rf_Sigma{1},'lower')};
        end
        function myirfs=irf(obj)
            if obj.is_sign_restriction && obj.exact_identification_flag==-1
                ndraws=obj.draws;
            else
                ndraws=1;
            end
            len=obj.irf_periods;
            n=obj.endo_nbr;
            myirfs=struct();
            for idraw=1:ndraws
                obj = find_structural_form(obj);
                for istate=1:obj.nstates
                    T_st=companion_form(obj.rf_B{istate},obj.constant,obj.nlags);
                    R_st=obj.R{istate};
                    for ishock=1:obj.exo_nbr
                        tmp=zeros(n,len);
                        y0=[R_st(:,ishock)
                            zeros((obj.nlags-1)*n,1)];
                        tmp(:,1)=y0(1:n,1);
                        for t=2:len
                            y0=T_st*y0;
                            tmp(:,t)=y0(1:n,1);
                        end
                        myirfs.(obj.exo_names{ishock})(:,:,idraw,istate)=tmp;
                    end
                end
            end
            for ishock=1:obj.exo_nbr
                data_=myirfs.(obj.exo_names{ishock});
                myirfs.(obj.exo_names{ishock})=rise_time_series(1,permute(data_,[2,1,3]),obj.endo_names);
            end
        end
        function obj = find_structural_form(obj,B,C)
            obj=estimate_reduced_form(obj);
            if nargin<3
                C=obj.rf_C;
                if nargin<2
                    B=obj.rf_B;
                end
            end
            for st=1:obj.nstates
                % generate a draw for the impact matrix
                %--------------------------------------
                CC_st=generate_draw(C{st});
                % find a rotation matrix of C consistent with the
                % restrictions
                %-------------------------------------------------
                P=compute_rotation(CC_st);
                % Compute corresponding short-run impact
                %---------------------------------------
                obj.R{st} = CC_st*P;
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
                    beta_temp = B{st}(:,2:end);
                else
                    beta_temp = B{st};
                end
                restrictions=obj.restrictions_structure;
                nrest=numel(restrictions);
                nrows=obj.endo_nbr*nrest;
                F=zeros(nrows,obj.exo_nbr);
                % process short-run restrictions
                %-------------------------------
                short=restrictions(~isinf(restrictions));
                if ~isempty(short)
                    max_order=short(end);
                    F(1:obj.endo_nbr,:)=L0;
                    T=companion_form(B{st},obj.constant,obj.nlags);
                    G=[L0;zeros(obj.endo_nbr*(obj.nlags-1),obj.exo_nbr)];
                    y0=G;
                    for iorder=1:max_order
                        y1=T*y0;
                        F((1:obj.endo_nbr)+obj.endo_nbr*iorder,:)=y1(1:obj.endo_nbr,:);
                        y0=y1;
                    end
                end   
                % prepare for the long run
                %-------------------------
                if isinf(restrictions(end))
                    beta = zeros(obj.endo_nbr,obj.endo_nbr);
                    for ii = 1:obj.nlags
                        beta_i=beta_temp(:,(ii-1)*obj.endo_nbr+1:ii*obj.endo_nbr);
                        beta = beta + beta_i;
                    end
                    Linf = (eye(obj.endo_nbr)-beta)\L0;
                    F((nrest-1)*obj.endo_nbr+1:end,:) = Linf;
                end
                
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
%             r.lag_structure_restrictions={};
            r.sign_restrictions={};
            r.draws = 1000; % Number of draws
            r.irf_periods=40;

            r.prior=struct('overall_tightness',1,...
                'relative_tightness_lags',1,...
                'relative_tightness_constant',0.1,...
                'tightness_on_lag_decay',1.2,...
                'weight_on_nvar_sum_coef',1,...
                'weight_on_single_dummy_initial',1);
            r.markov_chains=struct('name',{},'nstates',{},'controled_parameters',{},'transition_matrix',{});
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

function T=companion_form(B,constant,nlags)
if constant
    B = B(:,2:end);
end
n=size(B,1);
dd=n*(nlags-1);
T=[B;
    eye(dd),zeros(dd,n)];
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
