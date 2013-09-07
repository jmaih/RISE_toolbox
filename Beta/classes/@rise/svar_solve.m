function [obj,retcode]=svar_solve(obj,varargin)

% obj=rise(svar_struct)
% svar_truct=struct('model_type','svar','varnames',{{'v1','v2','v3'}});
% these varnames will automatically be varobs as well. The number of lags
% should be entered as construction for proper specification of the
% variable names. The restrictions can be entered at a later stage

Defaults=struct('svar_lags',4,'svar_restrictions','choleski','svar_order',[],'svar_restarts',1);

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    obj=Defaults;
    return
end

nobj=numel(obj);
if nobj>1
    retcode=nan(1,nobj);
    for iobj=1:nobj
        [obj(iobj),retcode(iobj)]=svar_solve(obj(iobj),varargin{:});
    end
    return
end
big_penalty=1e+8;

obj=set_options(obj,varargin{:});

svar_order=obj.options.svar_order;
if isempty(svar_order)
    svar_order={obj.orig_varendo.name};
end
debug=obj.options.debug;
svar_restarts=obj.options.svar_restarts;

% observations
obj=load_data(obj);
data=vertcat(obj.varobs.value);
varexo_data=vertcat(obj.varobs_exo.value);
% correct for the nans...
bad=any(isnan(data));
if any(bad)
    data=data(:,~bad);
%     warning('the observations above contain nans and have been removed from the dataset')
end
[nvar,smpl]=size(data);
% we need to get the number of lags...
nlags=obj.options.svar_lags;

y=data(:,nlags+1:end);
nobs=smpl-nlags;
k=nlags*nvar+1;
X=nan(k,nobs);
for ilag=1:nlags
    X((ilag-1)*nvar+1:ilag*nvar,:)=data(:,(nlags+1:end)-ilag);
end
% add the constant
X(end,:)=1;
% nexo=1;
constant_loc=size(X,1);
% add the other exogenous
if ~isempty(varexo_data)
    % check there is no constant in varexo_data
    low=min(varexo_data,[],2);
    high=max(varexo_data,[],2);
    fixed=low==high;
    if any(fixed)
        error('Constant not allowed in the list of observable variables')
    end
% %     varexo_data(fixed,:)=[];
    X=[X;
        varexo_data(:,nlags+1:end)];
%     nexo=nexo+size(varexo_data,1);
end

% solve the reduced-form
T0=y/X;

obj.T=[T0(:,1:constant_loc-1)
    eye((nlags-1)*nvar),zeros((nlags-1)*nvar,nvar)];

% the residuals
Resid=y-T0*X;
% Covariance of reduced-form errors or Short-run Covariance of endogenous
Resid_Resid=Resid*Resid';
SIGu=Resid_Resid/nobs; % <---SIGu=Resid*Resid'/(nobs-k);
LogLiktmp=my_svar_likelihood(SIGu);

% compute the steady state and store it right away
T_L_1=eye(nvar);
for ilag=1:nlags
    T_L_1=T_L_1-T0(:,(ilag-1)*nvar+1:ilag*nvar);
end
inv_T_L_1=T_L_1\eye(nvar);
ss_and_bgp_start_vals=inv_T_L_1*T0(:,constant_loc);
idx=(1:nvar)'*ones(1,nlags);
ss_and_bgp_start_vals=ss_and_bgp_start_vals(idx(:),1);

retcode=0;
orig_names={obj.orig_varendo.name};
orig_texnames={obj.orig_varendo.tex_name};
if numel(svar_order)~=numel(unique(svar_order))
    error('variables repeated in svar_order')
end
if numel(svar_order)~=numel(orig_names)
    error('number of variables in svar_order must match the number of observable variables in the SVAR')
end
% push the structural form
svar_restrictions=obj.options.svar_restrictions;
if ischar(svar_restrictions) && ismember(svar_restrictions,...
        {'choleski','cholesky','short_choleski','short_cholesky',...
        'long_choleski','long_cholesky'})
    svar_order=locate_variables(svar_order,orig_names); % SHOULD THIS BE VAROBS INSTEAD? COZ OF THE LAGS
    inv_order(svar_order)=1:nvar;
    switch svar_restrictions
        case {'choleski','cholesky','short_choleski','short_cholesky'}
            [R0,p] = chol(SIGu(svar_order,svar_order),'lower');
            if p
                retcode=22;
                return
            end
        case {'long_choleski','long_cholesky'}
            iT_L_1=T_L_1\eye(nvar);
            OMGu=iT_L_1*SIGu*iT_L_1';
            [F00,p] = chol(OMGu(svar_order,svar_order),'lower');
            if p
                retcode=22;
                return
            end
            R0=T_L_1(svar_order,svar_order)*F00;
        otherwise
        error('for svar_restrictions, did you mean cholesky or perhaps long_cholesky?')
    end
    R0=R0(inv_order,inv_order);
    penalty=0;
    fval=-my_svar_likelihood(R0*R0');
elseif iscell(svar_restrictions)
    [R0,fval,penalty]=get_restrictions();
else
    svar_restrictions=fcnchk(svar_restrictions,1);
    [varnames_,lb_,ub_]=svar_restrictions();
    [R0,fval,penalty]=get_restrictions({svar_restrictions,varnames_,lb_,ub_});
%     error('svar_restrictions only accepts strings and cellstr')
end
if debug
    disp('========= SVAR estimation output =========')
    disp(['nonlinear restrictions violation: ',num2str(penalty)])
    disp(['Likelihood of SVAR: ',num2str(-fval+penalty)])
    disp(['Likelihood of unrestricted VAR: ',num2str(LogLiktmp)])
    test=[' ',orig_names
        orig_names',num2cell(R0)
        ];
    disp('Estimated impact matrix (R)')
    disp(test)
    disp(['maximum divergence of covariances :',num2str(max(max(abs(SIGu-R0*R0'))))])
end

obj.log_lik=-fval+penalty;
exo_names={obj.varexo.name};
det_exo_names={obj.varobs_exo.name};
nexo=numel(exo_names);
det_flag=false(1,nexo);
for idet=1:numel(det_exo_names)
    det_flag(strcmp(det_exo_names{idet},exo_names))=true;
end
RR=nan(nvar,nexo);
RR(:,det_flag)=T0(:,constant_loc+1:end);
RR(:,~det_flag)=R0;
obj.R=[RR
    zeros((nlags-1)*nvar,nexo)];
% structural form
A00=-R0\eye(nvar); % note the minus sign
A0=eye(nvar*nlags);
A0(1:nvar,1:nvar)=A00;
Aminus=[A00*T0(:,1:constant_loc-1)
    -eye(nvar*(nlags-1)),zeros(nvar*(nlags-1),nvar)];
B=[eye(nvar)
    zeros(nvar*(nlags-1),nvar)];
new_names=orig_names;
for ilag=1:nlags-1
    new_names=[new_names,strcat('ZZZ_',orig_names,'_',int2str(ilag))]; %#ok<AGROW>
end
new_texnames=[orig_texnames,new_names(nvar+1:end)];
obj.NumberOfEndogenous(2)=nvar*nlags;
[new_names,tags]=sort(new_names);
new_texnames=new_texnames(tags);
Aplus=zeros(nvar*nlags); %#ok<NASGU>
A0=A0(tags,tags); %#ok<NASGU>
Aminus=Aminus(tags,tags); %#ok<NASGU>
B=B(tags,:); %#ok<NASGU>
ss_and_bgp_start_vals=ss_and_bgp_start_vals(tags);
obj.varendo=rise_variable.empty(0,1);
for ii=1:obj.NumberOfEndogenous(2)
    obj.varendo(ii,1)=rise_variable(new_names{ii},'tex_name',new_texnames{ii},...
        'id',ii,'det_steady_state',ss_and_bgp_start_vals(ii),...
        'balanced_growth',0);
end
% add the balance growth and push
ss_and_bgp_start_vals=[ss_and_bgp_start_vals
    zeros(size(ss_and_bgp_start_vals))];
obj.steady_state_and_balanced_growth_path=ss_and_bgp_start_vals;

quick_fill={'A0','Aminus','Aplus','B'};
for ii=1:numel(quick_fill)
    obj.(quick_fill{ii})=eval(quick_fill{ii});
end

    function [R,fval,penalty]=get_restrictions(restrict)
        F00=nan(nvar);
        R00=F00;
        fij_loc=[];
        rest_type='';
        old_param=[];
        old_nlc=0;
        old_R=[];
        old_loglik=[];
        if nargin<1
            nrest=size(svar_restrictions,1);
            nshort=0;
            nlong=0;
            for irest=1:nrest
                eqtn=strcmp(svar_restrictions{irest,1},orig_names);
                shock=strcmp(svar_restrictions{irest,2},orig_names);
                type=lower(svar_restrictions{irest,3});
                switch type
                    case {'short','short_run'}
                        R00(eqtn,shock)=svar_restrictions{irest,4};
                        nshort=nshort+1;
                    case {'long','long_run'}
                        F00(eqtn,shock)=svar_restrictions{irest,4};
                        nlong=nlong+1;
                    case 'sign'
                        error('sign restrictions not implemented at this point')
                    otherwise
                        error(['unknown type of restriction ',svar_restrictions{irest,3}])
                end
            end
            % test the order condition
            if nrest<.5*nvar*(nvar-1)
                error('Order condition:: the model is under-identified')
            elseif nrest==.5*nvar*(nvar-1)
                disp('Order condition:: the model is just identified')
            elseif nrest>.5*nvar*(nvar-1)
                disp('Order condition:: the model is over-identified')
            end
            if nlong==nrest
                param_loc=find(isnan(F00));
                rest_type='long';
            elseif nshort==nrest
                param_loc=find(isnan(R00));
                rest_type='short';
            else
                param_loc=find(isnan(R00));
                rest_type='mixed';
                fij_loc=find(~isnan(F00));
                % in case there is a mix of long and short run restrictions,
                % then we can impose the long-run restrictions using penalties
            end
            nparam=numel(param_loc);
            bounds=[-inf*ones(nparam,1),inf*ones(nparam,1)];
            % the elements on the diagonal should be positive as a
            % normalization assumption
            diag_id=find(eye(nvar));
            bounds(ismember(param_loc,diag_id),1)=0;
%             mainfunc=@(x)svar_residuals(x,param_loc,SIGu,T_L_1,inv_T_L_1,R00,F00,fij_loc,rest_type);
        else
            rest_func=restrict{1};
            reposition=locate_variables(orig_names,restrict{2});
            bounds=[restrict{3},restrict{4}];
            param_loc={rest_func,reposition};
%             mainfunc=@(x)svar_residuals(x,param_loc,SIGu);
            nparam=size(bounds,1);
        end
        [xout,fval]=restart('fmincon',svar_restarts);
        [~,R]=svar_lik(xout);
        penalty=norm(nonlinconstr(xout));
        
        function [xout,fval]=restart(optimizer,n)
            % get the bounds and remove the infinities to generate
            % meaningful start values
            lb=bounds(:,1);lb(isinf(lb))=50*sign(lb(isinf(lb)));
            ub=bounds(:,2);ub(isinf(ub))=50*sign(ub(isinf(ub)));
            xparam=bsxfun(@plus,lb,bsxfun(@times,ub-lb,rand(nparam,n)));
            xout=xparam;
            fval=nan(1,n);
            H=nan(nparam,nparam,n);
            for istart=1:n
                %-------------------------------------------------
                options=obj.options.optimset;
                PROBLEM=struct('solver',optimizer,'objective',@svar_lik,'x0',xparam(:,istart),'lb',bounds(:,1),'ub',bounds(:,2),'options',options,'nonlcon',@nonlinconstr);
                [xout(:,istart),fval(istart),H(:,:,istart)]=estimation_engine(PROBLEM,obj.options.hessian_type,obj.options.estim_blocks);
                %-------------------------------------------------
            end
            [~,tag]=min(fval);
            fval=fval(tag);
            xout=xout(:,tag);
        end
        function [minusloglik,R]=svar_lik(param)
            R=get_impact(param,param_loc,R00,F00,T_L_1,rest_type);
            SIG=R*R';
            LogLik=my_svar_likelihood(SIG);
            minusloglik=-LogLik;
            old_R=R;
            old_loglik=abs(minusloglik);
        end
        function [nlineqc,nleqc]=nonlinconstr(param)
            nlineqc=[];
            if isequal(param,old_param)||~strcmp(rest_type,'mixed')
                nleqc=old_nlc;
            else
                F=inv_T_L_1*old_R;
                nleqc=abs(F(fij_loc));
                nleqc=old_loglik*nleqc(:);
            end
            old_param=param;
        end
    end
    function loglik=my_svar_likelihood(SIG)
        loglik=-big_penalty;
        detSIG=det(SIG);
        if detSIG>0 && rcond(SIG)>1e-12
            iSIG=SIG\eye(nvar);
            if ~(any(any(isnan(iSIG))))
                loglik=-0.5*(nobs*log(detSIG)+nobs*nvar*log(2*pi)+trace(iSIG*Resid_Resid));
            end
        end
% %         if loglik<LogLiktmp
% %             keyboard
% %             loglik=-big_penalty;
% %         end
    end
end

function resvec=svar_residuals(param,param_loc,SIGu,...
    T_L_1,inv_T_L_1,R0,F0,fij_loc,rest_type)
if nargin<9
    rest_type=[];
    if nargin<8
        fij_loc=[];
        if nargin<7
            F0=[];
            if nargin<6
                R0=[];
                if nargin<5
                    inv_T_L_1=[];
                    if nargin<4
                        T_L_1=[];
                    end
                end
            end
        end
    end
end
R=get_impact(param,param_loc,R0,F0,T_L_1,rest_type);
Res2=[];
if strcmp(rest_type,'mixed')
    Res2=get_mix_discrepancy(R,inv_T_L_1,F0,fij_loc);
end

SIGR=R*R';
Res1=SIGR(:)-SIGu(:);
resvec=[Res1
    Res2];
end
function d=get_mix_discrepancy(R,inv_T_L_1,F0,fij_loc)
F=inv_T_L_1*R;
d=F0(fij_loc)-F(fij_loc);
d=d(:);
end
function R=get_impact(param,param_loc,varargin)
outside_restrictions=iscell(param_loc);
if outside_restrictions
    rest_func=param_loc{1};
    reposition=param_loc{2};
    R=rest_func(param);
    R=R(reposition,reposition);
else
    R0=varargin{1};
    F0=varargin{2};
    T_L_1=varargin{3};
    rest_type=varargin{4};
    switch rest_type
        case 'long'
            F=F0;
            F(param_loc)=param;
            R=T_L_1*F;
        case {'short','mixed'}
            R=R0;
            R(param_loc)=param;
    end
end
end

