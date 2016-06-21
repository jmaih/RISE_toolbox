classdef vstar % < rise_generic
    
    properties
        nlags=4
        thresholds
        parameters
        endogenous
        exogenous
        deterministic=[]
        constant=true
        data
        debug=false
        solution
    end
    
    properties(Hidden=true)%Access=protected,
        is_estimated
        lower_bound
        upper_bound
        num_regessors
%         data_variables
    end
    
    methods
        
        % Constructor
        %------------
        function obj=vstar(endo_names,varargin)
            
%             obj=obj@rise_generic();
            
            if nargin
                
                obj.endogenous=make_names(endo_names);
                
                xonames=strcat('EPS_',obj.endogenous.name);
                
                obj.exogenous=make_names(xonames);
                
                n=length(varargin);
                
                if rem(n,2)
                    
                    error('arguments must come in pairs')
                    
                end
                
                obj=process_options(obj,varargin{:});
                
                obj=set_baseline_parameters(obj);
                
            end
            
        end
        
        % abstract methods
        %-----------------
        varargout=set_solution_to_companion(varargin)
        
        varargout=problem_reduction(varargin)
        
        varargout=conclude_estimation(varargin)
        
        function obj=solve(obj)
            obj.solution=get_parameters(obj);
        end
        
        % Getting information about VSTAR objects
        %----------------------------------------
        varargout=addparam(varargin)
        
        varargout=comment(varargin)
        
        varargout=companion(varargin)
        
        varargout=eig(varargin)
        
        varargout=get(varargin)
        
        varargout=isexplosive(varargin)
        
        varargout=ispanel(varargin)
        
        varargout=isstationary(varargin)
        
        varargout=length(varargin)
        
        varargout=mean(varargin)
        
        varargout=nfitted(varargin)
        
        varargout=rngcmp(varargin)
        
        varargout=userdata(varargin)
        
        % Referencing VSTAR objects
        %--------------------------
        varargout=group(varargin)
        
%         varargout=subsasgn(varargin)
%         
%         varargout=subsref(varargin)
        
        % Simulation, forecasting and filtering
        %--------------------------------------
        varargout=ferf(varargin)
        
        varargout=filter(varargin)
        
        varargout=forecast(varargin)
        
        varargout=randsample(varargin)
        
        varargout=simulate(varargin)
        
        % Manipulation of VSTARs
        %-----------------------
        varargout=assign(varargin)
        
        varargout=alter(varargin)
        
        varargout=demean(varargin)
        
        % stochastic properties
        %----------------------
        varargout=acf(varargin)
        
        varargout=fmse(varargin)
        
        varargout=vma(varargin)
        
        % Estimation, identification and statistical tests
        %-------------------------------------------------
        function [obj,x1]=estimate(obj,options,varargin)
            
            p=obj.parameters;
            
            x0=p(obj.is_estimated);
            
            lb=obj.lower_bound(obj.is_estimated);
            
            ub=obj.upper_bound(obj.is_estimated);
            
%             x1=bee_gate(@objective,x0,lb,ub,options);
            
            x1=fmincon(@objective,x0,[],[],[],[],lb,ub,[],options);
            
            p(obj.is_estimated)=x1;
            
            obj.parameters=p;
            
            obj=solve(obj);
            
            function [minusLogLik,Incr,retcode,obj2]=objective(x)
                
                p(obj.is_estimated)=x;
                
                [LogLik,Incr,retcode,obj2]=...
                    likelihood_smooth_transition_var(p,obj);
                
                minusLogLik=-LogLik;
                
            end

        end
        
        varargout=infocrit(varargin)
        
        varargout=lrtest(varargin)
        
        varargout=portest(varargin)
        
        varargout=schur(varargin)
        
    end
    
    methods%(Access=private)
        function s=parameters2struct(obj)
            s=get_parameters(obj);
        end
    end
    
end

function [LogLik,Incr,retcode,obj]=likelihood_smooth_transition_var(params,obj)

% parameters are as follows
% p*(nlags*p+1+nx)*nregs unrestricted
% p standard deviations >0
% p*(p-1)/2 correlations [-1,1]
% p*(nregs-1) g parameters >= 0
% threshold parameters that must lie in the boundaries of the threshold
% variable and sorted...

debug=obj.debug;

obj.parameters=params;

obj=solve(obj);

s=obj.solution;

nlags=obj.nlags;

data=double(obj.data);

bad=any(isnan(data),2);

data(bad,:)=[];

offset=0;

endo_data=data(:,1:obj.endogenous.number).';

offset=offset+obj.endogenous.number;

exo_data=data(:,offset+(1:obj.deterministic.number)).'; % dummies and so on

offset=offset+obj.deterministic.number;

thresholds=obj.thresholds;

nthresh=numel(thresholds);

thresholds_data=data(:,offset+(1:nthresh)).';

nregs=nthresh+1;

[p,T0]=size(endo_data);

[y,X,T]=set_regressand_and_regressors();

[G,retcode]=get_transitions();

u=zeros(p,T); % = y;

if debug
    
    Incr0=nan(T,1);
    
end

SIG=s.C*s.C.';


LogLik=-1e+8;

Incr=[];

if retcode
    
    return
    
end

det_SIG=det(SIG);

is_failed=det_SIG<=1e-10;% ~all(isfinite(iSIG(:)));

if is_failed
    
    retcode=30002;
    
    return
    
end

iSIG=SIG\eye(p);

ldet_pl2pi=log(det_SIG)+p*log(2*pi);

for t=1:T
    
    Bt=s.B(:,:,1);
    
    for ireg=2:nregs
        
        Gt=G(ireg-1,t);
        
        Bt=Bt+Gt*s.B(:,:,ireg);
        
    end
    
    ut=y(:,t)-Bt*X(:,t);
    
    u(:,t)=ut;
    
    if debug
        
        Incr0(t)=ldet_pl2pi+ut.'*iSIG*ut;
        
    end
    
end

Incr=-.5*(diag((u.'*iSIG*u))+ldet_pl2pi);

if debug
    
    Incr0=-.5*Incr0;
    
    max(abs(Incr-Incr0))
    
    keyboard
    
end

LogLik=sum(Incr);

if ~utils.error.valid(LogLik)
    
    % nans or inf in likelihood
    retcode=300002;
    
end

    function [y,X,T,num_regessors]=set_regressand_and_regressors()
        
        [y,X,~,T,num_regessors]=vartools.set_y_and_x(...
            endo_data,exo_data,nlags,obj.constant);
        
    end

    function [G,retcode]=get_transitions()
        
        retcode=0;
        
        G=zeros(nregs-1,T0);
        
        for icol=1:nregs-1
            
            sij=thresholds_data(icol,:);
            
            g=obj.solution.thresholds{icol}.g;
            
            c=num2cell(obj.solution.thresholds{icol}.c);
            
            [G(icol,:),retcode]=thresholds(icol).func(sij,g,c{:});
            
            if retcode
                
                break
                
            end
            
        end
                    
        G=G(:,nlags+1:end);
        
    end

end

function [endog]=make_names(endo_names)


endog=cell(2,size(endo_names,2));

iter=0;

while ~isempty(endo_names)
    
    if endo_names{1}(1)=='"'
        
        if iter==0 ||...
                ~isempty(endog{2,iter})||...
                isempty(endog{1,iter})
            error(['description coming before model ',...
                'variable or two consecutive descriptions'])
            
        end
        
        endog{2,iter}=strrep(endo_names{1},'"','');
        
    else
        
        iter=iter+1;
        
        endog{1,iter}=endo_names{1};
        
    end
    
    endo_names=endo_names(2:end);
    
end

endog=endog(:,1:iter);

for ii=1:iter
    
    if isempty(endog{2,ii})
        
        endog{2,ii}=endog{1,ii};
        
    end
    
end

[~,tags]=sort(endog(1,:));

endog=endog(:,tags);

endog=struct('name',{endog(1,:)},...
    'tex_name',{endog(2,:)},'number',numel(tags));

end

function s=get_parameters(obj)

params=obj.parameters;

s=struct();

p=obj.endogenous.number;

thresholds=obj.thresholds;

nthresh=numel(thresholds);

nregs=nthresh+1;

s.B=zeros(p,obj.num_regessors,nregs);

offset=0;

for ireg_=1:nregs
    
    s.B(:,:,ireg_)=reshape(pull_params(1:p*obj.num_regessors),p,obj.num_regessors);
    
end

nchol=1*p*(p+1)/2;

sigpars=pull_params(1:p);

corr_pars=pull_params(1:nchol-p);

s.C=eye(p);

for icol=1:p-1
    
    s.C(icol+1:end,icol)=corr_pars(1:p-icol);
    
    corr_pars(1:p-icol)=[];
    
end

s.C=s.C*diag(sigpars);

s.thresholds=cell(1,nthresh);

for ii=1:nthresh
    
    np=obj.thresholds(ii).np;
    
    myparams=pull_params(1:np);
    
%         'data',thresh_data.data,...
    s.thresholds{ii}=struct('name',[obj.thresholds(ii).name,...
        '{',obj.thresholds(ii).lag,'}'],...
        'func',obj.thresholds(ii).func,...
        'g',myparams(1),'c',myparams(2:end));
    
end

    function pp=pull_params(span)
        
        stretch=offset+span;
        
        pp=params(stretch);
        
        offset=stretch(end);
        
    end

end

function obj=set_baseline_parameters(obj)

nlags=obj.nlags;

p=obj.endogenous.number;

nx=obj.deterministic.number;

obj.num_regessors=p*nlags+obj.constant+nx;

thresholds=obj.thresholds;

nthresh=numel(thresholds);

nregs=nthresh+1;

offset=0;

params=nan(10000,1);

lower_bound=params;

upper_bound=params;

is_estimated=true(size(params));

low=-10; high=10;

nx=obj.deterministic.number;

for ireg_=1:nregs
    
    pp_estim=true(p,obj.num_regessors);
    
    if ireg_==1
        
        pp=[1*eye(p),zeros(p,p*(nlags-1)+obj.constant+nx)];
        
    else
        thresh=thresholds(ireg_-1);
        
        bingo=locate_variables(thresh.controlled_vars,obj.endogenous.name);
        
        pp=zeros(p,obj.num_regessors);
        
        pp_estim(:)=false;
        
        pp_estim(bingo,:)=true;
        
    end
    
    add_params(1:p*obj.num_regessors,low,high,pp,pp_estim);
    
end

nchol=1*p*(p+1)/2;

% standard deviations
add_params(1:p,0,10,[],true(1,p));

% correlations
add_params(1:nchol-p,-1,1,[],true(1,nchol-p));

% transition functions
%---------------------
low_g=0;

g_start=0.5;

high_g=50;

if isempty(obj.data)
    
    error('data are needed to initialize the parameters to estimate')
    
end

for ireg_=2:nregs
    
    thresh=thresholds(ireg_-1);
    
    nameLag=[thresh.name,'{',thresh.lag,'}'];
    
    d=obj.data(nameLag);
    
    y=double(d);
    
    np=thresh.np;
    
    yl=nanmin(y); yu=nanmax(y);
    
    y_=linspace(yl,yu,np-1+2);
    
    low=[low_g;yl*ones(np-1,1)];
    
    high=[high_g;yu*ones(np-1,1)];
    
    g_c=[g_start;y_(2:end-1).'];
    
    add_params(1:np,low,high,g_c,true(1,np));
    
end

params(offset+1:end)=[];

lower_bound(offset+1:end)=[];

upper_bound(offset+1:end)=[];

is_estimated(offset+1:end)=[];

obj.is_estimated=is_estimated;

obj.parameters=params;

obj.lower_bound=lower_bound;

obj.upper_bound=upper_bound;


    function add_params(span,low,high,pp,pp_estim)
        
        stretch=offset+span;
        
        nn=numel(stretch);
        
        lower_bound(stretch)=low;
        
        upper_bound(stretch)=high;
        
        if nargin<4 || isempty(pp)
            
            pp=lower_bound(stretch)+...
                rand(nn,1).*(upper_bound(stretch)-lower_bound(stretch));
            
        end
        
        params(stretch)=pp;
        
        is_estimated(stretch)=pp_estim(:);
        
        offset=stretch(end);
        
    end

end

function obj=process_options(obj,varargin)

% fake deterministic
%--------------------
obj.deterministic=make_names({});
    

% make sure the data are processed last

while ~isempty(varargin) %for iarg=1:2:n
    
    if strcmp(varargin{1},'data') && length(varargin)>2
        
        varargin=[varargin(3:end),varargin(1:2)];
        
        continue
        
    end
    
    switch varargin{1}
        
        case 'thresholds'
            % we need to know:
            % - the name of the threshold variable
            % - its lags
            % - the endogenous variables affected by the threshold.
            % - If the list of endogenous variables affected is empty, we
            % assume the threshold affects all endogenous variables in the
            % system
            % e.g. {'loansRates{-3}',{'growth','inflation','interest rate'})

            thresh=varargin{2};
            
            pattern='(?<name>\w+)(\{|\()?(?<lag>[^})]+)?(\)|\})?';
            
            for ithresh=1:size(thresh,1)
                % 1-transition variable
                %---------------------
                tv0=thresh(ithresh,:);
                
                if iscell(tv0{1})
                    
                    if numel(tv0{1})~=2
                        
                        error('When the transition variable is in a cell, it most be followed by its description')
                    
                    end
                    
                    descript=strrep(tv0{1}{2},'"','');
                    
                    tv0{1}=tv0{1}{1};
                    
                else
                    
                    descript=tv0{1};
                    
                    lp=find(descript=='{');
                    
                    if isempty(lp)
                        
                        lp=find(descript=='(');
                        
                    end
                    
                    descript=descript(1:lp-1);
                    
                end
                
                tv=regexp(tv0{1},pattern,'names');
                
                tv.tex=descript;
                % 2-type of transition
                %---------------------
                tv=process_transition(tv0{2},tv);
                % 3-variables influenced by the transition
                %-----------------------------------------
                % in order to have different parameters estimated for a
                % transition, just declare several transitions with
                % different variables they control...
                controlled_vars=obj.endogenous.name;
                
                if numel(tv0)==3 && ~isempty(tv0{3})
                    
                    controlled_vars=tv0{3};
                    
                end
                
                if ~all(ismember(controlled_vars,obj.endogenous.name))
                    
                    error('all controlled variables must be endogenous')
                    
                end
                
                tv.controlled_vars=controlled_vars;
                
                if ithresh==1
                    
                    obj.thresholds=tv;
                    
                else
                    
                    obj.thresholds(ithresh)=tv;
                    
                end
                
                if isempty(tv.lag)
                    
                    tv.lag='0';
                    
                end
                
            end
            
        case 'constant'
            
            obj.constant=varargin{2};
            
            if ~islogical(obj.constant)
                
                error('constant must be a logical')
                
            end
            
        case 'nlags'
            
            obj.nlags=varargin{2};
            
            if ~(isreal(obj.nlags) && ...
                    isnumeric(obj.nlags) && ...
                    isscalar(obj.nlags) && ...
                    obj.nlags>0 && ...
                    floor(obj.nlags)==ceil(obj.nlags))
                
                error('nlags must be a positive constant')
                
            end
            
        case 'deterministic'
            
            obj.deterministic=make_names(varargin{2});
            
            if any(ismember(obj.deterministic.name,obj.endogenous.name))
                
                error('deterministic variables cannot be endogenous')
                
            end
            
        case 'data'
            % check that all the endogenous, deterministic and threshold
            % variables are provided for...
            mydata=varargin{2};
            
            if isa(mydata,'ts')
                
                mydata=pages2struct(mydata);
                
            end
            
            % collect data
            %-------------
            d=mydata.(obj.endogenous.name{1});
            
            myvars=obj.endogenous.name;
            
            for iv=2:obj.endogenous.number
                
                d=[d,mydata.(obj.endogenous.name{iv})];
                
            end
            
            for iv=1:obj.deterministic.number
                
                d=[d,mydata.(obj.deterministic.name{iv})];
                
                myvars=[myvars,obj.deterministic.name{iv}];
                
            end
            
            for iv=1:numel(obj.thresholds)
                
                thresh_data=mydata.(obj.thresholds(iv).name);
                
                d=[d,thresh_data{str2double(obj.thresholds(iv).lag)}];
                
                myvars=[myvars,[obj.thresholds(iv).name,'{',obj.thresholds(iv).lag,'}']];
            
            end
            
            d.varnames=myvars;
            
            obj.data=d;
            
        otherwise
            
            error(['unknown property ',parser.any2str(varargin{1})])
    end
    
    varargin=varargin(3:end);
    
end

end

function tv=process_transition(transfun,tv)

np=1;

if iscell(transfun)
    
    if ~strcmpi(transfun{1},'logisticn')
        
        error('transfun{1} must be logisticn')
        
    end
    
    np=transfun{2}+np;
    
    transfun=lower(transfun{1});
    
elseif ischar(transfun)
    
    transfun=lower(transfun);
    
    switch lower(transfun)
        
        case {'exponential','logistic'}
            
            np=np+1;
            
        case 'logistic2'
            
            np=np+2;
            
        otherwise
            
            error('unknown transition function')
            
    end
    
else
    
    error('transfun must be a char or a cell with two elements')
    
end

tv.np=np;

tv.func=str2func(transfun);

end