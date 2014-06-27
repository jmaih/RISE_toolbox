function FinalResults=automatic_model_selection(tsdata,endo_name,exo_names,options)%[final_mod,x]
% automatic model selection for single-equation models
% tsdata: Time series represented as a table
% endo_name: list of endogenous variables, provided they all share the same
% exogenous variables. If not, a loop has to be written. This could be
% improved
% exo_names: list of the exogenous variables
% options
% - lags : lag length: 1:lags for the endogenous and 0:lags for the
% exogenous
% - debug:  false by default
% - start_date
% - end_date
% - alpha : 0.05 by default
% - fixed: list of variables that cannot go out : 'constant' is the default

% TODO:
% 1- blocking when k>>n. This should be implemented as a function calling
% this one...
% 2- clearly define significance and validity: maybe add ARCH tests, etc...
% in order to ensure congruent reductions...
% 3- choice of the number of lags instead of letting the user painfully
% build the x matrix. One should first check that x does not contain a
% constant.

% this is inspired from the following papers:
% - G. Sucarrat and A. Escribano (2011): "Automated Model Selection in
% Finance: General-to-Specific Modelling of the Mean and Volatility
% Specifications"
% - D. F. Hendry and H.-M. Krolzig (2005): "The Properties of Automatic
% Gets Modelling." Economic Journal 115, C32–C61.
% - D. F. Hendry and H.-M. Krolzig (2001): "Automatic Econometric Model
% Selection using PcGets". London: Timberlake Consultants Press.

if nargin<4
    options=[];
end

switch tsdata.frequency
    case {'','D'}
        default_maxlags=1;
    case 'H'
        default_maxlags=2;
    case 'Q'
        default_maxlags=4;
    case 'M'
        default_maxlags=12;
    case 'W'
        default_maxlags=4;
end

diagnostics={'normality','autocorrelation','arch'};

default_options=struct('lags',1:default_maxlags,...
    'fixed',{'constant'},...
    'debug',false,...
    'start_date',[],...
    'end_date',[],...
    'alpha',0.05);
fresh_options=mysetfield(default_options,options);
% hard-coded for the moment
alpha=fresh_options.alpha;%*ones(3,1)
lags=fresh_options.lags;
fixed=fresh_options.fixed;
debug=fresh_options.debug;
start_date=fresh_options.start_date;
end_date=fresh_options.end_date;

if isempty(start_date)
    start_date=tsdata.TimeInfo(1).date;
end
if isempty(end_date)
    end_date=tsdata.TimeInfo(end).date;
end

first=find(strcmp(start_date,{tsdata.TimeInfo.date}));
last=find(strcmp(end_date,{tsdata.TimeInfo.date}));

if ischar(endo_name)
    endo_name=cellstr(endo_name);
end
if ischar(exo_names)
    exo_names=cellstr(exo_names);
end
exo_names=transpose(exo_names(:));
number_of_loops=numel(endo_name);

if ismember('constant',union(endo_name,exo_names))
    error([mfilename,':: variable name ''constant'' is reserved'])
end

yloc=locate_variables(endo_name,tsdata.varnames);
xloc=locate_variables(exo_names,tsdata.varnames);
data=double(tsdata);
first_row=data(1,:);
test=sum(abs(bsxfun(@minus,data,first_row)),1);
if any(test==0)
    error([mfilename,':: One variable is constant, please remove it'])
end
lags=sort(lags);
if any(lags<0) % should I allow for leads?
    error([mfilename,':: No leads (negative lags) allowed for the moment'])
end
maxlag=max(lags);

FinalResults=cell(1,number_of_loops);
for iloop=1:number_of_loops
    yname=endo_name{iloop};
    yloc_i=yloc(iloop);
    discard=xloc==yloc_i;
    xloc_i=xloc(~discard);
    xnames=exo_names(~discard);
    y0=data(first:last,yloc_i);
    x0=data(first:last,xloc_i);
    T=numel(y0);
    select0=maxlag+1:T;
    y=y0(select0);
    rhs=cell(1,0);
    x=[];
    for ilag=1:numel(lags)
        lag_i=lags(ilag);
        lag_string=['{-',num2str(lag_i),'}'];
        if lag_i~=0 % this is in case I add leads later on...
            % then add the endogenous variable's lag first
            x=[x,y0(select0-lag_i)]; %#ok<*AGROW>
            rhs=[rhs,{[yname,lag_string]}];
        end
        % add the lags of the other exogenous variables
        x=[x,x0(select0-lag_i,:)];
        rhs=[rhs,strcat(xnames,lag_string)];
    end
    % first pass at trimming
    it_nan=0;
    while it_nan<numel(y) && all(isnan([y(it_nan+1),x(it_nan+1,:)]))
        it_nan=it_nan+1;
    end
    % Now adjust sample according to the number of lags and it_nan
    y=y(it_nan+1:end);
    x=x(it_nan+1:end,:);
    real_start=maxlag+1+it_nan+(first-1);
    % finally add the constant term at the end
    rhs=[rhs,{'constant'}];
    fixed_indexes=[];
    if ~isempty(fixed)
        fixed_indexes=locate_variables(fixed,rhs);
    end
    x=[x,ones(size(x,1),1)];
    [FinalResults{iloop},disabled]=block_autometrics(fixed_indexes);
    FinalResults{iloop}.var_list=rhs(FinalResults{iloop}.list);
    recursive=[];
    try %#ok<TRYNC>
        fields={'b','bstd','tstat','R2','R2adj'};
        final_list=FinalResults{iloop}.list;
        [ymm,xmm]=remove_nan_rows(y,x(:,final_list));
        smpl=numel(ymm);
        lowest=numel(final_list);
        if smpl>lowest
            for iter=lowest+1:smpl
                res_iter=ols.ordinary_least_squares(ymm(1:iter),xmm(1:iter,:));
                for ifield=1:numel(fields)
                    if iter==lowest+1
                        recursive.(fields{ifield})=[];
                    end
                    recursive.(fields{ifield})=[recursive.(fields{ifield}),res_iter.(fields{ifield})];
                end
            end
        end
    end
    FinalResults{iloop}.Results.recursive=recursive;
    display_results(FinalResults{iloop},yname)
end

    function display_results(M,yname)
        nvar=numel(M.list);
        disp('=================================================================')
        disp('=                                                               =')
        disp([' Modeling [',yname,'] by OLS using Automagic Model Selection'])
        disp([' The GUM has ',int2str(numel(rhs)),' variables'])
        disp([' The estimation sample is: ',tsdata.TimeInfo(real_start).date,' -- ',tsdata.TimeInfo(last).date])
        if any(disabled)
            disp('=                                                               =')
            disp(upper(' The following tests failed at the GUM and were disabled'))
            disp(diagnostics(disabled))
        end
        disp('=                                                               =')
        disp('=================================================================')
        disp(' ')
        CC=cell(nvar,6);
        CC(1:end,1)=transpose(M.var_list);
        CC(1:end,2)=num2cell(M.Results.b); % coeffs
        CC(1:end,3)=num2cell(M.Results.bstd); % sd
        CC(1:end,4)=num2cell(M.Results.tstat); % tvak
        CC(1:end,5)=num2cell(M.Results.pval); % tprob
        CC(1:end,6)=num2cell(M.Results.partial_r2); % partial_r2
        CC=[{'','Coefficient','Std. Error','t-value','t-prob','partial_r2'};CC];
        disp(CC)
        disp(' ')
        k_1=M.Results.k-1;
        T_k=M.Results.T-M.Results.k;
        Ftest=M.Results.Ftest;
        CC={
            'sigma',sqrt(M.Results.sig2),'RSS',M.Results.RSS
            'R^2',M.Results.R2,['F(',int2str(k_1),',',int2str(T_k),')=',num2str(Ftest)],1-cdf('f',Ftest,k_1,T_k)
            'Adj.R^2',M.Results.R2adj,'log-likelihood',M.Results.loglik
            'no. of observations',M.Results.T,'no. of parameters',M.Results.k
            ['mean(',yname,')'],M.Results.ybar,['se(',yname,')'],M.Results.ystd
            };
        disp(CC)
        disp(' ')
        %                0.0304733                 0.156937493
        %                   0.415052   =   9.224 [0.000]**
        %
        % AR 1-7 test:      F(7,162)  =   1.6791 [0.1175]
        % ARCH 1-7 test:    F(7,169)  =   1.4519 [0.1878]
        % Normality test:   Chi^2(2)  =   8.8468 [0.0120]*
        % Hetero test:      F(26,156) =   1.8702 [0.0105]*
        % Hetero-X test:    F(104,78) =   1.8633 [0.0021]**
        % RESET23 test:     F(2,167)  =   3.7542 [0.0254]*
    end

    function [final_model,disabled]=block_autometrics(prefixed)
        if nargin<1
            prefixed=[];
        end
        pa=alpha;
        
        [T,k]=size(x);
        S0=prefixed; % candidate set
        kf=numel(prefixed);% number of terms that are forced in all models
        
        %% decide whether to block it or not. Allow for 10 degrees of freedom...
        use_blocks=T-k<10; % k/T>0.9;

        %%
        last_fixed=prefixed;
        if use_blocks
            BBAR_S=setdiff(1:k,S0);
            ii=0;
            stage='A';
            not_converged=true;
            pa_1=pa;
            pa_2=pa;
            shrinkage=1; % s in paper
            Jsteps=10;
            while not_converged
                Omitted=expansion_step(BBAR_S);
                
                if ii==0 && numel(Omitted)/numel(BBAR_S)<0.1
                    shrinkage=8*shrinkage;
                end
                
                % find S_{i+1} from model selection on union(S_{i},Oi)
                S0Oi=union(S0,Omitted);
                S1=reduction_step(S0Oi,S0,pa_1);
                
                chk_conv=change_stage();
                
                if chk_conv
                    not_converged=~isequal(S0,S1);
                end
                ii=ii+1;
                S0=S1;
                BBAR_S=setdiff(1:k,S0);
            end
            for ii=1:numel(prefixed)
                last_fixed(ii)=find(prefixed(ii)==S1);
            end
        else
            S1=1:k;
        end
        [final_model,disabled]=myautometrics(y,x(:,S1),pa,last_fixed);
        final_model.list=S1(final_model.list);
        
        function [chk_conv]=change_stage()
            chk_conv=false;
            switch stage
                case 'A'
                    stage='B';
                    pa_1=pa;
                    pa_2=pa_1;
                case 'C'
                    stage='D';
                    pa_1=4*pa;
                    pa_2=pa;
                otherwise
                    if isequal(S1,S0)
                        if strcmp(stage,'B')
                            stage='C';
                            pa_1=2*pa;
                            pa_2=pa_1;
                        else
                            chk_conv=true;
                        end
                    end
            end
        end
        
        function Omitted=expansion_step(BBAR_S)
            jj=0;
            not_converged_expansion=jj<Jsteps;
            while not_converged_expansion
                [BBAR_S,c0]=block_partitions(BBAR_S);
                % 1- Run B0 reductions union(B01,S1),...,union(B0B0,S1) keeping S1
                % fixed
                restricted_indexes=S0;
                B1=[];
                for ib=1:numel(BBAR_S)
                    if jj==0
                        BBAR_S{ib}=reduction_step(union(BBAR_S{ib},S0),restricted_indexes,shrinkage*pa_1);
                    else
                        BBAR_S{ib}=reduction_step(union(BBAR_S{ib},S0),restricted_indexes,shrinkage*pa_2);
                    end
                    % 2- let B1 be the selected regressors from all B0k,k=1,2,...,B0
                    B1=union(B1,BBAR_S{ib});
                end
                if small_enough(B1)
                    % 3- stop if dim(B1) small enough
                    jj=Jsteps;
                else
                    % restart using B1 blocks B1 (if necessary shrinking p-value)
                    jj=jj+1;
                    if k0j>=4*c0
                        shrinkage=min(shrinkage/16,1);
                    elseif k0j>=3*c0
                        shrinkage=min(shrinkage/8,1);
                    elseif k0j>=2*c0
                        shrinkage=min(shrinkage/4,1);
                    else
                        shrinkage=min(shrinkage/2,1);
                    end
                end
                BBAR_S=B1;
                not_converged_expansion=jj<Jsteps;
            end
            % 4- upon convergence after J steps, Oi=BJ
            Omitted=BBAR_S;
            
            function flag=small_enough(BB)
                k0j=numel(BB);
                flag=k0j<=c0; % the paper says ci, is it a typo?
            end
            function [Bout,c0]=block_partitions(Bin)
                ki=numel(S0);
                c0=min(128,round(0.4*(T-kf)));
                kb=round(max([c0-ki,c0/4,2]));
                nn=numel(Bin);
                Bin=Bin(randperm(nn));
                nbins=ceil(nn/kb);
                Bout=cell(1,nbins);
                for ic=1:nbins
                    Bout{ic}=Bin((ic-1)*kb+1:min(nn,ic*kb));
                end
            end
        end
        
        function new_indexes=reduction_step(indexes,fixed,pval)
            fixity=nan(size(fixed));
            for id=1:numel(fixity)
                fixity(id)=find(fixed(id)==indexes);
            end
            MM=myautometrics(y,x(:,indexes),pval,fixity);
            new_indexes=indexes(MM.list);
        end       
    end

    function [final_model,disabled]=myautometrics(y,x,alpha,fixed,isGUM)
        if nargin<5
            isGUM=true;
            if nargin<4
                fixed=[];
            end
        end
        [yc,xc]=remove_nan_rows(y,x);

        bad_variables=find(any(xc-real(xc)));
        if ~isempty(bad_variables)
            disp(bad_variables)
            error([mfilename,':: the variables above are not real'])
        end
        
        k=size(xc,2);
        list=(1:k);
        restricted=fixed;
        if ~isempty(restricted) && ~all(ismember(restricted,list))
            error([mfilename,':: restricted variables not in the list'])
        end
        
        % 1- check whether the unrestricted GUM is valid. Let k be the number of
        % regressors including the constant. Then k=k0+k1+1, where k1 is the number
        % of insignificant regressors. The constant is never removed...
        % [tstat,insignificant,re_ordered_list,isvalid,Results]=ols(list);
        GUM=struct('list',list,'active',true,'duplicated',false,'restricted',restricted,'Results',[]);
        [GUM,isvalid,failed_tests]=ols(GUM);
        
        disabled=failed_tests;% normality, autocorrelation, hetero
        
        switch isvalid
            case {0,1} % it is either good or some tests have failed
                k1=GUM.Results.insignificant;
                % 2- if the GUM is valid, define the number of insignificant variables
                % to be the number of search paths, that is k1.
                if k1
                    Models=struct('list',{},'active',{},'duplicated',{},'restricted',{},'Results',{});
                    for ii=1:k1
                        leader=list(GUM.Results.rank==ii);
                        Models(ii).list=setdiff(list,leader);
                        Models(ii).active=true;
                        Models(ii).duplicated=false;
                        Models(ii).restricted=restricted;
                    end
                    % 3- after removal of the first variable in a path, subsequent
                    % simplication in each path is undertaken using "single-path"
                    % GETS search, where the regressor with the highest p-value is
                    % sought deleted at  each simplification.
                    while any([Models.active])
                        for ii=1:numel(Models)
                            list_i=Models(ii).list;
                            for jj=ii+1:numel(Models)
                                if isequal(Models(jj).list,list_i)
                                    Models(ii).duplicated=true;
                                    break
                                end
                            end
                            if Models(ii).active && ~Models(ii).duplicated
                                Models(ii)=evolve_single_path(Models(ii));
                            end
                        end
                        duplicated=find([Models.duplicated]);
                        if ~isempty(duplicated)
%                             if debug
%                                 disp(['number of competing models BEFORE deletion is :',int2str(numel(Models))])
%                                 if numel(Models)==1
%                                     disp('about to delete the only remaining model')
%                                 end
%                             end
                            Models(duplicated)=[];
                            if debug
                                disp(['number of competing models after deletion is :',int2str(numel(Models))])
                            end
                        end
                    end
                    
                    if numel(Models)==0 % no reduction path was successful
                        final_model=GUM;
                    elseif numel(Models)==1
                        final_model=Models;
                    else
                        if isGUM
                            % create a general GUM of all final models and run
                            % again
                            newgum=[];
                            for ii=1:numel(Models)
                                newgum=union(newgum,Models(ii).list);
                            end
                            second_fixed=fixed;
                            for ii=1:numel(second_fixed)
                                second_fixed(ii)=find(fixed(ii)==newgum);
                            end
                            second_round=myautometrics(yc,xc(:,newgum),alpha,second_fixed,false);
                            % adjust for the variable list
                            final_model=second_round;
                            final_model.list=newgum(second_round.list);
                        else
                            % Now select by information criterion
                            AIC=Models(1).Results.criterion.AIC;
                            final_model=Models(1);
                            for ii=2:numel(Models)
                                if Models(ii).Results.criterion.AIC<AIC
                                    AIC=Models(ii).Results.criterion.AIC;
                                    final_model=Models(ii);
                                end
                            end
                        end
                    end
                else
                    if debug
                        disp(upper([mfilename,':: All variables are significant, terminal model is the GUM']))
                    end
                    final_model=GUM;
                end
            otherwise
                error([mfilename,':: too few observations even after blocking'])
        end
        
        function Model=rank_variables(Model)
            tmp=Model.Results.tstat;
            params_pval=Model.Results.pval;
            tmp=abs(tmp);
            % take care of the elite right here right now
            tmp(ismember(Model.list,Model.restricted))=inf;
            params_pval(ismember(Model.list,Model.restricted))=0;
            [~,tag]=sort(tmp);
            Model.Results.rank(tag)=1:numel(tmp);
            Model.Results.insignificant=sum(params_pval>=alpha);
            if debug && ~isempty(Model.restricted)
                if ~all(ismember(Model.restricted,Model.list))
                    keyboard
                end
            end
        end
        
        function M=evolve_single_path(M)
            % we attempt to delete one variable
            if isempty(M.Results)
                [M0,is_still_valid]=ols(M,disabled);
                if is_still_valid
                    M=M0;
                else
                    % Mark the path for deletion
                    M.duplicated=true;
                    return
                    %                     error([mfilename,':: this case is tricky, please send this example to junior.maih@gmail.com'])
                end
            end            
            next=next_candidate();
            M.active=~isempty(next);
            if M.active
                newlist=setdiff(M.list,next);
                M2=M; M2.list=newlist;
                [M2,is_still_valid]=ols(M2,disabled);
                if is_still_valid
                    % check whether the next candidate is part of the
                    % restriced list. It might be the case that because of
                    % some collinearity problem it was not possible to
                    % delete it earlier but that now it makes sense to do
                    % so.
                    loc=find(M2.restricted==next);
                    if ~isempty(loc)
                        M2.restricted(loc)=[];
                        % alternatively, one could just say that the guys
                        % that are in the restricted list can never be
                        % deleted. In that case, the search would stop when
                        % the next candidate is in the restricted list.
                        % But just how does that happen in the first place?
                        % another thing to do is to differentiate between
                        % the variables that have been imposed by the user
                        % and cannot be deleted and the variables that can
                        % be deleted because the algorithm at some step,
                        % found them deletable even if at some earlier
                        % stage this was not the case.
                        if debug
                            disp(['Variable ',int2str(next),' removed from the restricted list and from the model'])
                        end
                    end
                    M=M2;
                else
                    if ismember(next,M.restricted)
                        if debug
                            disp(['attempting to put variable ',int2str(next),' back into the list. Path will end here'])
                        end
                        M.active=false;
                    else
                        M.restricted=union(M.restricted,next);
                    end
                end
                M=rank_variables(M);
            end
            function next=next_candidate()
                % chooses the next candidate for deletion
                next=[];
                position=M.Results.rank==1;
                if M.Results.pval(position)>alpha && ~ismember(M.list(position),M.restricted)
                    next=M.list(position);
                end
            end
        end
        
        function [Model,isvalid,failed_tests]=ols(Model,disabled)
            if nargin<2
                disabled=[];
            end
            [ym,xm]=remove_nan_rows(yc,xc(:,Model.list));
            if isempty(ym)
                error([mfilename,':: all rows are nan'])
            end
            ki=numel(Model.list);
            isvalid=ki>0;
            failed_tests=[];
            if isvalid
                first_check=numel(ym)>=ki;
                isvalid = first_check-~first_check;
                if debug && isvalid==-1
                    disp([mfilename,':: too short sample '])
                    keyboard
                end
                if isvalid==1
                    Model.Results=ols.ordinary_least_squares(ym,xm);
                    Model=rank_variables(Model);
                    
                    tests_pval=[Model.Results.diagnostics.normality.pval
                        Model.Results.diagnostics.autocorrelation.pval
                        Model.Results.diagnostics.arch.pval]';
                    failed_tests=false(size(tests_pval));
                    if isempty(disabled)
                        disabled=false(size(tests_pval));
                    end
                    for it=1:numel(tests_pval)
                        failed_tests(it)=tests_pval(it)<alpha;
                    end
                    isvalid=~any(failed_tests(~disabled));
                end
            end
        end
    end

end

function [ym,xm]=remove_nan_rows(y,x)
nanrows=isnan(y)|any(isnan(x),2);
ym=y(~nanrows);
xm=x(~nanrows,:);
end
