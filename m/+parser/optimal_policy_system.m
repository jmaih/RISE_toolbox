function [modelBlock,dic,jac_toc]=optimal_policy_system(Utility,constraints,dic)

% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


[Policy_equations,sstate_Policy_equations,Policy_vars,...
    discount_name,discount]=policy_warmup(Utility,dic.add_welfare);

% select only the equations that are constraints, leaving aside all others
%--------------------------------------------------------------------------
good=strcmp(constraints(:,4),'normal');

other_equations=constraints(~good,:);

constraints=constraints(good,:);

% force discount to be in the parameter list
%--------------------------------------------
dic.parameters(end+1)=parser.listing('name',discount_name);

numberOfAdditionalEquations=numel(Policy_vars);

dic.auxiliary_variables.model=[dic.auxiliary_variables.model,Policy_vars];
% add utility and welfare
%-------------------------
nc=size(constraints,1);

for iv=1:numberOfAdditionalEquations
    
    dic.endogenous(end+1)=parser.listing('name',Policy_vars{iv});
    %%%%%dic.auxiliary_variables{end+1}=Policy_vars{iv};
    
    [constraints(nc+iv,:),dic]=parser.capture_equations(dic,...
        Policy_equations(iv,:),'model');
    
end

neqtns=size(constraints,1);

endo_names={dic.endogenous.name};

exo_names={dic.exogenous.name};

param_names={dic.parameters.name};

% symbols list
%--------------
symb_list=[dic.definitions(:).',endo_names,exo_names,param_names];

nvars=numel(endo_names);

[eqtns,wrt,occur,nleads,ncurrent,numVars]=get_equations_to_differentiate();

% add the steady states to the symbols list
%-------------------------------------------
sstate_names=regexp(eqtns,'steady_state_\w+','match');

sstate_names=unique([sstate_names{:}]);

symb_list=[sstate_names(:).',symb_list];

order=dic.is_optimal_policy_model+2*dic.is_optimal_simple_rule_model;

verbose=false;
tic

[derivs,args]=take_derivatives(eqtns,symb_list,wrt,order,verbose);

dic.planner_system.utils_derivs=[];

Lost_equations=[];

do_utility_differentiation();

if ~dic.is_optimal_simple_rule_model
    
    proto=splanar();
    % create map of derivatives according to... the map
    %---------------------------------------------------
    fx=proto; fx.func=0;
    
    fx=repmat(fx,neqtns+1,numVars);
    
    for ieqtn=1:neqtns+1
        
        this=derivs.derivatives(ieqtn);
        
        locs=this.location;
        
        locs=locs{2};
        
        for icol=1:numel(locs)
            
            obj=intercept_column(this,icol);
            
            xxx=locs(icol);
            
            fx(ieqtn,xxx)=obj;
            
        end
        
    end
    % collect the derivatives of the utility
    %---------------------------------------
    wx=fx(end,nleads+(1:ncurrent));
    
    fx(end,:)=[];
    
    % create multipliers
    %--------------------
    [MULT_LEAD,MULT,MULT_LAG]=create_multipliers();
    
    % pick up the last symbol
    %-------------------------
    pos = strcmp(symb_list,discount_name);
    
    discount_arg=args{pos};
    
    lead_equation=@(x)leader(x,[endo_names,exo_names]);
    
    lag_equation=@(x)lagger(x,[endo_names,exo_names]);
    
    % provision of steady-state calculation
    %---------------------------------------
    f_bf_1bf_wx=[num2cell(fx(:,nleads+(1:nvars)))
        num2cell(wx)];
    
    wxfx=wx;
    
    Lost_equations=cell(nvars,1);
    
    future=0;
    
    past=nleads+nvars;
    
    for ivar=1:nvars
        
        mfx=splanar(0);
        
        current=nleads+ivar;
        
        is_future=occur(1,ivar);
        
        if is_future
            
            future=future+1;
            
        end
        
        is_past=occur(3,ivar);
        
        if is_past
            
            past=past+1;
            
        end
        
        for ieqtn=1:neqtns
            % future
            %--------
            if is_future
                
                mfx=mfx+1/discount_arg*MULT_LAG(ieqtn)*lag_equation(fx(ieqtn,future));
                
                f_bf_1bf_wx{ieqtn,ivar}=f_bf_1bf_wx{ieqtn,ivar}+1/discount_arg*fx(ieqtn,future);
            
            end
            
            % current
            %---------
            if ~occur(2,ivar)
                
                error(['variable ',endo_names{ivar},' does not appear as current']);
           
            end
            
            mfx=mfx+MULT(ieqtn)*fx(ieqtn,current);
            
            % lag
            %-----
            if is_past
                
                mfx=mfx+discount_arg*MULT_LEAD(ieqtn)*lead_equation(fx(ieqtn,past));
                
                f_bf_1bf_wx{ieqtn,ivar}=f_bf_1bf_wx{ieqtn,ivar}+discount_arg*fx(ieqtn,past);
            
            end
            
            f_bf_1bf_wx{ieqtn,ivar}=char(f_bf_1bf_wx{ieqtn,ivar});
        
        end
        
        wxfx(ivar)=wxfx(ivar)-mfx;
        
        % rebuild the lost equations
        %----------------------------
        Lost_equations{ivar}=char(wxfx(ivar));
        
        % rebuild steady-state equations
        %--------------------------------
        f_bf_1bf_wx{ieqtn+1,ivar}=char(f_bf_1bf_wx{ieqtn+1,ivar});
    
    end
    
    old_size=size(f_bf_1bf_wx);
    
    f_bf_1bf_wx=f_bf_1bf_wx(:);
    
    good_fx=~strcmp(f_bf_1bf_wx,'0');
    
    f_bf_1bf_wx(~good_fx)=[];
    
    % clean up the various systems
    %-------------------------------
    Lost_equations=cleanup(Lost_equations);
    
    % add to the planner system
    %-----------------------------
    [f_bf_1bf_wx]=cleanup(f_bf_1bf_wx,true);
    
    dic.planner_system.static_mult_equations={f_bf_1bf_wx,old_size,good_fx};
    
    dic.planner_system.wrt=wrt(nleads+(1:nvars));

end

% remove the discount parameter created above
%---------------------------------------------
dic.parameters(end)=[];

% replace discount by its value also in the welfare equation
%-----------------------------------------------------------
constraints{end,1}(1,:)=strrep(constraints{end,1}(1,:),discount_name,discount);

modelBlock=[other_equations;constraints;Lost_equations];

% create auxiliary equations if needed: this may create issues for loose
% commitment later on
%-------------------------------------------------------------------------
[dic,modelBlock]=parser.create_auxiliary_equations(dic,modelBlock);

% add auxiliary equations to help solve the steady state
%--------------------------------------------------------
for iv=1:numberOfAdditionalEquations
    
    [tmp,dic]=parser.capture_equations(dic,sstate_Policy_equations(iv,:),'steady_state_model');
    
    dic.auxiliary_equations=[dic.auxiliary_equations
        tmp];

end

jac_toc=toc;

    function []=do_utility_differentiation()
        
        wrt_new=wrt;
        
        if dic.is_optimal_simple_rule_model
            
            utils_derivs=derivs;
            
        else
            % At this stage, the last equation is the utility. But to be on
            % the safe side, let's collect it from the model itself.
            % Let's differentiate it separately twice...
            %--------------------------------------------------------------
            % make sure the differentation is taken with respect to the
            % contemporaneous endogenous variables only
            wrt_new=dic.endogenous_list;
            
            [utils_derivs]=take_derivatives(eqtns(end),symb_list,wrt_new,2,verbose);
            
        end
        
        utils_derivs=splanar.print(utils_derivs);
            
        utils_derivs=utils_derivs(2);
        
        [utils_derivs.derivatives]=cleanup(utils_derivs.derivatives,true);
        
        % parse as steady state
        %------------------------
        utils_derivs.wrt=wrt_new;
        
        dic.planner_system.utils_derivs=utils_derivs;
        
        
    end

    function [eqtns,wrt,occur,nleads,ncurrent,numVars]=get_equations_to_differentiate()
        
        wrt_lags=cell(1,nvars);
        
        wrt_leads=wrt_lags;
        
        occur=false(3,nvars); %[leads,current,lags]'
        
        occur(2,:)=true;
        
        n_inner_eqtns=neqtns*(1-dic.is_optimal_simple_rule_model);
        
        eqtns=cell(n_inner_eqtns+1,1);
        % build equations
        %----------------
        for ieq=1:n_inner_eqtns+1
            
            if ieq<=n_inner_eqtns
                
                eqtn=constraints{ieq,1};
            
            else
                
                eqtn=Utility{1,1};
            
            end
            
            for icol_=1:size(eqtn,2)
                
                thyme=eqtn{2,icol_};
                
                if ~isempty(thyme) && thyme~=0
                    
                    vname=eqtn{1,icol_};
                    
                    vpos=strcmp(vname,endo_names);
                    
                    is_endo=any(vpos);
                    
                    if ~is_endo
                        
                        error('leads and lags in parameters not yet accommodated')
                    
                    end
                    
                    if thyme<0
                        
                        type='XLAG';
                        
                        occur(3,vpos)=true;
                    
                    else
                        
                        type='XLEAD';
                        
                        occur(1,vpos)=true;
                    
                    end
                    
                    new_var=[vname,sprintf('_%s_%0.0f',type,abs(thyme))];
                    
                    if ~any(strcmp(new_var,symb_list))
                        
                        symb_list=[symb_list,new_var];
                    
                    end
                    
                    if thyme<0
                        
                        wrt_lags{vpos}=new_var;
                    
                    else
                        
                        wrt_leads{vpos}=new_var;
                    
                    end
                    
                    eqtn{1,icol_}=new_var;
                
                end
                
            end
            
            eqtns{ieq}=cell2mat(eqtn(1,:));
            
            if strcmp(eqtns{ieq}(end),';')
                
                eqtns{ieq}(end)=[];
            
            end
            
        end
        
        % replace calls to steady states
        %--------------------------------
        eqtns=regexprep(eqtns,'(steady_state|\$)\((\w+)\)','steady_state_$2');
        
        % construct the list of the variables to differentiate
        %-----------------------------------------------------
        wrt=[wrt_leads(occur(1,:)),endo_names,wrt_lags(occur(3,:))];
        
        howmany=sum(occur,2);
        
        nleads=howmany(1);
        
        ncurrent=howmany(2);
        
        nlags=howmany(3);
        
        numVars=nleads+ncurrent+nlags;%numel(wrt);
        
    end

    function [MULT_LEAD,MULT,MULT_LAG]=create_multipliers()
        
        mult='MULT_';
        
        MULT_LEAD=repmat(proto,1,neqtns);
        
        MULT=MULT_LEAD;
        
        MULT_LAG=MULT_LEAD;
        
        mult_names=cell(1,neqtns);
        
        % ensure good representation for alphabetization
        %------------------------------------------------
        ndigits=size(num2str(neqtns),2);
        
        ooo=repmat('0',1,ndigits-1);
        
        xxxxx_=parser.listing('name','xxxxx_');
        
        xxxxx_=xxxxx_(1,ones(1,neqtns));
        
        for ieqtn_=1:neqtns
            
            this_digit=sprintf('%s%0.0f',ooo,ieqtn_);
            
            mult_name=sprintf('%s%s',mult,this_digit(end-ndigits+1:end));

            mult_names{ieqtn_}=mult_name;
            
            % add to the list of endogenous
            %-------------------------------
            xxxxx_(ieqtn_).name=mult_name;
            
            xxxxx_(ieqtn_).current_name=mult_name;
            
            % create arguments for splanar
            %------------------------------
            MULT(ieqtn_).func=mult_name;
            
            MULT_LEAD(ieqtn_).func=sprintf('%s_XLEAD_1',mult_name);
            
            MULT_LAG(ieqtn_).func=sprintf('%s_XLAG_1',mult_name);
        
        end
        
        dic.endogenous=[dic.endogenous,xxxxx_];
        
        dic.auxiliary_variables.model=[dic.auxiliary_variables.model,mult_names];
    
    end

    function [x,old_size]=cleanup(x,remove_time)
        
        if nargin<2
            
            remove_time=false;
            
        end
        
        x=strrep(x,'.*','*');
        
        x=strrep(x,'.^','^');
        
        x=strrep(x,'./','/');
        
        if remove_time
            
            x=regexprep(x,'(\w+)(_XLEAD_|_XLAG_)(\d+)','$1');
            
        else
            
            try_again=@do_it_again; %#ok<NASGU>
            
            x=regexprep(x,'(\w+)(_XLEAD_|_XLAG_)(\d+)','$1${try_again($2,$3)}');
            
        end
        
        % add ;
        %-------
        x=strcat(x,';');
        
        % replace discount by its value
        %-------------------------------
        x=strrep(x,discount_name,discount);
        
        % replace steady states (steady_state_C --> $(C))
        %-----------------------------------------------------------
        x=regexprep(x,'\<steady_state_(\w+)','\$($1)');
        
        % add a column of nan at the beginning and a column of opt pol at the end
        %-------------------------------------------------------------------------
        old_size=size(x);
        
        x=x(:);
        
        nlost=size(x,1);
        
        x=[repmat({nan},nlost,1),...
            x,...
            repmat({'optimal policy eqtns'},nlost,1)];
        
        [x,dic]=parser.capture_equations(dic,x,'model');
        
    end

end

function [Policy_equations,sstate_Policy_equations,Policy_vars,...
    discount_name,discount]=policy_warmup(Utility,add_welfare)
% welfare does not behave well in the system. It tends to create a
% near-nonstationarity, inducing numerical inaccuracies. In particular, the
% associated multiplier is close to a unit root, owing to the discount
% factor...
numberOfAdditionalEquations=1+add_welfare;

% prepare the discount
%----------------------
discount=cell2mat(Utility{3,1}(1,:));

leftpart=find(discount=='(',1,'first');

rightpart=find(discount==')',1,'last');

discount=discount(leftpart+1:rightpart-1);

% fake name for the discount factor
%------------------------------------
discount_name='discount___';

% additional information
%------------------------
Policy_equations={
    nan,sprintf('UTIL=%s',cell2mat(Utility{1,1}(1,:))),'Utility'
    nan,sprintf('WELF=(1-(%s))*UTIL+%s*WELF{+1};',...
    discount_name,discount_name),'Welfare'
    };

sstate_Policy_equations={
    nan,sprintf('UTIL=%s',cell2mat(Utility{1,1}(1,:))),'Utility'
    nan,'WELF=UTIL;','Welfare'
%     nan,sprintf('WELF=1/(1-(%s))*UTIL;',discount),'Welfare'
    };

Policy_vars={'UTIL','WELF'};

Policy_equations=Policy_equations(1:numberOfAdditionalEquations,:);

sstate_Policy_equations=sstate_Policy_equations(1:numberOfAdditionalEquations,:);

Policy_vars=Policy_vars(1:numberOfAdditionalEquations);

end

function [derivs,args,total_time]=take_derivatives(eqtns,symb_list,wrt,order,verbose)

if nargin<5
    
    verbose=false;
    
    if nargin<4
        
        order=1;
    
    end
    
end

args=splanar.initialize(symb_list,wrt);

original_funcs=eqtns;

for ifunc=1:numel(eqtns)
    
    [ioccur,eqtns{ifunc}]=parser.find_occurrences(eqtns{ifunc},symb_list);
    
    % re-create the function
    var_occur=symb_list(ioccur);
    
    argfun=cell2mat(strcat(var_occur,','));
    
    eqtns{ifunc}=str2func(['@(',argfun(1:end-1),')',eqtns{ifunc}]);
    
    original_funcs{ifunc}=eqtns{ifunc};
    
    arg_occur=args(ioccur);
    
    eqtns{ifunc}=eqtns{ifunc}(arg_occur{:});
    
end

tic

derivs=splanar.differentiate(eqtns,numel(wrt),order,verbose);

total_time=toc;

end

function x=lagger(x,varnames)

x=leader(x,varnames,false);

end

function x=leader(x,varnames,lead_it)

if nargin<3
    
    lead_it=true;
    
end

add=lead_it*1+(1-lead_it)*(-1);

nargs=numel(x.args);

xfunc=x.func;

if nargs
    
    for iarg=1:numel(x.args)
        
        x.args{iarg}=leader(x.args{iarg},varnames,lead_it);
        
    end
    
elseif ischar(xfunc)
    % apply only to variable names
    %-----------------------------
    if any(strcmp(xfunc,varnames))
        
        redo_it=@do_it; %#ok<NASGU>
        
        if isempty(strfind(xfunc,'_XLEAD_')) && isempty(strfind(xfunc,'_XLAG_'))
            
            if add<0
                
                xfunc=sprintf('%s_XLAG_%0.0f',xfunc,-add);
                
            else
                
                xfunc=sprintf('%s_XLEAD_%0.0f',xfunc,add);
                
            end
            
        else
            
            xfunc=regexprep(xfunc,'(\w+)(_XLEAD_|_XLAG_)(\d+)',...
                sprintf('$1${redo_it($2,$3,%0.0f)}',add));
            
        end
        
        x.func=xfunc;
        
    end
    
end

end

function string=do_it_again(lead_or_lag,d)

islag=strcmp(lead_or_lag,'_XLAG_');

if islag
    
    d=['-',d];
    
else
    
    d=['+',d];
    
end

string=['{',d,'}'];

end

function string=do_it(lead_or_lag,d,add)

string='';

islag=strcmp(lead_or_lag,'_XLAG_');

if islag
    
    digit=add-str2double(d);
    
else
    
    digit=str2double(d)+add;
    
end

if digit<0
    
    string=sprintf('_XLAG_%0.0f',abs(digit));
    
elseif digit>0
    
    string=sprintf('_XLEAD_%0.0f',abs(digit));
    
end

end