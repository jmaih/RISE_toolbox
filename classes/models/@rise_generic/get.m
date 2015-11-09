function [Reply,retcode]=get(obj,PropertyName)
% GET -- fetches information from rise_generic objects
%
% Syntax
% -------
% ::
%
%   [Reply,retcode]=get(obj,PropertyName)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|rfvar|svar]: model object
%
% - **PropertyName** [char]: name of the property or element desired. This
% includes
%
%   - **'structure'** [char]: derivatives + transition matrices + other
%   information need to solve the model
%
%   - **'definitions'** [char]: definitions. This requires the model to be
%   solved. If there are no definitions or if the definitions have been
%   substituted, the result will be an empty structure.
%
%   - **'solution'** [char]: solution of the model
%
%   - **'trend'|'growth'|'bgp'** [char]: balanced growth path. N.B. For
%   linear variables the BGP is x_t-x_{t-1}, whereas for log-linear
%   variables the BGP is x_t/x_{t-1}
%
%   - **'sstate'|'steadystate'|'steady_state'** [char]: steady state
%
%   - **'parameters'|'par_vals'** [char]: parameter values. It is also
%   possible to further taylor the output:
%       - 'parameters'|'parameters(default)'|'par_vals'|'par_vals(default)'
%       will give the same result (default) result
%       - '...(struct)' will return the parameters in vector of structures,
%       where each structure is a separate regime 
%       - '...(cell)' will return the parameters in a cell array in which
%       the first column is the list of the parameters and the subsequent
%       columns are the different regimes. 
%
%   - **'par_list'** [char]: list of parameters. Instead of the full list,
%   a sub-list or its complement (using a "~" sign in from of the
%   attribute) can also be queried: 
%       - '...(switching)' : list of parameters that are switching
%       - '...(trans_prob)' : list of transition probability parameters 
%       - '...(measurement_error)' : list of measurement-error parameters
%       - '...(in_use)' : list of parameters that are in use
%
%   - **'par_tex'** [char]: description of the parameters
%
%   - **'endo_list'** [char]: list of the endogenous variables. Instead of
%   the full list, a sub-list or its complement (using a "~" sign in from
%   of the attribute) can also be queried: 
%       - '...(lagrange_multiplier)' : Lagrange multipliers for optimal
%       policy models
%       - '...(static)': static variables or variables appearing only
%       contemporaneously in the model
%       - '...(predetermined)': predetermined variables i.e. variables
%       appearing with lags and not with leads
%       - '...(pred_frwrd_looking)': variables appearing with both a lead
%       and a lag
%       - '...(state)': endogenous state variables i.e. all variables
%       appearing with a lag
%       - '...(frwrd_looking)': variables appearing with leads but not lags
%       - '...(log_var)': variables for which a log-linear approximation is
%       declared from within the model file.
%       - '...(log_expanded)': variables for which a log-linear
%       approximation is requested after the model object is built.
%       - '...(auxiliary)': auxiliary variables automatically created by
%       RISE for support. Leads and lags greater than 1, lags or leads in
%       parameters or exogenous variables.
%       - '...(original)': endogenous variables declared in the model file
%       i.e. excluding the auxiliary variables
%       - '...(affect_trans_probs)': variables entering the calculation of
%       endogenous probabilities
%       - '...(hybrid_expect)': variables for which a hybrid expectation is
%       taken.
%       - '...(stationary)': variables that are not trending over time.
%
%   - **'endo_tex'** [char]: description of endogenous variables
%
%   - **'exo_list'** [char]: list of exogenous variables. Instead of the
%   full list, a sub-list or its complement (using a "~" sign in from of
%   the attribute) can also be queried: 
%       - '...(observed)' : list of exogenous variables that are observed
%       - '...(in_user)' : list of exogenous variables that appear in the
%       model block
%
%   - **'exo_tex'** [char]: description of exogenous variables
%
%   - **'obs_list'** [char]: list of observable variables. Instead of the
%   full list, a sub-list or its complement (using a "~" sign in from of
%   the attribute) can also be queried: 
%       - '...(endogenous)' : list of observable variables that are
%       endogenous.
%
%   - **'obs_tex'** [char]: description of observable variables
%
%   - **'chain_list'** [char]: list of markov chains
%
%   - **'chain_tex'** [char]: description of markov chains
%
%   - ***'regime_list'* [char]: list of regimes (i.e. composites of
%   states from different chains)
%
%   - **'regime_tex'** [char]: description of regimes
%
%   - **'state_list'** [char]: list of states of all the markov chains
%
%   - **'state_tex'** [char]: description of the states of all the markov
%   chains
%
%   - **'tex'|'description'** [char]: description for all the atoms in the
%   system.
%
%   - **'state_vars'** [char]: variables and their lag structure as
%   required for forecasting.
%
%   - **'mode'** [char]: parameters maximizing the posterior distribution
%
%   - **'prior_mean'** [char]: prior mean of the parameters
%
%   - **'start'** [char]: initial values for estimation (maximization of
%   the posterior)
%
%   - **pname** [char]: name of a particular parameter in the model
%
% Outputs
% --------
%
% - **Reply** []: value for the queried property/information
%
% - **retcode** [numeric]: 0 if an error is not encounted
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

% TODO: 
% create a separate get function for dsge
% ,'mean','median','post_sim_mode',,

if isempty(obj)
    Reply=struct();
    return
end
PropertyName(isspace(PropertyName))=[];

retcode=0;
if strcmpi(PropertyName,'structure')
    [~,retcode,Reply]=solve(obj);
    
elseif strcmp(PropertyName,'definitions')
    Reply=load_definitions();
    
elseif strcmpi(PropertyName,'solution')
    error('solution not gettable yet')
    
elseif ismember(lower(PropertyName),{'trend','growth','bgp'})
    Reply=load_steady_state_or_balanced_growth(obj.solution.bgp);
    
elseif ismember(lower(PropertyName),{'sstate','steadystate','steady_state'})
    if isempty(obj.solution)
        error('The model has not been solved')
    end
    Reply=load_steady_state_or_balanced_growth(obj.solution.ss);
    
elseif strncmpi(PropertyName,'parameters',length('parameters'))||...
        strncmpi(PropertyName,'par_vals',length('par_vals'))
    Reply=load_par_vals();
    
elseif strncmpi(PropertyName,'par_list',length('par_list'))
    Reply=load_list('parameters',length('par_list'));
    
elseif strncmpi(PropertyName,'par_tex',length('par_tex'))
    Reply=load_list('parameters',length('par_tex'),'tex_name');
    
elseif strncmpi(PropertyName,'endo_list',length('endo_list'))
    Reply=load_list('endogenous',length('endo_list'));
    
elseif strncmpi(PropertyName,'endo_tex',length('endo_tex'))
    Reply=load_list('endogenous',length('endo_tex'),'tex_name');
    
elseif strncmpi(PropertyName,'exo_list',length('exo_list'))
    Reply=load_list('exogenous',length('exo_list'));
    
elseif strncmpi(PropertyName,'exo_tex',length('exo_tex'))
    Reply=load_list('exogenous',length('exo_tex'),'tex_name');
    
elseif strncmpi(PropertyName,'obs_list',length('obs_list'))
    Reply=load_list('observables',length('obs_list'));
    
elseif strncmpi(PropertyName,'obs_tex',length('obs_tex'))
    Reply=load_list('observables',length('obs_tex'),'tex_name');
    
elseif any(strcmpi(PropertyName,{'chain_list','regime_list','state_list'}))
    item_=strrep(PropertyName,'list','names');
    Reply=obj.markov_chains.(item_);
    
elseif any(strcmpi(PropertyName,{'chain_tex','regime_tex','state_tex'}))
    tex_item=strcat(PropertyName,'_names');
    tex_item=obj.markov_chains.(tex_item);
    name_item=strrep(PropertyName,'tex','names');
    name_item=obj.markov_chains.(name_item);
    Reply=struct();
    for iname_=1:numel(name_item)
        Reply.(name_item{iname_})=tex_item{iname_};
    end
    
elseif any(strcmpi(PropertyName,{'tex','description'}))
    Reply=struct();
    items={'par_tex','endo_tex','exo_tex','chain_tex','regime_tex','state_tex'};
    for iii=1:numel(items)
        [Replyi]=get(obj,items{iii});
        Reply=utils.miscellaneous.mergestructures(Reply,Replyi);
    end
    
elseif strcmpi(PropertyName,'state_vars')
    Reply=load_state_variables();
    
elseif ismember(lower(PropertyName),{'mode','prior_mean','start'})%,'mean','median','post_sim_mode'
    Reply=load_parameters();
    
else
    par_list=get(obj,'par_list');
    ploc=locate_variables(PropertyName,par_list,true);
    if ~isnan(ploc)
        Reply=obj.parameter_values(ploc,:);
    else
        error(['unknown gettable property ',PropertyName])
    end
end


    function Reply=load_state_variables()
        if isa(obj,'dsge')
            s1=get(obj,'endo_list(state)');
        elseif isa(obj,'svar');
            % vars and svars
            s1=obj.endogenous.name;
        end
        s2={};
        if ~isempty(obj.observables)
            s2=get(obj,'obs_list(~endogenous)');
        end
        Reply=union(s1,s2);
        Reply=parser.lead_names(Reply,0);
        % remove lags
        %------------
        Relags=regexp(Reply,'\w+{.+','match');
        Relags=[Relags{:}];
        if ~isempty(Relags)
            Reply=setdiff(Reply,Relags);
        end
        ncols=numel(Reply);
        lags=zeros(1,ncols);
        if ~isempty(Relags)
            for ilag=1:numel(Relags)
                name=Relags{ilag};
                loc=strfind(name,'{');
                prefix=name(1:loc-1);
                lag=str2double(name(loc+1:end-1));
                pos=find(strcmp(prefix,Reply));
                if isempty(pos)
                    error(['"',prefix,'" is not in the list of state variables but surprisingly "',name,'" is'])
                end
                lags(pos)=max(lags(pos),abs(lag));
            end
        end
        Reply=[Reply;num2cell(lags)];
    end


    function Reply=load_parameters()
        type=lower(PropertyName);
        switch type
            case 'mode'
                xparam=obj.estimation.posterior_maximization.mode;
            case {'mean','median'}
                xparam=obj.estimation.posterior_simulation.(type);
            case 'post_sim_mode'
                xparam=obj.estimation.posterior_simulation.mode;
            case {'prior_mean','start'}
                xparam=vertcat(obj.estimation.priors.(type));
            otherwise
                error([mfilename,':: unrecognized type ',type])
        end
        param_names={obj.estimation.priors.name};
        Reply=struct();
        for iname=1:numel(param_names)
            valid_name=parser.param_texname_to_param_name(param_names{iname});
            Reply.(valid_name)=xparam(iname);
        end
    end


    function Reply=load_definitions()
        if isempty(obj.solution)
            error('The model has not been solved')
        end
        formulae=obj.definitions.dynamic;
        def_vals=obj.solution.definitions;
        ndefs=numel(formulae);
        Reply=struct();
        h=numel(def_vals);
        prototype=nan(1,h);
        for idef=1:ndefs
            for ireg=1:h
                prototype(ireg)=def_vals{ireg}(idef);
            end
            eqloc=strfind(formulae{idef},'=');
            Reply.(formulae{idef}(1:eqloc-1))=prototype;
        end
    end


    function Reply=load_steady_state_or_balanced_growth(item)
        
        if isempty(item{1})
            error('the steady state of the model has not been solved')
        end
        Reply=struct();
        endo_names=obj.endogenous.name;
        h=numel(item);
        tmp=nan(1,h);
        for ii=1:obj.endogenous.number
            for ireg=1:h
                tmp(ireg)=item{ireg}(ii);
            end
            Reply.(endo_names{ii})=tmp;
        end
    end


    function Reply=load_par_vals()
        type='default';
        left_par=PropertyName=='(';
        if any(left_par)
            left_par=find(left_par);
            if ~strcmp(PropertyName(end),')')
                error(['right parenthesis expected at the end of string, but found ',PropertyName(end)])
            end
            type=PropertyName(left_par+1:end-1);
        end
        Reply=struct();
        pnames=obj.parameters.name;
        pvals=obj.parameter_values;
        [np,nregs]=size(pvals);
        switch type
            case {'','default'}
                for ipar=1:np
                    Reply.(pnames{ipar})=pvals(ipar,:);
                end
            case 'struct'
                for ireg=1:nregs
                    for ipar=1:np
                        Reply(ireg).(pnames{ipar})=pvals(ipar,ireg);
                    end
                end
            case 'cell'
                Reply=[pnames(:),num2cell(pvals)];
            otherwise
                error(['unknown type ',type])
        end
    end


    function Reply=load_list(type,proplength,type2)
        if nargin<3
            type2='name';
        end
        Reply=obj.(type).(type2);
        if length(PropertyName)>proplength
            if ~strcmp(PropertyName(proplength+1),'(')
                error(['left parenthesis expected after list, but found ',PropertyName(proplength+1)])
            end
            if ~strcmp(PropertyName(end),')')
                error(['right parenthesis expected at the end of string, but found ',PropertyName(end)])
            end
            item=lower(PropertyName(proplength+2:end-1));
            ff=fieldnames(obj.(type));
            ff=regexp(ff,'(?<!\w+)is_\w+','match');
            ff=[ff{:}];
            negative=strcmp(item(1),'~');
            if negative
                item=item(2:end);
            end
            if ismember(['is_',item],ff)
                mygovern=obj.(type).(['is_',item]);
            else
                if strcmp(item,'stationary') && strcmpi(type,'endogenous')
                    log_vars=obj.endogenous.is_log_var;
                    bgp=cell2mat(obj.solution.bgp);
                    too_small=1e-10;
                    checklog=@(x)~any(abs(bgp(x,:)-1)>too_small);
                    checklev=@(x)~any(abs(bgp(x,:))>too_small);
                    mygovern=false(1,obj.endogenous.number);
                    for ivar=1:obj.endogenous.number
                        mygovern(ivar)=if_then_else(log_vars(ivar),...
                            checklog(ivar),checklev(ivar));
                    end
                else
                    disp(strrep(ff,'is_',''))
                    error(['"',item,'" is not a valid property. The valid ',...
                        ' properties are listed above'])
                end
            end
            if negative
                mygovern=~mygovern;
            end
            Reply=Reply(mygovern);
        elseif strcmp(type2,'tex_name')
            tmp=Reply;
            vnames=obj.(type).name;
            Reply=struct();
            for ivar=1:numel(vnames)
                Reply.(vnames{ivar})=tmp{ivar};
            end
        end
    end
end