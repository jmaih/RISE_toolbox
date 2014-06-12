function [Reply,retcode]=get(obj,PropertyName)

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
    Reply=load_steady_state_or_balanced_growth(obj.solution.ss);
elseif ismember(lower(PropertyName),{'par_vals','parameters'})% <--strcmpi(PropertyName,'par_vals')
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
elseif strcmpi(PropertyName,'chain_list')
    Reply=obj.markov_chains.chain_names;
elseif strcmpi(PropertyName,'regime_list')
    Reply=obj.markov_chains.regime_names;
elseif strcmpi(PropertyName,'state_list')
    Reply=obj.markov_chains.state_names;
elseif strcmpi(PropertyName,'state_vars')
    Reply=load_state_variables();
elseif ismember(lower(PropertyName),{'mode','mean','median','post_sim_mode','prior_mean','start'})
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
            valid_name=parser.param_name_to_valid_param_name(param_names{iname});
            Reply.(valid_name)=xparam(iname);
        end
    end
    function Reply=load_definitions()
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
        for ii=1:obj.endogenous.number(end)
            for ireg=1:h
                tmp(ireg)=item{ireg}(ii);
            end
            Reply.(endo_names{ii})=tmp;
        end
    end
    function Reply=load_par_vals()
        Reply=struct();
        pnames=get(obj,'par_list');
        for ipar=1:numel(pnames)
            Reply.(pnames{ipar})=obj.parameter_values(ipar,:);
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
            if ~ismember(['is_',item],ff)
                disp(strrep(ff,'is_',''))
                error(['"',item,'" is not a valid property. The valid properties are listed above'])
            end
            mygovern=obj.(type).(['is_',item]);
            if negative
                mygovern=~mygovern;
            end
            Reply=Reply(mygovern);
        end
    end
end