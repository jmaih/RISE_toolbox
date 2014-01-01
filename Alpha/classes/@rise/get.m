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
elseif strcmpi(PropertyName,'par_texnames')
    Reply=obj.parameters.tex_name;
elseif strncmpi(PropertyName,'par_list',8)
    Reply=load_list('parameters');
elseif strncmpi(PropertyName,'endo_list',9)
    Reply=load_list('endogenous');
elseif strncmpi(PropertyName,'exo_list',8)
    Reply=load_list('exogenous');
elseif strncmpi(PropertyName,'obs_list',8)
    Reply=load_list('observables');
elseif strcmpi(PropertyName,'chain_list')
    Reply=obj.markov_chains.chain_names;
elseif strcmpi(PropertyName,'regime_list')
    Reply=obj.markov_chains.regime_names;
elseif strcmpi(PropertyName,'state_list')
    Reply=obj.markov_chains.state_names;
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
            Reply.(param_names{iname})=xparam(iname);
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
        for ii=1:obj.endogenous.number(2)
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
    function Reply=load_list(type)
        Reply=obj.(type).name;
        proplength=length('endo_list');% =9
        if length(PropertyName)>proplength
            if ~strcmp(PropertyName(proplength+1),'(')
                error(['left parenthesis expected after endo_list, but found ',PropertyName(proplength+1)])
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