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
elseif strncmpi(PropertyName,'par_list',8)
    Reply=load_par_list();
elseif strcmpi(PropertyName,'par_texnames')
    Reply=obj.parameters.tex_name;
elseif strncmpi(PropertyName,'endo_list',9)
    Reply=load_endo_list();
elseif strncmpi(PropertyName,'exo_list',8)
    Reply=load_exo_list();
elseif strncmpi(PropertyName,'obs_list',8)
    Reply=load_obs_list();
elseif strcmpi(PropertyName,'chain_list')
    Reply=obj.markov_chains.chain_names;
elseif strcmpi(PropertyName,'regime_list')
    Reply=obj.markov_chains.regime_names;
elseif strcmpi(PropertyName,'state_list')
    Reply=obj.markov_chains.state_names;
else
    par_list=get(obj,'par_list');
    ploc=locate_variables(PropertyName,par_list,true);
    if ~isnan(ploc)
        Reply=obj.parameter_values(ploc,:);
    else
        error(['unknown gettable property ',PropertyName])
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
    function  Reply=load_par_list()
        Reply=obj.parameters.name;
        proplength=length('par_list');% =8
        if length(PropertyName)>proplength
            if ~strcmp(PropertyName(proplength+1),'(')
                error(['left parenthesis expected after par_list, but found ',PropertyName(proplength+1)])
            end
            if ~strcmp(PropertyName(end),')')
                error(['right parenthesis expected at the end of string, but found ',PropertyName(end)])
            end
            item=PropertyName(proplength+2:end-1);
            chain_loc=find(strcmpi(item,obj.markov_chains.chain_names));
            if ~isempty(chain_loc)
                mygovern=obj.parameters.governing_chain==chain_loc;
            else
                switch lower(item)
                    case 'switch'
                        mygovern=obj.parameters.is_switching;
                    case '~switch'
                        mygovern=~obj.parameters.is_switching;
                    case 'in_use'
                        mygovern=obj.parameters.is_in_use;
                    case '~in_use'
                        mygovern=~obj.parameters.is_in_use;
                    otherwise
                        error(['',item,''' is not recognized as the name of a markov chain'])
                end
            end
            Reply=Reply(mygovern);
        end
    end
    function Reply=load_obs_list()
        Reply=obj.observables.name;
        proplength=length('obs_list');% =8
        if length(PropertyName)>proplength
            if ~strcmp(PropertyName(proplength+1),'(')
                error(['left parenthesis expected after exo_list, but found ',PropertyName(proplength+1)])
            end
            if ~strcmp(PropertyName(end),')')
                error(['right parenthesis expected at the end of string, but found ',PropertyName(end)])
            end
            item=PropertyName(proplength+2:end-1);
            switch lower(item)
                case 'endo'
                    mygovern=obj.observables.is_endogenous;
                case '~endo'
                    mygovern=~obj.observables.is_endogenous;
                otherwise
                    error(['',item,''' is not recognized as a property of observable variables'])
            end
            Reply=Reply(mygovern);
        end
    end
    function  Reply=load_exo_list()
        Reply=obj.exogenous.name;
        proplength=length('exo_list');% =8
        if length(PropertyName)>proplength
            if ~strcmp(PropertyName(proplength+1),'(')
                error(['left parenthesis expected after exo_list, but found ',PropertyName(proplength+1)])
            end
            if ~strcmp(PropertyName(end),')')
                error(['right parenthesis expected at the end of string, but found ',PropertyName(end)])
            end
            item=PropertyName(proplength+2:end-1);
            switch lower(item)
                case 'obs'
                    mygovern=obj.exogenous.is_observed;
                case '~obs'
                    mygovern=~obj.exogenous.is_observed;
                case 'in_use'
                    mygovern=obj.exogenous.is_original;
                case '~in_use'
                    mygovern=~obj.exogenous.is_original;
                otherwise
                    error(['',item,''' is not recognized as a property of exogenous variables'])
            end
            Reply=Reply(mygovern);
        end
    end
    function Reply=load_endo_list()
        
        Reply=obj.endogenous.name;
        proplength=length('endo_list');% =9
        if length(PropertyName)>proplength
            if ~strcmp(PropertyName(proplength+1),'(')
                error(['left parenthesis expected after endo_list, but found ',PropertyName(proplength+1)])
            end
            if ~strcmp(PropertyName(end),')')
                error(['right parenthesis expected at the end of string, but found ',PropertyName(end)])
            end
            item=PropertyName(proplength+2:end-1);
            switch lower(item)
                case 'lagrange_multiplier'
                    mygovern=obj.endogenous.is_lagrange_multiplier;
                case '~lagrange_multiplier'
                    mygovern=~obj.endogenous.is_lagrange_multiplier;
                case 'static'
                    mygovern=obj.endogenous.is_static;
                case '~static'
                    mygovern=~obj.endogenous.is_static;
                case 'predetermined'
                    mygovern=obj.endogenous.is_predetermined;
                case '~predetermined'
                    mygovern=~obj.endogenous.is_predetermined;
                case 'pred_frwrd_looking'
                    mygovern=obj.endogenous.is_pred_frwrd_looking;
                case '~pred_frwrd_looking'
                    mygovern=~obj.endogenous.is_pred_frwrd_looking;
                case 'frwrd_looking'
                    mygovern=obj.endogenous.is_frwrd_looking;
                case '~frwrd_looking'
                    mygovern=~obj.endogenous.is_frwrd_looking;
                case 'original'
                    mygovern=obj.endogenous.is_original;
                case '~original'
                    mygovern=~obj.endogenous.is_original;
                case 'log_var'
                    mygovern=obj.endogenous.is_log_var;
                case '~log_var'
                    mygovern=~obj.endogenous.is_log_var; %
                case {'auxiliary'}
                    mygovern=obj.endogenous.is_auxiliary;
                case {'~auxiliary'}
                    mygovern=~obj.endogenous.is_auxiliary; %
                otherwise
                    error(['',item,''' is not recognized as a property of endogenous variables'])
            end
            Reply=Reply(mygovern);
        end
    end
end