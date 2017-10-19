function obj=process_thresholds(obj,thresh)

% we need to know:
% - the name of the threshold variable
% - its lags
% - the endogenous variables affected by the threshold.
% - If the list of endogenous variables affected is empty, we
% assume the threshold affects all endogenous variables in the
% system
% e.g. {'loansRates{-3}',{'growth','inflation','interest rate'})

pattern='(?<name>\w+)(\{|\()?(?<lag>[^})]+)?(\)|\})?';

for ithresh=1:numel(thresh.items)
    % 1-transition variable
    %---------------------
    tv0=thresh.items(ithresh);
    
    descript=regexprep(tv0.transition_description,'(\(|\{)(-)?\d+(\}|\))','');
        
    tv=regexp(tv0.transition_variable,pattern,'names');
    
    tv.tex=strrep(descript,'"','');
    % 2-type of transition
    %---------------------
    nc=size(tv0.threshold_priors,1);
    
    tv=process_transition(tv0.type,tv,nc);
    % 3-variables influenced by the transition
    %-----------------------------------------
    % in order to have different parameters estimated for a
    % transition, just declare several transitions with
    % different variables they control...
    controlled_vars=tv0.controlled;
    
    if isempty(controlled_vars)
        
        controlled_vars=obj.endogenous.name;
        
    end
        
    if ~all(ismember(controlled_vars,obj.endogenous.name))
        
        error('all controlled variables must be endogenous')
        
    end
    
    tv.controlled_vars=controlled_vars;
    
    tv.threshold_priors=tv0.threshold_priors;
    
    if isempty(tv.lag)
        
        tv.lag='0';
        
    end
    
    if ithresh==1
        
        obj.thresholds=tv;
        
    else
        
        obj.thresholds(ithresh)=tv;
        
    end
    
end

end

function tv=process_transition(transfun,tv,nc)

np=1+nc;
    
tv.np=np;

tv.func=str2func(transfun);

end