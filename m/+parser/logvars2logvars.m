function [dictionary,blocks,old_endo_names]=logvars2logvars(dictionary,blocks)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

logvarnames={dictionary.log_vars.name};
% old endo names will be useful if there is a problem in the model and has
% to be constructed prior to changing the names of the variables.
old_endo_names={dictionary.endogenous.name};

if ~isempty(logvarnames)
    if dictionary.is_linear_model
        error([mfilename,':: with ''log_vars'', the model cannot be declared as ''linear'''])
    end
    target_blocks={'model','steady_state_model','planner_objective'};
    for itarg=1:numel(target_blocks)
        target_blocks{itarg}=find(strcmp(target_blocks{itarg},{blocks.name}));
    end
    
    [exprLog,replaceLog,convLog]=regexp_logsetup(); %#ok<ASGLU>
    [expr,replace,convcoef]=regexp_setup(); %#ok<ASGLU>
    
    for iblk=1:numel(target_blocks)
        blocks(target_blocks{iblk}).listing(:,2)=regexprep(blocks(target_blocks{iblk}).listing(:,2),...
            exprLog,replaceLog);
        blocks(target_blocks{iblk}).listing(:,2)=regexprep(blocks(target_blocks{iblk}).listing(:,2),...
            expr,replace);
    end
    
    for ivar=1:numel(logvarnames)
        loc=strcmp(logvarnames{ivar},{dictionary.endogenous.name});
        dictionary.endogenous(loc).name=['LOG_',logvarnames{ivar}];
        dictionary.endogenous(loc).is_log_var=true;
    end
    % re-write the steady state model if we have equations such as exp(LOG_)=
    ssblock=blocks(target_blocks{2}).listing(:,2);
    for irow=1:size(ssblock,1)
        ssblock{irow}(isspace(ssblock{irow}))=[];
        equal_loc=strfind(ssblock{irow},'=');
        if isempty(equal_loc)||~strncmp(ssblock{irow},'exp(',4)
            continue
        end
        % remove the semi-colon at the end
        if strcmp(ssblock{irow}(end),';')
            ssblock{irow}=ssblock{irow}(1:end-1);
        end
        % find the matching closing parenthesis, which should be precisely
        % the first one on target, right before the equality sign
        right_par=equal_loc-1;
        if ~strcmp(ssblock{irow}(right_par),')')
            error(['unidentified problem in equation ',ssblock{irow},' please contact junior.maih@gmail.com'])
        end
        % now invert the equation
        ssblock{irow}=[ssblock{irow}(5:right_par-1),'=log(',ssblock{irow}(equal_loc+1:end),');'];
    end
    blocks(target_blocks{2}).listing(:,2)=ssblock;
end


function [exprLog,replaceLog,convLog]=regexp_logsetup()

nclogl='(?:log\()';
nclogr='(?:\))';
sstatel='(steady_state\()?';
sstater='(\))?';
ncl='(?<!\w+)';
vlist=parser.cell2matize(logvarnames);
left='(\(|\{)?';
pm='(\+|\-)?';
digits='(\d+)?';
right='(\)|\})?';
ncr='(?!\w+)';
% let nc = no capture
replaceLog='${convLog($1,$2,$3,$4,$5,$6,$7)}';
exprLog=[nclogl,sstatel,ncl,vlist,left,pm,digits,right,ncr,sstater,nclogr];
convLog=@converter;
    function out=converter(sstatel,vname,left,pm,digits,right,sstater)
        out=[sstatel,'LOG_',vname,left,pm,digits,right,sstater];
    end
end


function [expr,replace,convcoef]=regexp_setup()

sstatel='(steady_state\()?';
sstater='(\))?';
ncl='(?<!\w+)';
vlist=parser.cell2matize(logvarnames);
left='(\(|\{)?';
pm='(\+|\-)?';
digits='(\d+)?';
right='(\)|\})?';
ncr='(?!\w+)';
% let nc = no capture
replace='${convcoef($1,$2,$3,$4,$5,$6,$7)}';
expr=[sstatel,ncl,vlist,left,pm,digits,right,ncr,sstater];
convcoef=@coef_converter;
    function out=coef_converter(sstatel,vname,left,pm,digits,right,sstater)
        out=['exp(',sstatel,'LOG_',vname,left,pm,digits,right,sstater,')'];
    end
end

end
