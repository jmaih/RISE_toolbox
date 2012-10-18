function [strout,varlist]=analytical_symbolic_form(str,validNames,flag)
% turns x(1),y(1),ss(1),p(1) into x_1,y_1,ss_1,p_1 
varlist=[];
if ischar(validNames)
	validNames=cellstr(validNames);
end
var_types=validNames{1};
for ivar=2:numel(validNames)
    var_types=[var_types,'|',validNames{ivar}]; %#ok<AGROW>
end
switch flag
    case 'symbolic'
        strout=regexprep(str,['(',var_types,')\((\d*)\)'],'$1_$2');
        if nargout>1
            varlist = regexp(strout,['(',var_types,')_\d*'],'match');
            varlist=unique(varlist);
        end
    case {'analytic','analytical'}
        strout=regexprep(str,['(',var_types,')(\_)(\d*)'],'$1($3)');
    otherwise
        error([mfilename,':: unknown flag ',flag])
end


