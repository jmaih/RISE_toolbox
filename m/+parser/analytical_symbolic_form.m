function [strout,varlist]=analytical_symbolic_form(str,validNames,flag)
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

% turns x(1),y(1),ss(1),p(1) into x_1,y_1,ss_1,p_1
% this might not be very robust, I should add that the variables cannot
% be preceeded by another alphanumeric character or followed by a
% alphanumeric character.
varlist=[];
if ischar(validNames)
	validNames=cellstr(validNames);
end
validNames=validNames(:)';
var_types=cell2mat(strcat(validNames,'|'));
var_types=var_types(1:end-1);
NoANb='(?<![\w])';% no alphanumeric before
NoAf='(?![a-zA-Z])';% no alphabetic after
switch flag
    case 'symbolic'
        strout=regexprep(str,[NoANb,'(',var_types,')(\s*)\((\s*)(\d*)(\s*)\)'],'$1_$4');
        if nargout>1
            varlist = regexp(strout,[NoANb,'(',var_types,')_\d*',NoAf],'match');
            varlist=unique(varlist);
        end
    case {'analytic','analytical'}
        strout=regexprep(str,[NoANb,'(',var_types,')(\_)(\d*)'],'$1($3)');
    otherwise
        error([mfilename,':: unknown flag ',flag])
end


