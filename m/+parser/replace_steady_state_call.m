function string=replace_steady_state_call(string,flag)
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
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

if nargin<2
    flag='';
end

pattern='(?<!\w+)steady_state\((\d+)\)';
switch flag
    case 'symbolic'
        new_string='ss_$1';
    otherwise
        new_string='ss($1)';
end
string= regexprep(string,pattern,new_string);
    
end
