function [a_loc,eqtn_loc,var_loc]=decompose_parameter(x,headers)
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

underscores=find(x=='_');
aa=x(1:underscores(1)-1);
a_loc=[];
if nargin>1
    a_loc= strcmp(aa,headers);
end
eqtn_loc=str2double(x(underscores(1)+1:underscores(2)-1));
var_loc=str2double(x(underscores(2)+1:end));