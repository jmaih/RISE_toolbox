function [flag,msg]=is_forbidden_name(name)
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


fn={
    'regime','used for setting up simulations '
    };

start=strcmp(name,fn(:,1));
flag=any(start);
msg='';
if flag
    msg=sprintf('%s cannot be a declaration coz %s',name,fn{start,2});
end

end