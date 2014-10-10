function x=create_state_list(name,number)
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

x=cellstr(strcat(name,'_',num2str((1:number)')));
x=cellfun(@(z)z(~isspace(z)),x,'uniformOutput',false);
x=x(:)';
end
