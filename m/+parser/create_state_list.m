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
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if isscalar(number)
    
    number=(1:number).';
    
else
    
    number=number(:);
    
end

x=cellstr(strcat(name,'_',num2str(number)));

x=cellfun(@(z)z(~isspace(z)),x,'uniformOutput',false);

x=x(:).';

end
