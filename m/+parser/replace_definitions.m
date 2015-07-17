function batch=replace_definitions(batch,Definitions)
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


if ~isempty(Definitions)
    Definitions=regexp(Definitions,'=','split');
    Definitions=vertcat(Definitions{:});
    Definitions=strrep(Definitions,';','');
    for id=size(Definitions,1):-1:1
        % the parentheses here are not very efficient
        batch=strrep(batch,Definitions{id,1},['(',Definitions{id,2},')']);
    end
end
