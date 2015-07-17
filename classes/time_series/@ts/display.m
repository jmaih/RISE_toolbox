function display(self)
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

if numel(self)==1
    if self.NumberOfVariables
        disp(main_frame(self)) % does not display the name of the variable
    end
    index(self);
else
    disp(self)
end
end
