function n=numel(this,varargin)
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


if isempty(varargin)
    n=builtin('numel',this);
else
    n=1;
end

end