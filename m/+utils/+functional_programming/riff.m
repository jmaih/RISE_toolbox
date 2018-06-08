function out=riff(cond,f,varargin)
% INTERNAL FUNCTION
%

% variant of if_elseif where only the relevant function is evaluated.
% cond is a function whose evaluation returns a 1 x n vector of booleans
% f is a 1 x n cell array of functions

ping=find(cond(varargin{:}),1,'first');

out=f{ping}(varargin{:});

end