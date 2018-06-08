function out = if_elseif(varargin)
% INTERNAL FUNCTION
%


out=varargin{2*find([varargin{1:2:end}],1,'first')};
end

