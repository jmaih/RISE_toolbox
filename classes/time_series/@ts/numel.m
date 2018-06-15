function n=numel(this,varargin)
% INTERNAL FUNCTION
%

if isempty(varargin)
    n=builtin('numel',this);
else
    n=1;
end

end