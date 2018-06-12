function obj=prototypize(f)
% INTERNAL FUNCTION
%

persistent proto

if isempty(proto)

	proto=splanar(0);

end

obj=proto;

if nargin

	obj.func=f;

end

end
