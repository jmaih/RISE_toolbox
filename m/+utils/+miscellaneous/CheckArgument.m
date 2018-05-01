function A=CheckArgument(A,expected_size,id)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if expected_size(end)==1 && numel(expected_size)>2
    expected_size=expected_size(1:end-1);
end
if isempty(A)
   A=zeros(expected_size);
else
    if ~isequal(size(A),expected_size)
        error([mfilename,':: wrong size of argument ',int2str(id)])
    end
end
