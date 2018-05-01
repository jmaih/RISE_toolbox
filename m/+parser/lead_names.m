function out=lead_names(names,n)
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

if nargin<2
    n=1;
end
n=abs(n);

out=parser.concatenate_names_number(names,n);