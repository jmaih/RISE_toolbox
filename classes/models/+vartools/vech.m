function v=vech(A)
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

[nr,nc]=size(A);
if ~isequal(nr,nc)
    error([mfilename,':: matrix should be square'])
end
v=nan(.5*nr*(nr+1),1);
iter=0;
for ii=1:nr
    v(iter+(1:nr-ii+1))=A(ii:end,ii);
    iter=iter+nr-ii+1;
end