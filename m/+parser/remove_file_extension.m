function fname=remove_file_extension(fname)
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

thedot=find(fname=='.');
if ~isempty(thedot)
    fname=fname(1:thedot-1);
end
end