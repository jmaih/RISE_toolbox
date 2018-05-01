function parfor_save(filename,x) %#ok<INUSD>
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

save(filename,'x')

% function parfor_save(filename,varargin)
% save(filename,varargin{:})
%{
parfor ii=1:4
x=rand(10,10);
utils.parallel.parfor_save(sprintf('output%d.mat',ii),x)
end
%}