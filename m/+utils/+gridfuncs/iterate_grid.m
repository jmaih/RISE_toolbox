function R=iterate_grid(R)
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

% see also mygrid

[rr,cc]=size(R);
R=mat2cell(R,ones(rr,1),cc);
for irow=1:rr
    item=R{irow};
    stretch=1:item(end);
    n=numel(stretch);
    R{irow}=[item(ones(n,1),:),stretch(:)];
end
R=cell2mat(R);
