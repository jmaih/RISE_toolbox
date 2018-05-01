function c=minus(a,b)
% minus -- overload minus for cell
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

if ischar(a)
    a=cellstr(a);
end

if ischar(b)
    b=cellstr(b);
end

[c,ia]=setdiff(a,b);

[~,isort]=sort(ia);

c=c(isort);

end
