function c=compare_individuals(a,b)
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

c=1;
if (b.violstrength<a.violstrength)||...
        ((b.violstrength==a.violstrength) && b.f<a.f)
    c=2;
end
end