function xrotate(angle)
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

if nargin==0
    angle=90;
end
tmp=gcf;
kids=get(tmp,'children');
if ischar(angle)
    angle=str2double(angle);
end
for ii=1:numel(kids)
    if isprop(kids(ii),'xtick')
        tmp=get(kids(ii),'xtick');
        if ~is_serial(tmp(1));
            continue
        end
        rotateXLabels(kids(ii),angle)
    end
end

