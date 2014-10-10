function xrotate(angle)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

if nargin==0
    angle=90;
end
tmp=gcf;
kids=get(tmp,'children');
if ischar(angle)
    angle=str2double(angle);
end
for ii=1:numel(kids)
    tmp=get(kids(ii),'xtick');
    if ~is_serial(tmp(1));
        continue
    end
    rotateXLabels(kids(ii),angle)
end

