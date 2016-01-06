function []=dim(start,finish,colorstr)
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
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

% function []=shade(start,finish,colorstr);
%
%  start and finish are Nx1 vectors of starting and ending years.
%  The function shades between the start and finish pairs using colorstr

if nargin<3 || isempty(colorstr)
    colorstr=[159 182 205]/256; 
end;  % default is yellow

fig=gcf;
kids=get(fig,'children');
for ii=1:numel(kids)
    ylim=get(kids(ii),'ylim');
    if isnumeric(ylim)
        y=[ylim(1) ylim(2) ylim(2) ylim(1)];
        axes(kids(ii))
        hold on;
        for i=1:length(start);
            x=[start(i) start(i) finish(i) finish(i)];
            fill(x,y,colorstr);
        end;
        % Now, prevent the shading from covering up the lines in the plot.
        h = findobj(kids(ii),'Type','line');
        % set(h,'EraseMode','xor');
        
        h = findobj(kids(ii),'Type','patch');
        set(h,'EdgeColor','none');
        
        % This last one makes the tick marks visible
        set(kids(ii), 'Layer', 'top')
    end
end
