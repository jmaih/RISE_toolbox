function []=dim(start,finish,colorstr,fig)
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
%  start and finish are N x 1 vectors of starting and ending years.
%  The function shades between the start and finish pairs using colorstr

if nargin<4
    
    fig=[];
    
    if nargin<3 || isempty(colorstr)
        
        colorstr=[159 182 205]/256;
        
    end  % default is yellow
    
end

if isempty(fig)
    
    fig=gcf;
    
end

nfig=numel(fig);

if nfig > 1
    
    for ifig=1:nfig
        
        dim(start,finish,colorstr,fig(ifig))
        
    end
    
    return
    
end

kids=get(fig,'children');

for ii=1:numel(kids)
    
    ylim=get(kids(ii),'ylim');
    
    if isnumeric(ylim)
        
        y=[ylim(1) ylim(2) ylim(2) ylim(1)];
        
        axes(kids(ii))
        
        hold on
        
        for i=1:length(start)
            
            x=[start(i) start(i) finish(i) finish(i)];
            
            fill(x,y,colorstr);
            
        end
        % Now, prevent the shading from covering up the lines in the plot.
%         h = findobj(kids(ii),'Type','line');
        % set(h,'EraseMode','xor');
        
        h = findobj(kids(ii),'Type','patch');
        set(h,'EdgeColor','none');
        
        % This last one makes the tick marks visible
        set(kids(ii), 'Layer', 'top')
        
    end
    
end
