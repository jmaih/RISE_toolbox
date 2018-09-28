function []=dim(start,finish,colorstr,fig)
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

x=[start;start;finish;finish];

nx=size(x,2);

for ii=1:numel(kids)
    
    ylim=get(kids(ii),'ylim');
    
    if isnumeric(ylim)
        
        y=[ylim(1) ylim(2) ylim(2) ylim(1)].';
        
        y=y(:,ones(1,nx));
        
        axes(kids(ii))
        
        hold on
        
        mypatch=patch(x,y,colorstr,'EdgeColor','none');
        
        grand_kids=get(kids(ii),'children');
        
        grand_kids(grand_kids==mypatch)=[];
        
        grand_kids(end+1)=mypatch; %#ok<AGROW>
        
        set(kids(ii),'children',grand_kids);
                
        % This last one makes the tick marks visible
        set(kids(ii), 'Layer', 'top')
        
    end
    
end
