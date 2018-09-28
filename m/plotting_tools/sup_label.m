function [ax,h]=sup_label(text,whichLabel,supAxes)
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

% This modifies Ben Barrowes' suplabel (see Matlab Central).
% It places text as a title, xlabel, or ylabel on a group of subplots.
% Returns a handle to the label and a handle to the axis.
%  [ax,h]=sup_label(text,whichLabel,supAxes)
% returns handles to both the axis and the label.
%  ax=sup_label(text,whichLabel,supAxes)
% returns a handle to the axis only.
%  sup_label(text) with one input argument assumes whichLabel='x'
%
% whichLabel is any of 'x', 'y', 'yy', or 't', specifying whether the
% text is to be the xlabel, ylabel, right side y-label, or title
% respectively.
%
% supAxes is an optional argument specifying the Position of the
%  "super" axes surrounding the subplots.
%  supAxes defaults to [.08 .08 .84 .84]
%  specify supAxes if labels get chopped or overlay subplots
%
% EXAMPLE:
%  subplot(2,2,1);ylabel('ylabel1');title('title1')
%  subplot(2,2,2);ylabel('ylabel2');title('title2')
%  subplot(2,2,3);ylabel('ylabel3');xlabel('xlabel3')
%  subplot(2,2,4);ylabel('ylabel4');xlabel('xlabel4')
%  [ax1,h1]=sup_label('super X label');
%  [ax2,h2]=sup_label('super Y label','y');
%  [ax3,h2]=sup_label('super Y label (right)','yy');
%  [ax4,h3]=sup_label('super Title'  ,'t');
%  set(h3,'FontSize',30)
%
% SEE ALSO: text, title, xlabel, ylabel, zlabel, subplot,
%           suptitle (Matlab Central)

currax=findobj(gcf,'type','axes','-not','tag','sup_label');
currax=[currax;findobj(gcf,'tag','legend')];
if nargin < 3
    supAxes=[.08 .08 .84 .84];
    ah=findall(gcf,'type','axes');
    if ~isempty(ah)
        leftMin=inf;  bottomMin=inf;  leftMax=0;  bottomMax=0;
        axBuf=.04;
        set(ah,'units','normalized')
        ah=findall(gcf,'type','axes');
        for ii=1:length(ah)
            if strcmp(get(ah(ii),'Visible'),'on')
                thisPos=get(ah(ii),'Position');
                leftMin=min(leftMin,thisPos(1));
                bottomMin=min(bottomMin,thisPos(2));
                leftMax=max(leftMax,thisPos(1)+thisPos(3));
                bottomMax=max(bottomMax,thisPos(2)+thisPos(4));
            end
        end
        supAxes=[leftMin-axBuf,bottomMin-axBuf,...
            leftMax-leftMin+axBuf*2,bottomMax-bottomMin+axBuf*2];
    end
    if nargin < 2
        whichLabel = 'x';
        if nargin < 1
            help(mfilename);
            return;
        end
    end
end

if ~ischar(text) || ~ischar(whichLabel)
    error('text and whichLabel must be strings')
end
whichLabel=lower(whichLabel);

ax=axes('Units','Normal','Position',supAxes,'Visible','off','tag','sup_label');
if strcmp('t',whichLabel)
    set(get(ax,'Title'),'Visible','on')
    title(text,'interpreter','none');
elseif strcmp('x',whichLabel)
    set(get(ax,'XLabel'),'Visible','on')
    xlabel(text,'interpreter','none');
elseif strcmp('y',whichLabel)
    set(get(ax,'YLabel'),'Visible','on')
    ylabel(text,'interpreter','none');
elseif strcmp('yy',whichLabel)
    set(get(ax,'YLabel'),'Visible','on')
    ylabel(text,'interpreter','none');
    set(ax,'YAxisLocation','right')
end

% % restore all other axes
% for k=1:length(currax)
%     try
%         axes(currax(k)); %#ok<LAXES>
%     end
% end

if (nargout < 2)
    return
end
if strcmp('t',whichLabel)
    h=get(ax,'Title');
    set(h,'VerticalAlignment','middle')
elseif strcmp('x',whichLabel)
    h=get(ax,'XLabel');
elseif strcmp('y',whichLabel) || strcmp('yy',whichLabel)
    h=get(ax,'YLabel');
end