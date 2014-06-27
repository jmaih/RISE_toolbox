function [plot_handle,axis_handle]=plotyy(this1,this2,varargin)
this=[this1,this2];
datta=double(this);
if size(datta,2)>1 && size(datta,3)>1
    error([mfilename,':: cannot handle many variables and many pages simultaneously'])
else
    datta=squeeze(datta);
end

pp=plot_specs(this.date_number);

test=true;
if test
    [ax,h1,h2]=plotyy(pp.xdatenums,datta(:,1),pp.xdatenums,datta(:,2));
    if ~isempty(varargin)
        set([h1,h2],varargin{:});
    end
    set(ax,'xlim',pp.xlim,'XTick',pp.tickLocs,'XtickLabel',pp.xtick_labels)
else
    high=max(datta,[],1);
    low=min(datta,[],1);
    datta=bsxfun(@rdivide,bsxfun(@minus,datta,low),high-low);
    ax=plot(pp.xdatenums,datta,varargin{:});
    set(gca,'xlim',pp.xlim,'ylim',[min(datta(:)),max(datta(:))],'XTick',pp.tickLocs,'XtickLabel',pp.xtick_labels)
    axis tight
end

grid on

% linkaxes(ax,'x')

% rotateXLabels(gca,90);

if nargout>0
    plot_handle=ax;
    if nargout>1
        axis_handle=get(gca);
    end
end
end
