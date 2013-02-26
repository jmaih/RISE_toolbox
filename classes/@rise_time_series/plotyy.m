function [plot_handle,axis_handle]=plotyy(this1,this2,varargin)
this=[this1,this2];
datta=double(this);
if size(datta,2)>1 && size(datta,3)>1
    error([mfilename,':: cannot handle many variables and many pages simultaneously'])
else
    datta=squeeze(datta);
end

pp=plot_specs(this.TimeInfo,5);

% tmp=plotyy(pp.xdatenums,datta(:,1),pp.xdatenums,datta(:,2),varargin{:});

high=max(datta,[],1);
low=min(datta,[],1);
datta=bsxfun(@rdivide,bsxfun(@minus,datta,low),high-low);
tmp=plot(pp.xdatenums,datta,varargin{:});

set(gca,'xlim',pp.xlim,'XTick',pp.tickLocs,'XtickLabel',pp.xtick_labels)

grid on

if nargout>0
    plot_handle=tmp;
    if nargout>1
        axis_handle=get(gca);
    end
end
end
