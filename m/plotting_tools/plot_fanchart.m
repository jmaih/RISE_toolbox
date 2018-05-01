function hh=plot_fanchart(data,MainColor,nticks)
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

if nargin<3
    nticks=[];
    if nargin<2
        MainColor=[];
    end
end
if isempty(MainColor)
    MainColor='nb';
end

Colors=InterpolateColors(data.ci,MainColor);

edge='w';%[1 1 1]
transparency=1; % 1 opaque:,....,0 (invisible)
for jj=1:numel(data.ci)
    ydata=data.quantiles(:,jj);
    zdata=data.quantiles(:,end-jj+1);
    handle=plot_fill(data.date_numbers,ydata,zdata,Colors(end-jj+1,:));
    set(handle,'EdgeColor',edge,'FaceAlpha',transparency,'EdgeAlpha',transparency,'linestyle','none');
    hold on
end
plot(data.date_numbers,data.mean,'--k','linewidth',2)
pp=plot_specs(data.date_numbers,nticks);

set(gca,'xlim',pp.xlim,'XTick',pp.tickLocs,'XtickLabel',pp.xtick_labels) %...
%     'ylim',[ymin,ymax])


grid on
if nargout
    hh=gca();
end

end

function hh=plot_fill(x,ydown,yup,color)
good=~isnan(ydown) & ~isnan(yup);
filled=[yup(good);flipud(ydown(good))];
x=x(good);
x=[x(:);flipud(x(:))];
hh=fill(x,filled,color);
end

function BaseColors=InterpolateColors(ci,maincolor)
ci=ci(:);
n=numel(ci);

if isnumeric(maincolor) && size(maincolor,2)==3
    
    if size(maincolor,1)==1
        maincolor=ci*maincolor;
    elseif size(maincolor,1)~=n
        error('number of elements in confidence region different from the number of rows of the color map')
    end
    BaseColors=maincolor;
elseif ischar(maincolor)
    switch maincolor
        case {'y',[1 1 0]} % yellow
            BaseColors=ci*[1 1 0];
        case {'m',[1 0 1]} % magenta
            BaseColors=ci*[1 0 1];
        case {'c',[0 1 1]} % cyan
            BaseColors=ci*[0 1 1];
        case {'r',[1 0 0]} % red
            BaseColors=ci*[1 0 0];
        case {'g',[0 1 0]} %green
            BaseColors=ci*[0 1 0];
        case {'b',[0 0 1]} %blue
            BaseColors=ci*[0 0 1];
        case {'w',[1 1 1]} %white
            BaseColors=ci*[1 1 1];
        case {'nb'}
            BaseColors=[
                1 123 182
                65 156 200
                128 189 219
                192 222 237]/255;
            if n~=4
                error('number of elements in confidence region expected to be 4')
            end
        otherwise
            error('wrong color')
    end
else
            error('wrong color')
end

end