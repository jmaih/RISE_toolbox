function hh=plot_fill(x,ydown,yup,color)
good=~isnan(ydown) & ~isnan(yup);
filled=[yup(good);flipud(ydown(good))];
x=x(good);
x=[x(:);flipud(x(:))];
hh=fill(x,filled,color);
