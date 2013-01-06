function h=plot_window(db,start,finish,plotfunc,varargin)

db=window(db,start,finish);

if ischar(plotfunc)
    plotfunc=str2func(plotfunc);
end

h=plotfunc(db,varargin{:});