function hout=plot_window(db,start,finish,plotfunc,varargin)

db=window(db,start,finish);

if isempty(plotfunc)
    plotfunc=@plot;
end
if ischar(plotfunc)
    tmp=plotfunc;
    plotfunc=str2func(plotfunc);
else
    tmp=func2str(plotfunc);
end
plotyy_flag=strcmp(tmp,'plotyy');

if plotyy_flag
    vnames=db.varnames;
    % db(vnames{1}) would not work since we are inside a method of the
    % class. This is why we have to explicitly use subsref
    this1=subsref(db,struct('type','()','subs',{vnames(1)})); 
    this2=subsref(db,struct('type','()','subs',{vnames(2)}));
    h=plotyy(this1,this2,varargin{:});
else
    h=plotfunc(db,varargin{:});
end
if nargout
    hout=h;
end