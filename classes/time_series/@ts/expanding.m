function db=expanding(db,func,varargin)

window=[];

db=ts_roll_or_expand(db,func,window,varargin{:});

end