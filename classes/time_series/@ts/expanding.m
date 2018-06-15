function db=expanding(db,func,varargin)
% INTERNAL FUNCTION
%

window=[];

db=ts_roll_or_expand(db,func,window,varargin{:});

end