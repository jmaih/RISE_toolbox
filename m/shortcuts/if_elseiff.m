function varargout=if_elseiff(varargin)

% abs(if_elseiff({1},true,@cos,false,@sin)-0.54030230586814)<1e-15
% abs(if_elseiff({1},false,@cos,true,@sin)-0.841470984807897)<1e-15

[varargout{1:nargout}]=utils.functional_programming.if_elseiff(varargin{:});

end

