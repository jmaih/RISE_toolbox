function out = if_elseif(varargin)

out=varargin{2*find([varargin{1:2:end}],1,'first')};                                                               
end

