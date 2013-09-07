function c=if_elseif(varargin)

% check whether all the 'second' elements are zero
zero=true;
for ii=2:2:length(varargin)
    arg=varargin{ii};
    zero=zero && is_zero(arg);
    if ~zero
        break
    end
end

if zero
    c=planar(0);
else
    c=planar.multinary_operation('if_elseif',varargin{:});
end

end
