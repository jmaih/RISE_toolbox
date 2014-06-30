function this=horzcat(varargin)

for ii=1:length(varargin)
    if isempty(varargin{ii})
        if ii==1
            error('the first argument cannot be empty')
        end
    else
        if ~isequal(class(varargin{ii}),'ts')
            error([mfilename,':: input ',int2str(ii),' must be from class ts'])
        end
        if ii==1
            this=varargin{ii};
        else
            this=this&varargin{ii};
        end
    end
end
end
