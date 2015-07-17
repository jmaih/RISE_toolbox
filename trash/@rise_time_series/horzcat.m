function this=horzcat(varargin)
if length(varargin)<2
    error([mfilename,':: number of arguments should be at least 2'])
end
for ii=1:length(varargin)
    if ~isequal(class(varargin{ii}),'rise_time_series')
        error([mfilename,':: input ',int2str(ii),' must be from class rise_time_series'])
    end
    if ii==1
        this=varargin{ii};
    else
        this=this&varargin{ii};
    end
end
end
