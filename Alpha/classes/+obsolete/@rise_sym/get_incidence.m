function occur=get_incidence(varargin)
occur=[];
for iv=1:length(varargin)
    if isa(varargin{iv},'rise_sym') && ~isempty(varargin{iv}.incidence)
        if isempty(occur)
            occur=varargin{iv}.incidence;
        else
            occur=occur|varargin{iv}.incidence;
        end
    end
end
end