function out=irf(obj,varargin)
if isempty(obj)
    out=irf@rise_generic(obj);
    out=utils.miscellaneous.mergestructures(out,...
        struct('irf_anticipate',true));%,'irf_risk',true
    return
end

out=irf@rise_generic(obj,varargin{:});
end