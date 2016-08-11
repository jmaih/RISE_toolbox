function out=irf(obj,varargin)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if isempty(obj)
    
    out=irf@generic_switch(obj);
    
    out=utils.miscellaneous.mergestructures(out,...
        struct('irf_anticipate',true));%,'irf_risk',true
    
    return
    
end

obj=set(obj,varargin{:});

for iobj=1:numel(obj)
    
    obj(iobj).options.simul_anticipate_zero=~obj(iobj).options.irf_anticipate;
    
end

out=irf@generic_switch(obj);

end