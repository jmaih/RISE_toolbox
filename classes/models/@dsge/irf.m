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
    out=irf@rise_generic(obj);
    out=utils.miscellaneous.mergestructures(out,...
        struct('irf_anticipate',true));%,'irf_risk',true
    return
end

obj.options.simul_anticipate_zero=~obj.options.irf_anticipate;

out=irf@rise_generic(obj,varargin{:});
end