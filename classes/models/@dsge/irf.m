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
% Description
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

out=irf@rise_generic(obj,varargin{:});
end