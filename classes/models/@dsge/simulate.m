function varargout=simulate(obj,varargin)
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
    out=simulate@rise_generic(obj);
    [varargout{1:nargout}]=utils.miscellaneous.mergestructures(out,...
        struct('simul_sig',1,'simul_pruned',false,'simul_order',[]));
    return
end

% simul_order is solution order to use for simulation. Will most likely
% create problems with optimal policy if the user tries to simulate an optimal policy model
% with order of approximation higher than 1  

[varargout{1:nargout}]=simulate@rise_generic(obj,varargin{:});
end