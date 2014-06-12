function [D01,retcode]=compute_automatic_derivatives(derivatives,solve_order,varargin)
D01=utils.code.evaluate_automatic_derivatives(derivatives,solve_order,varargin{:});
good=all(cellfun(@(x)utils.error.valid(x),D01));
retcode=0;
if ~good
    retcode=2; % nans in jacobian
end
end
