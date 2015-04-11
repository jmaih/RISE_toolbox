function [H,issue]=numerical(fh,xbest,hessian_type)
% numerical -- compute the hessian numerically
%
% Syntax
% -------
% ::
%
%   [H,issue]=numerical(fh,xbest,hessian_type)
%
% Inputs
% -------
%
% - **fh** [char|function handle]: one-dimensional objective function
%
% - **xbest** [vector]: point at which the hessian has to be computed
%
% - **hessian_type** [{'fd'}|'opg']: type of hessian computed : finite
% differences or outer-product-gradient
%
% Outputs
% --------
%
% - **H** [d x d matrix]: hessian
%
% - **issue** [''|char]: description of any problem encountered during the
% calculation of the hessian.
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if nargin==0
    if nargout>1
        error([mfilename,':: when there are no inputs, nargout must be at most 1'])
    end
    H=struct('hessian_type','fd');
    return
end

if ischar(fh)
    fh=str2func(fh);
end

issue='';
switch lower(hessian_type)
    case 'fd'
        H(:,:,2) = utils.hessian.finite_differences(fh,xbest);
    case 'opg'
        H(:,:,2) = utils.hessian.outer_product(fh,xbest);
        if any(any(isnan(H(:,:,2))))||any(any(isinf(H(:,:,2))))
            issue='OPG unstable and inaccurate for calculation of Hessian, switched to finite differences';
            warning([mfilename,':: ',issue]) %#ok<WNTAG>
            warning([mfilename,':: OPG unstable for calculation of Hessian, switching to finite differences']) %#ok<WNTAG>
            H(:,:,2) = finite_difference_hessian(fh,xbest);
        end
    otherwise
        issue=['unknow hessian option ',hessian_type,' using finite differences'];
        warning([mfilename,':: ',issue]) %#ok<WNTAG>
        H(:,:,2) = utils.hessian.finite_differences(fh,xbest);
end

end
