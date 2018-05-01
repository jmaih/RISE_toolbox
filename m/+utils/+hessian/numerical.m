function [H,issue]=numerical(fh,xbest,hessian_type)
% numerical -- compute the hessian numerically
%
% ::
%
%
%   [H,issue]=numerical(fh,xbest,hessian_type)
%
% Args:
%
%    - **fh** [char|function handle]: one-dimensional objective function
%
%    - **xbest** [vector]: point at which the hessian has to be computed
%
%    - **hessian_type** [{'fd'}|'opg']: type of hessian computed : finite
%    differences or outer-product-gradient
%
% Returns:
%    :
%
%    - **H** [d x d matrix]: hessian
%
%    - **issue** [''|char]: description of any problem encountered during the
%    calculation of the hessian.
%
% Note:
%
% Example:
%
%    See also:

if nargin==0
       
    mydefaults=the_defaults();
    
    if nargout
        
        H=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return

end

if ischar(fh)
    
    fh=str2func(fh);
    
end

issue='';

switch lower(hessian_type)
    
    case 'fd'
        
        H = utils.hessian.finite_differences(fh,xbest);
        
    case 'opg'
        
        H = utils.hessian.outer_product(fh,xbest);
        
        if any(isnan(H(:)))
            
            issue='OPG unstable and inaccurate for calculation of Hessian, switched to finite differences';
            
            warning([mfilename,':: ',issue]) %#ok<WNTAG>
            
            warning([mfilename,':: OPG unstable for calculation of Hessian, switching to finite differences']) %#ok<WNTAG>
            
            H= finite_difference_hessian(fh,xbest);
        
        end
        
    otherwise
        
        issue=['unknow hessian option ',hessian_type,' using finite differences'];
        
        warning([mfilename,':: ',issue]) %#ok<WNTAG>
        
        H = utils.hessian.finite_differences(fh,xbest);

end

end

function d=the_defaults()

d={
    
'hessian_type','fd',@(x)ismember(x,{'fd','opg'}),...
' hessian_type must be ''fd'' or ''opg'''

};

end

    

