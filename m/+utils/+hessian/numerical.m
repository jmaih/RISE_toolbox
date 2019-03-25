%  Computes the hessian numerically
% 
%  ::
% 
%    [H,issue]=numerical(fh,xbest,hessian_type)
% 
%  Args:
% 
%     - **fh** [char|function handle]: one-dimensional objective function
%     - **xbest** [vector]: point at which the hessian has to be computed
%     - **hessian_type** [{'fd'}|'opg']: type of hessian computed : finite
%       differences or outer-product-gradient
% 
%  Returns:
%     :
% 
%     - **H** [d x d matrix]: hessian
%     - **issue** [''|char]: description of any problem encountered during the
%       calculation of the hessian.
% 
%