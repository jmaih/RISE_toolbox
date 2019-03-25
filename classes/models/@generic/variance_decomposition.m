%--- help for generic/variance_decomposition ---
%
%  Computes the variance decomposition for a DSGE model
% 
%  ::
% 
%     [Vardec,obj]=variance_decomposition(obj)
% 
%     [Vardec,obj]=variance_decomposition(obj,varargin)
% 
%  Args:
% 
%     obj (dsge | rise ): model object
% 
%     varargin : Pairwise list of extra arguments
% 
%       - **vardec_shocks** [empty|char|cellstr]: list of shocks
% 
%       - **vardec_periods** [vector|{[1 4 8 16 40 100 200 400]}]: periods
%         to consider for the decomposition
% 
%       - **vardec_theoretical** [false|{true}]: if true, the theoretical
%         variance decomposition is computed to the extent that this is
%         possible. If false, the decomposition is based on simulated data.
% 
%       - **vardec_ergodic** [{false}|true]: if true, compute the ergodic
%         variance decomposition for a regime-switching model. If false, the
%         computations are regime specific
% 
%  Returns:
%     :
% 
%     - **Vardec** [struct]: structure containing the variance decomposition
%       both for the
% 
%        - infinite horizon (**infinity**)
%        - for each period (**conditional**)
% 
%
%    Other functions named variance_decomposition
%
%       abstvar/variance_decomposition
%