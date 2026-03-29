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
%     - **obj** (dsge | rise ): model object
% 
%     - **varargin** : Pairwise list of extra arguments
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
% 
%     - **Vardec** [struct]: structure containing the variance decomposition
%       both for the
% 
%        - infinite horizon (**infinity**)
%        - for each period (**conditional**)
% 
%     - **retcode** [numeric]: = 0 if there is no problem
%  
%  N.B: The calculation of the variance is done by solving a Lyapunov
%  equation. There are different Lyapunov algorithms which can be set
%  through the option "lyapunov_algo". The available algorithms are:
% 
%    - 'doubling'  OR @lyap_solvers.doubling (the default)
%    - 'bicgstab'  OR @lyap_solvers.bicgstab
%    - 'bicgstabl'  OR @lyap_solvers.bicgstabl
%    - 'cgs'  OR @lyap_solvers.cgs
%    - 'fevdcov'  OR @lyap_solvers.fevdcov
%    - 'direct'  OR @lyap_solvers.direct
%    - 'fix_point'  OR @lyap_solvers.fix_point
%    - 'robust'  OR @lyap_solvers.robust
%    - 'sandwich'  OR @lyap_solvers.sandwich
%    - 'schur'  OR @lyap_solvers.schur
%    - 'bartels_stewart'  OR @lyap_solvers.bartels_stewart
% 
%  Notes: 
% 
%  - In the presence of a nonstationary model, it is not a good idea to
%    use the doubling algorithm : it won't converge. A different algorithm
%    such as "bartels_stewart" may work better. 
% 
%  - The algorithms listed above are invoked either as char/string or
%    directly as function handles prefixed with "lyap_solvers.". e.g.
%    'fevdcov' or @lyap_solvers.fevdcov
% 
%  - In the presence of a nonstationary model, it is not a good idea to
%    use the doubling algorithm : it won't converge. A different algorithm
%    such as "bartels_stewart" may work better.
%
%    Other uses of variance_decomposition
%
%       abstvar/variance_decomposition
%