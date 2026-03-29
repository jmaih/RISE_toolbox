%--- help for dsge/initial_sstate ---
%
%  initial_sstate : provides a template for initial values in the
%  calculation of the steady state of a dsge model
% 
%  Syntax ::
% 
%    db=initial_sstate(obj)
% 
%  Args:
% 
%     obj (rise | dsge): scalar or vector of model objects.
% 
%  Returns:
% 
%     - **db** [struct | cell]: structure or cell array of structures if
%        several models are given as input. the fields of a structure are
%        the names of the endogenous variables of the model. The number of
%        concatenated structures is the number of regimes of the model. each
%        field is a 1 x 2 cell array. e.g bounds.C={sstateInfo,bgpInfo}
%        where sstateInfo and bgpInfo are 1 x 3 vectors organized as
%        [start_value,lower_bound,upper_bound].
% 
%  Note:
% 
%     - The information in each structure is given in original units/levels
%       of the variables before any potential log transformation.
% 
%  See also dsge.sstate
%