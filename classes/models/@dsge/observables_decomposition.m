%--- help for dsge/observables_decomposition ---
%
%  Decomposes all variables of a DSGE model in terms of observables.
% 
%  ::
% 
%    weights=observables_decomposition(obj,select,xrange,db)
%    [weights,dec1,...,decn]=observables_decomposition(obj,select,xrange,db1,...,dbn)
% 
%  Args:
% 
%     obj (rise | dsge): scalar model object
%     select ({[]} | 'a' | 'att' | 'alpha' | 'v' | 'r' | 'epsilon' | 'eta'): type of
%       decomposition to perform
% 
%       - **a** : filter
%       - **att** : update
%       - **alpha** : smooth
%       - **v** : forecast errors
%       - **r** : variables calculated during the smoothing process
%       - **epsilon** : measurement errors
%       - **eta** : structural shocks
%       - **[]** : all the above
% 
%     xrange ({[]} | vector | cell array): start date and end date for the
%       decomposition
%     db (ts): database with the data to be used in the decomposition
% 
%  Returns:
%     :
% 
%     - **weights** [struct]: weights for the different elements requested from
%       the variable **select**
%     - **dec1** [struct]: hyper structure containing the time series for the
%       various decomposition types and variables
% 
%  Note:
% 
%     - if **select** is empty, all the decompositions are performed
%     - After doing the decomposition of different variables, one can take the
%       differences in the decomposition to perform, e.g., the decomposition of
%       forecast errors in terms of observables.
% 
%  See also:
%     - historical_decomposition
% 
%