%--- help for rdico.revectorize_function ---
%
%  vectorizes functions
% 
%    - In the presence of incidence (dynamic model), the variables are
%    vectorized following the indices present in the incidence matrix
% 
%    - Without incidence (static model), the variables are
%    vectorized in the order of their appearance. This order is recorded in
%    the second entry of the output with the third entry of the output being
%    the new variable name
%