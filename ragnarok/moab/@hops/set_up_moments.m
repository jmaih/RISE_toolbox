%--- help for set_up_moments ---
%
%  mom_vector can be
%  - a function handle: it returns a vector containing the moments at each
%  order or a matrix, in which case the columns represent the orders and the
%  rows the different shocks
%  - a vector: the number of elements corresponds to the order of the
%  moments.
%  - a matrix, in which case the columns represent the orders and the
%  rows the different shocks 
%  mom_vector is a vector or produces a vector, it is assumed that all
%  shocks have the same distribution
%