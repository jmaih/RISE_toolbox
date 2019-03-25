%  INTERNAL FUNCTION: computes A*OMEGA_i
% 
%  ::
% 
%    oo = dv_vz_omega(dv_vz,nz,code)
%    oo = dv_vz_omega(dv_vz,nz,P1,P2,...,Pm)
% 
%  Args:
% 
%     - **dv_vz** [matrix]: matrix for which the sum of permutations has to be
%       calculated
%     - **nz** [integer]: the number of columns of **dv_vz** has to be a
%       multiple of **nz**
%     - **code** [1|2|3|4|5|6|7|8|9]: pre-specified codes for permutations
%     - **P1,...,Pm** [vectors]: user-defined permutations
% 
%  Returns:
%     :
% 
%     - **oo** [matrix]: sum of permutations of **dv_vz**
% 
%