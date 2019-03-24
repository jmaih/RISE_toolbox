%  Compute the ABCD test statistics of Fernandez-Villaverde, Rubio-Ramirez, Sargent, and Watston (2007)
%  Uses the state-space system
% 
%  .. math:
% 
%     x(+1)=Ax+Bw(+1)
%     y(+1)=Cx+Dw(+1)
% 
%  and computes the eigenvalues of
% 
%  .. math:
% 
%     A-BC^(-1)D  (Condition 1)
% 
%  If all eigenvalues are smaller than 1, the poor man's invertibility
%  condition is satisfied and the structural shocks can be recovered from
%  the observables
% 
%  ::
% 
%     [test, A, B, C, D] = abcd_rise(m);
% 
%  Args:
%     - m (model object): model object
% 
%  Returns:
%     :
% 
%        - test : check above
%        - A : check above
%        - B : check above
%        - C : check above
%        - D : check above
% 
%  Reference:
%     Fernandez-Villaverde, Rubio-Ramirez,Sargent,
%     and Watson (2007), "ABCs (and Ds) of Understanding VARs", American
%     Economic Review, 97(3), 1021-1026
% 
%