%--- help for ts/chowlin ---
%
%  Temporal disaggregation using the Chow-Lin method
% 
%  ::
% 
%    [yh,res]=chowlin(y0,Xh0)
%    [yh,res]=chowlin(y0,Xh0,aggreg_type)
%    [yh,res]=chowlin(y0,Xh0,aggreg_type,estim_method)
%    [yh,res]=chowlin(y0,Xh0,aggreg_type,estim_method,ngrid)
% 
%  Args:
% 
%     y0 (ts object): low-frequency left-hand-side variable
% 
%     Xh0 (ts | struct): high-frequency right-hand-side explanatory variables
% 
%     aggreg_type ('flow' | {'average'} | 'index' | 'last' | 'first'): type of
%       aggregation:
% 
%       - 'flow' (or 1) is the sum,
%       - 'average','index' (or 2) is the average,
%       - 'last' (or 3) is the last element,
%       - 'first' (or 4) is the first element,
% 
%     estim_method ({0} | 1 | 2) : estimation method
% 
%       - 0 : Generalized/Weighted Least squares (with grid)
%       - 1 : Maximum Likelihood (with grid),
%       - 2 : Maximum Likelihood (without grid) using fmincon as optimizer
% 
%     ngrid (integer | {250}) : number of grid points
% 
%  Returns:
%     :
% 
%     - **yh** [ts] : (disaggregated) high-frequency left-hand-side variable
% 
%     - **res** [struct] : structure with further details on the computations
% 
%  REFERENCES:
%     - Chow, G. and Lin, A.L. (1971) "Best linear unbiased
%       distribution and extrapolation of economic time series by related
%       series", Review of Economic and Statistics, vol. 53, n. 4, p. 372-375.
%     - Bournay, J. y Laroque, G. (1979) "Reflexions sur la methode
%       d'elaboration  des comptes trimestriels", Annales de l'INSEE, n. 36, p.
%       3-30.
% 
%