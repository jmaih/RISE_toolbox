%  `IMMN` or IMM(N) filtering engine.
% 
%  Syntax:
%    [LogLik,Incr,retcode,Filters] = immn(syst,y,U,z,options)
%    [LogLik,Incr,retcode,Filters] = immn(syst,y,U,z,options,k)
%    [LogLik,Incr,retcode,Filters] = immn(...,hkkm_smoother)
% 
%  Inputs:
%    - `syst`: Internals of RISE
%    - `y`: Internals of RISE
%    - `U`: Internals of RISE
%    - `z`: Internals of RISE
%    - `options`: Internals of RISE
%    - `k`: k = 1, 2, 3, ... Order of the GPBN algorithm with  
%    - `hkkm_smoother` ({true}|false): do the Hashimzade,
%      Kirsanov, Kirsanova and Maih (2024) smoother. Otherwise do the naive
%      smoother.  
% 
%  Outputs:
%    - `LogLik`: Log likelihood
%    - `Incr`: Vector of log densities for each time t, t = 1, 2, ..., n. 
%      The sum of these gives the likelihood.
%    - `retcode`: return code (=0 if there is no problem)
%    - `Filters`: Structure with filtering information
% 
%  Examples of calls from a dsge model:
% ----------------------------------------
% 
%  imm2 + hkkm smoothing : 
%                  filter(m1,'kf_user_algo','immn')
% 
%  same as above : 
%                  filter(m1,'kf_user_algo',{'immn',2})
% 
%  imm2 + naive smoothing : 
%                  filter(m1,'kf_user_algo',{'immn',2,false})
% 
%  imm2 + hkkm smoothing : 
%                  filter(m1,'kf_user_algo',{'immn',2,true})
% 
%  See also: dsge.filter
%