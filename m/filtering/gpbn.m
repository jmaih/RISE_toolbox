%  `GPBN` or GPB(N) filtering engine.
% 
%  Syntax:
%    [LogLik,Incr,retcode,Filters] = gpbn(syst,y,U,z,options)
%    [LogLik,Incr,retcode,Filters] = gpbn(syst,y,U,z,options,k)
%    [LogLik,Incr,retcode,Filters] = gpbn(syst,y,U,z,options,k+1i)
%    [LogLik,Incr,retcode,Filters] = gpbn(syst,y,U,z,options,{k,collapseType})
%    [LogLik,Incr,retcode,Filters] = gpbn(...,hkkm_smoother)
% 
%  Inputs:
%    - `syst`: Internals of RISE
%    - `y`: Internals of RISE
%    - `U`: Internals of RISE
%    - `z`: Internals of RISE
%    - `options`: Internals of RISE
%    - `k`: k = 1, 2, 3, ... Order of the GPBN algorithm with Hashimzade,
%      Kirsanov, Kirsanova and Maih (2024) smoothing, heretofore hkkm.
%    - `k+1i`: the imaginary part indicates the naive smoothing is requested
%    - `collapseType` ({'bs'}|'cns'): Bar-Shalom (bs) or  consistent
%      collapsing in the GPBN spirit (cns)
%    - `hkkm_smoother` ({true}|false): do the hkkm_smoother smoother
%      otherwise do the naive smoother. This overrides the k+1i syntax,
%      which is maintained mostly for backward compatibility. 
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
%  gpb2 + bs collapse + hkkm smoothing : 
%                  filter(m1,'kf_user_algo','gpbn')
% 
%  same as above : 
%                  filter(m1,'kf_user_algo',{'gpbn',2})
% 
%  gpb2 + bs collapse + naive smoothing : 
%                  filter(m1,'kf_user_algo',{'gpbn',2+1i}) 
% 
%  same as above : 
%                  filter(m1,'kf_user_algo',{'gpbn',{2+1i,'bs'}})
% 
%  gpb2 + cns collapse + naive smoothing : 
%                  filter(m1,'kf_user_algo',{'gpbn',{2+1i,'cns'}})
% 
%  gpb2 + cns collapse + hkkm smoothing (override of naive) : 
%                  filter(m1,'kf_user_algo',{'gpbn',{2+1i,'cns'},true})
% 
%  gpb2 + cns collapse + naive smoothing (override of naive) : 
%                  filter(m1,'kf_user_algo',{'gpbn',{2+1i,'cns'},false})
% 
%  See also: dsge.filter
%