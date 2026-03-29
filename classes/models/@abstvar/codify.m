%--- help for abstvar/codify ---
%
%  codify: transforms expressions such as a3(v1,v2,... into a3(2,3,... This
%  can be used to facilitate the setting of priors on individual parameters
%  of the VAR to estimate.
% 
%  ::
% 
%     str=codify(self,str)
% 
%     str=codify(self,str,doCode)
% 
%  Args:
% 
%     - **self** (var object): var object
% 
%     - **str** (cellstr|string): string to transform. The string has one of
%       the following forms:
% 
%        - a3(v1,v2,...: svar and proxy_svar parameters
%        - b3(v1,v2,...: rfvar and prfvar parameters
%        - c(v1,v2,... : deterministic terms
%        - s(v1,v2,... : variance and covariance terms
% 
%     - **doCode** (true|{false}): 
% 
%        - if false, return expression with parentheses e.g. a3(2,3)
%        - if true, return expression with underscores e.g. a3_2_3
% 
%  Returns:
% 
%     str : transformed string
% 
%  Example: ::
% 
%    str=codify(psv,'a2(FFR,ygap)=0')
% 
%    str=codify(psv,'c(FFR,ygap,mpcoef,2)=0',true)
% 
%    str=codify(psv,'s(FFR,ygap,mpcoef,2)=0',true)
% 
%  .. note:: 
% 
%    - The function does not check whether the resulting parameter is
%      actually a parameter in the model
% 
%    - The function also codifies expressions and not just individual
%      parameters
%