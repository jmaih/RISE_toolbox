%--- help for rdico.shadowize_steady_state ---
%
%  at this stage, there are some expressions like :
%  - steady_state(v(r,3)) or  steady_state(v(r,t+1))
%  - steady_state(v(r,2)) or  steady_state(v(r,t))
%  - steady_state(v(r,1)) or  steady_state(v(r,t-1))
% 
%  But more generally we need to be able to handle situations such as
%  - steady_state(f(v(*))), or even
%  - steady_state(f(v(r1,c1),v(r2,c2),...,v(rn,cn))) 
%  i.e. the steady state is not applied to the
%  variable directly but rather to a function of the variable. The pain here
%  is that it requires us to go line by line and match by match
%