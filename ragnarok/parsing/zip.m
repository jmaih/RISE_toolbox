%  Used in the preparsing of for loops in the presence of multiple indexes
%  example
%  @for (aa,bb)=zip({'S','I','R'},{'Susceptible','Infected','Recovered'})
%  	BIG_U@{aa} "welfare @{bb}", C@{aa} "consumption (@{bb})"
%  @end
% 
%  produces
% 
%  BIG_US "welfare Susceptible", CS "consumption (Susceptible)"
%  BIG_UI "welfare Infected", CI "consumption (Infected)"
%  BIG_UR "welfare Recovered", CR "consumption (Recovered)"
%