%%Declarations: names are separated by a comma, a space or both
var	 X, PAI, R,  ZS, ZD

varexo ES, ED, ER

parameters nk_tp_1_2, nk_tp_2_1
% constant transition probabilities from a chain we call "nk"
%N.B: The chain name could be anything as long as it does not contain "_" (underscores)
% a valid transition probability name contains exactly 3 underscores
% N.B: In RISE, we declare only the off-diagonal elements of the transition matrix

% the remaining parameters switch depending on the specification of the model and so
% are included in the "specific" specifications.

model
	% N.B: time is denoted by () as in dynare or by {}. Below, we use the {} notation

   % Main equations
   X   = X{+1}-tau*(R-PAI{+1})+ZD;
   
   PAI = beta*PAI{+1}+kappa*X+ZS;
   
   R   = rhor*R{-1}+(1-rhor)*(gamma_1*PAI+gamma_2*X)+sigr*ER;

   % Shock processes
   ZD = rhod*ZD{-1}+sigd*ED;
   
   ZS = rhos*ZS{-1}+sigs*ES;

@#if indx_model==0
	@# include "param_original.rs"
@#elseif indx_model==1
	@# include "param_policyChangeOnly.rs"
@#elseif indx_model==2
	@# include "param_volatilityChangeOnly.rs"
@#elseif indx_model==3
	@# include "param_privateChangeOnly.rs"
@#end

parameterization
	% constant transition probabilities
	% those parameters do not change across the
	% specifications above and so they are included in this file
	nk_tp_1_2,   0.0128;
	nk_tp_2_1,   0;	
