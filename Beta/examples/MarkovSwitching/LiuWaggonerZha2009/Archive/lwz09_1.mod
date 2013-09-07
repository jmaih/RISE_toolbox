// Zheng Liu, Daniel F. Waggoner and Tao Zha (2009):
//Asymmetric expectation effects of regime shifts in monetary policy
// Review of Economic Dynamics 12 (2009) 284-303

// this model: page 284
// parameterization not given in the paper...

var A PAI;

varexo EPS_A;

parameters q_1_2 q_2_1 phi rhoa siga;

model
	phi*PAI=PAI{+1}+(1-rhoa)*A;

	A=rhoa*A{-1}+siga*EPS_A;
end

parameterization
	q_1_2, 1-0.95;
	q_2_1, 1-0.95;
	phi(q,1), .8;
	phi(q,2), 1.5;
	rhoa, .9;
	siga, .1;
end