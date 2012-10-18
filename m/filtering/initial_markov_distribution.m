function PAI00=initial_markov_distribution(q_now_lead,flag)
nreg=size(q_now_lead,1);
switch flag
case{0,'ergodic'}
	A=[eye(nreg)-transpose(q_now_lead);ones(1,nreg)];
	PAI00=A\[zeros(1,nreg),1]';
case{1,'diffuse'}
    PAI00=1/nreg*ones(nreg,1);
otherwise
	error([mfilename,':: unkown '])
end
