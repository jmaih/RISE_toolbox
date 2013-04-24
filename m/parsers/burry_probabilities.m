% 0-) Initialize differentiation session

% 1-) create all the functions including the probabilities, constant or not
eqtns=rise_sym_.system([probability_equations;model_equations],varlist,wrt);

% 2-) separate probabilities from model equations
nprobs=numel(probability_equations);
probability_equations=eqtns(1:nprobs);

eqtns=eqtns(nprobs+1:end);

% 3-) create kroneckers of transition matrices
% 3-a)
for chain=1:nchains
	Q=rise_sym_.empty(0);
	cumul=rise_sym_(0);
	Journal=cell(nstates)
	for s0=1:nstates
		for s1=1:nstate
			Journal{s0,s1}=[int2str(s0),';',int2str(s1)];
			if s0==s1
				continue
			end
			cumul=cumul+a_tp_i_j;
		end
		Q(s0,s0)=1-cumul;
	end
	if chain==1
		bigQ=Q;
		bigJournal=Journal;
	else
		bigQ=kron(bigQ,Q);
		bigJournal=cellkron(bigJournal,Journal);
	end
end

% 4-) create the function p(s0,s1)
% 4a-)create new rise_sym_ objects s0 and s1
s0=rise_sym('s0');
s1=rise_sym('s1');
bigp=0;
for irow=1:ncols
	for icol=1:ncols
		bigp=bigp+if_then_else(s0==irow & s1==icol,bigQ(irow,icol),0);
	end
end

% 5-) multiply the functions with the probabilities
funcvec=kron(eqtns,bigp);

% 6-) Take the derivatives of the system
