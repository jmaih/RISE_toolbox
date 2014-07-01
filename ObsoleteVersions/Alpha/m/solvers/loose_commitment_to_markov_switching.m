function [T,R,SS,BGP]=loose_commitment_to_markov_switching(obj,T,R,SS,BGP)
% this function puts the loose commitment solution into a markov switching
% form for estimation and simulation.

try
    narginchk(5,5)
catch me
    % for backward compatibility
%%    warning(me.message)
    error(nargchk(5,5,nargin,'string')) %#ok<NCHKN>
end

% should write a function called  that does
% what I do here, so that I can be used for simulation and for the
% calculation of moments too.
    nregs=size(R,4);
if obj.is_optimal_policy_model && ~obj.options.lc_reconvexify
    T=repmat(T,[1,1,nregs]);
    % if there exists a loosecommit chain, then find when the chain is in
    % discretion
    lc_chain_loc=find(strcmp('loosecommit',obj.markov_chains.chain_names));
    if ~isempty(lc_chain_loc)
        discretionary_regimes=cell2mat(obj.markov_chains.regimes(2:end,lc_chain_loc+1))==2;
        % here is the switch: in the alternative regime, we ignore past
        % promises
        T(:,obj.endogenous.is_lagrange_multiplier,discretionary_regimes)=0;
    end
end
% put everything into a cell format
%----------------------------------
tmpT=T;T=cell(1,nregs);
tmpR=R;R=cell(1,nregs);
for ireg=1:nregs
    if obj.options.lc_reconvexify||size(tmpT,3)==1
        T{ireg}=tmpT;
    else
        T{ireg}=tmpT(:,:,ireg);
    end
    R{ireg}=tmpR(:,:,:,ireg);
end
