function [Regimes,journal]=chain_grid(v)
% CHAIN_GRID -- computes the regimes from one or multiple markov chains
%
% ::
%
%
%   Regimes=CHAIN_GRID(a); : one chain with "a" states
%
%   Regimes=CHAIN_GRID([a,b]); : two chains with "a" states for the first
%   and "b" states for the second
%
%   Regimes=CHAIN_GRID([a1,a2,...,an]); : n chains with "a1" states for the
%   first, "a2" states for the second,...,"an" states for the nth
%
%   [Regimes,Journal]=CHAIN_GRID(...);
%
% Args:
%
%    - **a** [scalar|vector]: number of states in each one of the independent
%              Markov chains
%
% Returns:
%    :
%
%    - **Regimes** [matrix]: Combinations of all states. The matrix is of size
%      prod([a1,a2,...,an]) x n
%
%    - **Journal** [cell array]: Description of the transition matrix. Each
%      cell contains a matrix in which the first row is the regime today and
%      the second one is the regime tomorrow. The cell array is of size
%      prod([a1,a2,...,an]) x prod([a1,a2,...,an])
%
% Note:
%
% Example:
%
%    Regimes=chain_grid(10); : one chain with 10 regimes
%
%    Regimes=chain_grid([10,5]) : 2 chains with 5 regimes for the second one
%
%    Regimes=chain_grid([10,5,3]): 3 chains with 3 regimes on the third one
%
%    See also:

% v is a vector of the number of states in each markov chain
% Regimes is an prod(v) x NumberOfParameters matrix describing the different
% composite regimes
% example
% journal gives information about the transition matrix, which is the
% tensor product of individual chains' transition matrices.

if ~isequal(v,floor(v))
    
    error([mfilename,':: elements in v must be integers'])
    
end

Regimes=utils.gridfuncs.mygrid(v);

if nargout>1
    
    pv=prod(v);
    
    journal=cell(pv);
    
    for ii=1:pv
        
        for jj=1:pv
            
            journal{ii,jj}=Regimes([ii,jj],:);
            
        end
        
    end
    
end

end

%  THE SCRIPT BELOW IS NOT TO BE DELETED. I MOVE THE FUNCTION UNDER
%  UTILITIES SO THAT I CAN USE IT FOR OTHER STUFF. BUT I KEEP A COPY HERE
%  AS A REFERENCE IN CASE I EVER DELETE THE OTHER ONE OR MODIFY IT.
% % % % function Regimes=mygrid(v)
% % % % 	n=numel(v);	
% % % % 	% first construct the indexes of the intervals
% % % % 	Regimes=transpose(1:v(1)); % first to v(1)th interval
% % % % 	
% % % % 	for p=2:n
% % % % 	    [rg,cg]=size(Regimes);
% % % % 	    vp=v(p);
% % % % 	    Ip=repmat(transpose(1:vp),1,rg);
% % % % 	    G0=nan(rg*vp,cg);
% % % % 	    for ii=1:rg
% % % % 	        G0((ii-1)*vp+1:ii*vp,:)=Regimes(ii*ones(vp,1),:);
% % % % 	    end
% % % % 	    Regimes=[G0,Ip(:)];
% % % % 	end
% % % % end