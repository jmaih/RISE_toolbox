function map_=map_regimes(obj,big2small)

% MAP_REGIMES -- map regimes from small ones to big ones and vice versa.
% Useful when dealing with loose commitment.
%
% ::
%
%
%   map_=MAP_REGIMES(obj)
%
%   map_=MAP_REGIMES(obj,big2small)
%
% Args:
%
%    - **obj** [dsge|rise]: model object
%
%    - **big2small** [{true}|false]: decides the direction of the mapping
%
% Returns:
%    :
%
%    - **map_** [1xh|1xbigh]: vector assigning the regimes
%
% Note:
%
% Example:
%
%    See also:

% 

if nargin<2
    
    big2small=true;
    
end

bigh=obj.markov_chains.regimes_number;

h=obj.markov_chains.small_markov_chain_info.regimes_number;

if h==bigh
    
    map_=1:h;
    
    return
    
end
    
loose_com_col= strcmp(obj.markov_chains.chain_names,parser.loose_commit);

big_regimes=cell2mat(obj.markov_chains.regimes(2:end,2:end));

small_regimes=cell2mat(obj.markov_chains.small_markov_chain_info.regimes(2:end,2:end));

big_regimes(:,loose_com_col)=[];

% map small regimes into big ones. NB: regimes are decided by rows and not
% by columns 
%-------------------------------------------------------------------------

if big2small
    
    map_=nan(1,h);
    
    for ireg=1:h
        
        map_(ireg)=find(all(bsxfun(@minus,small_regimes(ireg,:),...
            big_regimes)==0,2),1,'first');
        
    end
    
else
    
    map_=nan(1,bigh);
    
    for ireg=1:bigh
        
        map_(ireg)=find(all(bsxfun(@minus,big_regimes(ireg,:),...
            small_regimes)==0,2));
        
    end
    
end

end