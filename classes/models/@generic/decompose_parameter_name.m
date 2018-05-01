function [position,regime_states,pname,chain,state]=...
decompose_parameter_name(pname,markov_chains,param_names,...
governing_chain)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


regimes=cell2mat(markov_chains.regimes(2:end,2:end));

chain_names=markov_chains.regimes(1,2:end);

position=find(strcmp(pname,param_names));

if ~isempty(position)
    
    state=1;
    
    chain='const';
    
    chain_id=find(strcmp(chain,chain_names));
    
else
    
    ptex=parser.param_name_to_param_texname(pname,chain_names);
    
    left_par=strfind(ptex,'(');
    
    right_par=strfind(ptex,')');
    
    comma=strfind(ptex,',');
    
    if isempty(left_par)||isempty(right_par)||isempty(comma)
        
        error(['"',ptex,'" or "',pname,'" is not recognized as a parameter name'])
    
    end
    
    pname=ptex(1:left_par-1);
    
    position=find(strcmp(pname,param_names));
    
    if isempty(position)
        
        error(['"',pname,'" is not recognized as a parameter name'])
        
    end
    
    chain=ptex(left_par+1:comma-1);
    
    state=str2double(ptex(comma+1:right_par-1));
    
    chain_id=find(strcmp(chain,chain_names));

end

if ~(chain_id==governing_chain(position))
    
    error(['parameter ',pname,' is not controlled by ',chain_names(chain_id)])

end

% locate the state in the regimes
regime_states=find(regimes(:,chain_id)==state);

if isempty(regime_states)
    
    error([sprintf('%0.0f',state),' is not a valid state for parameter ',pname])

end

end
