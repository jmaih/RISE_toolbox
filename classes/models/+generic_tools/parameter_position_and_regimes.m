function [position,regime_states]=parameter_position_and_regimes(pnames,...
param_names,governing_chain,chain_names,...
grand_chains_to_small,regimes)
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

if ischar(pnames)
    pnames=cellstr(pnames);
end
npar=numel(pnames);

position=nan(npar,1);
regime_states=cell(1,npar);
for ipar=1:npar
    [position(ipar),regime_states{ipar}]=retrieval_engine(pnames{ipar});
end
    function [position,regime_states]=retrieval_engine(pname)
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
            not_good=isempty(left_par)||isempty(right_par)||isempty(comma);
            if ~not_good
                pname=ptex(1:left_par-1);
                position=find(strcmp(pname,param_names));
                not_good=isempty(position);
            end
            if not_good
                error(['"',pname,'" not recognized as a parameter name'])
            end
            chain=ptex(left_par+1:comma-1);
            state=str2double(ptex(comma+1:right_par-1));
            chain_id=find(strcmp(chain,chain_names));
        end
        govChain=grand_chains_to_small(governing_chain(position));
        if isnan(govChain)
            error(['parameter ',pname,' is not controlled by the specified chain'])
        elseif ~(chain_id==govChain)
            error(['parameter ',pname,' is not controlled by ',chain_names{chain_id}])
        end
        % locate the state in the regimes
        regime_states=find(regimes(:,chain_id)==state);
        if isempty(regime_states)
            error([sprintf('%0.0f',state),' is not a valid state for parameter ',pname])
        end
    end
end