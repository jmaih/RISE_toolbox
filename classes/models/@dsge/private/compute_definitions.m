function [obj,retcode]=compute_definitions(obj,pp)
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


definitions=obj.routines.definitions;
number_of_regimes=obj.markov_chains.small_markov_chain_info.regimes_number;
if nargin<2
    pp=obj.parameter_values;
end
% evaluate definitions
retcode=0;
defs=cell(1,number_of_regimes);
for ii=number_of_regimes:-1:1 % backward to optimize speed
    tmp=utils.code.evaluate_functions(definitions,pp(:,ii)); 
    if any(isnan(tmp))||any(isinf(tmp))||any(~isreal(tmp))
        retcode=5;
        if obj.options.debug
            utils.error.decipher(retcode)
        end
        break
    end
    defs{ii}=tmp; %#ok<*EVLC>
end
obj.solution.definitions=defs;
