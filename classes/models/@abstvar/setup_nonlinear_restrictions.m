function self=setup_nonlinear_restrictions(self,express,rplfunc)
% setup_nonlinear_restrictions - sets nonlinear restrictions
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
%    - uses estim_nonlinear_restrictions, which should be a cell array. Each
%    item of the array is a string of the form
%      - 'f(p1,p2,...,pn)>=h(p1,p2,...,pn)'
%      - 'f(p1,p2,...,pn)>h(p1,p2,...,pn)'
%      - 'f(p1,p2,...,pn)<=h(p1,p2,...,pn)'
%      - 'f(p1,p2,...,pn)<h(p1,p2,...,pn)'
%      - 'pj=h(p1,p2,...,pn)'
%
%    - In the statements above,
%      - eqtn [digits|variable name]
%      - vbl [digits|variable name]
%      - lag [digits]
%      - chain [char]
%      - state [digits]
%
% Example:
%
%    See also:

RestrictionsBlock=self.estim_.linear_restrictions;

nc=numel(RestrictionsBlock);

is_inequality=false(1,nc);

for ii=1:nc
    
    is_inequality(ii)=any(RestrictionsBlock{ii}=='<')||...
        any(RestrictionsBlock{ii}=='>');
    
end

RestrictionsBlock=RestrictionsBlock(is_inequality);

self.estim_.linear_restrictions(is_inequality)=[];

if isempty(RestrictionsBlock)
    
    return
    
end

if isstruct(RestrictionsBlock)
    
    RestrictionsBlock=RestrictionsBlock.original;
    
end

RestrictionsBlock=cellfun(@(x)x(~isspace(x)),RestrictionsBlock,...
    'uniformOutput',false);

% put in symbolic form right here right now
%------------------------------------------
RestrictionsBlock=regexprep(RestrictionsBlock,express,'${rplfunc($1,$2,$3,$4,$5,$6)}');

param_names=self.parameters;

nparams=numel(param_names);

governing_chain=nan(1,nparams);

for ii=1:numel(self.markov_chains)
    
    these_params=self.markov_chains(ii).param_list;
    
    governing_chain(locate_variables(these_params,param_names))=ii;
    
end

% chain_names=self.markov_chain_info.small_markov_chain_info.chain_names;
chain_names=self.markov_chain_info.small_markov_chain_info.regimes(1,2:end);

regimes=cell2mat(self.markov_chain_info.small_markov_chain_info.regimes(2:end,2:end));

endo_names=self.endogenous;

RestrictionsBlock=parameterize(RestrictionsBlock);

[self.estim_.nonlinres,~,derived_parameters]=...
    generic_tools.nonlinear_restrictions_engine(...
    param_names,regimes,chain_names,governing_chain,...
    RestrictionsBlock);

if ~isempty(derived_parameters)
    
    error('derived parameters should not appear here')
    
end

    function restr=parameterize(restr)
        
        express=['\<',...
            '(a|b)',...1
            '(\d+)',...2
            '(?:\()',...
            '(\d+)',...3
            '(?:,)',...
            '(\w+)',...4
            '(?:,)?',...
            '(\w+)?',...5
            '(?:,)?',...
            '(\d+)?',...6
            '(?:\))\>'];
        
        replace=@engine; %#ok<NASGU>
        
        restr=regexprep(restr,express,'${replace($1,$2,$3,$4,$5,$6)}');
        
        function str=engine(a,lag,eqtn,vbl,chain,state)
            
            vbl=vbl2vbl(vbl);
            
            str=[a,lag,'_',eqtn,'_',vbl];
            
            if ~isempty(chain)
                
                str=[str,'(',chain,',',state,')'];
                
            end
            
            function v=vbl2vbl(v)
                
                v=locate_variables(v,endo_names);
                
                v=int2str(v);
                
            end
            
        end
        
    end

end