function nmc=initialize_markov_chain(name,nstates,varargin)
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

% parse the information in the name/value pairs
pnames_defaults = {
    'is_endogenous' ,nan,@(x)islogical(x) && isscalar(x)
    'param_list',{},@(x)ischar(x)||iscellstr(x)
    'param_list_tex',{},@(x)ischar(x)||iscellstr(x)
    'is_switching',[],@(x)islogical(x)
    'duration', [],@(x)isa(x,'double') && all(real(x)>0)
    };
if nargin==0
    if nargout
        nmc=parser.initialize_markov_chain('LouisPergaudBahoyaMaih',6);
        nmc=nmc(1:0);
    else
        disp(pnames_defaults(:,1:3))
    end
    return
end

assert(isvarname(name) && ~any(name=='_'),...
    ['The name for a markov chain must be a valid variable name and not ',...
    ' contain any "_". ',name,' is not valid']);

[is_endogenous,param_list,param_list_tex,is_switching,duration]=...
    parse_arguments(pnames_defaults,varargin{:});

if ~isempty(param_list) && isempty(param_list_tex)
    param_list_tex=param_list;
end

nparam=numel(param_list);
assert(numel(param_list_tex)==nparam,'number of parameters inconsistent')

state_names=parser.create_state_list(name,nstates);

nmc=struct('name',name,...
    'number_of_states',nstates,...
    'state_names',{state_names(:).'},...
    'state_tex_names',{state_names(:).'},...
    'is_endogenous' ,is_endogenous,...
    'param_list',{param_list(:).'},...
    'param_list_tex',{param_list_tex(:).'},...
    'is_switching',is_switching,...
    'duration',duration);

end
