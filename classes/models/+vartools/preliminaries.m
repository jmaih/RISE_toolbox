function out=preliminaries(varargin)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

out=struct.empty(0);
if nargin==0
    return
end
out=struct();
r=varargin{1};
varargin=varargin(2:end);
rtmp=eval([r.model_class,'.template()']);
rfields=fieldnames(r);
test=rfields(~ismember(rfields,fieldnames(rtmp)));
if ~isempty(test)
    disp(test)
    error(['the fields above should not appear in build up of a "',r.model_class,'" object'])
end
out.construction_data=r;
endo_names=collect_names(r.endogenous);
out.exogenous=strcat('EXO_',endo_names);
new_obs=endo_names(~ismember(endo_names,r.observables));
r.observables=[r.observables,new_obs];
observs=collect_names(r.observables);

for ifield=1:numel(rfields)
    ff=rfields{ifield};
    if isempty(strfind(ff,'restriction'))
        out.(ff)=r.(ff);
    else
        out.restrictions.(ff)=r.(ff);
    end
end
det_vars_loc=~ismember(observs,endo_names);
det_vars=observs(det_vars_loc);
out.exogenous=union(out.exogenous,det_vars);
% exog=exog(~det_vars_loc);
n=numel(endo_names);
out.nx=numel(det_vars)+out.constant;
nlags=out.nlags;
allien=out.block_exogenous(~ismember(out.block_exogenous,endo_names));
if ~isempty(allien)
    disp(allien)
    error('the variables above should be endogenous for them to be block exogenous')
end
out.is_block_exogenous=ismember(endo_names,out.block_exogenous);

param_template=vartools.set_structural_matrices_structure(out.model_class,n,nlags,out.nx,out.is_block_exogenous);

% create the parameters
%-----------------------
out.markov_chains=vartools.reset_markov_chains(param_template,out.markov_chains,endo_names);

% finally
%--------
out.remains=varargin;
out.param_template=param_template;

    function c=collect_names(x)
        tmpx=regexp(x,'".+"','match');
        tmpx=[tmpx{:}];
        c=sort(setdiff(x,tmpx));
    end
end