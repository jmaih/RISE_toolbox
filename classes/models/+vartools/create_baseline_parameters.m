function [nan_list,fixed]=create_baseline_parameters(param_template,param_names)
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


nan_list=[];
fixed=struct();
nout=nargout;
for ii=1:size(param_template,2)
    [jj,kk]=find(isnan(param_template{2,ii}));
    if isempty(jj)
        continue
    end
    nan_list=[nan_list;cellstr(strcat(param_template{1,ii},'_',int2str(jj),'_',int2str(kk)))]; %#ok<*AGROW>
    if nout==2
        [jj,kk]=find(~isnan(param_template{2,ii}));
        thisList=cellstr(strcat(param_template{1,ii},'_',int2str(jj),'_',int2str(kk)));
        if ~isempty(jj)
            for ll=1:numel(thisList)
                fixed.(thisList{ll})=param_template{2,ii}(jj(ll),kk(ll));
            end
        end
    end
end
if nargin>1 && ~isempty(param_names)
    fix_list=fieldnames(fixed);
    bogus=~ismember(fix_list,param_names);
    fixed=rmfield(fixed,fix_list(bogus));
end
