function obj=setup_general_restrictions(obj)
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
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

% this can be useful also for checking restrictions such as explosive VARs,
% etc. and not just restrictions on DSGEs
nobj=numel(obj);
for ii=1:nobj
    if ~isempty(obj(ii).options.estim_general_restrictions)
        if iscell(obj(ii).options.estim_general_restrictions)
            obj(ii).general_restrictions_data=...
                @(x)obj(ii).options.estim_general_restrictions{1}(x,...
                obj(ii).options.estim_general_restrictions{2:end});
        else
            obj(ii).general_restrictions_data=obj(ii).options.estim_general_restrictions;
        end
    end
end
% estim_general_restrictions=cell(1,nobj);

