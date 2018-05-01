function obj=setup_measurement_errors(obj)
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


param_names=obj.parameters.name;
obj.measurement_errors_restrictions=[];

for ii=1:sum(obj.observables.number) 
    % pick only the endogenous observables
    if ~obj.observables.is_endogenous(ii)
        continue
    end
    vi=obj.observables.name{ii};
    loc=find(strcmp(['stderr_',vi],param_names));
    % the line above is a bit inefficient. I have to change it. No big
    % deal, it is done only once
    if ~isempty(loc)
        obj.measurement_errors_restrictions=...
            [obj.measurement_errors_restrictions;ii,loc];
    end
end
