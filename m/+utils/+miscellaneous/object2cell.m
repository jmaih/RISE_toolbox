function [C,hyper_class]=object2cell(obj)
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

% nobj=numel(obj);
if isempty(obj)
    C=obj;
    hyper_class='';
    return
end
hyper_class=class(obj(1));
fields=fieldnames(obj);
nfields=numel(fields);
C=cell(3,nfields);
for ifield=1:nfields
    C{1,ifield}=fields{ifield};
    myclass=class(obj(1).(fields{ifield}));
    C{3,ifield}=myclass;
    switch myclass
        case 'char'
            C{2,ifield}=transpose({obj.(fields{ifield})});
        case {'double','logical'}
            C{2,ifield}=vertcat(obj.(fields{ifield}));
        otherwise
            error([mfilename,':: unknown class ',myclass])
    end
end
