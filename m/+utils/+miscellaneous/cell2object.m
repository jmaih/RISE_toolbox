function obj=cell2object(C,the_class)
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

fields=C(1,:);
vals=C(2,:);
types=C(3,:); %#ok<NASGU>
nobj=size(vals{1},1);
obj=eval([the_class,'.empty(0,1)']);
nfields=numel(fields);
for iobj=1:nobj
    the_string='';
    for ifield=1:nfields
        eval(['vv_',int2str(ifield),'=extract(vals{ifield},types{ifield},iobj);'])
        the_string=[the_string,'fields{',int2str(ifield),'},vv_',int2str(ifield),',']; %#ok<AGROW>
    end
    obj(iobj,1)=eval([the_class,'(',the_string(1:end-1),')']);
end

function v=extract(vals,type,id) %#ok<DEFNU>
switch type
    case {'double','logical'}
        v=vals(id,:,:);
    case 'char'
        v=vals{id};
    otherwise
        error([mfilename,':: cannot extract type ',type])
end