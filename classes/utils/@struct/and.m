function h=and(varargin)
% and -- overload and for structures
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

nargs=length(varargin);
h=struct();
[allfields,specfields,typical]=get_all_fields_names();
nfields=numel(allfields);
for ifield=1:nfields
    batch=cell(1,nargs);
    name=allfields{ifield};
    for ii=1:nargs
        batch{ii}=typical;
        if any(strcmp(name,specfields{ii}))
            batch{ii}=varargin{ii}.(name);
        end
    end
    try
        h.(name)=cat(2,batch{:});
    catch
        warning(['field ',name,' failed'])
    end
end

    function [ff,fi,typical]=get_all_fields_names()
        fi=cell(1,nargs);
        for iarg=1:nargs
            hi=varargin{iarg};
            if numel(hi)~=1
                error('cannot concatenate vectors of structures')
            end
            fi{iarg}=fieldnames(hi);
            if iarg==1
                typical=nan(size(varargin{iarg}.(fi{iarg}{1})));
                ff=fi{iarg}(:).';
            end
            ff=[ff,fi{iarg}(:).']; %#ok<AGROW>
        end
    end

end
