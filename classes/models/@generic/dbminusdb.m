function db=dbminusdb(obj,odb,cdb)
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

if isempty(obj)

    db=cell(0,4);

    return

end


db=struct();

endo_names=obj.endogenous.name;

is_log_var=obj.log_vars;

fields=fieldnames(odb);

for ifield=1:numel(fields)
    
    name=fields{ifield};
    
    if ~isfield(cdb,name)
        
        continue
        
    end
    
    lv=false;
    
    pos=find(strcmp(name,endo_names));
    
    if ~isempty(pos)
        
        lv=is_log_var(pos);
        
    end
    
    if lv
        
        db.(name)=odb.(name)/cdb.(name);
        
    else
        
        db.(name)=odb.(name)-cdb.(name);
        
    end
    
end

end