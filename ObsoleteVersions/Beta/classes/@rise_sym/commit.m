function obj=commit(obj)
if isempty(obj.args)||~isempty(obj.ref)
    return
end

the_map=rise_sym.push('get_map');

[eqtn,obj]=compact_form(obj);
if ~isKey(the_map.fid,eqtn);
    Count=the_map.fid.Count+1;
    ref=['ref_',the_map.tag,'_',sprintf('%.0f',Count)];
    obj.ref=ref; % handle will push this information automatically.
    obj.key=eqtn; % handle will push this information automatically.
    % handle will give a handle rather than writing the full object
    % into newline, hopefully this will lead to efficiency gains
    newline=struct('diff_refs',{[]},'obj',obj,'rank',Count);%...'ref',ref,'key',eqtn,
    % this is a new equation and has not been differentiated yet, so there
    % is not need to initialize the element at position diff_loc
    the_map.fid(eqtn)=newline;
else
    newline=the_map.fid(eqtn);
    % handle will push this information automatically, even inside the
    % container. At least I hope so.
    newline.obj.ncalls=newline.obj.ncalls+1;
    obj=newline.obj;
%%    obj.ref=newline.obj;
%%    obj.key=newline.obj;
end

rise_sym.push('set_map',the_map);

end

function [eqtn,obj]=compact_form(obj)
Args=obj.args;
for iarg=1:numel(Args)
    if ~isa(obj.args{iarg},'rise_sym')
        % handle will push this information automatically.
        obj.args{iarg}=rise_sym(obj.args{iarg});
    end
    Args{iarg}=get_reference(obj.args{iarg});
end
eqtn=rise_sym.recompose(obj.func,Args{:});
end
