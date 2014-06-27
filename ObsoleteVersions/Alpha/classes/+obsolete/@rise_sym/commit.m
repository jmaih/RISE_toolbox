function obj=commit(obj)
global rise_sym_main_map
if isempty(obj.args)||~isempty(obj.ref)
    return
end

[eqtn,obj]=recompose(obj);
if ~isKey(rise_sym_main_map.fid,eqtn);
    Count=double(rise_sym_main_map.fid.Count)+1;
    ref=['ref_',rise_sym_main_map.tag,'_',sprintf('%.0f',Count),'_']; % adding an underscore latter so that we can use strrep, which should be faster, rather than regexp
    obj.ref=ref; % handle will push this information automatically.
    obj.key=eqtn; % handle will push this information automatically.
    % handle will give a handle rather than writing the full object
    % into newline, hopefully this will lead to efficiency gains
    newline=struct('diff_refs',{[]},'obj',obj,'rank',Count);%...'ref',ref,'key',eqtn,
    % this is a new equation and has not been differentiated yet, so there
    % is not need to initialize the element at position diff_loc
    rise_sym_main_map.fid(eqtn)=newline;
else
    newline=rise_sym_main_map.fid(eqtn);
    % handle will push this information automatically, even inside the
    % container. At least I hope so.
%    newline.obj.ncalls=newline.obj.ncalls+1;
    obj=newline.obj;
end

end
