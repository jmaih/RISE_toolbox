function eqtn=compact_form(obj)
Args=obj.args;
for iarg=1:numel(Args)
    if ~isa(obj.args{iarg},'rise_sym')
        % handle will push this information automatically.
        obj.args{iarg}=rise_sym(obj.args{iarg});
    end
    Args{iarg}=get_reference(obj.args{iarg});
end
eqtn=rise_sym.recompose(obj.func,Args);
end