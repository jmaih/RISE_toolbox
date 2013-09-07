function Reply=push(Command,arg2,arg3)
persistent rise_sym_main_map

switch Command
    case 'tag'
        rise_sym_main_map.tag=arg2;
        Reply=[];
    case 'fetch_value'
        Reply=rise_sym_main_map.fid(arg2.key);
    case 'get_map'
        Reply=rise_sym_main_map;
    case 'get_derivative' %
        obj=arg2;
        wrt=arg3;
        Reply=fetch_derivative(obj,wrt);
    case 'get_count'
        Reply=rise_sym_main_map.fid.Count;
    case 'set_map' %
        rise_sym_main_map=arg2;
        Reply=[];
    case 'get_zero'
        Reply=rise_sym_main_map.zero;
    case 'get_one'
        Reply=rise_sym_main_map.one;
    case 'initialize'
        nvar=numel(arg2);
        rise_sym_main_map=struct(...
            'fid',containers.Map(),...
            'nwrt',nvar,...
            'wrt',{arg2},...
            'zero',rise_sym(0),...
            'one',rise_sym(1));
        Reply=[];
    otherwise
        error('unknown command %s',Command)
end
    function d=fetch_derivative(obj,wrt)
        if isnumeric(obj)
            d=rise_sym(0);
        else
            old_line=rise_sym_main_map.fid(obj.key);
            batch=old_line.diff_refs;
            id=wrt.id;
            if isempty(batch)||isempty(batch{id})
                d=diff_private(obj,wrt);
                if ~isnumeric(d)
                    if isempty(batch)
                        batch=cell(1,rise_sym_main_map.nwrt);
                    end
                    batch{id}=d.key;
                    old_line.diff_refs=batch;
                    rise_sym_main_map.fid(obj.key)=old_line;
                end
            else
                deriv_line=rise_sym_main_map.fid(batch{id});
                d=deriv_line.obj; %<-- d=old_line.obj;
            end
        end
        
        function dx=diff_private(x,wrt)
            nargs=numel(x.args);
            d_args=cell(1,nargs);
            odd=true;
            for iarg=1:nargs
                if odd && strcmp(x.func,'if_elseif')
                    d_args{iarg}=x.args{iarg};
                else
                    d_args{iarg}=diff(x.args{iarg},wrt);
                end
                odd=~odd;
            end
            dx=compose_derivatives(x,d_args);
        end
    end

end
