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
    case 'set_map' %
        rise_sym_main_map=arg2;
        Reply=[];
    case 'get_zero' 
        Reply=rise_sym_main_map.zero;
    case 'get_one' 
        Reply=rise_sym_main_map.one;
    case 'initialize' %
        wrt=struct('name',{},'id',{});
        if ischar(arg2)
            arg2=cellstr(arg2);
        end
        nvar=numel(arg2);
        for ivar=1:nvar
            wrt(ivar).name=arg2{ivar};
            wrt(ivar).id=ivar;
        end
        rise_sym_main_map=struct(...
            'fid',containers.Map(),...
            'nwrt',nvar,...
            'wrt',{wrt},...
            'zero',rise_sym(0),...
            'one',rise_sym(1));
        Reply=wrt;
    otherwise
        error('unknown command %s',Command)
end
    function d=fetch_derivative(obj,wrt)        
        old_line=rise_sym_main_map.fid(obj.key);
        batch=old_line.diff_refs;
        id=wrt.id;
        if isempty(batch)||isempty(batch{id})
            d=diff_private(obj,wrt);
            if isempty(batch)
                batch=cell(1,rise_sym_main_map.nwrt);
            end
            if isnumeric(d)
                d=rise_sym(d);
                batch{id}='constant';
            else
                batch{id}=d.key;
            end
            old_line.diff_refs=batch;
            rise_sym_main_map.fid(obj.key)=old_line;
        else
            if strcmp(batch{id},'constant')
                d=rise_sym(0);
            else
                d=old_line.obj;
            end
        end
        
        function dx=diff_private(x,wrt)
            nargs=numel(x.args);
            d_args=cell(1,nargs);
            for iarg=1:nargs
                d_args{iarg}=diff(x.args{iarg},wrt);
            end
            dx=compose_derivatives(x,d_args);
        end
    end

end
