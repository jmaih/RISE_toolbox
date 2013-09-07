function d=diff(x,wrt)
global rise_sym_main_map

nrows=numel(x);
ncols=numel(wrt);
if nrows*ncols>1
    d=multi_algorithm(x,wrt);
else
    d=uni_algorithm(x,wrt);
end
    function d=multi_algorithm(x,wrt)
        d=rise_sym.empty(0);
        for irow=1:nrows
            for jcol=1:ncols
                d(irow,jcol)=uni_algorithm(x(irow),wrt(jcol));
            end
        end
    end

    function d=uni_algorithm(x,wrt)
        if isempty(x.incidence)||~x.incidence(wrt.id)
            d=rise_sym_main_map.zero; % rise_sym.push('get_zero');%d=rise_sym(0);
        elseif isempty(x.args)
            if x.incidence(wrt.id)
                d=rise_sym_main_map.one; % rise_sym(1);%d=rise_sym.push('get_one');
            else
                d=rise_sym_main_map.zero; %rise_sym(0);%d=rise_sym.push('get_zero');
            end
        else
            if isnumeric(x)
                d=rise_sym(0);
            else
                old_line=rise_sym_main_map.fid(x.key);
                batch=old_line.diff_refs;
                id=wrt.id;
                if isempty(batch)||isempty(batch{id})
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
                    d=compose_derivatives(x,d_args);
                    if ~isnumeric(d)
                        if isempty(batch)
                            batch=cell(1,rise_sym_main_map.nwrt);
                        end
                        batch{id}=d.key;
                        old_line.diff_refs=batch;
                        rise_sym_main_map.fid(x.key)=old_line;
                    end
                else
                    deriv_line=rise_sym_main_map.fid(batch{id});
                    d=deriv_line.obj; %<-- d=old_line.obj;
                end
            end
        end
    end
end


