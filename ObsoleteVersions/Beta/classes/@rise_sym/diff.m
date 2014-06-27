function d=diff(x,wrt)

nrows=numel(x);
ncols=numel(wrt);
if nrows*ncols>1
    d=multi_algorithm(x,wrt);
else
    d=uni_algorithm(x,wrt);
end

    function d=uni_algorithm(x,wrt)
        if isempty(x.incidence)||~x.incidence(wrt.id)
            d=rise_sym(0);%d=rise_sym.push('get_zero');
        elseif isempty(x.args)
            if x.incidence(wrt.id)
                d=rise_sym(1);%d=rise_sym.push('get_one');
            else
                d=rise_sym(0);%d=rise_sym.push('get_zero');
            end
        else
            d=rise_sym.push('get_derivative',x,wrt);
        end
    end
    function d=multi_algorithm(x,wrt)
        d=rise_sym.empty(0);
        for irow=1:nrows
            for jcol=1:ncols
                d(irow,jcol)=uni_algorithm(x(irow),wrt(jcol));
            end
        end
    end
end


