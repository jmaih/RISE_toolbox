function d=diff(x,wrt)

[nrows,ncols]=size(x);
if nrows>1||ncols>1
    d=cell(nrows,ncols);
    for irow=1:nrows
        for icol=1:ncols
            d{irow,icol}=diff(x(irow,icol),wrt);
        end
    end
    return
end

if isempty(x.incidence)
    % numbers/vectors and variables which are not part of differentiation
    d=0;
elseif isempty(x.args)
    % variables that are part of differentiation
    d=x.incidence(wrt);
    if any(d)
        d=double(d);
    else
        % no need to carry around a vector
        d=0;
    end
else
    % functions of variables
    if_elseif_flag=strcmp(x.func,'if_elseif');
    if_then_else_flag=strcmp(x.func,'if_then_else');
    nargs=numel(x.args);
    d_args=x.args;
    odd=true;
    for iarg=1:nargs
        if (odd && if_elseif_flag)||(iarg==1 && if_then_else_flag)
            continue
        end
        d_args{iarg}=diff(x.args{iarg},wrt);
        odd=~odd;
    end
    d=compose_derivatives(x,d_args);
end
% make sure we remain planar
d=planar(d);

end


