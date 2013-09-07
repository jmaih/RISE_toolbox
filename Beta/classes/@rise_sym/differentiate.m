function derivs=differentiate(model_equations,transition_matrices,varlist,wrt,order)
if nargin<4
    order=3;
end
if ischar(wrt)
    wrt=cellstr(wrt);
end

wrt=rise_sym.initialize_differentiation_session(wrt);
rise_sym.push('tag',sprintf('%.0f',0));

eqtns=rise_sym.equation2rise_sym(model_equations,varlist,{wrt.name});

if ~isempty(transition_matrices)
    TM=rise_sym.transition_matrices2transition_matrix(transition_matrices,varlist,wrt);    
    eqtns=burry_probabilities(eqtns,TM);
end

nwrt=numel(wrt);
derivs=struct('derivs',{},'map',{},'order',{});
last_guy=[];
for istep=1:order+1
    derivs(istep).order=istep-1;
    G=['G',sprintf('%.0f',derivs(istep).order)];
    derivs(istep).map={};
    if istep==1
        derivs(istep).derivs=eqtns;
        [nrows,ncols]=size(eqtns);
        for irow=1:nrows
            for jcol=1:ncols
                new_ref=[G,'(',sprintf('%.0f',irow),',',sprintf('%.0f',jcol),')'];
                derivs(istep).derivs(irow,jcol)=rise_sym.swap_references(...
                    derivs(istep).derivs(irow,jcol),new_ref,derivs(istep).order);
            end
        end
    else
        rise_sym.push('tag',sprintf('%.0f',derivs(istep).order));
        [nrows,ncols]=size(derivs(istep-1).derivs);
        derivs(istep).derivs=rise_sym.empty(0);
        old_map=derivs(istep-1).map;
        is_needed=false(1,ncols*nwrt);
        for ifunc=1:nrows
            iter=0;
            jter=0;
            funcstr=sprintf('%.0f',ifunc);
            for icol=1:ncols
                for ivar=1:nwrt
                    iter=iter+1;
                    if ifunc==1
                        if ~isempty(old_map)
                            last_guy=old_map{icol};
                        end
                        if isempty(last_guy)||last_guy(end)<=ivar
                            derivs(istep).map=[derivs(istep).map,{[last_guy,ivar]}];
                            is_needed(iter)=true;
                        end
                    end
                    if ~is_needed(iter)
                        continue
                    end
                    jter=jter+1;
                    derivs(istep).derivs(ifunc,jter)=get_derivative(...
                        derivs(istep-1).derivs(ifunc,icol),ivar,jter);
                end
            end
        end
    end
end
% now write the maps in a more convenient way
for istep=1:numel(derivs)
    derivs(istep).map=cell2mat(derivs(istep).map')';
end

% rise_sym.close_differentiation_session();

    function d=get_derivative(obj,ivar,icol)
        d=diff(obj,wrt(ivar));
        new_ref_=[G,'(',funcstr,',',sprintf('%.0f',icol),')'];
        d=rise_sym.swap_references(d,new_ref_,derivs(istep).order);
    end

end

function eqtns=burry_probabilities(eqtns,TM)
% Burrying transition probabilities
if ~isempty(TM)
    [nrows,ncols]=size(TM);
    mycell=cell(2,nrows*ncols);
    iter=0;
    for irow=1:nrows
        rowstr=sprintf('%0.f',irow);
        for jcol=1:ncols
            iter=iter+1;
            mycell(:,iter)={['s0==',rowstr,' & s1==',sprintf('%0.f',jcol),];TM(irow,jcol)};
        end
    end
    mycell=mycell(:)';
    TM=if_elseif(mycell{:});
    eqtns=kron(TM,eqtns);
end
end