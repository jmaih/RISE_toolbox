function [TM,Journal,Regimes,derivsForm]=transition_matrices2transition_matrix(transition_matrices,varlist,wrt)
if nargout>3
    % this will be important for writing the function
    wrt=rise_sym.initialize_differentiation_session(wrt);
    rise_sym.push('tag',sprintf('%.0f',0));
end
% first add the transition matrix indexes to the varlist so that they are
% recognized as variables
test=union(varlist,{'s0','s1'});
if numel(test)<numel(varlist)+2
    error('s0 and s1 are reserved words')
end
varlist=test;
chain_names=fieldnames(transition_matrices);
nchains=numel(chain_names);
Journal={};
TM=[];
nstates=nan(1,nchains);
for ichain=1:nchains
    tm=transition_matrices.(chain_names{ichain});
    tm=rise_sym.equation2rise_sym(tm,varlist,{wrt.name});
    [nrows,ncols]=size(tm);
    nstates(ichain)=nrows;
    JJ=cell(nrows,ncols);
    for irow=1:nrows
        for jcol=1:ncols
            JJ{irow,jcol}=[sprintf('%0.f',irow),';',sprintf('%0.f',jcol)];
        end
    end
    update_transition(JJ,tm)
end

if nargout>1
    [Regimes,Journal_check]=chain_grid(nstates);
    
    if ~isequal(Journal_check,Journal)
        error('Journal does not check')
    end
    
    if nargout>3
        derivsForm=struct('derivs',TM,'map',{{}},'order',0);
        G=['G',sprintf('%.0f',derivsForm.order)];
        [nrows,ncols]=size(TM);
        for irow=1:nrows
            for jcol=1:ncols
                new_ref=[G,'(',sprintf('%.0f',irow),',',sprintf('%.0f',jcol),')'];
                derivsForm.derivs(irow,jcol)=rise_sym.swap_references(...
                    derivsForm.derivs(irow,jcol),new_ref,derivsForm.order);
            end
        end
    end
end

    function update_transition(jnl,tm)
        if isempty(Journal)
            Journal=jnl;
            TM=tm;
        else
            TM=kron(TM,tm);
            siz0=size(Journal);
            siz1=size(jnl);
            siz=siz0.*siz1;
            Jtmp=Journal;
            Journal=cell(siz);
            for io=1:siz0(1)
                for jo=1:siz0(2)
                    semicol=find(Jtmp{io,jo}==';');
                    TODAY_=Jtmp{io,jo}(1:semicol-1);
                    TOMORROW_=Jtmp{io,jo}(semicol+1:end);
                    for i1=1:siz1(1)
                        rr=(io-1)*siz1(1)+i1;
                        for j1=1:siz1(2)
                            cc=(jo-1)*siz1(2)+j1;
                            semicol=find(jnl{i1,j1}==';');
                            today_=jnl{i1,j1}(1:semicol-1);
                            tomorrow_=jnl{i1,j1}(semicol+1:end);
                            Journal{rr,cc}=[TODAY_,',',today_,';',TOMORROW_,',',tomorrow_];
                        end
                    end
                end
            end
        end
    end
end