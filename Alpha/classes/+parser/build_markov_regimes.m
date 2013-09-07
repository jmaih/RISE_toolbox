function [TM,markov_chain_info,myifelseif]=build_markov_regimes(transition_matrices,markov_chains)

chain_names=sort(fieldnames(transition_matrices));
markov_chain_info=struct();
markov_chain_info.chains_number=numel(chain_names);
Journal={};
TM=[];
nstates=nan(1,markov_chain_info.chains_number);
state_names={};
for ichain=1:markov_chain_info.chains_number
    tm=transition_matrices.(chain_names{ichain});
    [nrows,ncols]=size(tm);
    nstates(ichain)=nrows;
    JJ=cell(nrows,ncols);
    for irow=1:nrows
        state_names=[state_names,[chain_names{ichain},'_',sprintf('%0.0f',irow)]]; %#ok<AGROW>
        for jcol=1:ncols
            JJ{irow,jcol}=[sprintf('%0.f',irow),';',sprintf('%0.f',jcol)];
        end
    end
    update_transition(JJ,tm)
end

markov_chain_info.journal=Journal;
markov_chain_info.chain_names=chain_names(:)';
markov_chain_info.chain_tex_names=markov_chain_info.chain_names;
markov_chain_info.state_names=sort(state_names);
markov_chain_info.state_tex_names=markov_chain_info.state_names;
markov_chain_info.chain_is_endogenous=sparse([markov_chains.is_endogenous]);

if nargout>1
    [Regimes,Journal_check]=chain_grid(nstates);
    markov_chain_info.regimes_number=size(Regimes,1);
    tmp=cell(markov_chain_info.regimes_number+1,markov_chain_info.chains_number+1);
    tmp(2:end,2:end)=num2cell(Regimes);
    tmp(1,2:end)=chain_names;
    for ireg=1:markov_chain_info.regimes_number
        tmp{ireg+1,1}=['regime_',sprintf('%0.0f',ireg)];
    end
    markov_chain_info.regimes=tmp;
    markov_chain_info.regime_names=transpose(tmp(2:end,1));
    markov_chain_info.regime_tex_names=markov_chain_info.regime_names;

    if ~isequal(Journal_check,Journal)
        error('Journal does not check')
    end
    if nargout>2
        myifelseif='if_elseif(';
        for irow=1:size(TM,1)
            s0=sprintf('%0.0f',irow);
            for icol=1:size(TM,2)
                s1=sprintf('%0.0f',icol);
                myifelseif=[myifelseif,'s0==',s0,' & s1==',s1,',',TM{irow,icol},',']; %#ok<AGROW>
            end
        end
        myifelseif=[myifelseif(1:end-1),')'];
    end
end
markov_chain_info=orderfields(markov_chain_info);

    function update_transition(jnl,tm)
        if isempty(Journal)
            Journal=jnl;
            TM=tm;
        else
            TM=parser.kron(TM,tm);
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
