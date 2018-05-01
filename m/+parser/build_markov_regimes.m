function [TM,markov_chain_info,myifelseif]=build_markov_regimes(...
transition_matrices,markov_chains,issorted)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if nargin<3
    
    issorted=true;
    
end

recreate_transition_matrix()

chain_names=fieldnames(transition_matrices);

if issorted
    
    chain_names=sort(chain_names);
    
end

markov_chain_info=struct();

chains_number=numel(chain_names);

offending_chain_names={parser.loose_commit()};

Journal={};

TM=[];

nstates=nan(1,chains_number);

for ichain=1:chains_number
    
    tm=transition_matrices.(chain_names{ichain});
    
    [nrows,ncols]=size(tm);
    
    nstates(ichain)=nrows;
    
    JJ=cell(nrows,ncols);
    
    for irow=1:nrows
        
        for jcol=1:ncols
            
            JJ{irow,jcol}=[irow;jcol];
            
        end
        
    end
    
    update_transition(JJ,tm,chain_names{ichain})
    
end

markov_chain_info.journal=Journal;

markov_chain_info.chain_is_endogenous=sparse([markov_chains.is_endogenous]);

if nargout>1
    
    [markov_chain_info,Journal_check]=update_markov_chains_info_add_regimes(...
        markov_chain_info,nstates,chain_names);
    
    markov_chain_info.small_markov_chain_info=markov_chain_info;
    
    offending_chain_loc=ismember(chain_names,offending_chain_names);
    
    if any(offending_chain_loc)
        
        small_nstates=nstates(~offending_chain_loc);
        
        small_chain_names=chain_names(~offending_chain_loc);
        
        markov_chain_info.small_markov_chain_info=...
            update_markov_chains_info_add_regimes(...
            markov_chain_info.small_markov_chain_info,...
            small_nstates,small_chain_names);
        
    end
    
    markov_chain_info.grand_chains_to_small=locate_variables(markov_chain_info.chain_names,...
        markov_chain_info.small_markov_chain_info.chain_names,true);
    
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

    function update_transition(jnl,tm,chain_name)
        
        if isempty(Journal)
            
            Journal=jnl;
            
            TM=tm;
            
        else
            
            if ~any(strcmp(chain_name,offending_chain_names))
                
                TM=parser.kron(TM,tm);
                
            end
            
            siz0=size(Journal);
            
            siz1=size(jnl);
            
            siz=siz0.*siz1;
            
            Jtmp=Journal;
            
            Journal=cell(siz);
            
            for io=1:siz0(1)
                
                for jo=1:siz0(2)
                    
                    TODAY_=Jtmp{io,jo}(1,:);
                    
                    TOMORROW_=Jtmp{io,jo}(2,:);
                    
                    for i1=1:siz1(1)
                        
                        rr=(io-1)*siz1(1)+i1;
                        
                        for j1=1:siz1(2)
                            
                            cc=(jo-1)*siz1(2)+j1;
                            
                            today_=jnl{i1,j1}(1,:);
                            
                            tomorrow_=jnl{i1,j1}(2,:);
                            
                            Journal{rr,cc}=[[TODAY_,today_];[TOMORROW_,tomorrow_]];
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end

    function [reg_info,Journal]=update_markov_chains_info_add_regimes(reg_info,nstates,chain_names)
        
        [tmp,Journal]=utils.gridfuncs.chain_grid(nstates);
        
        reg_info.regimes_number=size(tmp,1);
        
        reg_info.chains_number=numel(chain_names);
        
        reg_info.regime_names=parser.create_state_list('regime',reg_info.regimes_number);
        
        reg_info.regimes=[
            {''},chain_names(:)'
            reg_info.regime_names(:),num2cell(tmp)
            ];
        
        reg_info.regime_tex_names=reg_info.regime_names;
        
        reg_info.chain_names=sort(chain_names(:))';
        
        reg_info.chain_tex_names=reg_info.chain_names;
        %-------------
        
        state_names={};
        
        state_tex_names={};
        
        for ichain__=1:reg_info.chains_number
            
            cn=reg_info.chain_names{ichain__};
            
            loc=strcmp(cn,{markov_chains.name});
            
            new_states=markov_chains(loc).state_names;
            
            new_state_tex_names=markov_chains(loc).state_tex_names;
            
            state_names=[state_names,new_states]; %#ok<AGROW>
            
            state_tex_names=[state_tex_names,new_state_tex_names]; %#ok<AGROW>
            
        end
        
        reg_info.state_names=state_names;
        
        reg_info.state_tex_names=state_tex_names;
        
    end

    function recreate_transition_matrix()
        
        if ~isempty(transition_matrices)
            
            return
            
        end
        
        transition_matrices=struct();
        
        for ii=1:numel(markov_chains)
            
            cname=markov_chains(ii).name;
            
            ns=markov_chains(ii).number_of_states;
            
            % fake news: just to be able to use some algorithm later on
            %-----------------------------------------------------------
            transition_matrices.(cname)=repmat({sprintf('1/%0.0f',ns)},ns,ns);
            
        end
        
    end

end