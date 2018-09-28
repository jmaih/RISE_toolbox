function retcode=code2file(xcell,fname,directory)
% INTERNAL FUNCTION
%

retcode=0;

if isempty(xcell)
    
    retcode=603; % empty cell
    
    return

end

default_output_name='out_';

output_list=default_output_name;

input_list='';

if iscell(xcell)
    
    [xout]=main_engine();
    
    xout=[{[default_output_name,'=[']};xout(:);{'];'}];

elseif isstruct(xcell)
    
    derivative_fields={'size','derivatives'};
    
    eval_fields={'code','argins','argouts'};
    
    if all(isfield(xcell,derivative_fields))
        
        [xout]=derivative_engine();
    
    elseif all(isfield(xcell,eval_fields))
        
        [xout]=eval_engine();
    
    else
        % then it must be a transition matrix
        [xout]=transition_matrix_engine();
    
    end
    
elseif isa(xcell,'function_handle')
    
    retcode=utils.code.code2file({xcell},fname,directory);
    
    return

else
    
    error('first input must be a cell, a structure or a function handle')

end

if ~retcode
    % write to file
    %--------------
    file_writer(fname,xout,input_list,output_list,directory);

end

    function [xout]=transition_matrix_engine()
        
        xstruct=xcell;
        
        chain_names=fieldnames(xstruct);
        
        is_loose_commit=any(strcmp(chain_names,parser.loose_commit()));
        
        xout={[default_output_name,'=struct();'];'Q=1;'};
        
        if is_loose_commit
           
            xout=[
                xout
                {'Qinit=1;'}
                ];
        
        end
        
        for iname=1:numel(chain_names)
            
            chain=chain_names{iname};
            
            xcell={xstruct.(chain)};
            
            [xout_i]=main_engine();
            
            xout=[
                xout
                [default_output_name,'.',chain,'=',xout_i{1},';']
                {['Q=kron(Q,',default_output_name,'.',chain,');']}
                ]; %#ok<AGROW>
            
            if is_loose_commit && ~strcmp(chain,parser.loose_commit())
                
                xout=[
                    xout
                    {['Qinit=kron(Qinit,',default_output_name,'.',chain,');']}
                    ]; %#ok<AGROW> 
            
            end
            
        end
        
        xout=[
            xout
            {['[',default_output_name,'.Q,retcode]=utils.code.validate_transition_matrix(Q);']}
            {'if ~retcode'}
            ];
        
        if ~is_loose_commit
            
            xout=[
                xout
                {'Qinit=Q;'}
                ];
        
        end
        
        xout=[
            xout
            {['[',default_output_name,'.Qinit,retcode]=utils.code.validate_transition_matrix(Qinit);']}
            {'end'}
            ];
        
        output_list=['[',output_list,',retcode]'];
    
    end

    function [xout]=derivative_engine()
        
        tmp=xcell;
        
        order=numel(tmp);
       
        xout={};
        
        this_output_name=cellfun(@(x)x(~isspace(x)),...
            strcat({default_output_name},...
            num2str((1:order)')),'uniformOutput',false);
        
        prologue={};
        
        different_orders={'first','second','third','fourth','fifth',...
            'sixth','seventh','eigth','nineth','tenth'};
        
        end_prologue=0;
        
        for io=1:order
            
            [xout_io]=do_one_order(io);
            
            subfunc_name=sprintf('do_%s_order',different_orders{io});
            
            add_on={sprintf('if nargout >%0.0f',io-1)};
                        
            end_prologue=end_prologue+1;
           
            prologue=[
                prologue
                add_on
                {[this_output_name{io},'=',subfunc_name,...
                input_list,';']}
                ]; %#ok<AGROW>
            
            % write subfunction to file
            %---------------------------
            file_writer(subfunc_name,xout_io,input_list,output_list,...
                [directory,'/private']);

        end
        
        xout=[
            prologue
            repmat({'end'},end_prologue,1)% closing the if nargout > n-1
            xout
            {'end'}
            ];
        
        output_list=cell2mat(strcat(this_output_name(:)',','));
        
        output_list=['[',output_list(1:end-1),']'];
        
        function [xout]=do_one_order(oo)
            
            xcell=tmp(oo).derivatives(:,1);
            
            locs=tmp(oo).derivatives(:,2);
            
            % output initialization
            tmp_size=full(tmp(oo).size);
            
            xxx=sprintf('%s=spalloc(%0.0f,%0.0f,%0.0f);',...
                default_output_name,tmp_size(1),tmp_size(2),...
                tmp(oo).nnz_derivs);
            
            [xout]=main_engine();
            
            incr=1000;
            
            siz_xout=incr;
            
            xout_=cell(siz_xout,1);
            
            iter_out=0;
            
            if ~isempty(xcell)
                
                for irows=1:size(xcell,1)
                    
                    if ~isempty(xout{irows})
                        
                        cols=locs{irows}(2:end);
                        
                        strcols=stringify_indexes(cols);
                        
                        assign(locs{irows}{1},strcols,xout{irows})
                        
                        % reassess
                        %---------
                        if oo>1
                            
                            for ii=1:numel(cols)
                                
                                cols_i=cols{ii};
                                
                                if numel(cols_i)==1
                                    
                                    continue
                                    
                                end
                                
                                stud=cols_i(1);
                                
                                strcols=stringify_indexes(cols_i(2:end));
                                
                                string=sprintf('%s(%0.0f,%0.0f)',...
                                    default_output_name,locs{irows}{1},stud);
                                
                                assign(locs{irows}{1},strcols,string)
                                
                            end
                            
                        end
                
                    end
                    
                end
                
                % add initialization
                xout=[{xxx};xout_(1:iter_out)];
            
            else
                
                xout={xxx};
        
            end
            
            xout=[xout;{'end'}];
            
            function assign(row,strcols,string)
                
                if iter_out==siz_xout
                    
                    xout_=[xout_;cell(incr,1)];
                    
                    siz_xout=siz_xout+incr;
                    
                end
                
                iter_out=iter_out+1;
                
                xout_{iter_out}=sprintf('%s(%0.0f,%s)=%s;',...
                    default_output_name,row,...
                    strcols,string);
                
            end
            
        end
        
    end

    function strcols=stringify_indexes(indexes)
        
        n=numel(indexes);
        
        cell_type=iscell(indexes);
        
        ind=1;
        
        strcols=sprintf('%0.0f',pull_out_one());
        
        for ind=2:n
        
            strcols=[strcols,',',sprintf('%0.0f',pull_out_one())]; %#ok<AGROW>
        
        end
        
        if n>1
        
            strcols=['[',strcols,']'];
    
        end
        
        function n=pull_out_one()
            
            if cell_type
                
                n=indexes{ind}(1);
                
            else
                
                n=indexes(ind);
                
            end
            
        end
        
    end

    function [xout]=main_engine()
        
        n=numel(xcell);
        
        xout=xcell;
        
        entry_gate='';
        
        for item=1:n
            
            if ~isempty(xcell{item})
                
                if ~isa(xcell{item},'function_handle')
                    
                    retcode=604;%<-- error('all elements in xcell should be function handles')
                    
                    return
                                
                end
                
                xout{item}=func2str(xcell{item});
                
                if isempty(entry_gate)
                
                    right_parenth=find(xout{item}==')',1,'first');
                    
                    entry_gate=xout{item}(1:right_parenth);
            
                end
                % remove the first occurrence only
                %----------------------------------
                xout{item}=xout{item}(right_parenth+1:end);
                
            end
            
        end
                
        if isempty(input_list)
        
            input_list=strrep(entry_gate,'@','');
    
        end
        
    end

    function [code]=eval_engine()
        
        if isempty(input_list)
            
            input_list=cell2mat(strcat(xcell.argins,','));
        
            input_list=['(',input_list(1:end-1),')'];
        
        end
        
        output_list=cell2mat(strcat(xcell.argouts,','));
        
        output_list=['[',output_list(1:end-1),']'];
        
        if isempty(xcell.code)
            
            code={};
        
            retcode=603;
       
        else
            
            code=regexp(xcell.code,';','split');
            % replace the nargout_ which is used in evaluation with nargout, which is
            % used in the normal function
            code=regexprep(code,'(?<!\w+)narg(out|in)_(?!\w+)','narg$1');
            
            code=code(:);
            
            code=code(cell2mat(cellfun(@(x)~isempty(x),code,'uniformOutput',false)));
        
            code=strcat(code,';');
    
        end
        
    end

end

function file_writer(fname,xout,input_list,output_list,directory)

if isempty(input_list)
    
    input_list='(~,~,~,~,~,~,~,~,~,~)';
    
end

fid=fopen([fname,'.m'],'w');

fprintf(fid,'%s\n',['function ',output_list,'=',fname,input_list,...
    '%#ok<*INUSL,*INUSD>']);

fprintf(fid,'%s\n',['% Code automagically generated by RISE on ',...
    datestr(now)]);

fprintf(fid,'%s\n\n',['%--------------------------------------------------',...
    '------------']);

for icod=1:numel(xout)
    
    if ~isempty(xout{icod})
        
        fprintf(fid,'%s\n\n',[xout{icod}]);
        
    end
    
end

fclose(fid);

if ~isdir(directory)
    
    mkdir(directory)
    
end

movefile([fname,'.m'],directory,'f')

end
