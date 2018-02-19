function [express,engine]=parameter_parsing_tool(variables,pnames)

mkvc=parser.cell2matize(variables.markov_chains);

endovars=parser.cell2matize(variables.endogenous);

exovars=parser.cell2matize(variables.exogenous);

koef='a|b|c|s';

lag='\d+';

third=['\d+','|',mkvc];

second=[third,'|',endovars,'|',exovars];

fourth='\d+';

eqtn='\d+|\w+';

express=['\<(',koef,')',... one letter
    '(',lag,')?'... possibly but not necessarily followed by a digit
    '\(',... opening parenthesis
    '(',eqtn,')',... digits or endogenous variable
    ',?',... possibly followed by comma (standard deviations may not be)
    '(',second,')?',... possibly followed by digits or word character
    ',?',... possibly followed by a comma
    '(',third,')?',... possibly followed by a word character
    ',?',... possibly followed by a comma
    '(',fourth,')?',... possibly followed by digits
    '\)' % closing parenthesis
    ];

engine=@replacer;

    function out=replacer(koef,lag,eqtn,second,third,fourth)
        
        if ~isempty(eqtn)
            
            pos=find(strcmp(eqtn,variables.endogenous));
            
            if ~isempty(pos)
                
                eqtn=int2str(pos);
                
            end
            
        end
        
        switch koef
            
            case {'a','b'}
                
                out=do_symbolization(true);
                
            case {'c','s'}
                
                out=do_symbolization(false);
                
        end
        
        if ~any(strcmp(out,pnames))
            
            error('parameter not found in the list of estimated params')
            
        end
        
        function out=do_symbolization(islag)
            
            out=koef;
            
            if islag
                % a2(1,gdp), a2(1,gdp,mc,3)
                if isempty(lag)
                    
                    error('lag missing for a lag term')
                    
                end
                
                out=[out,lag];
                
            else
                % c(1,0), c(1,xv), c(1,xv,mc,3)
                if ~isempty(lag)
                    
                    error('deterministic terms cannot have lags')
                    
                end
                
            end
            
            is_ready=false;
            
            if islag||strcmp(koef,'c')
                
                vbl=check_variable(second,'endogenous');
                
            else
                
                if all(isstrprop(second,'digit'))
                    % proceed as normal
                    vbl=second;
                    
                elseif any(strcmp(second,variables.markov_chains))
                    
                    out=[out,'_',eqtn,'_',second];
                    
                    if ~isempty(third)||isempty(fourth)
                        
                        error('wrong specification of a standard deviation term')
                        
                    end
                    
                    state=fourth;
                    
                    out=[out,'_',state];
                    
                    is_ready=true;
                    
                else
                    
                    error('covariance and standard deviation terms cannot contain variables')
                    
                end
                
            end
            
            if ~is_ready
                
                out=[out,'_',eqtn,'_',vbl];
                
                if ~isempty(third)
                    
                    mc=check_variable(third,'markov_chains');
                    
                    if isempty(fourth)
                        
                        error('state missing')
                        
                    end
                    
                    state=fourth;
                    
                    out=[out,'_',third,'_',state];
                    
                end
                
            end
            
        end
        
        function vbl=check_variable(vbl,group)
            
            if isempty(vbl)
                
                error('variable is empty')
                
            end
            
            group0=group;
            
            try
                
                group=variables.(group);
                
            catch
                
                error(['group ',group,' does not exist'])
                
            end
            
            if ~all(isstrprop(vbl,'digit'))
                
                loc=find(strcmp(vbl,group));
                
                if isempty(loc)
                    
                    error(['variable ',vbl,...
                        ' does not exist among ',group0])
                    
                end
                
                vbl=int2str(loc);
                
            end
            
        end
        
    end

end