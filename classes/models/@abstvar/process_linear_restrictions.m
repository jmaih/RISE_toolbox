function self=process_linear_restrictions(self)

self=setup_nonlinear_restrictions(self);

restr=union(self.linear_restrictions,self.linear_restrictions_prime_time);

nc=numel(restr);

if nc==0
    
    return
    
end

pnames=self.mapping.estimList;

mypositions=struct();

for iii=1:numel(pnames)
    
    mypositions.(pnames{iii})=iii;
    
end

np=numel(pnames);

variables=struct('endogenous',{self.endogenous},...
    'exogenous',{self.exogenous},...
    'markov_chains',{self.mapping.regimes(1,2:end)});

[express,rplfunc]=opening_remarks(variables,pnames); %#ok<ASGLU>

A=cell(nc,1);

b=zeros(nc,1);

for ic=1:nc
    
    [A{ic},b(ic)]=process_one_restriction(restr{ic});
    
end

A=cell2mat(A);

% chop only if it is a reduced-form VAR with constant parameters
%----------------------------------------------------------------
isConstant=numel(variables.markov_chains)==1 && ~self.optimize;

if isConstant && ~ isempty(self.nonlinres)
    
    error('nonlinear restrictions not allowed in the constant-parameter VARs')
    
end

% number of coefficients per equation
%------------------------------------
K=(self.nx+self.nlags*self.nvars)*self.ng;

neqtns=self.nvars*self.ng;

npc=K*neqtns;

cutoff=npc*isConstant+np*(1-isConstant);

self.linres=vartools.linear_restrictions(A(:,1:cutoff),b);

    function [a,b]=process_one_restriction(restri)
        
        orig_restri=restri;
        
        restri(isspace(restri))=[];
        
        restri=strrep(restri,';','');
        
        try
            % put in symbolic form
            restri=regexprep(restri,express,'${rplfunc($1,$2,$3,$4,$5,$6)}');
            
            eq=find(restri=='=');
            
            if ~isempty(eq)
                % flush everything to the left-hand side
                restri=[restri(1:eq-1),'-(',restri(eq+1:end),')'];
                
            end
            
            [fn,fnr,args,positions]=functionalize(restri);
            
            [a,b]=differentiate_first_order();
            
        catch
            
            error(['Problem with restriction ',orig_restri])
            
        end
        
        function [a,b]=differentiate_first_order()
            
            % to be replaced with
            
            a=zeros(1,np);
            
            step=sqrt(eps);
            
            nargs=numel(args);
            
            argins=num2cell(zeros(1,nargs));
            
            f0=fn(argins{:});
            
            b=-f0;
            
            for iarg=1:nargs
                
                argi=argins;
                
                argi{iarg}=argi{iarg}+step;
                
                f1=fn(argi{:});
                
                a(positions(iarg))=(f1-f0)/step;
                
            end
            
            if self.debug
                
                test=utils.numdiff.jacobian(fnr,ones(nargs,1));
                
                disp(orig_restri)
                
                disp([test;a(positions)])
                
                keyboard
                
            end
            
        end
        
    end

    function [fn,fnr,args,positions]=functionalize(restri)
        
        expr=parser.cell2matize(pnames(:).');
        
        items=regexp(restri,['\<',expr,'\>'],'match');
        
        args=items;
        
        nargs=numel(args);
        
        items=cell2mat(strcat(items,','));
        
        fn=['@(',items(1:end-1),')',restri];
        
        fnr=re_encode_function(fn);
        
        fn=str2func(fn);
        
        fnr=str2func(fnr);
        
        positions=find_positions();
        
        function positions=find_positions()
            
            positions=locate_variables(args,pnames);
            
        end
        
        function fn=re_encode_function(fn)
            
            rp=find(fn==')',1,'first');
            
            xxx_='xxx_';
            
            fn=[fn(1:2),xxx_,fn(rp:end)];
            
            for ii=1:nargs
                
                fn=regexprep(fn,['\<',args{ii},'\>'],[xxx_,'(',int2str(ii),',:)']);
                
            end
            
        end
        
    end

end


function [express,engine]=opening_remarks(variables,pnames)

mkvc=parser.cell2matize(variables.markov_chains);

endovars=parser.cell2matize(variables.endogenous);

exovars=parser.cell2matize(variables.exogenous);

koef='a|b|c|s';

lag='\d+';

third=['\d+','|',mkvc];

second=[third,'|',endovars,'|',exovars];

fourth='\d+';

eqtn='\d+';

% express=['\<(',koef,')',... one letter
%     '(',lag,')?'... possibly but not necessarily followed by a digit
%     '(?:\(|_)',... opening parenthesis
%     '(',eqtn,')',... digits
%     '(?:,|_)?',... possibly followed by comma (standard deviations may not be)
%     '(',second,')?',... possibly followed by digits or word character
%     '(?:,|_)?',... possibly followed by a comma
%     '(',third,')?',... possibly followed by a word character
%     '(?:,|_)?',... possibly followed by a comma
%     '(',fourth,')?',... possibly followed by digits
%     '\)?' % closing parenthesis
%     ];

express=['\<(',koef,')',... one letter
    '(',lag,')?'... possibly but not necessarily followed by a digit
    '\(',... opening parenthesis
    '(',eqtn,')',... digits
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