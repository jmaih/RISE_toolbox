function self=process_linear_restrictions(self)

variables=struct('endogenous',{self.endogenous},...
    'exogenous',{self.exogenous},...
    'markov_chains',{self.mapping.regimes(1,2:end)});

pnames=self.estim_.links.estimList;

[express,rplfunc]=abstvar.parameter_parsing_tool(variables,pnames); 

self=setup_nonlinear_restrictions(self,express,rplfunc);

restr=union(self.estim_.linear_restrictions,...
    self.linear_restrictions_prime_time);

nc=numel(restr);

if nc==0
    
    return
    
end

mypositions=struct();

for iii=1:numel(pnames)
    
    mypositions.(pnames{iii})=iii;
    
end

np=numel(pnames);

A=cell(nc,1);

b=zeros(nc,1);

for ic=1:nc
    
    [A{ic},b(ic)]=process_one_restriction(restr{ic});
    
end

A=cell2mat(A);

% chop only if it is a reduced-form VAR with constant parameters
%----------------------------------------------------------------
isConstant=numel(variables.markov_chains)==1 && ~self.optimize;

if isConstant && ~ isempty(self.estim_.nonlinres)
    
    error('nonlinear restrictions not allowed in the constant-parameter VARs')
    
end

% number of coefficients per equation
%------------------------------------
K=(self.nx+self.nlags*self.nvars)*self.ng;

neqtns=self.nvars*self.ng;

npc=K*neqtns;

cutoff=npc*isConstant+np*(1-isConstant);

self.estim_.linres=vartools.linear_restrictions(A(:,1:cutoff),b);

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