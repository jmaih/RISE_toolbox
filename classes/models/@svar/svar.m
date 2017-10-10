classdef svar < abstvar
    
    properties(Constant,Hidden)
        
        optimize = true
        
    end
    
    methods(Access=protected,Hidden)
        
        function self=prime_time(self)
            % do the main prime time from the superclass. then add some
            % more stuff...
            self=prime_time@abstvar(self);
            
            self=prime_time_engine(self);
            
        end
        
    end
    
    methods
        
        function self=svar(varargin)
            
            self=self@abstvar(varargin{:});
            
            if nargin>0
                
                self=abstvar.recreate_parameters(self,0);
                
                % set prime-time restrictions i.e. a0_i_i =1
                self=prime_time(self);
                
            end
            
        end
                        
        function varargout=irf(self,shock_names,irf_periods,params)
            
            n=nargin;
            
            set_defaults()
            
            Rfunc=identification(self);
            
            [varargout{1:nargout}]=irf@abstvar(self,shock_names,irf_periods,params,Rfunc);
            
            function set_defaults()
                
                if n < 4
                    
                    params=[];
                    
                    if n< 3
                        
                        irf_periods=[];
                        
                        if n<2
                            
                            shock_names=[];
                            
                        end
                        
                    end
                    
                end
                
                params=solve(self,params);
                
                if isempty(irf_periods),irf_periods=40; end
                
            end
            
        end
        
        function Rfunc=identification(varargin)
            
            Rfunc=@engine;
            
            function [R,retcode]=engine(p)
                
                retcode=0;
                
                nregs=size(p.A,3);
                
                for ireg=1:nregs
                    
                    A0=p.A(:,p.nx+(1:p.nvars),ireg);
                    
                    Ri=A0\diag(p.S0(:,:,ireg));
                    
                    if ireg==1
                        
                        R=Ri(:,:,ones(nregs,1));
                        
                    end
                    
                    R(:,:,ireg)=Ri;
                    
                end
                
            end
            
        end
        
        varargout=print_structural_form(varargin)
        
    end
    
    methods(Sealed)
                
    end
    
end

function self=prime_time_engine(self)

N=1000;

mylinres=cell(1,N);

iter=0;

% scan all markov chains for a0_i_i
for ic=1:numel(self.markov_chains)
    
    mc=self.markov_chains(ic);
    
    myparams=regexp(mc.param_list,'a0_\d+_\d+','match');
    
    myparams=[myparams{:}];
    
    if isempty(myparams)
        
        continue
        
    end
    
    d=find_diag_terms();
    
    myparams=abstvar.reinflate(myparams(d),mc.name,mc.number_of_states);  
    
    for jj=1:numel(myparams)
        
        iter=iter+1;
        
        mylinres{iter}=[myparams{jj},'=1'];
    
    end
    
end

self.linear_restrictions_prime_time=mylinres(1:iter);

    function d=find_diag_terms()
        
        n=numel(myparams);
        
        d=false(1,n);
        
        for ii=1:n
            
            a0=myparams{ii};
            
            unds=find(a0=='_');
            
            d(ii)=strcmp(a0(unds(1)+1:unds(2)-1),...
                a0(unds(2)+1:end));
            
        end
        
    end

end
