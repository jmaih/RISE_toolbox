function [sol,M]=solve(self,varargin)

[sol,M]=solve@abstvar(self,varargin{:});

resolve();

    function resolve()
        
        if ~self.is_panel||...
                self.nregs>1||...
                strcmp(self.homogeneity,'unrestricted')
            
            return
            
        end
        
        nx=self.nx; ng=self.ng; h=size(sol.Q.Q,1); nvars=self.nvars; nlags=self.nlags;
        
        [kng,T0]=size(self.estim_.X);
        
        % number of coefficients per equation if we had a plain VAR
        k=kng/ng;
        
        batch_rows=cell(1,ng);
        
        batch_cols=cell(1,ng);
        
        for ig=1:ng
            
            [batch_rows{ig},batch_cols{ig}]=abstvar.map_panel(nvars,nx,nlags,ig,ng);
            
        end
        
        switch self.homogeneity
            
            case 'pooled'
                
                pooled_solution()
                
            case 'meanGroup'
                
                mean_group_solution()
                
            case 'independent'
                
                independent_solution()
                
            case {'static','dynamic'}
                % - static: static homogeneity: deterministic coefficients common
                
                % - dynamic: dynamic homogeneity: lag coefficients common,
                %   different constants
                
                separate_covariances()
                
            otherwise
                
                error(['unrecognized homogeneity "',self.homogeneity,'"'])
                
        end
        
        function separate_covariances()
            
            sol.S(:)=0;
            
            for ih=1:h
                
                % residuals
                ri=self.estim_.Y-sol.B(:,:,ih)*self.estim_.X;
                                
                for g=1:self.ng
                    
                    Sig=the_covariance(ri(batch_rows{g},:));
                    
                   sol.S(:,:,ih)=replace(sol.S(:,:,ih),Sig,batch_rows);
                    
                end
                
            end
            
        end
                
        function independent_solution()
            
            S=sol.S;
            
            sol.S(:)=0;
            
            for g=1:ng
                
                sol.S(batch_rows{g},batch_rows{g},:)=S(batch_rows{g},batch_rows{g},:);
                
            end
            % the covariance was computed with respect to kng, which is
            % incorrect under independence...
            sol.S=sol.S*(T0-kng)/(T0-k);
            
        end
        
        function pooled_solution()
            
            sol.S(:)=0;
            
            for ih=1:h
                
                % residuals
                ri=self.estim_.Y-sol.B(:,:,ih)*self.estim_.X;
                
                S=cell(1,ng);
                
                for g=1:self.ng
                    
                    S{g}=ri(batch_rows{g},:);
                    
                end
                
                S=the_covariance(cell2mat(S));
                
                sol.S(:,:,ih)=replace(sol.S(:,:,ih),S,batch_rows);
                
            end
            
        end
        
        function mean_group_solution()
            
            sol.S(:)=0;
            
            for ih=1:h
                
                B=0;
                
                for g=1:ng
                    
                    B=B+sol.B(batch_rows{g},batch_cols{g},ih);
                    
                end
                
                B=B/ng;
                
                sol.B(:,:,ih)=replace(sol.B(:,:,ih),B,batch_rows,batch_cols);
                
                % residuals
                ri=self.estim_.Y-sol.B(:,:,ih)*self.estim_.X;
                
                S=0;
                
                for g=1:self.ng
                    
                    Sig=the_covariance(ri(batch_rows{g},:));
                    
                    S=S+Sig;
                    
                end
                
                S=S/ng;
                
                sol.S(:,:,ih)=replace(sol.S(:,:,ih),S,batch_rows);
                
            end
                        
        end
            
        function X=replace(X,x,rows,cols)
            
            if nargin<4
                
                cols=rows;
                
            end
            
            for gg=1:self.ng
                
                X(rows{gg},cols{gg})=x;
                
            end
            
        end
        
        function c=the_covariance(d)
            
            T=size(d,2);
            
            c=d*d.'/(T-k);
            
        end
        
    end

end
