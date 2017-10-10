function self=prime_time(self)

if ~self.is_panel||strcmp(self.homogeneity,'unrestricted')
    
    return
    
end

ng=self.ng;

[dynamic_own_terms,dynamic_cov_terms]=dynamic_terms();

[static_own_terms,static_cov_terms]=static_terms();

% HERE WE DO THE RESTRICTIONS CORRESPONDING TO E.G. THE PANEL VAR. OTHER
% RESTRICTIONS COULD BE ADDED WITH THE STRUCTURAL VAR BUT AS FAR AS WE ARE
% CONCERNED AT THIS STAGE THEY DON'T MATTER

%  IN ALL CASES WE LET THE FIRST COUNTRY BY THE REFERENCE COUNTRY

N=10000;

mylinres=cell(1,N);

iter=0;

switch self.homogeneity
    
    case 'pooled'
        % - pooled: everything pooled: = dynamic + static, all cross unit
        % reactions are zero
        same_dynamics()
        
        same_statics()
        
        same_variance()
        
    case 'meanGroup'
        % - meanGroup: average across groups
        if ~isempty(self.prior)
            
            error(['the Mean Group estimator for panel data is a classical ',...
                'estimator and does not require priors'])
            
        end
        
        no_covariation()
        
        same_variance()        
        
    case 'independent'
        % - meanGroup: average across groups
        if strcmp(self.homogeneity,'meanGroup') && ~isempty(self.prior)
            
            error(['the Mean Group estimator for panel data is a classical ',...
                'estimator and does not require priors'])
            
        end
        
        no_covariation()
        
        different_variance()
        
    case 'static'
        % - static: static homogeneity: deterministic coefficients common
        same_statics()
        
        no_dynamic_covariation()
        
        different_variance()
        
    case 'dynamic'
        % - dynamic: dynamic homogeneity: lag coefficients common,
        %   different constants
        same_dynamics()
        
        no_static_covariation()
        
        different_variance()
        
    otherwise
        
        error(['unrecognized homogeneity "',self.homogeneity,'"'])
        
end

mylinres=mylinres(1:iter);

self.linear_restrictions_prime_time=mylinres;

    function different_variance()
        
        if self.nregs==1
            % in the case of one regime, the estimation is done
            % analytically... and so the restrictions cannot be applied
            % beforehand to the covariance matrix, unless one decide to do
            % estimation through optimization, something that I have not
            % yet sorted out. Moreover, this will not apply to SVAR and
            % PSVAR systems.
            return
            
        end
        
        run_no_covariation_engine('s',dynamic_cov_terms,dynamic_own_terms,[],true)
        
    end

    function same_variance()
        
        if self.nregs==1
            % in the case of one regime, the estimation is done
            % analytically... and so the restrictions cannot be applied
            % beforehand to the covariance matrix, unless one decide to do
            % estimation through optimization, something that I have not
            % yet sorted out. Moreover, this will not apply to SVAR and
            % PSVAR systems.
            return
            
        end
        
        % S_i_i = S_1_1
        
        run_same_engine('s',dynamic_own_terms,dynamic_own_terms,[],true)
        
        run_no_covariation_engine('s',dynamic_cov_terms,dynamic_own_terms,[],true)
        
    end

    function same_dynamics()
        
        no_dynamic_covariation()
        
        % A_i_i = A_1_1
        
        run_same_engine('b',dynamic_own_terms,dynamic_own_terms,1:self.nlags)
        
    end

    function same_statics()
        
        no_static_covariation()
        
        % C_i = C_j
        
        run_same_engine('c',static_own_terms,dynamic_own_terms)
        
    end

    function no_covariation()
        
        no_static_covariation()
        
        no_dynamic_covariation()
        
    end

    function no_dynamic_covariation()
        
        % A_i_j=0
        
        run_no_covariation_engine('b',dynamic_cov_terms,dynamic_own_terms,...
            1:self.nlags)
        
    end

    function no_static_covariation()
        
        % C_i_j=0
        
        run_no_covariation_engine('c',static_cov_terms,dynamic_own_terms)
        
    end

    function [own_terms,cov_terms]=static_terms()
        
        % the rows are the same as those for the dynamic terms. Here we do
        % the columns
        
        own_terms=cell(1,ng);
        
        cov_terms=cell(1,ng);
        
        nd=numel(self.exogenous)+self.constant;
        
        for ii=1:ng
            
            own_terms{ii}=ii:ng:ng*nd;
            
            cov_terms{ii}=setdiff(1:ng*nd,own_terms{ii});
            
        end
        
    end

    function [own_terms,cov_terms]=dynamic_terms()
        
        nvar=numel(self.endogenous);
        
        own_terms=cell(1,ng);
        
        cov_terms=cell(1,ng);
        
        for ii=1:ng
            
            own_terms{ii}=ii:ng:ng*nvar;
            
            cov_terms{ii}=setdiff(1:ng*nvar,own_terms{ii});
            
        end
        
    end

    function pname=populate_name(pname0)
        
        for ii=1:numel(self.markov_chains)
            
            cname=self.markov_chains(ii).name;
            
            plist=self.markov_chains(ii).param_list;
            
            if any(strcmp(pname0,plist))
                
                nstates=self.markov_chains(ii).number_of_states;
                
                break
                
            end
            
        end
        
        pname=abstvar.reinflate(pname0,cname,nstates);
        
    end
        
    function run_same_engine(stud,column_terms,row_terms,lags,iscov)
        
        if nargin < 5
            
            iscov=false;
            
            if nargin<4
                
                lags=[];
                
            end
            
        end
        
        main_row_terms=row_terms{1};
        
        main_col_terms=column_terms{1};
        
        for g=2:ng
            
            my_row_terms=row_terms{g};
            
            my_col_terms=column_terms{g};
            
            for ii=1:numel(my_row_terms)
                
                irow=my_row_terms(ii);
                
                stud_i=main_row_terms(ii);
                
                for jj=1:numel(my_col_terms)
                    
                    icol=my_col_terms(jj);
                    
                    if iscov && irow<icol
                        
                        continue
                        
                    end
                    
                    stud_j=main_col_terms(jj);
                    
                    if isempty(lags)
                        
                        yourname=sprintf('%s_%0.0f_%0.0f',stud,irow,icol);
                        
                        myname=sprintf('%s_%0.0f_%0.0f',stud,stud_i,stud_j);
                        
                        do_it(yourname,myname)
                        
                    else
                        
                        for ilag=lags
                            
                            yourname=sprintf('%s%0.0f_%0.0f_%0.0f',stud,ilag,irow,icol);
                            
                            myname=sprintf('%s%0.0f_%0.0f_%0.0f',stud,ilag,stud_i,stud_j);
                            
                            do_it(yourname,myname)
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        function do_it(yourname,myname)
            
            yourname=populate_name(yourname);
            
            myname=populate_name(myname);
            
            for kk=1:numel(yourname)
                
                iter=iter+1;
                
                mylinres{iter}=sprintf('%s=%s',yourname{kk},myname{kk});
                
            end
            
        end
        
    end
        
    function run_no_covariation_engine(stud,column_terms,row_terms,lags,iscov)
        
        if nargin<5
            
            iscov=false;
            
            if nargin<4
                
                lags=[];
                
            end
            
        end
        
        if iscov && ~isempty(lags)
            
            error('covariance cannot have lags ')
            
        end
        
        for g=1:ng
            
            opponents=column_terms{g};
            
            myterms=row_terms{g};
                
            for irow=myterms
                
                for icol=opponents
                    
                    if iscov && irow<icol
                        
                        continue
                        
                    end
                    
                    if isempty(lags)
                        
                        pname=sprintf('%s_%0.0f_%0.0f',stud,irow,icol);
                        
                        do_it(pname)
                    else
                        
                        for ilag=1:self.nlags
                            
                            pname=sprintf('%s%0.0f_%0.0f_%0.0f',stud,ilag,irow,icol);
                            
                            do_it(pname)
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        function do_it(pname)
            
            % find the estimation name: inflate the markov chain
            %---------------------------------------------------
            pname=populate_name(pname);
            
            % avoid issues with the constant chain
            if ischar(pname),pname=cellstr(pname); end %
            
            for ii=1:numel(pname)
                
                iter=iter+1;
                
                mylinres{iter}=sprintf('%s=0',pname{ii});
                
            end
            
        end
        
    end

end