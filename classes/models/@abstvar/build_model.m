function [mc,switch_prior,restr]=build_model(m,vnames,options)

% this function creates inputs for model setting and estimation of a
% structural VAR identified by Choleski restrictions.
% the inputs are:
% - m [string] nvmc or nvABCmc, where n is the number of states for the
% variances, m the number of states for the coefficients. A, B, C are the
% equations identified by variable names.
% - vnames [cellstr]: list of the endogenous variables ordered according
% to their position in the VAR
% - options [struct|{}]: default options for beta and dirichlet
% distributions
% 
% the outputs are
% - mc [struct]: structure containing information about the markov chains
%
% - switch_prior [struct] structure containing the priors on the transition
% probabilities
%
% - restr [cell array]: list of (Choleski) restrictions. state-identifying
% restrictions are not added yet

% models={'1v|1c','2v|1c','3v|1c','1v|2c','2v|2c','3v|2c',...
%     '3v|S2c','3v|SC2c','3v|SCP2c','3v|SRM2c','3v|RM2c','3v|RMC2c'
% 'S2v|C2v|M2c'};

if nargin<3
    
    options=struct();
    
    options.beta=struct('weight_diag',0.85,'std_off_diag',0.1);
    
    options.dirichlet=struct('weight_diag',0.8,'std_diag',0.1);
    
end

mc=struct('name',{},...
    'number_of_states',{},...
    'controlled_parameters',{},...
    'endogenous_probabilities',{},...
    'probability_parameters',{});

switch_prior=struct();

restr=cell(0,1);

aControls=cell(3,0);

splits=split_model(m,vnames);

% remove chains with one state
%-----------------------------
States=[splits.nstates];

splits(States==1)=[];

if isempty(splits)
    
    do_choleski_restrictions()
    
    return
    
end

ndirich=0; coef_iter=0; vol_iter=0;

all_types={splits.type};

is_coeff=strcmp(all_types,'c');

ncoeff=sum(is_coeff);

nvol=sum(~is_coeff);

aControls=cell(3,ncoeff);

for isp=1:numel(splits)
    
    this=splits(isp);
    
    if strcmp(this.type,'c')
        
        coef_iter=coef_iter+1;
        
        mcname=set_chain_name(this.type,coef_iter,ncoeff);
        
        [c,a]=create_controls(this.eqtns,'c','a');
        
        cas={c,a};
        
        if isempty(this.eqtns)
            % by default you control all the equations
            this.eqtns=1:numel(vnames);
            
        end
        
        aControls{1,coef_iter}=this.eqtns;
        
        aControls{2,coef_iter}=mcname;
        
        aControls{3,coef_iter}=this.nstates;
        
%         typical_param='a0';
        
    else
        
        vol_iter=vol_iter+1;
        
        mcname=set_chain_name(this.type,vol_iter,nvol);
        
        cas={create_controls(this.eqtns,'s')};
        
    end
    
    mc(isp)=struct('name',mcname,...
        'number_of_states',this.nstates,...
        'controlled_parameters',{cas},...
        'endogenous_probabilities',[],...
        'probability_parameters',[]);
    
    do_priors(mcname,this.nstates);
    
%     do_state_identification_restrictions(mcname,this.nstates,typical_param)
    
end

do_choleski_restrictions()

    function do_state_identification_restrictions(mcname,nstates,pname)
        
        for is0=2:nstates
            
            restr=[restr
                {sprintf('%s(%s,%0.0f)>=%s(%s,%0.0f)',...
                pname,mcname,is0,...
                pname,mcname,is0-1)}]; %#ok<AGROW>
            
        end
        
    end

    function do_choleski_restrictions()
        
        % set the Choleski identification
        %--------------------------------
        neqtns=numel(vnames);
        
        for ieqtn=1:neqtns
            
            mcname_=''; nstates=[];
            
            for icol=1:size(aControls,2)
                
                if ismember(ieqtn,aControls{1,icol})
                    
                    mcname_=aControls{2,icol}; 
                    
                    nstates=aControls{3,icol};
                    
                    break
                    
                end
                
            end
            
            for jeqtn=ieqtn+1:neqtns
                
                left=sprintf('a0(%0.0f,%s',ieqtn,vnames{jeqtn});
                
                if isempty(mcname_)
                    
                    batch={[left,')=0']};
                    
                else
                    
                    batch=cell(0,1);
                    
                    for istate=1:nstates
                        
                        batch=[batch
                            {sprintf('%s,%s,%0.0f)=0',left,mcname_,istate)}]; %#ok<AGROW>
                        
                    end
                    
                end
                
                restr=[restr
                    batch]; %#ok<AGROW>
                
            end
            
        end
        
        % set the normalization : no need as the diagonal elements are
        % already unity
        %-------------------------------------------------------------
        
    end

    function do_priors(mc,nstates)
        
        is_dirichlet=nstates>2;

        if is_dirichlet
            
            weight_diag=options.dirichlet.weight_diag;
            
            std_diag=options.dirichlet.std_diag;
            
            common_mean=(1-weight_diag)/(nstates-1);
            
        else
            
            weight_diag=options.beta.weight_diag;
            
            std_off_diag=options.beta.std_off_diag;
            
        end
        
        for ii=1:nstates
            
            if is_dirichlet
                
                ndirich=ndirich+1;
                
                dname=sprintf('dirichlet_%0.0f',ndirich);
                
                space=cell(2,nstates-1);
                
            end
            
            diter=0;
            
            for jj=1:nstates
                
                if ii==jj,continue,end
                
                tp=sprintf('%s_tp_%0.0f_%0.0f',mc,ii,jj);
                
                if is_dirichlet
                    
                    diter=diter+1;
                    
                    space{1,diter}=tp;
                    
                    space{2,diter}=common_mean;
                    
                else
                    
                    dname=tp;
                    
                    space={1-weight_diag,1-weight_diag,std_off_diag,'beta'};
                    
                end
                
            end
            
            if is_dirichlet
                
                space=[{std_diag},space(:).'];
                
            end
            
            switch_prior.(dname)=space;
            
        end
        
    end

end
        
function varargout=create_controls(eqtns,varargin)

varargout=varargin;

if isempty(eqtns)
    
    return
    
end

n=numel(eqtns);

eqtns=num2cell(eqtns(:).');

eqtns=cellfun(@(x)int2str(x),eqtns,'uniformOutput',false);

eqtns=cell2mat(strcat(eqtns,','));

eqtns=eqtns(1:end-1);

if n>1
    
    eqtns=['[',eqtns,']'];
    
end

for ii=1:length(varargin)
    
    varargout{ii}=[varargin{ii},'(',eqtns,')'];
    
end

end

function cn=set_chain_name(type_,iter,n)

switch type_
    case 'c'
        cn='syncoef';
    case 'v'
        cn='synvol';
end

if n>1
    
    cn=sprintf('%s%0.0f',cn,iter);
    
end


end

function splits=split_model(m,vList)

split_signs={'&','|'};

splits=first_split();

splits=second_split(splits);

    function s=second_split(splits)
        
        s=struct('eqtns',{},'nstates',{},'type',{});
        
        for isp=1:numel(splits)
            
            this=splits{isp};
            
            d=find(isstrprop(this,'digit'));
            
            if isempty(d)
                
                error([this,':: code should have a digit'])
                
            elseif numel(d)>1 && ~all(d(2:end)-d(1:end-1)==1)
                
                error([this,':: digits should be consecutive'])
                
            end
            
            s(isp).eqtns=set_equations(this(1:d(1)-1));
            
            s(isp).nstates=str2double(this(d));
            
            if s(isp).nstates==1 && ~isempty(s(isp).eqtns)
                
                error([this,':: When the number of states is 1 one cannot have a separate markov chain for specified variables'])
                
            end
            
            s(isp).type=this(d(end)+1:end);
            
            if ~ismember(s(isp).type,{'c','v'})
                
                error([this,':: types should be "c" or "v"'])
                
            end
            
        end
        
        function e=set_equations(eqtns)
            
            eqtns0=eqtns;
            
            eqtns(isspace(eqtns))=[];
            
            e=[];
            
            batch='';
            
            while ~isempty(eqtns)
                
                e1=eqtns(1);
                
                if any(strcmp(e1,{'(',')'}))
                    
                    e=flush_(e);
                    
                else
                    
                    batch=[batch,e1]; %#ok<AGROW>
                    
                end
                
                eqtns=eqtns(2:end);
                
            end
            
            e=flush_(e);
            
            function e=flush_(e)
                
                if isempty(batch)
                    
                    return
                    
                end
                
                pos=find(strcmp(batch,vList));
                
                if isempty(pos)
                    
                    error(['variable "',batch,'" not listed among endogenous'])
                    
                end
                
                if ismember(pos,e)
                    
                    error(['variable "',batch,'" listed more than once in ',eqtns0])
                    
                end
                
                e=[e,pos];
                
                batch='';
                
            end
            
        end
        
    end

    function splits=first_split()
        
        splitters=[];
        
        for ii=1:numel(split_signs)
            
            this=find(m==split_signs{ii});
            
            splitters=[splitters,this(:).']; %#ok<AGROW>
                        
        end
        
        if any(splitters==1)||any(splitters==length(m))
            
            error(['splitters "',split_signs,'" cannot occur at the beginning or at the end of a model code'])
            
        end
        
        max_splits=100;
        
        splits=cell(1,max_splits);
        
        iter_split=1;
        
        while ~isempty(m)
            
            m1=m(1);
            
            if any(strcmp(m1,split_signs))
                
                iter_split=iter_split+1;
                
                if iter_split>max_splits
                    
                    error('Your model has too many switches')
                    
                end
                
            else
                
                splits{iter_split}=[splits{iter_split},m1];
                
            end
            
            m=m(2:end);
            
        end
        
        splits=splits(1:iter_split);
        
    end

end