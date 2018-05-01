function blocks=create_estimation_blocks(obj,blocks)
% H1 line
%
% ::
%
%
% Args:
%
%    -**obj** [rise|dsge]: model object
%
%    -**blocks** [cell array]: blocking information
%
% Returns:
%    :
%
% Note:
%
%    this function separates the parameters to estimate into blocks controlled
%    by one or several markov chains, depending on the information provided in
%    blocks
%    case 1: blocks={[1,3,5,7],[2,9,10],[6,20],...} The numbers represent the
%    order in which the estimated parameters are entered
%    case 2: blocks={{'alpha','beta'},'gam',{'delta','upsil','omicr'},...}. In
%            this case, the first block includes alpha and beta, which are
%            either parameter names or chain names. If they are parameters,
%            then alpha and beta will be estimated simultaneously. If they are
%            chain names, then all the parameters controlled by alpha and all
%            the parameters controlled by beta will be estimated
%            simultaneously
%    The estimated parameters that are controlled by a markov chain can be
%    entered either as pname(chain,state) or as pname_chain_state
%
% Example:
%
%    See also:

if isempty(obj)
    
    blocks=cell(0,4);
    
    return
    
end

if nargin<2
    
    blocks=[];

end

if isempty(obj.estimation.priors)
    
    warning('no estimated parameters')
    
    return

end

chain_names=obj.markov_chains.chain_names;

all_estim_chains={obj.estimation.priors.chain};

param_names={obj.estimation.priors.name};

% transform the names from name(chain,state) to the name_chain_state
%-------------------------------------------------------------------
param_names=parser.param_texname_to_param_name(param_names);

npar=numel(param_names);

nblks=numel(blocks);

if nblks==0
    
    random_blocks()

else
    
    ontheline=1:npar;
    
    reblocks=cell(1,nblks);
    
    for iblok=1:nblks
        
        chains=blocks{iblok};
        
        continue_flag=false;
        
        if iscell(chains)||ischar(chains)
            % the user has entered either the parameter names or the chain
            % names. We don't know yet
        elseif isnumeric(chains) % {[1,4,5],[3,2,10],...}
            % here the user has entered the ranks of the parameters
            if numel(unique(chains))~=numel(chains)
                
                error(['parameters duplicated in block ',int2str(iblok)])
            
            end
            
            % put the parameter locations into the internal bloc
            %---------------------------------------------------
            reblocks{iblok}=chains;
            
            ontheline=setdiff(ontheline,chains);
            
            continue_flag=true;
        
        else
            
            error('elements of blocks must be cellstr arrays or char')
        
        end
        
        if ischar(chains)
            
            chains=cellstr(chains);
        
        end
        
        if ~ continue_flag
            % make valid names since we now know we have cell arrays of
            % strings
            %-----------------------------------------------------------
            chains=parser.param_texname_to_param_name(chains);
            
            for ic=1:numel(chains)
                
                if ismember(chains{ic},chain_names)
                    
                    if continue_flag
                        
                        error('mixing parameter names and chain names is not allowed')
                    
                    end
                    
                elseif ismember(chains{ic},param_names)
                    
                    loc=find(strcmp(chains{ic},param_names));
                    
                    % put the parameter locations into the internal bloc
                    %---------------------------------------------------
                    reblocks{iblok}=[reblocks{iblok},loc];
                    
                    ontheline(ontheline==loc)=[];
                    
                    continue_flag=true;
                    
                    continue
                
                else
                    
                    error([chains{ic},' not recognized as a chain name or a parameter name'])
                
                end
                
            end
            
        end
        
        if continue_flag
            
            continue
        
        end
        
        % the user has entered the chain names as blocks
        %-----------------------------------------------
        discard=false(1,numel(ontheline));
        
        for ipar=1:numel(ontheline)% this is dynamic
            
            pchain=all_estim_chains{ontheline(ipar)};
            % there could be several chains lumped together
            if ismember(pchain,chains)
                
                reblocks{iblok}=[reblocks{iblok},ontheline(ipar)];
                
                discard(ipar)=true;
            
            end
            
        end
        
        ontheline(discard)=[];
    
    end
    
    if ~isempty(ontheline)
        
        disp(parser.param_name_to_param_texname(param_names(ontheline)))
        
        error('those parameters were not assigned a block')
    
    end
    
    blocks=reblocks;

end

for iblok=1:nblks
    
    tmp=unique(blocks{iblok});
    
    if numel(tmp)~=numel(blocks{iblok})
        
        error(['parameter repeated in block ',int2str(iblok)])
    
    end
    
end

    function random_blocks()
        
        maxblk=min(10,npar);
        
        nblks=ceil(npar/maxblk);
        
        blocks=cell(1,nblks);
        
        d1=randperm(npar);
        
        for iblk=1:nblks
            
            blocks{iblk}=d1((iblk-1)*maxblk+1:min(maxblk*iblk,npar));
        
        end
        
    end

end


