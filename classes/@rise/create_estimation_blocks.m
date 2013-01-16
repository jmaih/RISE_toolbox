function blocks=create_estimation_blocks(obj,blocks)
% this function separates the parameters to estimate into blocks controled
% by one or several markov chains, depending on the information provided in
% blocks

if isempty(obj)
    blocks=struct();
    return
end

if nargin<2
    blocks=[];
end
if isempty(obj.estimated_parameters)
    warning('no estimated parameters')
    return
end

param_names={obj.estimated_parameters.name};
npar=numel(param_names);

nblks=numel(blocks);
if nblks==0
    random_blocks()
else
    ontheline=1:npar;
    reblocks=cell(1,nblks);
    chain_names={obj.markov_chains.name};
    for iblok=1:nblks
        chains=blocks{iblok};
        if ~iscell(chains)&&~ischar(chains)
            error('elements of blocks must be cellstr arrays or char')
        end
        if ischar(chains)
            chains=cellstr(chains);
        end
        continue_flag=false;
        for ic=1:numel(chains)
            if ismember(chains{ic},chain_names)
                if continue_flag
                    error('mixing parameter names and chain names is not allowed')
                end
            elseif ismember(chains{ic},param_names) 
                reblocks{iblok}=[reblocks{iblok},find(strcmp(chains{ic},param_names))];
                continue_flag=true;
                continue
            else
                error([chains{ic},' not recognized as a chain name or a parameter name'])
            end
        end
        if continue_flag
            continue
        end
        discard=false(1,numel(ontheline));
        for ipar=1:numel(ontheline)% this is dynamic
            pname=param_names{ontheline(ipar)};
            [c,istp,cname]=detect_controling_chain(pname);
            feed_it=false;
            if istp
                if ismember(cname,chains)
                    feed_it=true;
                end
            else
                if ismember(c,chains)
                    feed_it=true;
                end
            end
            if feed_it
                reblocks{iblok}=[reblocks{iblok},ontheline(ipar)];
                discard(ipar)=true;
            end
        end
        ontheline(discard)=[];
    end
    if ~isempty(ontheline)
        disp(param_names(ontheline))
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

function [c,flag,cname]=detect_controling_chain(pname)
[flag,~,cname]=is_transition_probability(pname);
    c='const';
if ~flag
    left_par=strfind(pname,'(');
    if ~isempty(left_par)
        comma=strfind(pname,',');
        c=pname(left_par+1:comma-1);
    end
end
end
