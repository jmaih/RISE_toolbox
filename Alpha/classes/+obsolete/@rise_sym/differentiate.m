function [derivs,zeroth_order]=differentiate(model_equations,varlist,wrt,incidence,partitions,order,func_name)
if nargin<7
    func_name=[];
    if nargin<6
        order=[];
        if nargin<5
            partitions=[];
            if nargin<4
                incidence=[];
                if nargin<3
                    wrt=[];
                    if nargin<2
                        error('insufficient number of arguments')
                    end
                end
            end
        end
    end
end
global rise_sym_main_map

if ischar(wrt)
    wrt=cellstr(wrt);
end
nwrt=numel(wrt);
if isempty(order),order=2;end

if isempty(wrt),wrt=varlist;end

if isempty(incidence),incidence=true(1,nwrt);end

if isempty(partitions),partitions={'a';nwrt};end

if numel(incidence)~=nwrt
    error('incidence and wrt should have the same number of elements')
end
labels=partitions(1,:);
parts=cell2mat(partitions(2,:));
nparts=numel(parts);
if sum(parts)~=nwrt
    error('sum of partitions does not match the number of variables to differentiate with respect to')
end
long_labels='';
for ipart=1:nparts
    % multiplying letters with numbers give numbers. We need to take the
    % char to put it back to normal char
    long_labels=[long_labels,char(labels{ipart}*ones(1,parts(ipart)))]; %#ok<AGROW>
end
wrt_locs=cumsum(incidence);
% nwrt_eff=sum(incidence);
wrt_locs(~incidence)=nan;
orig_wrt=struct();
for ii=1:nwrt
    % current position/original position
    orig_wrt(ii).name=wrt{ii};
    orig_wrt(ii).id=wrt_locs(ii);
    % the variables with no incidence will not have an imaginary part
end
wrt=orig_wrt; clear orig_wrt
wrt_wise=true;
% rise_sym.initialize_differentiation_session(wrt(incidence));

rise_sym_main_map=struct(...
    'fid',containers.Map(),...
    'nwrt',numel(wrt(incidence)),...
    'wrt',{wrt(incidence)},...
    'zero',rise_sym(0),...
    'one',rise_sym(1),...
    'tag',sprintf('%.0f',0));

eqtns=rise_sym.equation2rise_sym(model_equations,varlist,{wrt(incidence).name});
[neqtns,ncols]=size(eqtns);

nwrt=numel(wrt);
derivs=struct('derivs',{},'map',{},'order',{},'pairings',{},'size',{},...
    'last_line',{},'tag',{});
myclasses=[];
for istep=1:order+1
    this_order=istep-1;
    derivs(istep).order=this_order;
    if istep==1
        derivs(istep).derivs=eqtns;
        derivs(istep).map={[];1;[]}; % index -> target -> location
        derivs(istep).size=[neqtns,ncols];
        derivs(istep).last_line=rise_sym_main_map.fid.Count;
        % no tag for the zero order
        % prepare for higher orders
        %--------------------------
        previous_indexes=derivs(istep).map;
        continue
    end
    derivs(istep).tag=['G',sprintf('%.0f',derivs(istep).order)];
    myclasses=build_grid(myclasses,nparts);
    [mypartitions,these_labels]=set_partitions();
    % do this only once and use the results for finding permutations of
    % arbitrary vectors.
    proto_permutation=cell2mat(mypermutation((1:this_order)));
    
    rise_sym_main_map.tag=sprintf('%.0f',derivs(istep).order);
    derivs(istep).derivs=rise_sym.empty(0);
    
    ncols=size(previous_indexes,2);
    current_indexes=cell(3,ncols*nwrt);  % index -> target -> location
    ci=0; % current index
    store_index=0;
    for irt=1:nwrt
        for icol=1:ncols
            last_index=previous_indexes{1,icol};
            % get target column to differentiate
            %-----------------------------------
            target=previous_indexes{2,icol};
            ci=ci+1;
            new_index=[last_index,irt];
            current_indexes{1,ci}=new_index;
            location=current_indexes{3,ci};
            if isempty(location)
                % set the location to nan so that the permutations of the
                % same index are not recomputed in future rounds
                location=nan;
                if ~(isempty(target)||~incidence(irt))
                    store_index=store_index+1;
                    % check if at least one derivative is different from zero,
                    % in which case the column qualifies for storage.
                    is_non_zero_column=false;
                    for ieq=1:neqtns
                        derivs(istep).derivs(ieq,store_index)=get_derivative(...
                            derivs(istep-1).derivs(ieq,target),irt);
                    end
                    if is_non_zero_column
                        % set the next target
                        location=store_index;
                        current_indexes{2,ci}=store_index;
                    else
                        derivs(istep).derivs(:,store_index)=[];
                        store_index=store_index-1;
                    end
                end
                % find all permutations
                perms=new_index(proto_permutation);
                % locate the permutations
                ypred=locate_permutation(perms,nwrt,wrt_wise);
                % set the location to all permutations
                current_indexes(3,ypred)={location};
            end
            pair_derivative(new_index,location);
        end
    end
    derivs(istep).map=current_indexes;
    derivs(istep).size=[neqtns,store_index];
    for jlab=1:numel(these_labels)
        oo_=these_labels{jlab};
        mypartitions.(oo_)=conclude(mypartitions.(oo_));
    end
    derivs(istep).pairings=mypartitions;
    derivs(istep).last_line= rise_sym_main_map.fid.Count;
    % prepare for next order
    %--------------------------
    previous_indexes=current_indexes;
end

% store some important fields: order and map
%-------------------------------------------
order_map=rmfield(derivs,{'derivs','map'});%'order','pairings');

% write the codes: the main derivatives and the zeroth order
%----------------
[derivs,zeroth_order]=rise_sym.print(derivs,func_name);%,chop_output

% format output
%--------------
derivs={derivs,order_map};

% % close the session
% %------------------
% rise_sym.close_differentiation_session();

    function [pp,this_labels]=set_partitions()
        nlabs=size(myclasses,1);
        pp=struct();
        this_labels=cell(1,nlabs);
        for ilab=1:nlabs
            ooo=[labels{myclasses(ilab,:)}];
            pp.(ooo)=rise_pair(ooo);
            this_labels{ilab}=ooo;
        end
    end

    function d=get_derivative(obj,ivar)
        d=diff(obj,wrt(ivar));
        is_non_zero_column=is_non_zero_column||~is_zero(d);
    end

    function pair_derivative(index,location)
        if isnan(location)
            location=[];
        end
        ooo=long_labels(index);
        update(mypartitions.(ooo),location);
    end
end

