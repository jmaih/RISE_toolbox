function [allOrders,zeroth_order]=print(derivs,...
    min_ncount) % minimum number of occurrences at which substitution occurs
global  rise_sym_main_map
if nargin<2
    min_ncount=[];
end
expansion_factor=300;
if isempty(min_ncount)
    min_ncount=3;
end

derivTags={derivs.tag}; % flag for elements not to be tagged
nitems=numel(derivs);
if numel(derivTags)~=nitems
    error('the number of untagged should be the same as the number of derivatives groupings')
end
% dividing lines in the unordered map: it comes from Java and it is a
% uint64 object. Although it does not matter, we double it.
partitions=double([derivs.last_line]); 

% now get the map and re-order it
%--------------------------------
Count=double(rise_sym_main_map.fid.Count);
rise_sym_main_map=values(rise_sym_main_map.fid);
unordered_map=cell(Count,3);
prototype=cell(1,3);
for ic=1:Count
    vi=rise_sym_main_map{ic};
    prototype{1}=vi.obj.ref;
    prototype{2}=vi.obj.key;
    prototype{3}=vi.obj.ncalls;
    unordered_map(vi.rank,:)=prototype;
end

% for fast evaluation, now can replace & with &&
%-----------------------------------------------
unordered_map(:,1:2)=strrep(unordered_map(:,1:2),'&','&&');

% dejaVu=cell(expansion_factor,2); % containers.Map()
% idv=0;
dejaVu=containers.Map();
code=cell(1,nitems);
nout=0;
for item=1:nitems
    prologue=cell(0,3);
    auxiliary_equations=cell(0,3);
    
    % tag only where necessary
    %-------------------------
    if ~isempty(derivTags{item})
        obj=derivs(item).derivs;
        order=derivs(item).order;
        if order>9
            error('string replacements will not execute propertly for orders greater than 9. Contact junior.maih if you believe order greater than 9 are important')
        end
        order_string=sprintf('%0.0f',derivs(item).order);
        [nrows,ncols]=size(obj);
        Gtag=derivTags{item};
        auxiliary_equations=retag_map(obj);
        % write the prologue
        %-------------------
        nout=nout+1;
        prologue=[{'if nargout',sprintf('%0.0f',nout-1),inf} % lhs, rhs, #calls
            {Gtag,['sparse(',sprintf('%0.0f',nrows),',',sprintf('%0.0f',ncols),')'],inf} % lhs, rhs, #calls
            ];
    end
    
    % find the appropriate batch
    %---------------------------
    lastrow=partitions(1);
    % if lastrow==0 this could be an indication that the derivatives for
    % the higher order are 0. This could have served as a basis for
    % detecting whether the model is linear had it not been for the
    % inclusion of parameter derivatives, which need not enter the model
    % linearly. In any case, this is just a necessary condition but not a
    % sufficent one. For instance the function exp(x) has all its
    % higher-order derivatives equal to exp(x), but only one such equation
    % will be stored in the map, while the objects corresponding to the
    % derivatives will just have a pointer to that value.
    
%     if lastrow==0,keyboard, end
    
    code{item}=[
        prologue
        unordered_map(1:lastrow,:)
        auxiliary_equations
        ];
    % inefficient?
    
    locs=strncmp('G',code{item}(:,1),1);
    code{item}(locs,3)={inf};
    
    partitions=partitions(2:end)-lastrow;
    unordered_map=unordered_map(lastrow+1:end,:);
end

% Trim the zeroth_order separately if it is tagged
%-------------------------------------------------
% note that trim returns only two columns...
zeroth_order=[];
if ~isempty(derivTags{1})
    zeroth_order=rise_sym.trim(code{1},min_ncount);
    zeroth_order=strcat(zeroth_order(:,1),'=',zeroth_order(:,2),';');
    % close
    zeroth_order=[zeroth_order;'end;']; % semicolon needed for later parsing
    % replace nargout= with nargout>
    zeroth_order=strrep(zeroth_order,'nargout=','nargout>');
end

% put all the code together and trim it
%--------------------------------------
if nitems>1
    item=nitems+1;
    while item>2
        item=item-1;
        code{item-1}=[code{item-1};code{item}];
        code(item)=[];
    end
    allOrders=rise_sym.trim(code{1},min_ncount);
    % need to put line breaks here before proceeding
    allOrders=strcat(allOrders(:,1),'=',allOrders(:,2),';');
    allOrders=[allOrders;repmat({'end;'},nout,1)];
    % replace nargout= with nargout>
    allOrders=strrep(allOrders,'nargout=','nargout>');
else
    allOrders=zeroth_order;
end

    function extra_equations=retag_map(thisobj)
        % retag only the selected location. If a location is not selected,
        % it won't have auxiliary equations.
        extra_equations=cell(expansion_factor,3);
        iextra=0;
        for cc=1:ncols
            for rr=1:nrows
                if is_zero(thisobj(rr,cc))
                    continue
                end
                G=[Gtag,'(',sprintf('%0.0f',rr),',',sprintf('%0.0f',cc),')'];
                if isnumeric(thisobj(rr,cc))
                    add_row_to_cell_extra_equations({G,sprintf('%0.0f',thisobj(rr,cc).func),inf}); % lhs, rhs, #calls
                elseif isempty(thisobj(rr,cc).args)
                    add_row_to_cell_extra_equations({G,thisobj(rr,cc).func,inf}); % lhs, rhs, #calls
                else
                    if ~isKey(dejaVu,thisobj(rr,cc).ref);
                        if strcmp(order_string,thisobj(rr,cc).ref(5))
                            % then we know that the reference will appear
                            % on the left hand side
                            % a reference of order i cannot be found in batch i-1
                            %----------------------------------------------
                            replace_occurrences_in_map(thisobj(rr,cc).ref,G);
                            dejaVu(thisobj(rr,cc).ref)=G;
                        else
                            % the reference was constructed in some earlier
                            % order
                            %----------------------------------------------
                            add_row_to_cell_extra_equations({G,thisobj(rr,cc).ref,inf});  % lhs, rhs, #calls
                        end
                    else
                        add_row_to_cell_extra_equations({G,dejaVu(thisobj(rr,cc).ref),inf});  % lhs, rhs, #calls
                    end
                end
            end
        end
        extra_equations=extra_equations(1:iextra,:);
        function replace_occurrences_in_map(ref,G)
            unordered_map(:,1:2)=regexprep(unordered_map(:,1:2),['(?<!\w+)',ref,'(?!\w+)'],G);
        end
        
        function add_row_to_cell_extra_equations(myrow)
            iextra=iextra+1;
            if size(extra_equations,1)==iextra
                equation_ncols=size(extra_equations,2);
                extra_equations(end+expansion_factor,:)=cell(1,equation_ncols);
            end
            extra_equations(iextra,:)=myrow;
        end
    end
end