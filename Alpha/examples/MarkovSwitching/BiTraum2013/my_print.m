function [allOrders,zeroth_order]=my_print(derivs,...
    derivTags,... flag for elements not to be tagged
    min_ncount) % minimum number of occurrences at which substitution occurs
% validNames={'y','x','ss','param','def'};
% 
%     zeroth_order=analytical_symbolic_form(zeroth_order,validNames,'analytic');
%     allOrders=analytical_symbolic_form(allOrders,validNames,'analytic');

if nargin<3
    min_ncount=[];
end
expansion_factor=300;
if isempty(min_ncount)
    min_ncount=3;
end

nitems=numel(derivs);
if numel(derivTags)~=nitems
    error('the number of untagged should be the same as the number of derivatives groupings')
end
% dividing lines in the unordered map
partitions=[derivs.last_line]; 

% now get the map and re-order it
%--------------------------------
vals=rise_sym.push('get_map');
Count=vals.fid.Count;
vals=values(vals.fid);
unordered_map=cell(Count,2);
prototype=cell(1,2);
for ic=1:Count
    vi=vals{ic};
    prototype{1}=vi.obj.ref;
    prototype{2}=vi.obj.key;
    unordered_map(vi.rank,:)=prototype;
end
clear vals 
% list_in=regexp(unordered_map(:,1),'(?!w)G\d+\(\d+,\d+\)(?<!w)','match');
% list_in=[list_in{:}];

dejaVu=cell(expansion_factor,2);
idv=0;
code=cell(1,nitems);
for item=1:nitems
    prologue=cell(0,2);
    auxiliary_equations=cell(0,2);
    
    % tag only where necessary
    %-------------------------
    if ~isempty(derivTags{item})
        obj=derivs{item};
        [nrows,ncols]=size(obj);
        Gtag=derivTags{item};
        auxiliary_equations=retag_map(obj);
        % write the prologue
        %-------------------
        prologue={Gtag,['sparse(',sprintf('%0.0f',nrows),',',sprintf('%0.0f',ncols),');']};
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
    if lastrow<=0
        keyboard
    end
    
    code{item}=[
        prologue
        unordered_map(1:lastrow,:)
        auxiliary_equations
        ];
    partitions=partitions(2:end)-lastrow;
    unordered_map=unordered_map(lastrow+1:end);
end

% Trim the zeroth_order separately if it is tagged
%-------------------------------------------------
zeroth_order=[];
if ~derivTags(1)
    zeroth_order=rise_sym.trim(code{1},min_ncount);
end

% put all the code together and trim it
%--------------------------------------
if nitems>1
    allOrders=rise_sym.trim(code,min_ncount);
else
    allOrders=zeroth_order;
end

    function extra_equations=retag_map(thisobj)
        % retag only the selected location. If a location is not selected,
        % it won't have auxiliary equations.
        extra_equations=cell(expansion_factor,2);
        iextra=0;
        for cc=1:size(thisobj,2)
            for rr=1:size(thisobj,1)
                if is_zero(thisobj(rr,cc))
                    continue
                end
                G=[Gtag,'(',sprintf('%0.0f',irow),',',sprintf('%0.0f',icol),')'];
                if isnumeric(thisobj(rr,cc))
                    add_row_to_cell_extra_equations({G,sprintf('%0.0f',thisobj(rr,cc).func)});
                elseif isempty(thisobj(rr,cc).args)
                    add_row_to_cell_extra_equations({G,thisobj(rr,cc).func});
                else
                    loc=find(strcmp(thisobj(rr,cc).ref,dejaVu(:,1)));
                    if isempty(loc)
                        % a reference of order i cannot be found in batch i-1
                        unordered_map=replace_occurrences_in_map(thisobj.ref,G);
                        add_row_to_cell_deja_vu({thisobj.ref,G});
                    else
                        add_row_to_cell_extra_equations({G,dejaVu{loc,2}});
                    end
                end
            end
        end
        extra_equations=extra_equations(1:iextra,:);
        function replace_occurrences_in_map(ref,G)
            unordered_map=regexprep(unordered_map,['(?<!\w+)',ref,'(?!\w+)'],G);
        end
        function add_row_to_cell_deja_vu(myrow)
            idv=idv+1;
            if size(dejaVu,1)==idv
                dejaVu(end+(1:expansion_factor),:)={};
            end
            dejaVu(idv,:)=myrow;
        end
        function add_row_to_cell_extra_equations(myrow)
            iextra=iextra+1;
            if size(extra_equations,1)==iextra
                extra_equations(end+(1:expansion_factor),:)={};
            end
            extra_equations(iextra,:)=myrow;
        end
    end
end