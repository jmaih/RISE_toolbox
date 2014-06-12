function d=sort(db,dim,mode,vname)
if nargin<4
    vname=[];
end
if nargin<3 || isempty(mode)
    mode='ascend';
end
if ~any(strcmp(mode,{'ascend','descend'}))
    error('3rd argument should be "ascend" or "descend"')
end
if nargin<2 || isempty(dim)
    dim=1;
end
if ~any(dim==1:3)
    error('dim must be either 1, 2 or 3')
end
nvars=db.NumberOfVariables;

d=main_frame(db,false);

switch dim
    case 1
        % pick one column
        col_choice=1;
        if nvars>1
            if isempty(vname)
                error(['the database has many variables, a variable has to',...
                    ' be chosen in the 4th argument (arguments 2 and 3 can be empty)'])
            end
            col_choice=locate_variables(vname,db.varnames);
            if numel(col_choice)>1
                error(['many variables with name ',vname,' in the database'])
            end
            % a column (single variable) has to be chosen, hence, variables must
            % have names
        end
        [~,tags]=sort(db.data(:,col_choice,:),dim,mode);
        for ipage=1:size(tags,3) %<-- db.NumberOfPages
            tags_i=[1;tags(:,:,ipage)+1];
            d(:,:,ipage)=d(tags_i,:,ipage);
        end
    case {2,3}
        if dim==2 && nvars>1 && ~all(strcmp(db.varnames{iname},db.varnames{1}))
            error('all variables must have the same name')
        end
        [~,tags]=sort(db.data,dim,mode);
        tmp=d(2:end,2:end,:);
        tmp=tmp(tags);
        d(2:end,2:end,:)=tmp;
end

end