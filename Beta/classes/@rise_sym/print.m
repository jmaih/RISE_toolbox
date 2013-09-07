function finalOutput=print(derivatives,func_name,chop_output,debug)
if nargin<4
    debug=false;
    if nargin<3
        chop_output=[];
        if nargin<2
            func_name=[];
        end
    end
end
if isempty(chop_output)
    chop_output=true;
end

vals=rise_sym.push('get_map');
Count=vals.fid.Count;
vals=values(vals.fid);
unordered_map=cell(Count,3);
prototype=cell(1,3);
for ic=1:Count
    vi=vals{ic};
    prototype{1}=vi.obj.ref;
    prototype{2}=vi.obj.key;
    prototype{3}=vi.obj.ncalls;
    unordered_map(vi.rank,:)=prototype;
end
clear vals
list_in=regexp(unordered_map(:,1),'(?!w)G\d+\(\d+,\d+\)(?<!w)','match');
list_in=[list_in{:}];

unordered_map=rise_sym.trim_derivatives(unordered_map(:,1:2),debug);

unordered_map=strcat(unordered_map(:,1),'=',unordered_map(:,2),';');
validNames={'y','x','ss','param','def'};
unordered_map=analytical_symbolic_form(unordered_map,validNames,'analytic');

G=cell(1,numel(derivatives));
for iorder=1:numel(derivatives)
    ord=iorder-1;
    G{iorder}=['G',int2str(ord)];
end
out=cell(1000,1);
itercount=0;
for iorder=1:numel(derivatives)
    ord=iorder-1;
    iter=0;
    % note the apparently un-necessary semicolon added. it will be useful
    % for separating equations later on
    if chop_output
        update_output(['if nargout_>',int2str(ord),';'])
    end
    this=derivatives(iorder).derivs;
    [nrows,ncols]=size(this);
    update_output([G{iorder},'=sparse(',int2str(nrows),',',int2str(ncols),');'])
    while iter<size(unordered_map,1)
        iter=iter+1;
        theRow=unordered_map{iter};
        if strcmp(theRow(1),'G')
            leftpar=strfind(theRow,'(');
            thisOrder=str2double(theRow(2:leftpar-1));
        elseif strcmp(theRow(1),'r')
            underscores=strfind(theRow,'_');
            thisOrder=str2double(theRow(underscores(1)+1:underscores(2)-1));
        else
            error('string unrecognized')
        end
        if ord~=thisOrder
            iter=iter-1;
            break
        end
        % write the row
        update_output(theRow)
    end
    unordered_map=unordered_map(iter+1:end);
    % write the remaining items
    for irow=1:nrows
        irow_=int2str(irow);
        for jcol=1:ncols
            jcol_=int2str(jcol);
            obj=this(irow,jcol);
            qij=[G{iorder},'(',irow_,',',jcol_,')'];
            if is_zero(this(irow,jcol))
                continue
            elseif isnumeric(obj.func)
                update_output([qij,'=',sprintf('%0.16g',obj.func),';'])
            elseif isempty(obj.args)
                update_output([qij,'=',obj.func,';'])
            else
                loc=find(strcmp(list_in,qij));
                if ~isempty(loc)
                    list_in(loc)=[];
                    continue
                else
                    update_output([qij,'=',obj.ref,';'])
                end
            end
        end
    end
    if chop_output
        % I put a semicolon, which will be useful later on when separating the
        % equations
        update_output('end;')
    end
end

out=out(1:itercount);
finalOutput=struct('code',cell2mat(out(:)'),'argins',...
    {[validNames,'s0','s1']},'argouts',{G});

if ~isempty(func_name)
    code2file(finalOutput,func_name)
end

    function update_output(in_string)
        itercount=itercount+1;
        if itercount==size(out,1),out{end+1000}=[];end
        out{itercount}=in_string;
    end
end