function [finalOutput,zeroth_order]=print(derivatives,func_name,chop_output,debug)
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
unordered_map=cell(Count,2);
prototype=cell(1,2);
for ic=1:Count
    vi=vals{ic};
    prototype{1}=vi.obj.ref;
    prototype{2}=vi.obj.key;
%     prototype{3}=vi.obj.ncalls;
    unordered_map(vi.rank,:)=prototype;
end
clear vals
list_in=regexp(unordered_map(:,1),'(?!w)G\d+\(\d+,\d+\)(?<!w)','match');
list_in=[list_in{:}];
validNames={'y','x','ss','param','def'};
correct_form=@(x)analytical_symbolic_form(x,validNames,'analytic');

% % locate the first occurrence of first-order derivatives G1
% %----------------------------------------------------------
% first_start=find(strncmp(unordered_map(:,1),'G1',2),1,'first');
% zero_order_map=unordered_map(1:first_start-1,:);
% % expand the zero-order map with information about the zero-order
% % derivatives
% objzero=derivatives(1).derivs;
% G0='G0';
% zero_order_map(end+(1:numel(objzero)))={};
% ziter=0;
% zero_iter=first_start-1;
% while ziter<numel(objzero)
%     ziter=ziter+1;
%     if is_zero(objzero(ziter))
%         continue
%     end
%     zero_iter=zero_iter+1;
%     if isnumeric(objzero(ziter).func)
%         thispush=[sprintf('%0.16g',objzero(ziter).func),';'];
%     elseif isempty(objzero(ziter).args)
%         thispush=[correct_form(objzero(ziter).func),';'];
%     else
%         thispush=[correct_form(objzero(ziter).ref),';'];
%     end
%     zero_order_map(zero_iter,:)={[G0,'(',ziter,',1)'],thispush};
% end
% 
% % trim the zero-th order
% %-----------------------
% zero_order_map=rise_sym.trim_derivatives(zero_order_map,debug);

zeroth_order=[];

% trim the zero-th and higher orders
%-----------------------------------
unordered_map=rise_sym.trim_derivatives(unordered_map,debug);

unordered_map=strcat(unordered_map(:,1),'=',unordered_map(:,2),';');
unordered_map=analytical_symbolic_form(unordered_map,validNames,'analytic');

nderiv=numel(derivatives);
G=cell(1,nderiv);
for ii=1:nderiv
    G{ii}=['G',int2str(ii-1)];
end
toDelete=false(1,nderiv);
out=cell(1000,1);
itercount=0;
write_code();
finalOutput=finalize(itercount,G(~toDelete));

if ~isempty(func_name)
    code2file(finalOutput,func_name)
end

    function write_code()
        for iorder=1:nderiv
            ord=iorder-1;
            iter=0;
            this=derivatives(iorder).derivs;
            [nrows,ncols]=size(this);
            % note the apparently un-necessary semicolon added. it will be useful
            % for separating equations later on
            % try and match G{iorder}
            test=regexp(unordered_map,G{iorder},'once');
            test=[test{:}];
            if isempty(test)
               toDelete(iorder)=true;
            end
            if ord>0 && chop_output
                update_output(['if nargout_>',int2str(ord-1),';'])
            end
            if ~toDelete(iorder)
                nrows_big=derivatives(iorder).size(1);
                ncols_big=derivatives(iorder).size(2);
                update_output([G{iorder},'=sparse(',int2str(nrows_big),',',int2str(ncols_big),');'])
            end
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
                %--------------
                update_output(theRow)
            end
            unordered_map=unordered_map(iter+1:end);
            % write the remaining items if ord>0
            if ord>0
                for irow=1:nrows
                    irow_=int2str(irow);
                    for jcol=1:ncols
                        jcol_=sprintf('%0.0f',jcol);%int2str(jcol);
                        obj=this(irow,jcol);
                        qij=[G{iorder},'(',irow_,',',jcol_,')'];
                        if is_zero(this(irow,jcol))
                            continue
                        elseif isnumeric(obj.func)
                            update_output([qij,'=',sprintf('%0.16g',obj.func),';'])
                        elseif isempty(obj.args)
                            update_output([qij,'=',correct_form(obj.func),';'])
                        else
                            loc=find(strcmp(list_in,qij));
                            if ~isempty(loc)
                                list_in(loc)=[];
                                continue
                            else
                                update_output([qij,'=',correct_form(obj.ref),';'])
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
        end
    end

    function update_output(in_string)
        itercount=itercount+1;
        if itercount==size(out,1),out{end+1000}=[];end
        out{itercount}=in_string;
    end

    function finalOutput=finalize(linecount,GG)
        xout=out(1:linecount);
        finalOutput=struct('code',cell2mat(xout(:)'),'argins',...
            {[validNames,'s0','s1']},'argouts',{GG});
    end
end