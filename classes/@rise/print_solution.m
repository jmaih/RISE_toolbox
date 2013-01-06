function print_solution(obj,varlist,file2save2,equation_format)

if nargin<4
    equation_format=[];
    if nargin<3
        file2save2=[];
        if nargin<2
            varlist=[];
        end
    end
end
if isempty(equation_format)
    equation_format=false;
end
if isempty(varlist)
    varlist=char({obj(1).varendo.name});
end
precision='%8.6f';

% get the location of the variables
ids=locate_variables(varlist,{obj(1).varendo.name});

if ~isempty(file2save2)
    fid=fopen(file2save2,'w');
else
    fid=1;
end

nobj=numel(obj);
string='';
for kk=1:numel(obj)
    if nobj>1
        string=int2str(kk);
    end
    fprintf(fid,'\n%s\n',['MODEL ',string,' SOLUTION']);
    
    T=obj(kk).T;
    if isempty(T)
        fprintf(fid,'\n%s\n','MODEL HAS NOT BEEN SOLVED');
    else
        R=obj(kk).R;
        SS=obj(kk).steady_state_and_balanced_growth_path(1:obj(kk).NumberOfEndogenous(2),:);
        BGP=obj(kk).steady_state_and_balanced_growth_path(obj(kk).NumberOfEndogenous(2)+1:end,:);
        Risk=vertcat(obj(kk).varendo.risk);
        solver=obj(kk).options.solver;
        if isnumeric(solver)
            solver=int2str(solver);
        end
        fprintf(fid,'%s :: %s\n','SOLVER',solver);
        endo_names={obj(kk).varendo.name};
        exo_names={obj(kk).varexo.name};
        [~,~,order,h]=size(R);
        T(abs(T)<1e-9)=0;
        R(abs(R)<1e-9)=0;
        SS(abs(SS)<1e-9)=0;
        BGP(abs(BGP)<1e-9)=0;
        SS=SS(ids,:);
        BGP=BGP(ids,:);
        Risk(abs(Risk)<1e-9)=0;
        if order>1
            fprintf(fid,'\n%s\n','PRINTING RESULTS ONLY FOR ORDER 1');
        end
        for ii=1:h
            % first select the matrices corresponding to the variables of
            % interest
            Tii=T(ids,:,ii);
            Rii=R(ids,:,1,ii);
            Risk_ii=Risk(ids,ii);
            % then keep only the state variables in T. It is important to call
            % function "any" with the dimension in order to avoid problems if
            % we are reporting results for only one variable.
            state_cols=any(Tii,1);
            Tii=Tii(:,state_cols);
            state_vars={};
            if any(state_cols)
                state_vars=strcat(endo_names(state_cols),'{-1}');
            end
            % do the same with the shocks
            shock_cols=any(Rii,1);
            Rii=Rii(:,shock_cols);
            shock_names=exo_names(shock_cols);
            var_names=endo_names(ids);
            
            GrandState=[{'steady state','risk','bal. growth'},state_vars,shock_names];
            GrandMat=[SS(:,ii),Risk_ii,BGP(:,ii),Tii,Rii];
            % trim some more
            mat_cols=any(GrandMat,1);
            GrandMat=GrandMat(:,mat_cols);
            GrandState=GrandState(mat_cols);
            
            data=[{'Endo Name'},var_names
                GrandState',num2cell(GrandMat')];
            %             data=cell(2,0);
            %             % first columns
            %             data=[data,{'Endo Name',[{'steady state'},{'risk'},{'bal. growth'},state_vars,shock_names]'}']; %#ok<AGROW>
            %             for vv=1:numel(var_names)
            %                 data=[data,{var_names{vv},[SS(vv,ii);Risk_ii(vv);BGP(vv,ii);Tii(vv,:)';Rii(vv,:)']}']; %#ok<AGROW>
            %             end
            
            B=concatenate(data,precision);
            body_format='\n';
            % start from the end
            for bb=size(B,2):-1:1
                body_format=['%',int2str(B{2,bb}),'s ',body_format]; %#ok<AGROW>
            end
            nrows=size(B{1,1},1);
            number_of_headers=size(B,2);
            fprintf(fid,'\n%s %4.0f\n','Regime',ii);
            if equation_format
                B0=B{1,1};
                for icols=2:number_of_headers
                    Bi=B{1,icols};
                    tmp=[Bi(1,:),'='];
                    for irows=2:size(Bi,1)
                        tmp=[tmp,Bi(irows,:),B0(irows,:),'+']; %#ok<AGROW>
                    end
                    tmp(isspace(tmp))=[];
                    tmp=strrep(tmp,'+-','-');
                    tmp=strrep(tmp,'-+','-');
                    tmp=strrep(tmp,'steadystate','');
                    tmp=tmp(1:end-1);
                    fprintf(fid,'%s \n',tmp);
                end
            else
                for rr=1:nrows
                    data_ii=cell(1,number_of_headers);
                    for jj=1:number_of_headers
                        data_ii{jj}=B{1,jj}(rr,:);
                    end
                    fprintf(fid,body_format,data_ii{:});
                end
            end
        end
    end
    
end
if ~isempty(file2save2)
    fclose(fid);
end

function B=concatenate(data,precision)

B=cell(2,0);
[span,nargs]=size(data);

for ii=1:nargs
    A='';
    for jrow=1:span
        add_on=data{jrow,ii};
        if isnumeric(add_on)
            add_on=num2str(add_on,precision);
        end
        A=char(A,add_on);
    end
    A=A(2:end,:);
    B=[B,{A,size(A,2)+2}']; %#ok<AGROW>
end

% function B=concatenate(data,precision)
%
% B=cell(2,0);
% nargs=size(data,2);
%
% for ii=1:nargs
%     if ii==1
%         span=numel(data{2,1});
%     end
%     A=data{1,ii};
%     add_on=data{2,ii};
%     if iscellstr(add_on)
%         add_on=char(add_on);
%     end
%     if isempty(add_on)
%         A=char(A,repmat('--',span,1));
%     elseif isnumeric(add_on)
%         for jj=1:span
%             A=char(A,num2str(add_on(jj),precision));
%         end
%     elseif ischar(A)
%         A=char(A,add_on);
%     else
%         error([mfilename,':: unknown type'])
%     end
%     B=[B,{A,size(A,2)+2}']; %#ok<AGROW>
%
% end