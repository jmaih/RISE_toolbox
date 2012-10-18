function eqn_obj=symbolic_derivator_up_to_3rd_order(endo_names,exo_names,par_names,lead_lag_incidence,equations,order)
% the dictionary contains all the names of the parameters, the endogenous,
% the exogenous

% 1- replace leads and lags by something else

if ischar(equations)
    equations=cellstr(equations);
end
number_of_equations=numel(equations);
% remove possible semi-colons in each equation
for eq=1:number_of_equations
    equations{eq}=strrep(equations{eq},';','');
end

if ischar(endo_names)
    endo_names=cellstr(endo_names);
end
if ischar(exo_names)
    exo_names=cellstr(exo_names);
end
if ischar(par_names)
    par_names=cellstr(par_names);
end
dictionary=union(union(endo_names,exo_names),par_names);
dictionary=dictionary(:);

% augment the dictionary with lags and leads using information in the
% lead_lag_incidence matrix. Replace (-1) with __L and (+1) with __F.
lagged_vars_id=find(lead_lag_incidence(:,1));
nlags=numel(lagged_vars_id);
auxiliary_lagged_vars=cell(nlags,1);
old_lagged_vars=cell(nlags,1);
for ii=1:nlags
    oldsubstr1=strcat(endo_names{lagged_vars_id(ii)},'(-1)');
    oldsubstr2=strcat(endo_names{lagged_vars_id(ii)},'{-1}');
    newsubstr=strcat(endo_names{lagged_vars_id(ii)},'__L');
    auxiliary_lagged_vars{ii}=newsubstr;
    old_lagged_vars{ii}=oldsubstr2;
    % replace the equations simultaneously
    for eq=1:number_of_equations
        equations{eq}=strrep(equations{eq},oldsubstr1,newsubstr);
        equations{eq}=strrep(equations{eq},oldsubstr2,newsubstr);
    end
end
lead_vars_id=find(lead_lag_incidence(:,3));
nleads=numel(lead_vars_id);
auxiliary_lead_vars=cell(nleads,1);
old_lead_vars=cell(nleads,1);
for ii=1:nleads
    oldsubstr1=strcat(endo_names{lead_vars_id(ii)},'(+1)');
    oldsubstr2=strcat(endo_names{lead_vars_id(ii)},'{+1}');
    oldsubstr3=strcat(endo_names{lead_vars_id(ii)},'(1)');
    oldsubstr4=strcat(endo_names{lead_vars_id(ii)},'{1}');
    newsubstr=strcat(endo_names{lead_vars_id(ii)},'__F');
    auxiliary_lead_vars{ii}=newsubstr;
    old_lead_vars{ii}=oldsubstr4;
    % replace the equations simultaneously
    for eq=1:number_of_equations
        equations{eq}=strrep(equations{eq},oldsubstr1,newsubstr);
        equations{eq}=strrep(equations{eq},oldsubstr2,newsubstr);
        equations{eq}=strrep(equations{eq},oldsubstr3,newsubstr);
        equations{eq}=strrep(equations{eq},oldsubstr4,newsubstr);
    end
end
% add the auxiliary variables to the dictionary
dictionary=[auxiliary_lagged_vars;dictionary;auxiliary_lead_vars];
% list the variables to differentiate with respect to, ordering them
% according to lead_lag_incidence and adding exogenous at the end
with_respect_to=[auxiliary_lagged_vars
    endo_names(:)
    auxiliary_lead_vars
    exo_names(:)];

auxiliary_vars=[auxiliary_lagged_vars
    auxiliary_lead_vars];
old_vars=[old_lagged_vars
    old_lead_vars];

% create symbolic expressions for all the elements in the dictionary
for ii=1:numel(dictionary)
    eval([dictionary{ii},'=sym(dictionary{ii});'])
end

% symbolize every equation
for eq=1:number_of_equations
    equations{eq}=sym(equations{eq});
end

% create equation objects
number_of_variables=numel(with_respect_to);
the_zero=sym('0');


% take one equation, find all the derivates and put them into an equation
% object
% % % % % % % eq=0;
% % % % % % % eq_nbr=numel(equations);
% % % % % % % eqn_obj=rise_equation.empty(0);
% % % % % % % Crixus='';
% % % % % % % while eq<eq_nbr
% % % % % % %     eq=eq+1;
% % % % % % %     Crixus=char(Crixus,['%equation ',int2str(eq)]);
% % % % % % %     eqn_obj(eq,1)=rise_equation(equations{eq},with_respect_to,par_names,eq);
% % % % % % %     M1=repmat({the_zero},1,number_of_variables);
% % % % % % %     M1_=repmat({'0'},1,number_of_variables);
% % % % % % %     Crixus=char(Crixus,['M1=zeros(1,',int2str(number_of_variables),');']);
% % % % % % %     
% % % % % % %     M2_={};
% % % % % % %     M3_={};
% % % % % % %     if order>1
% % % % % % %         M2=repmat({the_zero},number_of_variables*ones(1,2));
% % % % % % %         M2_=repmat({'0'},number_of_variables*ones(1,2));
% % % % % % %         Crixus=char(Crixus,['M2=zeros(',int2str(number_of_variables),'*ones(1,2));']);
% % % % % % %         if order >2
% % % % % % %             M3=repmat({the_zero},number_of_variables*ones(1,3));
% % % % % % %             M3_=repmat({'0'},number_of_variables*ones(1,3));
% % % % % % %             Crixus=char(Crixus,['M3=zeros(',int2str(number_of_variables),'*ones(1,3));']);
% % % % % % %         end
% % % % % % %     end
% % % % % % %     
% % % % % % %     if ~isequal(equations{eq},the_zero)
% % % % % % %         for i1=1:number_of_variables
% % % % % % %             M1{i1}=diff(equations{eq},with_respect_to{i1});
% % % % % % %             M1_{i1}=replace_back(char(M1{i1}),auxiliary_vars,old_vars);
% % % % % % %             if ~isequal(M1{i1},the_zero)
% % % % % % %                 Crixus=char(Crixus,['M1(',int2str(i1),')=',M1_{i1},';']);
% % % % % % %                 if order>1
% % % % % % %                     for i2=1:i1
% % % % % % %                         M2{i1,i2}=diff(M1{i1},with_respect_to{i2});
% % % % % % %                         if ~isequal(M2{i1,i2},the_zero)
% % % % % % %                             M2_{i1,i2}=replace_back(char(M2{i1,i2}),auxiliary_vars,old_vars);
% % % % % % %                             Crixus=char(Crixus,['M2(',int2str(i1),',',int2str(i2),')=',M2_{i1,i2},';']);
% % % % % % %                             if i2<i1
% % % % % % %                                 M2_{i2,i1}=M2_{i1,i2};
% % % % % % %                                 Crixus=char(Crixus,['M2(',int2str(i2),',',int2str(i1),')=M2(',int2str(i1),',',int2str(i2),');']);
% % % % % % %                             end
% % % % % % %                             if order>2
% % % % % % %                                 for i3=1:i2
% % % % % % %                                     M3{i1,i2,i3}=diff(M2{i1,i2},with_respect_to{i3});
% % % % % % %                                     if ~isequal(M3{i1,i2,i3},the_zero)
% % % % % % %                                         M3_{i1,i2,i3}=replace_back(char(M3{i1,i2,i3}),auxiliary_vars,old_vars);
% % % % % % %                                         Crixus=char(Crixus,['M3(',int2str(i1),',',int2str(i2),...
% % % % % % %                                             ',',int2str(i3),')=',M3_{i1,i2,i3},';']);
% % % % % % %                                         if i1<i2 && i2<i3
% % % % % % %                                             M3_{i1,i3,i2}=M3_{i1,i2,i3};
% % % % % % %                                             Crixus=char(Crixus,['M3(',int2str(i1),',',int2str(i3),...
% % % % % % %                                                 ',',int2str(i2),')=M3(',int2str(i1),',',int2str(i2),',',int2str(i3),');']);
% % % % % % %                                             M3_{i2,i1,i3}=M3_{i1,i2,i3};
% % % % % % %                                             Crixus=char(Crixus,['M3(',int2str(i2),',',int2str(i1),...
% % % % % % %                                                 ',',int2str(i3),')=M3(',int2str(i1),',',int2str(i2),',',int2str(i3),');']);
% % % % % % %                                             M3_{i2,i3,i1}=M3_{i1,i2,i3};
% % % % % % %                                             Crixus=char(Crixus,['M3(',int2str(i2),',',int2str(i3),...
% % % % % % %                                                 ',',int2str(i1),')=M3(',int2str(i1),',',int2str(i2),',',int2str(i3),');']);
% % % % % % %                                             M3_{i3,i1,i2}=M3_{i1,i2,i3};
% % % % % % %                                             Crixus=char(Crixus,['M3(',int2str(i3),',',int2str(i1),...
% % % % % % %                                                 ',',int2str(i2),')=M3(',int2str(i1),',',int2str(i2),',',int2str(i3),');']);
% % % % % % %                                             M3_{i3,i2,i1}=M3_{i1,i2,i3};
% % % % % % %                                             Crixus=char(Crixus,['M3(',int2str(i3),',',int2str(i2),...
% % % % % % %                                                 ',',int2str(i1),')=M3(',int2str(i1),',',int2str(i2),',',int2str(i3),');']);
% % % % % % %                                         end
% % % % % % %                                     end
% % % % % % %                                 end
% % % % % % %                             end
% % % % % % %                         end
% % % % % % %                     end
% % % % % % %                 end
% % % % % % %             end
% % % % % % %         end
% % % % % % %     end
% % % % % % %     eqn_obj(eq,1)=eqn_obj(eq,1).set_properties('first_order_derivatives',M1_,...
% % % % % % %         'second_order_derivatives',M2_,...
% % % % % % %         'third_order_derivatives',M3_);
% % % % % % % end

eq=0;
eq_nbr=numel(equations);
eqn_obj=rise_equation.empty(0);
nvar_str=int2str(number_of_variables);
Crixus=['M1=zeros(',int2str(eq_nbr),',',nvar_str,');'];
if order>1
    Crixus=char(Crixus,['M2=zeros([',int2str(eq_nbr),',',nvar_str,'*ones(1,2)]);']);
    if order>2
        Crixus=char(Crixus,['M3=zeros([',int2str(eq_nbr),',',nvar_str,'*ones(1,3)]);']);
    end
end
while eq<eq_nbr
    eq=eq+1;
    if ~isequal(equations{eq},the_zero)
        for ii=1:number_of_variables
            M_i=diff(equations{eq},with_respect_to{ii});
            M_i_store=char(M_i);
%             M_i_store=replace_back(char(M_i),auxiliary_vars,old_vars);
            if ~isequal(M_i,the_zero)
                Crixus=char(Crixus,['M1(',int2str(eq),',',int2str(ii),')=',M_i_store,';']);
                if order>1
                    for jj=1:ii
                        M_i_j=diff(M_i,with_respect_to{jj});
                        if ~isequal(M_i_j,the_zero)
                            M_i_j_store=char(M_i_j);
%                             M_i_j_store=replace_back(char(M_i_j),auxiliary_vars,old_vars);
                            M2_str=['M2(',int2str(eq),',',int2str(ii),',',int2str(jj),')'];
                            Crixus=char(Crixus,[M2_str,'=',M_i_j_store,';']);
                            if jj<ii
                                Crixus=char(Crixus,['M2(',int2str(eq),',',int2str(jj),',',int2str(ii),')=',M2_str,';']);
                            end
                            if order>2
                                for kk=1:jj
                                    M_i_j_k=diff(M_i_j,with_respect_to{kk});
                                    if ~isequal(M_i_j_k,the_zero)
                                        M_i_j_k_store=char(M_i_j_k);
%                                         M_i_j_k_store=replace_back(char(M_i_j_k),auxiliary_vars,old_vars);
                                        M3_str=['M3(',int2str(eq),',',int2str(ii),',',int2str(jj),',',int2str(kk),')'];
                                        Crixus=char(Crixus,[M3_str,'=',M_i_j_k_store,';']);
                                        if ii<jj && jj<kk
                                            Crixus=char(Crixus,['M3(',int2str(eq),',',int2str(ii),',',int2str(kk),...
                                                ',',int2str(jj),')=',M3_str,';']);
                                            Crixus=char(Crixus,['M3(',int2str(eq),',',int2str(jj),',',int2str(ii),...
                                                ',',int2str(kk),')=',M3_str,';']);
                                            Crixus=char(Crixus,['M3(',int2str(eq),',',int2str(jj),',',int2str(kk),...
                                                ',',int2str(ii),')=',M3_str,';']);
                                            Crixus=char(Crixus,['M3(',int2str(eq),',',int2str(kk),',',int2str(ii),...
                                                ',',int2str(jj),')=',M3_str,';']);
                                            Crixus=char(Crixus,['M3(',int2str(eq),',',int2str(kk),',',int2str(jj),...
                                                ',',int2str(ii),')=',M3_str,';']);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
disp(Crixus)

keyboard

function outstr=replace_back(instr,auxiliary_vars,old_vars)
outstr=instr;
if ~strcmp(outstr,'0')
    for aa=1:numel(auxiliary_vars)
        outstr=strrep(outstr,auxiliary_vars{aa},old_vars{aa});
    end
end
outstr(isspace(outstr))=[];

