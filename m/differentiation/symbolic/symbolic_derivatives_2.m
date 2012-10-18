function DerivativesPlan=symbolic_derivatives(equations,symb_list,...
    with_respect_to,order)
if nargin<4
    order=1;
end

% 1- replace leads and lags by something else

if ischar(equations)
    equations=cellstr(equations);
end
number_of_equations=numel(equations);
% remove possible semi-colons in each equation
for eq=1:number_of_equations
    equations{eq}=strrep(equations{eq},';','');
end
if ischar(symb_list)
    symb_list=cellstr(symb_list);
end
if ischar(with_respect_to)
    with_respect_to=cellstr(with_respect_to);
end

number_of_variables=numel(with_respect_to);
% construct a table with the occurrences of the variables to differentiate
equation_variable_map=false(number_of_equations,number_of_variables);
delimiters=[char([9:13,32]),'.(){};/*-+=^,[] '];
for eq=1:number_of_equations
    eqtn= equations{eq};
    while ~isempty(eqtn)
        [tok,eqtn]=strtok(eqtn,delimiters); %#ok<STTOK>
        loc=find(strcmp(tok,with_respect_to));
        if ~isempty(loc)
            equation_variable_map(eq,loc)=true;
        end
    end
end


% create symbolic expressions for all the elements in the dictionary
for ii=1:numel(symb_list)
    eval([symb_list{ii},'=sym(symb_list{ii});'])
end

% symbolize every equation
for eq=1:number_of_equations
    equations{eq}=sym(equations{eq});
end

DerivativesPlan=struct;
the_zeros=[sym('0'),sym('0.0')];
dimension_squeeze=number_of_equations==1;
for eq=1:number_of_equations
    % 0th order
    % Left is the number of times the equation is to be differentiated in
    % the next round
    Base0={1,[],equations{eq}};
    for oo=1:order
        if eq==1
            DerivativesPlan(oo).order=oo;
            if dimension_squeeze && oo>1
                DerivativesPlan(oo).size=number_of_variables*ones(1,oo);
            else
                DerivativesPlan(oo).size=[number_of_equations,...
                    number_of_variables*ones(1,oo)];
            end
            DerivativesPlan(oo).indices=[];
            DerivativesPlan(oo).derivs=cell(0,1);
            % those two lines will be filled in outside of this function
            % but they are initialized here.
            DerivativesPlan(oo).shadow=[];
            DerivativesPlan(oo).formatted=[];
        end
        derivoo=[];
        nbase=size(Base0,1);
        for bb=1:nbase
            indices=Base0{bb,2};
            iter=0;
            for ii=Base0{bb,1}:number_of_variables
                iter=iter+1;
                % differentiate only if the variable appears in the
                % equation in the first place. Ideally we should also check
                % this for all higher orders before taking the
                % derivs...
                if equation_variable_map(eq,ii) && ~ismember(Base0{bb,3},the_zeros)
                    derivative=diff(Base0{bb,3},with_respect_to{ii});
                    if ~ismember(derivative,the_zeros)
                        if dimension_squeeze
                            eq_index=[];
                        else
                            eq_index=eq;
                        end
                        derivoo=[derivoo;{Base0{bb,1}+iter-1,[indices,ii],derivative}]; %#ok<AGROW>
                        DerivativesPlan(oo).indices=[DerivativesPlan(oo).indices
                            eq_index,indices,ii];
                        DerivativesPlan(oo).derivs=[DerivativesPlan(oo).derivs;...
                            {char(derivative)}];
                    end
                end
            end
        end
        Base0=derivoo;
    end
end

