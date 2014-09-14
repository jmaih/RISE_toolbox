function v=load_start_values(names,db,date,v)
% output arguments
%   - v : vector of start values with the same number of elements as the
%   number of names
% input arguments
%   - date: string, integer(if annual dates), serial date
%   - db: ts object or struct with ts objects as fields
%   - names: string or cellstr
%   - v (optional): initial vector of starting values. It has to have the
%   same number of elements as names. if absent, it is initialized with zeros. 
% example
% db=ts('1990q1',rand(100,5),{'v1','v2','v3','v4','v5'})
% v=utils.forecast.load_start_values({'v4','v5','v5_AUX_L_8{+2}','v5_AUX_F_3'},db,'2011Q4')

if ischar(names)
    names=cellstr(names);
end

nv=numel(names);

if isa(db,'ts')
    db=pages2struct(db);
elseif ~isstruct(db)
    error('db must be a ts or a struct with ts fields')
end

% missing variables will start at zero
if nargin<4
    v=zeros(nv,1);
end
if numel(v)~=nv
    error('4th argument size does not match the number of names')
end
serial_date=date2serial(date);

dbnames=fieldnames(db);

for iv=1:nv
    vname1=names{iv};
    in_database=isfield(db,vname1);
    if in_database
        lead_lag=0;
        found_name=vname1;
    else
        [vname2,lead_lag,in_database]=get_name_and_lead_lag(vname1);
        in_database=in_database && lead_lag<0;
        if in_database
            found_name=vname2;
        end
    end
    if in_database
        data=db.(found_name);
        if ~isa(data,'ts')
            error(['field ',found_name,' in db must be a ts'])
        end
        target_date=serial_date+lead_lag;
        v(iv)=double(data(target_date));
    end
end

    function [vname_out,lead_lag,in_database]=get_name_and_lead_lag(vname_in)
        
        % check first for curly braces
        %-----------------------------
        [vname_out,lead_lag]=curly_braces_check(vname_in);
        
        in_database=any(strcmp(vname_out,dbnames));
        
        if ~in_database
            % if the variable is auxiliary, do further processing
            %----------------------------------------------------
            vname0=parser.translate_auxiliary_names(vname_out);
            if ~strcmp(vname0,vname_out)% then variable has curly braces
                [vname_out,lead_lag2]=curly_braces_check(vname0);
                lead_lag=lead_lag+lead_lag2;
                in_database=any(strcmp(vname_out,dbnames));
            end
        end
        
        function [name_out,lead_lag]=curly_braces_check(name_in)
            lead_lag=0;
            name_out=name_in;
            left_curly_brace=find(name_in=='{');
            if ~isempty(left_curly_brace)
                right_curly_brace=find(name_in=='}');
                if isempty(right_curly_brace)
                    error(['variable name ',name_in,' is incorrect'])
                end
                lead_lag=str2double(name_in(left_curly_brace+1:right_curly_brace-1));
                name_out=name_in(1:left_curly_brace-1);
            end
        end
    end

end

