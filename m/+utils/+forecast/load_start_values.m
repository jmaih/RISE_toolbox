function v=load_start_values(names,db,date,v)
% load_start_values - load the start values for forecasting
%
% ::
%
%
%   v=load_start_values(names,db,date)
%   v=load_start_values(names,db,date,v)
%
% Args:
%
%    - **names** [char|cellstr]: names of variables to load from the database
%
%    - **db** [ts|struct]: ts time series data
%
%    - **date** [char|numeric|serial date]: date for which the data are
%      requested
%
%    - **v** (optional)[vector,3-dimensional array|[]]: start values. If
%      provided, it should have the same number of rows as the number of names
%
% Returns:
%    :
%
%    - **v** [vector,3-dimensional array|[]]: start values
%
% Note:
%
%    - leads are ignored
%
%    - if a variable is not found in the database, it is initialized at 0 in
%      the first page and nan in the subsequent pages
%
% Example:
%
%    db=ts('1990q1',rand(100,5,3),{'v1','v2','v3','v4','v5'})
%    v=utils.forecast.load_start_values({'v4','v5','v5_AUX_L_8{+2}','v5_AUX_F_3'},db,'2011Q4')
%
%    See also: data_request

if ischar(names)
    names=cellstr(names);
end

nv=numel(names);
% missing variables will start at zero
if nargin<4
    v=[];
end
if isempty(v)
    v=zeros(nv,1);
end
if ~isa(db,'ts')&&~isstruct(db)
    error('db must be a ts or a struct with ts fields')
end

% first collect the data so as to make the pages consistent
%-----------------------------------------------------------
db=ts.collect(db);
npages=get(db,'NumberOfPages');

if numel(v)~=nv
    error('4th argument size does not match the number of names')
end
serial_date=date2serial(date);

dbnames=get(db,'varnames');

% the database may have several pages : conditioning information
% note that going from 
%----------------------------------------------------------------
v=full(v); % avoid problems with sparse
if npages>1
    if size(v,3)>1
        error('it is not allowed to have initial conditions with more than one page when the database has more than one page')
    end
    % only pick the first page of v in case v has many pages
    v=cat(3,v(:,1,1),nan(nv,1,npages-1));
end

for iv=1:nv
    vname1=names{iv};
    in_database=any(strcmp(vname1,dbnames));
    if in_database
        lead_lag=0;
        found_name=vname1;
    else
        [vname2,lead_lag,in_database]=get_name_and_lead_lag(vname1);
        % ignore the leads
        %------------------
        in_database=in_database && lead_lag<0;
        if in_database
            found_name=vname2;
        end
    end
    if in_database
        data=db(found_name);
        if ~isa(data,'ts')
            error(['field ',found_name,' in db must be a ts'])
        end
        target_date=serial_date+lead_lag;
        if all(ismember(target_date,data.date_numbers))
            tmp=double(data(target_date));
            v(iv,1,:)=permute(tmp,[2,1,3]);
        end
    end
end

% remove the nan pages starting from the end until we hit a page that does
% not contains nans
%--------------------------------------------------------------------------
for ipage=npages:-1:2
    vlast=v(:,:,end);
    if all(isnan(vlast))
        v(:,:,end)=[];
    else
        break
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

