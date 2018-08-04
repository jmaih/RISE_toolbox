function this=collect(varargin)
% collect - brings together several time series objects into a one time series
%
% ::
%
%
%   - this=collect(v1,v2,...,vn)
%
% Args:
%
%    - **Vi** [cell|ts|struct]: time series in ts format:
%      - cell: When Vi is a cell, then its format should be {vname,ts} i.e.
%      the first element is the name of the variable and the second is the
%      data for the variable. In this case, the data must be a single time
%      series
%      - ts:
%      - struct: the fields of the structure should be of the ts format.
%
% Returns:
%    :
%
%    - **this** [ts]: a time series with many columns and potentially many
%      pages
%
% Note:
%
% Example:
%
%    See also:

nn=length(varargin);

exitflag= nn==0||isempty(varargin{1})||...
    (nn==1 && isa(varargin{1},'ts'));

if exitflag
    
    if nn==0
        
        this=ts.empty(0);
        
    else
        
        this=varargin{1};
        
    end
    
    return
    
end

cellnames=cell(1,0);

celldata=cell(1,0);

celldescription=cell(1,0);

ii=0;

offset=0;

while ~isempty(varargin)
    
    ii=ii+1;
    
    db_i=varargin{1};
    
    if isa(db_i,'cell')
        
        if ~isa(db_i{2},'ts')
            
            error('element in the second cell must be an input')
            
        end
        
        if db_i{2}.NumberOfVariables>1
            
            error('when input is a cell, the number of variables should be 1')
        
        end
        
        cellnames=[cellnames,db_i(1)];
        
        celldata=[celldata,db_i(2)];
        
        celldescription=[celldescription,db_i{2}.description];
        
        offset=offset+1;
        
    elseif isa(db_i,'ts')
        % if there are many names then check that nargin==1
        if isempty(db_i.varnames{1})
            
            error([mfilename,':: database ',int2str(ii),' must have an internal name'])
        
        end
        
        cellnames=[cellnames,db_i.varnames(:).'];
        
        celldata=[celldata,cell(1,db_i.NumberOfVariables)];
        
        celldescription=[celldescription,db_i.description];
        
        if db_i.NumberOfVariables==1
            
            offset=offset+1;
            
            celldata{offset}=db_i;
        
        else
            
            for ivar=1:db_i.NumberOfVariables
                
                offset=offset+1;
                % calling subsref from within a method does not work. I has
                % to be done explicitly using the functional form
                % celldata{offset}=db_i(db_i.varnames{ivar});
                S=struct('type','()','subs',{db_i.varnames(ivar)});
                
                celldata{offset}=subsref(db_i,S);
            
            end
            
        end
        
    elseif isa(db_i,'struct')
        
        fields=fieldnames(db_i);
        
        n0=numel(fields);
        
        cellnames=[cellnames,fields(:).']; %#ok<*AGROW>
        
        celldata=[celldata,cell(1,n0)];
        
        for jj=1:n0
            
            offset=offset+1;
            
            celldata{offset}=db_i.(fields{jj});
            
            newdescrip=db_i.(fields{jj}).description;
            
            if isempty(newdescrip)
                
                newdescrip={''};
                
            end
            
            if ischar(newdescrip)
                
                newdescrip=cellstr(newdescrip);
                
            end
            
            celldescription=[celldescription,newdescrip(:).'];
            
        end
        
    else
        
        error([mfilename,':: input ',int2str(ii),' must be either instance of class ''ts'' or a cell as {''varname'',ts object} or a struct of ts objects'])
    
    end
    
    varargin=varargin(2:end);

end

unique_names=unique(cellnames);

n_unic=numel(unique_names);

duplicates=false(1,n_unic);

locs_unic=cell(1,n_unic);

for ii=1:numel(unique_names)
    
    locs_unic{ii}=find(strcmp(unique_names{ii},cellnames));
    
    if numel(locs_unic{ii})>1
        
        duplicates(ii)=true;
        
    end
    
end

if any(duplicates)
    
    disp(unique_names(duplicates))
    
    warning([mfilename,':: the variables above are duplicated'])
    
end

% find the first date, the last date, the highest frequency
first_date=celldata{1}.date_numbers(1);

last_date=celldata{1}.date_numbers(end);

frequencies_codes={'','H','Q','M'
    1,2,3,4};

loc= strcmp(celldata{1}.frequency,frequencies_codes(1,:));

frequency_code=frequencies_codes{2,loc};

highest_frequency=celldata{1}.frequency;

npages=nan(1,offset);

% also check the number of variables
for ii=1:offset
    % redo the time series if the cell contains variable "with many variables"
    if celldata{ii}.NumberOfVariables>1
        
        vnames=celldata{ii}.varnames;
        
        mtee=cellfun(@isempty,vnames,'uniformoutput',false);
        
        mtee=[mtee{:}];
        
        if ~all(mtee)
            
            error('sub-variables are not expected to have names when using collect')
        
        end
        
        data_=double(celldata{ii});
        
        celldata{ii}=ts(celldata{ii}.start,permute(data_,[1,3,2]),...
            celldata{ii}.description);
    
    end
    
    if ii>1
        
        loc= strcmp(celldata{ii}.frequency,frequencies_codes(1,:));
        
        if frequencies_codes{2,loc}>frequency_code
            
            highest_frequency=celldata{ii}.frequency;
            
            frequency_code=frequencies_codes{2,loc};
            
            first_date=convert_date(first_date,highest_frequency);
            
            last_date=convert_date(last_date,highest_frequency);
            
        end
        
        first_last=convert_date(celldata{ii}.date_numbers([1,end]),highest_frequency);
        
        if first_last(1)<first_date
            
            first_date=first_last(1);
            
        end
        
        if first_last(end)>last_date
            
            last_date=first_last(end);
            
        end
        
    end
    
    npages(ii)=celldata{ii}.NumberOfPages;
    
end

newdates=(first_date:last_date)';

newdata=nan(numel(newdates),n_unic,max(npages));

for kk=1:n_unic
    
    locs=locs_unic{kk};
    
    for iloc=1:numel(locs)
        
        ii=locs(iloc);
        
        % step 1, convert the dates to the new frequency
        dat_n_ii=convert_date(celldata{ii}.date_numbers,highest_frequency);
        
        oldata=celldata{ii}.data;
        
        locations=nan(size(dat_n_ii));
        
        for jj=1:numel(locations)
            
            loc=find(abs(newdates-dat_n_ii(jj))<1e-9);% <--loc=find(newdates==dat_n_ii(jj));
            
            if isempty(loc)
                
                error([mfilename,':: date not found in new vector'])
            
            end
            
            locations(jj)=loc;
        
        end
        
        template=newdata(locations,kk,1:npages(ii));
        
        % override the template only in places that are nan
        free=isnan(template);
        
        template(free)=oldata(free);
        
        if any(template(~free)) && any(oldata(~free))
            
            warning('due to duplications, some data were not copied...')
        
        end
        % find the date locations
        newdata(locations,kk,1:npages(ii))=template;
    
    end
    
end
% trim the dates and data in case there are further trailing
% nans. Do this using the fake data in case there are many pages
fakedata=newdata(:,:);

first_good=find(any(~isnan(fakedata),2),1,'first');

last_good=find(any(~isnan(fakedata),2),1,'last');

newdates=newdates(first_good:last_good);

newdata=newdata(first_good:last_good,:,:);
% I could always build a special type of ts instead of
% calling ts. But the reason I need ts is because I
% started working with it earlier

n_names=numel(cellnames);

discard=false(1,n_names);

for ii=1:numel(unique_names)
    
    v=unique_names{ii};
    
    pos=find(strcmp(v,cellnames));
    
    discard(pos(2:end))=true;
    
end

this=ts(newdates,newdata,unique_names,celldescription(~discard));

end

function this=convert_date(this0,newfreq)

dec=serial2dec(this0);

oldfreq=dec(1).frequency;

year=[dec.year];

period=[dec.period];

freqtable={'','H','Q','M'
    1,2,3,4};

oldfreq_id=find(strcmp(oldfreq,freqtable(1,:)));

if isempty(oldfreq_id)
    
    error([mfilename,':: unknown frequency ''',oldfreq,''])
    
end

newfreq_id=find(strcmp(newfreq,freqtable(1,:)));

if isempty(newfreq_id)
    
    error([mfilename,':: unknown frequency ''',newfreq,''])
    
end

if oldfreq_id==newfreq_id
    
    this=this0;
    
else
    
    if oldfreq_id>newfreq_id
        
        error([mfilename,':: date conversion from high to low frequency is undefined'])
        
    end
    
    if oldfreq_id==1 % Annual
        
        if newfreq_id==2 % 'H'
            
            period=2;
            
        elseif newfreq_id==3 % 'Q'
            
            period=4;
            
        elseif newfreq_id==4% 'M'
            
            period=12;
            
        end
        
    elseif oldfreq_id==2 % H
        
        if newfreq_id==3% 'Q'
            
            period=2*period;
            
        elseif newfreq_id==4% 'M'
            
            period=6*period;
            
        end
        
    elseif oldfreq_id==3 % Q
        
        if newfreq_id==4% 'M'
            
            period=3*period;
            
        end
        
    elseif oldfreq_id==4 % M
        % this case is not expected to happen
    end
    
    freq=frequency2num(newfreq);
    
    this=dec2serial(year,period,freq);
    
end

end
