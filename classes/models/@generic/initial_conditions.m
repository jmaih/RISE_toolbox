function db=initial_conditions(obj,endhist_date,growth_type,ncond,~)

if isempty(obj)
    
    db=cell(0,4);
    
    return
    
end

if nargin<4
    
    ncond=0;
    
end

end_sample_date=endhist_date;

db=growth_database(obj,endhist_date,end_sample_date,growth_type);

if ncond
    
    fields=fieldnames(db);
    
    for ifield=1:numel(fields)
        
        name=fields{ifield};
        
        data=db.(name).data;
        
        if ifield==1
            
            [nr,nc]=size(data);
            
            proto=cat(3,data,nan(nr,nc,ncond));
            
            proto_ts=ts(db.(name).start,proto);
            
        end
        
        proto(:,:,1)=data;
        
        db.(name)=reset_data(proto_ts,proto);
        
    end
    
end

end