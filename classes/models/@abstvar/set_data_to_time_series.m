function dout=set_data_to_time_series(self,din,vnames,start_date)

dout=struct();

if isstruct(din)
    
    fnames=fieldnames(din);
    
    ng=numel(fnames);
    
    for ii=1:ng
        
        dout.(fnames{ii})=set_data_to_time_series(self,din.(fnames{ii}),vnames,start_date);
        
    end
    
    return
    
end

rtmp=struct();

for iv=1:numel(vnames)
    
    d=permute(din(iv,:,:,:),[2,3,4,1]);
    
    if iv==1
        
        proto=ts(start_date,d);
        
    end
    
    rtmp.(vnames{iv})=reset_data(proto,d);
    
end

dout=rtmp;

end

