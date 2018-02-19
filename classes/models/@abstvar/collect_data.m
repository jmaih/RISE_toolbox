function [y,x,date_range]=collect_data(self,date_range)

ynames=self.endogenous;

xnames=self.exogenous;

yxnames=[ynames,xnames];

ng=max(1,numel(self.members));%self.kdata.ng;

if ng==1
    
    yx=collectdata(self.estim_.data,yxnames);
    
else
    
    yx=collectdata(self.estim_.data.(self.members{1}),yxnames);
    
    yx=yx(:,:,ones(1,ng));
    
    for ig=2:ng
        
        yx(:,:,ig)=collectdata(self.estim_.data.(self.members{ig}),yxnames);
        
    end
    
end

% split d and x
yloc=locate_variables(ynames,yxnames);

xloc=locate_variables(xnames,yxnames);

y=yx(yloc,:,:);

x=yx(xloc,:,:);

add_constants()

if ng>1
    % do the senguesse
    %-----------------
    y=vartools.xpand_panel(y);
    
    x=vartools.xpand_panel(x);
    
end

    function add_constants()
        
        sd=size(yx,2);
        
        if self.constant
            
            wan=ones(1,sd,ng);
            
            x=cat(1,wan,x);
            
        end
        
    end

    function d=collectdata(data,vnames)
        
        nvars_=numel(vnames);
        
        if nvars_==0
            
            d=[];
            
            return
            
        end
        
        tmp=ts.collect(data);
        
        d=double(tmp);
        
        where=locate_variables(vnames,tmp.varnames);
        
        dn=tmp.date_numbers;
        
        date_range=abstvar.reset_range(date_range,dn);
        
        start=find(dn==date_range(1));
        
        finish=find(dn==date_range(2));
        
        d=permute(d(start:finish,where,:,:),[2,1,3,4]);
                
    end

end
