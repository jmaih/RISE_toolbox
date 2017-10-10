function [d,x,date_range]=collect_data(self,date_range)

ynames=self.endogenous;

xnames=self.exogenous;

yxnames=[ynames,xnames];

ng=max(1,numel(self.members));%self.kdata.ng;

if ng==1
    
    dx=collectdata(self.data,yxnames);
    
else
    
    dx=collectdata(self.data.(self.members{1}),yxnames);
    
    dx=dx(:,:,ones(1,ng));
    
    for ig=2:ng
        
        dx(:,:,ig)=collectdata(self.data.(self.members{ig}),yxnames);
        
    end
    
end

% split d and x

ny=numel(ynames);

d=dx(1:ny,:,:);

x=dx(ny+1:end,:,:);

add_constants()

if ng>1
    % do the senguesse
    %-----------------
    d=xpand_panel(d);
    
    x=xpand_panel(x);
    
end

    function d=xpand_panel(d)
        % going from nvar x T x npages to nvar*npages x T
        d=permute(d,[2,3,1]);
        
        d=d(:,:)';
        
    end

    function add_constants()
        
        sd=size(d,2);
        
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
