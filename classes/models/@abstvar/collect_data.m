function [y,x,date_range,regimes]=collect_data(self,date_range,is_fixed_regime)
% INTERNAL FUNCTION: Parse data into a useable format
%
% Note:
%    It is assumed that data is set in self.estim_.data before the execution of this function
%

if nargin<3
    
    is_fixed_regime=[];
    
end

if isempty(is_fixed_regime)
    
    is_fixed_regime=false;
    
end

ynames=self.endogenous;

xnames=self.exogenous;

if is_fixed_regime
    
    xnames=[xnames(:).',{'hist_regimes'}];
    
end

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

regimes=[];

if is_fixed_regime
    
    regimes=x(end,:,:);
    
    x=x(1:end-1,:,:);
    
    for ig=2:ng
        
        if any(regimes(1,:,ig)-regimes(1,:,1))
            
            error('all regimes across panel members should be the same')
            
        end
        
    end
    
    regimes=regimes(1,:,1);
    
end

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
