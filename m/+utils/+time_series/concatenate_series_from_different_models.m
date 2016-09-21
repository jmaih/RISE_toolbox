function output=concatenate_series_from_different_models(dbcell)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

% dbcell is a cell array of structures of time series objects
% all elements should have the same dimensions else there will be a crash
% returns a structure whose fields are the concatenated time series objects

[varList,ncols,dn_start,dn_last]=do_intersection(dbcell);

output=struct();

dn0=dn_start:dn_last;

nrows=numel(dn0);

nobj=numel(dbcell);

tank=nan(nrows,nobj,ncols);

mod_names=strcat('model_',cellfun(@num2str,num2cell(1:nobj),'uniformOutput',false));

reg_names=strcat('regime_',cellfun(@num2str,num2cell(1:ncols),'uniformOutput',false));

for ivar=1:numel(varList)
    
    for imod=1:numel(dbcell)
        
        this=dbcell{imod}.(varList{ivar});
        
        dn=this.date_numbers;
        
        datta=double(this);
        
        ncols_imod=size(datta,2);
        
        if ncols_imod==1 && ncols > 1
            
            ncols_imod=ncols;
            
            datta=datta(:,ones(1,ncols));
            
        end
        
        % better start search from the beginning
        start=find(dn(1)==dn0,1,'first');
        
        % better start search from the end
        finish=find(dn(end)==dn0,1,'last');
        
        tank(start:finish,imod,1:ncols_imod)=datta;
        
    end
    
    if ncols > 1
        
        for ireg=1:ncols
            
            output.(reg_names{ireg}).(varList{ivar})=...
                ts(dn_start,tank(:,:,ireg),mod_names);
            
        end
        
    else
        
        output.(varList{ivar})=ts(dn_start,tank,mod_names);
        
    end
    
end

end


function [vlist,ncols,dn_start,dn_last]=do_intersection(dbcell)

vlist=fieldnames(dbcell{1});

this=dbcell{1}.(vlist{1});

dn1=this.date_numbers;

dn_start=dn1(1);

dn_last=dn1(end);

datta=double(this);

[~,ncols]=size(datta);

for imod2=2:numel(dbcell)
    
    vlist2=fieldnames(dbcell{imod2});
    
    vlist=intersect(vlist,vlist2);
    
    this=dbcell{imod2}.(vlist2{1});
    
    dn2=this.date_numbers;
    
    dn_start=min(dn_start,dn2(1));
    
    dn_last=max(dn_last,dn2(end));
    
    datta=double(this);
    
    [~,ncols2]=size(datta);
    
    ncols=max(ncols,ncols2);
    
end

end