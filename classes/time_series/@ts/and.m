function this3=and(this1,this2)
% this1= ts(1990,rand(10,1))
% this2= ts('1990Q1',rand(10,1))
% this3= this1 & this2

if nargin~=2
    error([mfilename,':: number of arguments should be 2'])
end
this3=cat(2,this1,this2);

end

% if isempty(this2)
%     this3=this1;
% elseif isempty(this1)
%     this3=this2;
% else
%     [dn1,dn1_max,dn1_min,nvar1,npages1,freq1]=decompose_series(this1);
%     [dn2,dn2_max,dn2_min,nvar2,npages2,freq2]=decompose_series(this2);
% 
%     if freq1~=freq2
%         error('data should have same frequency')
%     end
%     dn_max=max(dn1_max,dn2_max);
%     dn_min=min(dn1_min,dn2_min);
%     dn=dn_min:dn_max;
%     nvar=nvar1+nvar2;
%     npages=max(npages1,npages2);
%     nobs=numel(dn);
%     
%     datta=nan(nobs,nvar,npages);
%     rows1=ts.set_locations(dn1,dn);
%     rows2=ts.set_locations(dn2,dn);
%     
%     datta(rows1,1:nvar1,1:npages1)=this1.data;
%     datta(rows2,nvar1+(1:nvar2),1:npages2)=this2.data;
%     
%     this3=ts(dn,datta,[this1.varnames(:);this2.varnames(:)]);
% end