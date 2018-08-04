function [dn,dn_max,dn_min,nvar,npages,freq]=decompose_series(this)
% INTERNAL FUNCTION
%

dn=this.date_numbers;

dn_max=max(dn);

dn_min=min(dn);

nvar=this.NumberOfVariables;

npages=this.NumberOfPages;

freq=serial2frequency(dn);

if ~all(freq==freq(1))
    
    error('data should have same frequency')
    
end

freq=freq(1);

end
