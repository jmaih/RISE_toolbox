function [dn,dn_max,dn_min,nvar,npages,freq]=decompose_series(this)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

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
