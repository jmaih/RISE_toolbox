function p=period2period(p0,freq0,freq1,head)
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

% transform a period a particular frequency into a period of another
% frequency
% p0 is the initial period and freq0 is the initial frequency. freq1 is the
% new frequency, while head is a flag the denotes the beginning of a period
% or the end.
%            for annual data, p=1
%            for bi-annual data, p=[1,2]
%            for quarterly data, p=[1,2,3,4]
%            for mmonthly data, p=[1,2,3,4,5,6,7,8,9,10,11,12]
% e.g. 1990Q1 = 1990M1 if head and 1990M3 if tail
% so 1 = period2period(1,4,12,1) and 3 = period2period(1,4,12,0)

if nargin<4
    head=[];
end
freq0=frequency2num(freq0);
freq1=frequency2num(freq1);
if numel(p0)==2 && isempty(head)
    p=p0;
    p(1)=period2period(p0(1),freq0(1),freq1(1),true);
    p(2)=period2period(p0(2),freq0(1),freq1(1),false);
    return
elseif isempty(head)
    error('head or tail must be specified in the last argument')
end

if head
    p=floor((p0-1)*freq1(1)/freq0(1)+1);
else
    p=ceil(p0*freq1(1)/freq0(1));
end

end