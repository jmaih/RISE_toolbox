function [start,finish]= nber_dates(freq)
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

if nargin<1 || isempty(freq)
    freq='Q';
end
%  Official NBER recession dates, from
%   http://www.nber.org/cycles/cyclesmain.html  (see NBERcycles.html)

dates=[
% Peak Year,Month   Trough Year, Month
1857 6      1858 12;
1860 10     1861 6 ;
1865 4      1867 12;
1869 6      1870 12;
1873 10     1879 3 ;
1882 3      1885 5 ;
1887 3      1888 4 ;
1890 7      1891 5 ;
1893 1      1894 6 ;
1895 12     1897 6 ;
1899 6      1900 12;
1902 9      1904 8 ;
1907 5      1908 6 ;
1910 1      1912 1 ;
1913 1      1914 12;
1918 8      1919 3 ;
1920 1      1921 7 ;
1923 5      1924 7 ;
1926 10     1927 11;
1929 8      1933 3 ;
1937 5      1938 6 ;
1945 2      1945 10;
1948 11     1949 10;
1953 7      1954 5 ;
1957 8      1958 4 ;
1960 4      1961 2 ;
1969 12     1970 11;
1973 11     1975 3 ;
1980 1      1980 7 ;
1981 7      1982 11;
1990 7      1991 3 ;
2001 3      2001 11;
2007 12     2009 06;  % Starts Dec 2007, End=June 2009
];
switch upper(freq)
    case 'M'
        periods=dates(:,[2,4]);
    case 'Q'
        periods=ceil(dates(:,[2,4])/4);
    otherwise
end
start=strcat(num2str(dates(:,1)),freq,num2str(periods(:,1)));
finish=strcat(num2str(dates(:,3)),freq,num2str(periods(:,2)));

