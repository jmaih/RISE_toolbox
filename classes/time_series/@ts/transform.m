function data=transform(datam,type)
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


% type is one of the following:
% 1 or 'level': Untransformed data - level (default)
% 2 or 'pct_ch_ar': Percentage change at compound annual rate (% chg)
% 3 or 'pct_ch': Period to period percentage change (% change)
% 4 or 'yoy_pct_ch': Year over year percentage change (yr/yr % chg)
% 5 or 'diff': Period to period difference change
% 6 or 'yoy_diff': Year over Year difference change
% 7 or 'log_ch_ar': %Log Change - compound annual rate (ann l-chg)
% 8 or 'log_ch': %Log Change - period to period (log-change)
% 9 or 'yoy_log_ch': Year to year log change (yr-yr l-chg)
% this follows the definitions of Haver
        
if nargin<2
    type=[];
end
if isempty(type)
    type=1;
end
switch datam.frequency
    case 'Q'
        npy=4;
    case 'M'
        npy=12;
    case 'H'
        npy=2;
    case ''
        npy=1;
    otherwise
        error(['frequency ',datam.frequency,' not implemented'])
end
switch type
    case {1,'level'}% Untransformed data - level (default)
        data=datam;
    case {2,'pct_ch_ar'}%Percentage change at compound annual rate (% chg)
        data=100*((datam/datam{-1})^npy-1);
    case {3,'pct_ch'}% Period to period percentage change (% change)
        data=100*((datam/datam{-1})-1);
    case {4,'yoy_pct_ch'} %Year over year percentage change (yr/yr % chg)
        data=100*((datam/datam{-npy})-1);
    case {5,'diff'}%Period to period difference change
        data=datam-datam{-1};
    case {6,'yoy_diff'}%Year over Year difference change
        data=datam-datam{-npy};
    case {7,'log_ch_ar'}%Log Change - compound annual rate (ann l-chg)
        data=npy*100*log(datam/datam{-1});
    case {8,'log_ch'}%Log Change - period to period (log-change)
        data=100*log(datam/datam{-1});
    case {9,'yoy_log_ch'}%Year to year log change (yr-yr l-chg)
        data=100*log(datam/datam{-npy});
    otherwise
        error(['unknown transformation ',type])
end