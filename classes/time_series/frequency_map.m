function fmap=frequency_map()
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

fmap={
    'W',52
    'M',12
    'Q',4
    'H',2
    '',1
    };
fmap=struct('strings',{fmap(:,1)},'code',cell2mat(fmap(:,2)));
end