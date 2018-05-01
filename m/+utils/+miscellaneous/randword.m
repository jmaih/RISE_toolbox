function word=randword(length)
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

alphabet={'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};
numbers={'0','1','2','3','4','5','6','7','8','9'};
symbols={'_'}; %'§','@','£','#','$','&',
a_nbr=numel(alphabet);
n_nbr=numel(numbers);
s_nbr=numel(symbols);

total=a_nbr+n_nbr+s_nbr;
aw=a_nbr/total;
nw=n_nbr/total;
sw=1-aw-nw;

cs=[0,cumsum([aw,nw,sw])];
cs(end)=1;

word=alphabet{randsample(a_nbr,1)};
for j=2:length
    r=rand;
    if r>=cs(1) && r<cs(2)
        draw=alphabet{randsample(a_nbr,1)};
    elseif r>=cs(2) && r<cs(3)
        draw=numbers{randsample(n_nbr,1)};
    elseif r>=cs(3)
        draw=symbols{randsample(s_nbr,1)};
    end
    word=[word,draw]; %#ok<AGROW>
end
