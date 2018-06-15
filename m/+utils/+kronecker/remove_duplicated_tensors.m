function p=remove_duplicated_tensors(p)
% INTERNAL FUNCTION
%

ncases=size(p,1);

discard=false(1,ncases);

for ii=1:ncases-1

    if discard(ii)

        continue

    end

    this=p(ii,:);

    tmp=bsxfun(@minus,p(ii+1:end,:),this);

    doublons=find(all(tmp==0,2))+ii;

    discard(doublons)=true;

end

p=p(~discard,:);

end