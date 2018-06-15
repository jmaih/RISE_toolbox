function c=encode_map(m,siz)
% INTERNAL FUNCTION
%

ind=sub2ind(siz,m(:,2),m(:,3));

c=[m(:,1),ind];

end
