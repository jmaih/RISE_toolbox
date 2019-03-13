function n=norm(varargin)

fval=varargin{1};

n=fval.'*fval; % this used by LM so do not change to below...

% % n=norm(fval,inf);

end