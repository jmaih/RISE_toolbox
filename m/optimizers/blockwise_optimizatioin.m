function blockwise_optimization(objective,x0,lb,ub,blocks,optimfunc,varargin)

% check literature on large-scale global optimization

npar=numel(x0);
if isempty(blocks)
	blocks={1:npar};
end
if isempty(optimfunc)
	optimfunc=@local_optimize;
end
nblks=numel(blocks);

while ~done
	for blk=1:nblks
		if rand<0.001
			this=1:npar;
			disp('optimizing the whole vector')
		else
			this=blocks{blk};
		end
		x0(this)=optimfunc(wrapper,x0(this),lb(this),ub(this),x0);
	end
end

	function f=wrapper(z,x)
		x(this)=z;
		f=objective(x,varargin{:});
	end

end

