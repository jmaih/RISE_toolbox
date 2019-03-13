function opts=set_defaults(opts)

options=struct('Display','off',...% if this is set to zero, all screen output is suppressed
    'TolFun',sqrt(eps),...
    'MaxIter',1000);

if isempty(opts)
    
    opts=struct();
    
elseif ~isstruct(opts)
    
    error('opts must be empty or a structure')
    
end

fopts=fieldnames(options);

for iopt=1:numel(fopts)
    
    if ~isfield(opts,fopts{iopt})
        
        opts.(fopts{iopt})=options.(fopts{iopt});
        
    end
    
end

end