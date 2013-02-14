function obj=compute_definitions(obj,pp)

definitions=obj.func_handles.definitions;

number_of_regimes=obj.NumberOfRegimes;
if nargin<2
    parameters_image=obj.parameters_image;
    startval_loc= strcmp('startval',parameters_image(1,:));
    pp=parameters_image{2,startval_loc};
end
def=nan(numel(obj.definitions),number_of_regimes);
% evaluate definitions
for ii=number_of_regimes:-1:1 % backward to optimize speed
    def(:,ii)=online_function_evaluator(definitions,pp(:,ii)); %#ok<*EVLC>
end

if ~isempty(def)
    [nrows,ncols]=size(def);
    tmp=mat2cell(def,ones(nrows,1),ncols);
    [obj.definitions(1:nrows).value]=(tmp{:});
end
