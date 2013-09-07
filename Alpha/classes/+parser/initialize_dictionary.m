function dictionary = initialize_dictionary()
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

dictionary=struct();

dictionary.definitions={};
% unlike the declared list, definitions will never be sorted
dictionary.known_words={'steady_state','argzero','x0_','x1_','param_obj','commitment','discount',...
    'log','exp','cos','sin','normpdf','normcdf'};
% Never change the line below without updating determine_status (which has
% been expanded to allow for xx_ssmdef_ with as many digits as possible at
% the end.
dictionary.known_words=[dictionary.known_words,strcat({'xx_ssmdef_'},cellstr(int2str((1:9)'))')];
% dictionary.steady_state_parameters={};
dictionary.time_varying_probabilities={};
dictionary.symbols={'#','!',')','(','}','{',']','[',',',';','.','=','@'};
dictionary.add_operators={'+','-'};
dictionary.mult_operators={'*','^','/'};
dictionary.relational_operators={'<','>'};

distr_list=what('+distributions');
for ilist=1:numel(distr_list)
    if ~isempty(distr_list(ilist).m)
        mydistrlist=distr_list(ilist).m;
        break
    end
end
% for some reason I do not understand, this sometimes returns a 2 x 1
% structure instead of a 1 x 1. I have experienced it when using parallel
% computing, but not otherwise.
dictionary.Distributions=strrep(mydistrlist,'.m','_pdf');
% % % dictionary.Distributions={'uniform_pdf','unif_pdf','normal_pdf','norm_pdf',...
% % %     'gamma_pdf','gam_pdf','beta_pdf','invg_pdf','inv_gamma_pdf',...
% % %     'inv_gamma2_pdf','dirichlet_pdf'};
dictionary.syntax_special={'y]','param]',')]','x]','[(',']=','@f','(@','])',...
    '([','[n','n]',...
    '>=','<=','param>','>param',')>','>(','>f','(cn','cn)','cn,',',cn',',n',...
    '>n','<n'};
% last line added for the parameter restrictions block
dictionary.syntax_function={'y,','x,','param,','),','f,','n,',',y',',x',...
    ',param',',f',',n',',('};
dictionary.syntax_time={'y(','y{','x(','x{','n}','{+','}+','}*','{n','})','};','+}','}='};
dictionary.syntax_typical={'y)','y+','y*','y;','x)','x+','x*','x;','param)','def)',...
    'param+','def+','param*','def*','param;','def;','))',')+',')*',');','(y','(x',...
    '(param','(def','((','(+','(f','(n','+y','+x','+param','+def','+(','+f','+n',...
    '*y','*x','*param','*def','*f','*n','*(','f)','f(','n.','.n','f;','n)','n+',...
    'n*','nn','n;','=.','=param','=def','=x','=y','=f','=(','=+','=n',...
    'param=','def=','x=','y=','f=','n=',')=','tvp='};
% endogenous(y),exogenous(x),parameter(param), known word or function (f),
% number(n), definition(def)

dictionary.is_linear_model=false;
dictionary.input_list={'y','x','ss','param','def','s0','s1'};

end

