function dictionary = initialize_dictionary()
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

dictionary=struct();

dictionary.auxiliary_variables=struct('model',{{}},'ssmodel_solved',{{}});

dictionary.definitions={};
% unlike the declared list, definitions will never be sorted
dictionary.known_words={'steady_state','argzero','x0_','x1_','param_obj','commitment','discount',...
    'log','exp','cos','sin','normpdf','normcdf','$'};

% dictionary.steady_state_parameters={};
dictionary.time_varying_probabilities={};

dictionary.symbols={'#','!','?',')','(','}','{',']','[',',',';','.','=','@'};

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

dictionary.syntax_special={'y]','param]',')]','x]','[(',']=','@f','(@','])',...
    '([',...
    '>=','<=','param>','>param',')>','>(','>f','(cn','cn)','cn,',',cn',',n',...
    '>n','<n'};
% last line added for the parameter restrictions block
dictionary.syntax_function={'y,','x,','param,','def,','),','f,','n,',',y',',x',...
    ',param',',def',',f',',n',',('};

dictionary.syntax_time={'y(','y{','x(','x{','n}','{+','}+','}*','{n','})','};','+}','}='};

dictionary.syntax_typical={'y)','y+','y*','y;','x)','x+','x*','x;','param)','def)',...
    'param+','def+','param*','def*','param;','def;','))',')+',')*',');','(y','(x',...
    '(param','(def','((','(+','(f','(n','+y','+x','+param','+def','+(','+f','+n',...
    '*y','*x','*param','*def','*f','*n','*(','f)','f(','n.','.n','f;','n)','n+',...
    'n*','nn','n;','=.','=param','=def','=x','=y','=f','=(','=+','=n',...
    'param=','def=','x=','y=','f=','n=',')=','tvp=','*+','([','])','[n','n]',',+'};
% endogenous(y),exogenous(x),parameter(param), known word or function (f),
% number(n), definition(def)

dictionary.determine_status=@determine_status;

dictionary.input_list=parser.input_list();

end

function [status,loc]=determine_status(x,dict,neat)

if nargin<3

    neat=false;

end

if neat %nargout==1

    if nargout>1
    
        error('functional programming does not allow more than one output')
    
    end
    
    status=if_elseif(any(strcmp(x,{dict.endogenous.name})),'y',... %(1)
    any(strcmp(x,{dict.exogenous.name})),'x',... %(2)
    any(strcmp(x,{dict.parameters.name})),'param',... %(3)
    any(strcmp(x,dict.add_operators)),'+',... %(4)
    any(strcmp(x,dict.mult_operators)),'*',... %(5)
    any(strcmp(x,dict.symbols)),x,... %(6)
    ~isnan(str2double(x)),'n',... %(7)
    any(strcmp(x,dict.relational_operators)),'>',... %(8)
    any(strcmp(x,dict.chain_names)),'cn',... %(9)
    any(strcmp(x,dict.time_varying_probabilities)),'tvp',... %(10)
    any(strcmp(x,dict.definitions)),'def',... %(11)
    any(strcmp(x,dict.known_words))||...%(12)
    exist([x,'.m'],'file')... %(13)
    ||(strncmp(x,'xx_ssmdef_',10) && all(isstrprop(x(11:end),'digit'))),... %(14)
    'f',...
    1,'unknown');

else
    
    status='unknown';
    
    loc=find(strcmp(x,dict.endogenous_list),1,'first');
    
    if ~isempty(loc)
    
        status='y'; %(1)
    
    else
        
        loc=find(strcmp(x,dict.exogenous_list),1,'first');
        
        if ~isempty(loc)
        
            status='x'; %(2)
        
        else
            
            loc=find(strcmp(x,dict.parameters_list),1,'first');
            
            if ~isempty(loc)
            
                status='param'; %(3)
            
            else
                
                loc=find(strcmp(x,dict.add_operators),1,'first');
                
                if ~isempty(loc)
                
                    status='+'; %(4)
                
                else
                    
                    loc=find(strcmp(x,dict.mult_operators),1,'first');
                    
                    if ~isempty(loc)
                    
                        status='*'; %(5)
                    
                    else
                        
                        loc=find(strcmp(x,dict.symbols),1,'first');
                        
                        if ~isempty(loc)
                        
                            status=x; %(6)
                        
                        else
                            
                            flag=~isnan(str2double(x));
                           
                            if flag
                            
                                loc=flag;
                            
                            end
                            
                            if ~isempty(loc)
                            
                                status='n'; %(7)
                            
                            else
                                
                                loc=find(strcmp(x,dict.relational_operators),1,'first');
                                
                                if ~isempty(loc)
                                
                                    status='>'; %(8)
                                
                                else
                                    
                                    loc=find(strcmp(x,dict.chain_names),1,'first');
                                    
                                    if ~isempty(loc)
                                    
                                        status='cn'; %(9)
                                    
                                    else
                                        
                                        loc=find(strcmp(x,dict.time_varying_probabilities),1,'first');
                                        
                                        if ~isempty(loc)
                                        
                                            status='tvp'; %(10)
                                        
                                        else
                                            
                                            loc=find(strcmp(x,dict.definitions),1,'first');
                                            
                                            if ~isempty(loc)
                                            
                                                status='def'; %(11)
                                            
                                            else
                                                
                                                loc=find(strcmp(x,dict.known_words),1,'first');
                                                
                                                if ~isempty(loc)||... %(12)
                                                        exist([x,'.m'],'file')||... % no need for location  %(13)
                                                        (strncmp(x,'xx_ssmdef_',10) && all(isstrprop(x(11:end),'digit'))) % no need for location %(14)
                                                
                                                    status='f';
                                                
                                                end
                                                
                                            end
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

end

