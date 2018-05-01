function equation=greekify(equation)
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

% so far I will just greekify names but later on, I might also greekify
% equations
greek_letters={'alpha','beta','gamma','delta','epsilon','kappa',...
    'lambda','mu','nu','omega','phi','pi','chi','psi','rho',...
    'sigma','tau','upsilon','Sigma','Pi','Lambda','Omega','Gamma'};
% % DELIMITERS=[char([9:13,32]),'[]{}(),;=+-*/^@'];
DELIMITERS=parser.delimiters();%DEL=[char([9:13,32]),'[]{}(),;=+-*/^@><'];
% 1-\beta \frac{\left( 1-\frac{\kappa }{2}\left( \Pi _{t}-1\right) ^{2}\right)
% Y_{t}}{\left( 1-\frac{\kappa }{2}\left( \Pi _{t+1}-1\right) ^{2}\right)
% Y_{t+1}}\frac{1}{\exp \left( \mu _{t+1}\right) }\frac{R_{t}}{\Pi _{t+1}}
for ii=1:numel(equation)
    if isempty(equation(ii).tex_name)
        equation(ii).tex_name=greekify_intern(equation(ii).name);
    end
end
    function greek=greekify_intern(equation)
        greek='';
        while ~isempty(equation)
            [tokk,rest_]=strtok(equation,DELIMITERS);
            if isempty(tokk)
                greek=[greek,equation];
            else
                loc_=strfind(equation,tokk);
                loc_=loc_(1);
                greek=[greek,equation(1:loc_-1)];
                span=length(tokk);
                for i1=1:numel(greek_letters)
                    g=greek_letters{i1};
                    gspan=length(g);
                    if span>=gspan
                        if strcmp(tokk(1:span),g) && (span==gspan||strcmp(tokk(span+1),'_'))
                            right=tokk(span+2:end);
                            tokk=['$\',g];
                            if ~isempty(right)
                                tokk=[tokk,' _{',right,'}']; % NB: there is a space
                            end
                            tokk=[tokk,'$'];
                        end
                    end
                end
                greek=[greek,tokk];
            end
            equation=rest_;
        end
    end
end
