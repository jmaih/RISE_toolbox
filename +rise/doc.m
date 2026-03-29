%   DOC Reference page in Help browser.
%   
%      DOC opens the Help browser, if it is not already running, and 
%      otherwise brings the Help browser to the top.
%    
%      DOC FUNCTIONNAME displays the reference page for FUNCTIONNAME in
%      the Help browser. FUNCTIONNAME can be a function or block in an
%      installed MathWorks product.
%    
%      DOC METHODNAME displays the reference page for the method
%      METHODNAME. You may need to run DOC CLASSNAME and use links on the
%      CLASSNAME reference page to view the METHODNAME reference page.
%    
%      DOC CLASSNAME displays the reference page for the class CLASSNAME.
%      You may need to qualify CLASSNAME by including its package: DOC
%      PACKAGENAME.CLASSNAME.
%    
%      DOC CLASSNAME.METHODNAME displays the reference page for the method
%      METHODNAME in the class CLASSNAME. You may need to qualify
%      CLASSNAME by including its package: DOC PACKAGENAME.CLASSNAME.
%    
%      DOC FOLDERNAME/FUNCTIONNAME displays the reference page for the
%      FUNCTIONNAME that exists in FOLDERNAME. Use this syntax to display the
%      reference page for an overloaded function.
%    
%      DOC USERCREATEDCLASSNAME displays the help comments from the
%      user-created class definition file, UserCreatedClassName.m, in an
%      HTML format in the Help browser. UserCreatedClassName.m must have a
%      help comment following the classdef UserCreatedClassName statement
%      or following the constructor method for UserCreatedClassName. To
%      directly view the help for any method, property, or event of
%      UserCreatedClassName, use dot notation, as in DOC
%      USERCREATEDCLASSNAME.METHODNAME. 
% 
%      Examples:
%         doc abs
%         doc fixedpoint/abs  % ABS function in the Fixed-Point Designer Product
%         doc handle.findobj  % FINDOBJ method in the HANDLE class
%         doc handle          % HANDLE class
%         doc containers.Map  % Map class in the containers method
%         doc sads            % User-created class, sads
%         doc sads.steer      % steer method in the user-created class, sads
%
%    Documentation for doc
%       doc doc
%
%