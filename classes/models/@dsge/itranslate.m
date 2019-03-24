%--- help for dsge/itranslate ---
%
%  translate -- Translates comprehensible atoms (model variables) into RISE codes
% 
%  ::
% 
%    outList=translate(obj,inList)
%    outList=translate(obj,inList,order)
% 
%  Args:
% 
%     obj (rise | dsge): scalar model object.
% 
%     inList (char | cellstr): List of atoms to translate e.g. C,
%        X{+1}, lambda_x, EPS_A,
% 
%     order (integer|0|{1}): If order>1 the returned list of atoms is a list
%        of kroneckers in which the number of elements in each items is the
%        order of the kronecker. This is useful for instance to understand
%        what combination of variables makes up a column in a
%        differentiation. If order=0, the translation is with respect to the
%        static, rather than the dynamic model
% 
%  Returns:
%     :
% 
%     - **outList** [cellstr]: List of translated atoms
% 
%  Example:
%     :
% 
%     list=get(obj,'endo_list')
%     outList=itranslate(obj,list)
%     outList=itranslate(obj,list,0) % contemporaneous or steady state
% 
%     list=get(obj,'exo_list')
%     outList=itranslate(obj,list)
% 
%     list=get(obj,'param_list')
%     outList=itranslate(obj,list)
% 
%     list=get(obj,'def_list')
%     outList=itranslate(obj,list)
% 
%  See also : dsge/translate
%     :
% 
% 
%