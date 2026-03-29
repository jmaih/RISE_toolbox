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
%    Other uses of itranslate
%
%       rdico/itranslate
%