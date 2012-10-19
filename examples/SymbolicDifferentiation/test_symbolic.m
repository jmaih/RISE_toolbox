%% housekeeping
clc,clear classes
%% create some variables
fprintf(1,'%s\n\n','1- creating variables');
a=sad_tree('a');b=sad_tree('b');c=sad_tree('c');
%% create a function
fprintf(1,'%s\n\n','2- creating a function');
func=@(a,b,c)exp(a+2*log(b+c)-a*atan(b*c));
%% run the function with the variables to create a tree
fprintf(1,'%s\n\n','3- creating a tree');
tree=func(a,b,c);
%% differentiate wrt variable a
fprintf(1,'%s\n\n','4- differentiating with respect to a');
da=tree.diff(a)
%% print the derivative
fprintf(1,'%s\n\n','5- printing the derivative wrt a');
char(da)
%% reflag the tree 
fprintf(1,'%s\n\n','6- flagging the original tree to create auxiliary variables: none will be created if all operations are only done once');
tree.re_flag_tree
%% differentiate wrt b
fprintf(1,'%s\n\n','7- differentiating with respect to b');
db=tree.diff(b)
%% print the derivate
fprintf(1,'%s\n\n','8- printing the derivative wrt b');
char(db)
%% reflag the tree again
fprintf(1,'%s\n\n','9- re-flagging the original tree again: in this case, repeated operations lead to the creation of auxiliary variables');
tree.re_flag_tree
%% reprint the derivatives
fprintf(1,'%s\n\n','10- print the shortened derivatives');
char(da)
char(db)
%% print the expanded derivatives
fprintf(1,'%s\n\n','11- print the expanded derivatives')
char(da,1)
char(db,1)
%% print the tree
fprintf(1,'%s\n\n','12- viewing the operation tree in the original tree')
tree.print
%% get the references
fprintf(1,'%s\n\n','13- getting the auxiliary variables or references')
tree.collect_references
