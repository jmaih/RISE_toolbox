%% housekeeping
clc,clear classes
%% create some variables
disp('1- creating variables')
a=sad_tree('a');b=sad_tree('b');c=sad_tree('c');
%% create a function
disp('2- creating a function')
func=@(a,b,c)exp(a+2*log(b+c)-a*atan(b*c));
%% run the function with the variables to create a tree
disp('3- creating a tree')
tree=func(a,b,c);
%% differentiate wrt variable a
disp('4- differentiating with respect to a')
da=tree.diff(a)
%% print the derivative
disp('5- printing the derivative wrt a')
char(da)
%% reflag the tree 
disp('6- flagging the original tree to create auxiliary variables: none will be created if all operations are only done once')
tree.re_flag_tree
%% differentiate wrt b
disp('7- differentiating with respect to b')
db=tree.diff(b)
%% print the derivate
disp('8- printing the derivative wrt b')
char(db)
%% reflag the tree again
disp('9- re-flagging the original tree again: in this case, repeated operations lead to the creation of auxiliary variables')
tree.re_flag_tree
%% reprint the derivatives
disp('10- print the shortened derivatives')
char(da)
char(db)
%% print the expanded derivatives
disp('11- print the expanded derivatives')
char(da,1)
char(db,1)
%% print the tree
disp('12- viewing the operation tree in the original tree')
tree.print
%% get the references
disp('13- getting the auxiliary variables or references')
tree.collect_references
