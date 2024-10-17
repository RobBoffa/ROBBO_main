function [c,ceq,Dc,Dceq] = Scalarization(x,F,NL_const,arg,m,q,opt)
% f = [f1,f2,df1,df2] =  F(x,arg)        Problem objective functions
f = F(x,arg); 
if opt.SpecifyConstraintGradient && opt.SpecifyObjectiveGradient
    [c,ceqBOP,Dc,DceqBOP] = NL_const(x,arg);
    df1 = f(3:2+length(x));
    df2 = f(3+length(x):end);
    Dceq = [df2-m*df1, DceqBOP];

else
    [c,ceqBOP] = NL_const(x,arg);
    Dc  = [];
    Dceq = [];
end  

ceq = [f(2)-m*f(1)-q; ceqBOP];


end

