function [c,ceq,dc,dceq] = Scalarization(x,F,NL_const,arg,m,q)
%[f1,f2,df1,df2] =  F(x,arg)        Problem objective functions
%[c,ceq,dc,dceq] = NL_const(x,arg)  Problem constraint

outputs = cell(1, 4);
[outputs{:}]  = NL_const(x,arg); % problem constraints
f   = F(x,arg);
c   = outputs{1};
ceq = [f(2)-m*f(1)-q; outputs{2}];
dc  = outputs{3};
if isnan(outputs{4})
    dceq = outputs{4};
else
    dceq = [f(4)-m*f(3); outputs{4}];
end

end

