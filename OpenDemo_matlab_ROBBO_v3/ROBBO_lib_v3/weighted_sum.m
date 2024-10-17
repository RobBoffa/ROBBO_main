function [cf, Dcf]= weighted_sum(x,F,W,arg,opt)
    
f = F(x,arg); 
df1 = f(3:2+length(x));
df2 = f(3+length(x):end);

cf = W(1)*f(1) + W(2)*f(2);

if opt.SpecifyConstraintGradient && opt.SpecifyObjectiveGradient
    Dcf = W(1)*df1 + W(2)*df2;
else
    Dcf = [];
end


end

