function [cBOP,ceq,DcBOP,Dceq] = costr_ut(x,F,arg,NL_const,Utopia,anch_ind,opt)
f =  F(x,arg);

if opt.SpecifyConstraintGradient && opt.SpecifyObjectiveGradient
    [cBOP,ceqBOP,DcBOP,DceqBOP] = NL_const(x,arg);
    df = [f(3:2+length(x)), f(3+length(x):end)];
    Dceq = [DceqBOP, df(:,anch_ind)];
else
   [cBOP,ceqBOP] = NL_const(x,arg);
   Dceq  = [];
   DcBOP = [];
end    

    f = F(x,arg);
    ceq_add = f(anch_ind)-Utopia(anch_ind);
    ceq = [ceqBOP ; ceq_add];


end

