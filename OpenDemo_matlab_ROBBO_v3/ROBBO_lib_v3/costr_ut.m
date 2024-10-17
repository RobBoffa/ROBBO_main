function [c,ceq] = costr_ut(x,F,arg,NL_const,Utopia,anch_ind)
    [c,ceq_prev] = NL_const(x,arg);
    z = F(x,arg);
    ceq_add = z(anch_ind)-Utopia(anch_ind);
    ceq = [ceq_prev;ceq_add];
end

