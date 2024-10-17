function CF = weighted_sum(x,F,W,arg)
    
Obj = F(x,arg); 
CF = W(1)*Obj(1) + W(2)*Obj(2);

end

