function [Un] = compUn(uP,lP,m)
deltaX = lP(1)-uP(1);
deltaY = uP(2)-lP(2);

Uny = sqrt(deltaY^2+(deltaY/m)^2); 
Unx = sqrt(deltaX^2+(m*deltaX)^2); 
Un = min([Uny,Unx]);
end

