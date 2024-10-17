function [UB,LB,PF,UN] = Boud(p,mode)

% [~,Np] = size(p);

Nsamp = 10e4;
UB = zeros(2,Nsamp);
LB = zeros(2,Nsamp);

V  = linspace(p(1,1),p(1,end),Nsamp);

i = 1;
while i <= Nsamp
    v = V(i);
    index1 = find(p(1,:) <= v, 1, 'last');
    index2 = index1 + 1;

    p1  =  p(:,index1);
    p2  =  p(:,index2);
    L  = abs(p1-p2);
    l = (L(1)-L(2))/2;

    allv = find(V>=p1(1) & V<=p2(1));
    for j=1:length(allv)
        k = allv(j);
        vj = V( k);
        d1 = vj   - p1(1);
        d2 = p2(1)- vj;
        if d1 <= l
            UB(:, k) = [vj;  p1(2)+d1];
            LB(:, k) = [vj;  p1(2)-d1];
        elseif d1>l && d2>l
            if p1(2)<p2(2)
                UB(:, k) = [vj; p1(2)+d1];
                LB(:, k) = [vj; p2(2)-d2];
            else
                UB(:, k) = [vj; p2(2)+d2];
                LB(:, k) = [vj; p1(2)-d1];
            end
        elseif d1>l && d2<=l
            UB(:, k) = [vj; p2(2)+ d2];
            LB(:, k) = [vj; p2(2)-d2];
        else
            error('Not consistent interval')
        end
    end
    i = allv(end)+1;
end

if strcmp(mode,"central")
    PF = 0.5*(UB+LB);
    UN = 0.5*(UB(2,:)-LB(2,:));
elseif strcmp(mode,"linear")
    Lin = interp1(p(1,:),p(2,:),V);
    PF = [V;
          Lin];
    distUB = UB(2,:) - PF(2,:);
    distLB = PF(2,:) - LB(2,:); 
    UN = max(distUB,distLB);
else 
    error('Mode.est must be: central, linear')
end

end

