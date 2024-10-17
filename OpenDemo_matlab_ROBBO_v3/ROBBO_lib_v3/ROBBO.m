function [report] = ROBBO(F,BOP,tol,SolvOptions,Mode)

%% Optimization problem features
NL_const =  BOP.NL_const;
xLB      =  BOP.xLB;
xUB      =  BOP.xUB;
A        =  BOP.A;
b        =  BOP.b;
Aeq      =  BOP.Aeq;
beq      =  BOP.beq;
arg      =  BOP.arg;
x0       =  BOP.x0;

% RBBO default mode
if isempty(Mode)
    Mode = struct('tol','percentage','samp','uniform','est','central','iterNumber','auto','Nstop',0,'SpecifyAnchorPoints',[]);
end

if isempty(SolvOptions)
    SolvOptions = optimoptions('fmincon',...
        'SpecifyConstraintGradient',false,...
        'SpecifyObjectiveGradient',false);
end

%% Anchor Points
if isempty(Mode.SpecifyAnchorPoints) 
    % compute Utopia
    
    % compute weak upper anchor point
    W = [1, 0];
    sub_xanc1 = fmincon(@(x)weighted_sum(x,F,W,arg,SolvOptions),x0,A,b,Aeq,beq,xLB,xUB,@(x)NL_const(x,arg),SolvOptions); 
    fsub_anc1   = F(sub_xanc1,arg);
    
    % compute weak lower anchor point
    W = [0,1];
    sub_xanc2 = fmincon(@(x)weighted_sum(x,F,W,arg,SolvOptions),x0,A,b,Aeq,beq,xLB,xUB,@(x)NL_const(x,arg),SolvOptions); 
    fsub_anc2   = F(sub_xanc2,arg);
    
    Utopia = [fsub_anc1(1),fsub_anc2(2)]';
    
    % compute strong upper anchor point
    W = [0, 1];
    xanc1 = fmincon(@(x)weighted_sum(x,F,W,arg,SolvOptions),sub_xanc1,A,b,Aeq,beq,xLB,xUB,@(x)costr_ut(x,F,arg,NL_const,Utopia,1,SolvOptions),SolvOptions); 
    fanch1   = F(xanc1,arg);
    anc1     = fanch1(1:2);
    % compute strong lower anchor point
    W = [1, 0];
    xanc2 = fmincon(@(x)weighted_sum(x,F,W,arg,SolvOptions),sub_xanc2,A,b,Aeq,beq,xLB,xUB,@(x)costr_ut(x,F,arg,NL_const,Utopia,2,SolvOptions),SolvOptions); 
    fanch2   = F(xanc2,arg);
    anc2 = fanch2(1:2);
else
    anc1 =  Mode.SpecifyAnchorPoints.FirstAnchorPoint;
    anc2 =  Mode.SpecifyAnchorPoints.SecondAnchorPoint;
    xanc1 = Mode.SpecifyAnchorPoints.FirstAnchorSolution;
    xanc2 = Mode.SpecifyAnchorPoints.SecondAnchorSolution;

end


% compute distances along f1 and f2 axes (PF ranges)
dist = abs(anc2-anc1); 

%% Rotation Matrix
% Compute tolerances
if strcmp(Mode.tol,'percentage')  % tolerances expressed as percentage of the PF ranges
    d1 = dist(1)*tol(1);
    d2 = dist(2)*tol(2);
elseif strcmp(Mode.tol,'absolute') % tolerances as absolute values
    d1 = tol(1);
    d2 = tol(2);
else
     error('Mode.tol must be: percentage, absolute')
end
% Copute the rotation matrix based on the tolerances
R = sqrt(2)/2*[1 -1;1 1]*[1/d1 0; 0 1/d2];
% set the scalarization weights (of the weighted sum as W)
W = [d1,d2]; 

%% Plot anchor poits
% ancorpoints in (q,v) ref sys
Ranc1  =  R*anc1; % anchor points in (q,v) "rotated" ref sys
Ranc2  =  R*anc2;
% Compute point distances along q and v
Rdist  =  abs(Ranc2-Ranc1);

% Upper bound number of samples for Bisection
l = Rdist(1);
nsampleGreed = 2;
i = 0;
while l>2*sqrt(2)
    i = i+1;
    l = Rdist(1)/2^i;
    nsampleGreed= nsampleGreed+2^(i-1);
end

% Upper bound number of samples for uniform SetMembership sampling
nsampleUni = ceil(Rdist(1)/(2*sqrt(2)))+1;

% Uniform partitioning of V0
Vuni = linspace(Ranc1(1),Ranc2(1),nsampleUni);

% Generate Animation plot
figure(1), set(gcf, 'Color', 'w');  figWidth = 1600; % length in pixel
                        figHeight = 900; % high in pixel
                        set(gcf, 'Position', [100, 100, figWidth, figHeight]);
                        pos1 = [0.1, 0.5, 0.8, 0.45];
                        pos2 = [0.1, 0.1, 0.8, 0.25];
            
                        subplot('Position', pos1);title('Pareto front (q,v)','Interpreter', 'latex','FontSize',15)
                        hold on
                        if strcmp(Mode.samp,'uniform'), xline(Vuni,":"), end
                        scatter(Ranc1(1),Ranc1(2),'red','filled')
                        scatter(Ranc2(1),Ranc2(2),'red','filled')
                        hold off
                        title('ROBBO: P.front estimation via Set membership','Interpreter', 'latex')
                        xlabel('$v$','Interpreter', 'latex','FontSize',10)
                        ylabel('$q$','Interpreter', 'latex','FontSize',10)
                        ax = gca;
                        ax.LineWidth = 1; 
                        ax.FontSize = 18; 
                        subplot('Position', pos2);
                        title('Uncertainty','Interpreter', 'latex')
                        xlabel('$v$','Interpreter', 'latex','FontSize',10)
                        ylabel('$\lambda_{\mathcal{D}}(v)$','Interpreter', 'latex','FontSize',10)
                        ax = gca;
                        ax.LineWidth = 1; 
                        ax.FontSize = 18; 


% GIF - initialization
Nfig = 1;
drawnow
filename = fullfile('Animation.gif');
frame = getframe(Nfig);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
delayTime = 0.2; % Gif frame delay
imwrite(imind, cm, filename, 'gif',  'LoopCount',inf,  'DelayTime', delayTime);

% Initialize PF points collector
P = zeros(2,nsampleGreed);
P(:,1:2)=[Ranc1,Ranc2];

% Compute Bounds/Estimate/Uncertainty
nsamp = 2;
[UB,LB,PF,UN] = Boud(P(:,1:nsamp),Mode.est);

% Update Plot
figure(1),  
            subplot('Position', pos1);
                hold on
                if strcmp(Mode.samp,'uniform'), xline(Vuni,":"), end
                plot(UB(1,:),UB(2,:),'-g','LineWidth',1.5),plot(LB(1,:),LB(2,:),'-g','LineWidth',1.5),plot(PF(1,:),PF(2,:),'-r','LineWidth',1.5)
                hold off
                title('ROBBO: P.front estimation via Set membership','Interpreter', 'latex')
                xlabel('$v$','Interpreter', 'latex','FontSize',10)
                ylabel('$q$','Interpreter', 'latex','FontSize',10)
            ax = gca;
            ax.LineWidth = 1; 
            ax.FontSize = 18; 
            subplot('Position', pos2);
                hold on,plot(PF(1,:),UN,'-b','LineWidth',1.5),hold off,title('Uncertainty','Interpreter', 'latex');
                xlabel('$v$','Interpreter', 'latex','FontSize',10)
                ylabel('$\lambda_{\mathcal{D}}(v)$','Interpreter', 'latex','FontSize',10)
            ax = gca;
            ax.LineWidth = 1; 
            ax.FontSize = 18; 

%% Main cycle

a_cost   =   d2/d1;          % angular coefficient of constraint f2=a_cost*f1+b_cost
if strcmp(Mode.samp,'greedy') || strcmp(Mode.samp,'bisection')
    Nstep               =   nsampleGreed;             % max number of iterations
elseif strcmp(Mode.samp,'uniform')
    Nstep               =   nsampleUni;               % max number of iterations
else
    error('Mode.samp must be: greedy, uniform, bisection')
end

if strcmp(Mode.iterNumber,'manual')  
    Nsolutions = Nstep-2;
elseif strcmp(Mode.iterNumber,'auto') 
    Nsolutions = Mode.Nstop-2;
else
     error('Mode.iterNumber must be: manual, auto')
end

solutions = zeros(Nsolutions,length(x0));
solutions([1,2],:) = [xanc1;xanc2];

for i = 1:Nsolutions-2
    Un_max = 0;
    % Compute the next sample coordinate
    if (strcmp(Mode.samp,'greedy') && strcmp(Mode.est,'central'))   ||  strcmp(Mode.samp,'bisection')
        for j = 1:nsamp-1
            P1 = P(:,j);
            P2 = P(:,j+1);
            Rdist  = abs(P2-P1);
            if strcmp(Mode.est,'central')
                Un_interval   = 0.5*(Rdist(1)-Rdist(2)); % initial uncertainty
            elseif strcmp(Mode.est,'linear')
                Un_interval   = (Rdist(1)^2-Rdist(2)^2)/(2*Rdist(1)); % initial uncertainty
            else
                error('Mode.est must be: central, linear')
            end
            % find interval with max uncertainty
            if Un_interval>=Un_max
                max_int = j;
                Un_max = Un_interval;
            end
        end
        
        P1 = R\P(:,max_int);
        P2 = R\P(:,max_int+1);
        mp = 0.5*(P1+P2);

    elseif strcmp(Mode.samp,'greedy') && strcmp(Mode.est,'linear')
        [val,ind] = max(UN);
        vj = UB(1,ind);
        index1 = find(P(1,1:nsamp) <= vj, 1, 'last');
        index2 = index1 + 1;
        P1  =  P(:,index1);
        P2  =  P(:,index2);
        d1 = vj - P1(1);
        d2 = P2(1) - vj;
        ub = min(P1(2)+d1,P2(2)+d2);
        lb = max(P1(2)-d1,P2(2)-d2);
        Un_interval = val;
        if Un_interval>Un_max
                max_int = index1;
                Un_max = Un_interval;
                rmp     = [vj;0.5*(ub + lb)];
        elseif Un_interval==Un_max
                if abs(vj-0.5*(P1(1)+P2(1))) < abs(rmp(1)-0.5*(P1(1)+P2(1)))
                    max_int = index1;
                    Un_max = Un_interval;
                    rmp     = [vj;0.5*(ub + lb)];
                end
        end
        mp = R\rmp;

    elseif strcmp(Mode.samp,'uniform')
        rmp = Vuni(:,1);
        for j=1:nsampleUni-1
            vj = Vuni(1,j);
            index1 = find(P(1,1:nsamp) <= vj, 1, 'last');
            index2 = index1 + 1;
            P1  =  P(:,index1);
            P2  =  P(:,index2);
            d1 = vj - P1(1);
            d2 = P2(1) - vj;
            ub = min(P1(2)+d1,P2(2)+d2);
            lb = max(P1(2)-d1,P2(2)-d2);

            if strcmp(Mode.est,'central')
                Un_interval = 0.5*(ub - lb);
            elseif strcmp(Mode.est,'linear')
                Lin = interp1([P1(1,:),P2(1,:)],[P1(2,:),P2(2,:)],vj);
                distUB = ub - Lin;
                distLB = Lin - lb; 
                Un_interval = max(distUB,distLB);
            else 
                error('Mode.est must be: central, linear')
            end

            
            if Un_interval>Un_max
                max_int = index1;
                Un_max = Un_interval;
                rmp     = [vj;0.5*(ub + lb)];
            elseif Un_interval==Un_max
                if abs(vj-0.5*(P1(1)+P2(1))) < abs(rmp(1)-0.5*(P1(1)+P2(1)))
                    max_int = index1;
                    Un_max = Un_interval;
                    rmp     = [vj;0.5*(ub + lb)];
                end
            end
        end
        mp = R\rmp;
    else
        error("Error in Mode")
    
    end
 
    if Un_max<=sqrt(2)
        break
    end

    b_cost =   mp(2)-a_cost*mp(1);
    xstar  = fmincon(@(x)weighted_sum(x,F,W,arg,SolvOptions),x0,A,b,Aeq,beq,xLB,xUB,@(x)Scalarization(x,F,NL_const,arg,a_cost,b_cost,SolvOptions),SolvOptions);
    fpnew  = F(xstar,arg);
    pnew   = fpnew(1:2);
    Rpnew  = R*pnew;

    P(:,max_int+2:nsamp+1) = P(:,max_int+1:nsamp);
    P(:,max_int+1) = Rpnew;
    nsamp  = nsamp+1;

    solutions(max_int+2:nsamp+1,:) = solutions(max_int+1:nsamp,:);
    solutions(max_int+1,:) = xstar;


    [UB,LB,PF,UN] = Boud(P(:,1:nsamp),Mode.est);
    
    figure(1)
        subplot('Position', pos1);
            plot(PF(1,:),PF(2,:),'-r','LineWidth',1.5)
            hold on
            h_bound = plot(UB(1,:), UB(2,:), 'g', 'LineWidth', 1.5);
            scatter(P(1,1:nsamp),P(2,1:nsamp),'red','filled')
            plot(LB(1,:), LB(2,:), 'g', 'LineWidth', 1.5);
            hold off
            title('ROBBO: P.front estimation via Set membership','Interpreter', 'latex')
            xlabel('$v$', 'Interpreter', 'latex','FontSize',10)
            ylabel('$q$', 'Interpreter', 'latex','FontSize',10)
            if strcmp(Mode.samp,'uniform'), xline(Vuni,":"), end
            legend({'P. front estimate','Bounds','P. front samples'},'Location', 'best','FontSize',15,'Interpreter', 'latex')
            ax = gca;
            ax.LineWidth = 1; 
            ax.FontSize = 18; 
        subplot('Position', pos2);
            plot(PF(1,:),UN,'-b','LineWidth',1.5),hold off,title('Uncertainty','Interpreter', 'latex');
            xlabel('$v$','Interpreter', 'latex','FontSize',10)
            ylabel('$\lambda_{\mathcal{D}}(v)$','Interpreter', 'latex','FontSize',10)
            ax = gca;
            ax.LineWidth = 1; 
            ax.FontSize = 18; 

    

    frame = getframe(Nfig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','append');
    imwrite(imind, cm, filename, 'gif',  'WriteMode', 'append', 'DelayTime', delayTime);


end

riP   = R\P(:,1:nsamp);
riUB  = R\UB;
riLB  = R\LB;
riPF   = R\PF;

figure(2), set(gcf, 'Color', 'w');  
            set(gcf, 'Position', [100, 100, figWidth, figHeight]);
            subplot('Position', pos1)
            plot(riPF(1,:),riPF(2,:),'-r','LineWidth',1.5),
            title('P. front estimate validation','Interpreter', 'latex')
            hold on
            plot(riUB(1,:),riUB(2,:),'g','LineWidth',1.5)
            scatter(riP(1,:),riP(2,:),'red','filled'),
            plot(riLB(1,:),riLB(2,:),'g','LineWidth',1.5)
            hold off
            xlabel('$f_1$','Interpreter', 'latex','FontSize',10)
            ylabel('$f_2$','Interpreter', 'latex','FontSize',10)
            ax = gca;
            ax.LineWidth = 1; 
            ax.FontSize = 18; 
            legend({'P. front estimate','Bounds','Computed solutions'},'Location', 'best','FontSize',15,'Interpreter', 'latex')
        subplot('Position', pos2)
            plot(riPF(1,:),UN,'-b','LineWidth',1.5)
            title('Uncertainty','Interpreter', 'latex')
            xlabel('$f_1$','Interpreter', 'latex','FontSize',10)
            ylabel('$\lambda_{\mathcal{D}}(f_1)$','Interpreter', 'latex','FontSize',10)
            ax = gca;
            ax.LineWidth = 1; 
            ax.FontSize = 18; 

    

%% final report
report = struct();
report.nsampleGreed = nsampleGreed;
report.nsampleUni   = nsampleUni;
report.nsamp        = nsamp;
report.plotQV       = struct('P',P,'UB',UB,'LB',LB,'PF',PF);
report.plotf1f2     = struct('riP',riP,'riUB',riUB,'riLB',riLB,'riPF',riPF);
report.input        = struct('BOP',BOP,'tol',tol,'SolvOptions',SolvOptions,'Mode',Mode);
report.solutions    = solutions(1:nsamp,:);
end


