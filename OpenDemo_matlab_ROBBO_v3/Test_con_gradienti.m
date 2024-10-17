clc
close all
clearvars

% add ROBBO library path
addpath('ROBBO_lib_v3')

%% BOP Problem Setup
% Example 2
% rotated anchor point different q-coordinate
w1 = 1;  % move the second anchor point (linear P.front)
w2 = 2;  % move the first anchor point  (linear P.front)
c  = 10; % move both the anchor points  (linear P.front)
tol = [0.02, 0.04]; % User tolerances expressed as percentage (see "Mode")

%% ROBBO input: struct setup
arg      = [5,8,3];   % Objective functions argument 
% fmincon inputs, see "fmincon" documentation
xLB      = [0  0];           
xUB      = [5, 3];
A        = [];
b        = [];
Aeq      = [];
beq      = [];
x0       = [0 0];
SolvOptions = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', true, ...
    'Display', 'iter-detailed', ...
    'OptimalityTolerance', 1e-12, ...
    'StepTolerance', 1e-12, ... 
    'ConstraintTolerance', 1e-8);

% input struct1: Bi-objective opt problem info
BOP = struct('NL_const',@NL_const, 'xLB', xLB,'xUB',xUB, 'A', A, 'b', b,'Aeq',Aeq,'beq',beq, 'x0', x0, 'arg' ,arg);

% manual anchor points
P1 = [9.84e-13,50]';
P2 = [136,4]';
xanch1 = [4.06e-07,2.88e-07];
xanch2 = [5,2];
AnchorPointsStruct = struct('FirstAnchorPoint', P1, 'SecondAnchorPoint', P2, 'FirstAnchorSolution', xanch1, 'SecondAnchorSolution', xanch2);

% input struct2: ROBBO functioning mode
Mode = struct('tol','percentage','samp','uniform','est','central','iterNumber','manual','Nstop',10,'SpecifyAnchorPoints',[]);


%% Algorithm and exact PF
[report] = ROBBO(@F,BOP,tol,SolvOptions,Mode);

% Funzione obiettivo anonima che passa 'arg' come parametro
objFun = @(x) F(x, arg);

% Vincoli non lineari anonimi
nonlcon = @(x) NL_const(x, arg);
 SolvOptions = optimoptions('gamultiobj','PopulationSize',100);
% Chiamata a paretosearch
[x, fval] =  gamultiobj(objFun,2,A,b,Aeq,beq,xLB,xUB,nonlcon,SolvOptions);
% Plot del fronte di Pareto trovato
figure(2)
subplot('Position',[0.1, 0.5, 0.8, 0.45])
hold on
scatter(fval(:,1), fval(:,2), 'black','filled');
legend({'P. front estimate','Bounds','Computed solutions','','GA solutions'},'Location', 'best','FontSize',15,'Interpreter', 'latex')
grid on;

% 
% % P.front analytical expression (Plot)
% figure(2)
% subplot('Position',[0.1, 0.5, 0.8, 0.45])
% hold on
% filename = fullfile('Animation.gif');
% delayTime = 0.8; % Gif frame delay
% frame = getframe(2);
% im = frame2im(frame);
% [imind,cm] = rgb2ind(im,256);
% imwrite(imind,cm,filename,'gif','WriteMode','append');
% imwrite(imind, cm, filename, 'gif',  'WriteMode', 'append', 'DelayTime', delayTime);

%% Functions and contraints definition (BOP: Bi-objective Oprimization Problem)
% Vector-valued objective functions
function [f] = F(x, arg)
    f1 = 4*x(1)^2+ 4*x(2)^2;
    f2 = (x(1)-arg(1))^2 + (x(2)-arg(1))^2;

    df1 = [8*x(1),8*x(2)];
    df2 = [2*(x(1)-arg(1)), 2*(x(2)-arg(1))];

    f  = [f1,f2,df1,df2]';

end

% Nonlinear constraints
function [c, ceq,dc,dceq] = NL_const(x,arg)
    c   = [(x(1)-arg(1))^2 + x(2)^2-25; 7.7 - ((x(1)-arg(2))^2+(x(2)+arg(3))^2)];
    dc = [ [2*(x(1)-arg(1));2*x(2)], [2*(x(1)-arg(2));2*(x(2)+arg(3))] ];
    ceq = [];
    dceq  = []; 

end


