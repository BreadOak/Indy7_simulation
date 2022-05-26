%%%%%%%%%% Indy7 Simulation %%%%%%%%%%
clc
clear all
close all

% Load Indy7 information
addpath('MR/') 
[Slist, Mlist, Glist, M, w, p, robot] = load_urdf("UR5.urdf",6);

% Simulation
dt = 0.01;
endTime = 3;
thetalist_list = {};

% Initial value setting
thetalist   = [0, 0, 0, 0, 0, 0]';
dthetalist  = [0, 0, 0, 0, 0, 0]';
ddthetalist = [0, 0, 0, 0, 0, 0]';
taulist     = [0, 0, 0, 0, 0, 0]';
Ftip        = [0, 0, 0, 0, 0, 0]';
g           = [0; 0; -9.8];

for t = linspace(0,endTime,endTime/dt)

    % Inverse Dynamics for control
    mMat = MassMatrix(thetalist, Mlist, Glist, Slist);
    cMat = VelQuadraticForces(thetalist, dthetalist, Mlist, Glist, Slist);
    gMat = GravityForces(thetalist, g, Mlist, Glist, Slist);
    
    % Controller

    % Forward Dynamics for simulation
    ddthetalist = ForwardDynamics(thetalist, dthetalist, taulist, g, Ftip, Mlist, Glist, Slist);
    
    % Euler's Method
    [thetalistNext, dthetalistNext] = EulerStep(thetalist, dthetalist, ddthetalist, dt);
    thetalist  = thetalistNext;
    dthetalist = dthetalistNext;
    thetalist_list{end+1} = thetalist;
end

%% Draw robot
config = homeConfiguration(robot);
f1 = figure;
for i = 1 : 1 : length(thetalist_list)
    thetalist = thetalist_list{i};
    for j = 1 : 1 : length(thetalist)
        config(j).JointPosition = thetalist(j);
    end
    show(robot,config,"FastUpdate",1,"PreservePlot",0);
    %view(120,35)
    drawnow
end