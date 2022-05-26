
%%%%%%%%%% Indy7 Simulation %%%%%%%%%%
clc
clear all
close all

% Load Indy7 information
addpath('MR/') 
[Slist, Mlist, Glist, M, w, p, robot] = load_urdf("indy7.urdf",6);

% Simulation
dt = 0.01;
endTime = 15;
thetalist_list = {};
err_q_list = {};
err_t_list = {};

% Initial value setting
x_init      = [0, 0, 0, 0, 0, 0]';
thetalist   = [0, 0, 0, 0, 0, 0]';
dthetalist  = [0, 0, 0, 0, 0, 0]';
ddthetalist = [0, 0, 0, 0, 0, 0]';
taulist     = [0, 0, 0, 0, 0, 0]';
Ftip        = [0, 0, 0, 0, 0, 0]';
g           = [0; 0; -9.8];

% PID gain
% Kp = 10.*[1,1,1,1,1,1]';
% Ki = []';
% Kd = 5.*[1,1,1,1,1,1]';
Kp = 10.*[1,1,1,1,1,1]';
Ki = []';
Kd = 5.*[1,1,1,1,1,1]';

% Desired Pos&Vel
x_des     = [0, 0, 0, 0, 0, 0.8]';
x_dot_des = [0, 0, 0, 0, 0, 0]';

q_des     = [0, pi/6, 0, 0, pi/3, 0]';
q_dot_des = [0, 0, 0, 0, 0, 0]';

% Create desired trajectory
method = 5;
x_start = eul2rotm(x_init(1:3).');
x_start = RpToTrans(x_start, x_init(4:6));
x_end = eul2rotm(x_des(1:3).');
x_end = RpToTrans(x_end, x_des(4:6));
Desired_trajectory = CartesianTrajectory(x_start, x_end, endTime, endTime/dt, method);

count = 0;
for t = linspace(0,endTime,endTime/dt)

    % Trajectory
    count = 1 + count;
    traj_T = cell2mat(Desired_trajectory(count));
    traj_expmat = MatrixLog6(traj_T);
    x_des = se3ToVec(traj_expmat);

    % Inverse Dynamics for control
    mMat = MassMatrix(thetalist, Mlist, Glist, Slist);
    cMat = VelQuadraticForces(thetalist, dthetalist, Mlist, Glist, Slist);
    gMat = GravityForces(thetalist, g, Mlist, Glist, Slist);
    
    %Js = JacobianSpace(Slist, thetalist);
    Jb = JacobianBody(Slist, thetalist);
    Jb = [zeros(1,6); zeros(1,6); zeros(1,6); Jb(4,:); Jb(5,:); Jb(6,:)];
    %Jm = inv(mMat)*Js'*inv((Js*inv(mMat)*Js'));
    Jb = [Jb(4,:); Jb(5,:); Jb(6,:)];
    T = FKinSpace(M, Slist, thetalist);
    expmat = MatrixLog6(T);
    x_cur = se3ToVec(expmat);
    err_t  = x_des - x_cur;

    % Joint space error
%     err_q  = q_des - thetalist;
%     err_qd = q_dot_des - dthetalist;
    err_q  = Jm*err_t;
    err_qd = Jm*x_dot_des - dthetalist;

    % Controller
    taulist = mMat*(Kp.*err_q + Kd.*err_qd) + cMat + gMat; 

    % Forward Dynamics for simulation
    ddthetalist = ForwardDynamics(thetalist, dthetalist, taulist, g, Ftip, Mlist, Glist, Slist);
    
    % Euler's Method
    [thetalistNext, dthetalistNext] = EulerStep(thetalist, dthetalist, ddthetalist, dt);
    thetalist  = thetalistNext;
    dthetalist = dthetalistNext;
    thetalist_list{end+1} = thetalist;
    err_q_list{end+1} = err_q;
    err_t_list{end+1} = err_t;
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

%% Plot error
% joint space error
err_q1 = zeros(1,length(err_q_list));
err_q2 = zeros(1,length(err_q_list));
err_q3 = zeros(1,length(err_q_list));
err_q4 = zeros(1,length(err_q_list));
err_q5 = zeros(1,length(err_q_list));
err_q6 = zeros(1,length(err_q_list));
time  = zeros(1,length(err_q_list));
f2 = figure;
for i = 1 : 1 : length(err_q_list)
    errlist_q = err_q_list{i};
    err_q1(i) = errlist_q(1);
    err_q2(i) = errlist_q(2);
    err_q3(i) = errlist_q(3);
    err_q4(i) = errlist_q(4);
    err_q5(i) = errlist_q(5);
    err_q6(i) = errlist_q(6);
    time(i)  = i*dt; 
end
subplot(6,1,1);
plot(time, err_q1)
grid on
subplot(6,1,2);
plot(time, err_q2)
grid on
subplot(6,1,3);
plot(time, err_q3)
grid on
subplot(6,1,4);
plot(time, err_q4)
grid on
subplot(6,1,5);
plot(time, err_q5)
grid on
subplot(6,1,6);
plot(time, err_q6)
grid on

% task space error
err_tx = zeros(1,length(err_t_list));
err_ty = zeros(1,length(err_t_list));
err_tz = zeros(1,length(err_t_list));
time  = zeros(1,length(err_t_list));
f3 = figure;
for i = 1 : 1 : length(err_t_list)
    errlist_t = err_t_list{i};
    err_tx(i) = errlist_t(4);
    err_ty(i) = errlist_t(5);
    err_tz(i) = errlist_t(6);
    time(i)  = i*dt; 
end
subplot(3,1,1);
plot(time, err_tx)
grid on
subplot(3,1,2);
plot(time, err_ty)
grid on
subplot(3,1,3);
plot(time, err_tz)
grid on