
%%%%%%%%%% Indy7 Task-Space Motion Control Simulation %%%%%%%%%%
clc
clear all
close all

% Load Indy7 information
addpath('MR/') 
[Slist, Mlist, Glist, M, w, p, robot] = load_urdf("indy7.urdf",6);

% Simulation parameter
dt = 0.01;
endTime = 5;

% Initial value setting
x_init      = [0, 0, 0, -0.0000, -0.1865, 1.3275]'; % initial end-effector position
thetalist   = [0, 0, 0, 0, 0, 0]';
dthetalist  = [0, 0, 0, 0, 0, 0]';
prev_Xe     = [0, 0, 0, 0, 0, 0]';
prev_err_p  = [0, 0, 0]'; 
thetalist_list = {};
err_t_list     = {};
count = 0;

% PI gain
Kp = 10.*[1,1,1,1,1,1]';
Ki = 25.*[1,1,1,1,1,1]';

% Desired Position
x_des     = [0, 0, 0, 0.5, -0.5, 0.5]';

% Create desired trajectory
method = 5;
x_start = eul2rotm(x_init(1:3).');
x_start = RpToTrans(x_start, x_init(4:6));
x_end = eul2rotm(x_des(1:3).');
x_end = RpToTrans(x_end, x_des(4:6));
Desired_trajectory = CartesianTrajectory(x_start, x_end, endTime, endTime/dt, method);

for t = linspace(0,endTime,endTime/dt)

    % Desired trajectory
    count = 1 + count;
    traj_T = cell2mat(Desired_trajectory(count));
    [traj_R, traj_p] = TransToRp(traj_T);

    % Forward kinematics
    T = FKinSpace(M, Slist, thetalist);
    [R, p] = TransToRp(T);

    % Body jacobian
    Js = JacobianSpace(Slist, thetalist);
    invT = TransInv(T);
    Jb = Adjoint(invT)*Js;            % 6-DOF
    Jb = [Jb(4,:); Jb(5,:); Jb(6,:)]; % 3-DOF

    % Hybrid configuration representation of Task-Space error
    err_p = traj_p - p;
    err_R = so3ToVec(MatrixLog3(R'*traj_R));
    Xe = [err_R; err_p]; 

    % Task-Space motion controller
    % [6-DOF]
%     dthetalist = pinv(Jb)*(Kp.*Xe + Ki.*(prev_Xe + dt*Xe)); % 6-DOF
%     prev_Xe = Xe;

    % [3-DOF]
    dthetalist = pinv(Jb)*(Kp(4:6).*err_p + Ki(4:6).*(prev_err_p + dt*err_p)); % 3-DOF
    prev_err_p = err_p;

    % Euler's method
    thetalist  = thetalist + dt*dthetalist;
    thetalist_list{end+1} = thetalist;
    err_t_list{end+1} = x_des(4:6) - p;

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
% task space error
err_tx = zeros(1,length(err_t_list));
err_ty = zeros(1,length(err_t_list));
err_tz = zeros(1,length(err_t_list));
time   = zeros(1,length(err_t_list));

for i = 1 : 1 : length(err_t_list)
    errlist_t = err_t_list{i};
    err_tx(i) = errlist_t(1);
    err_ty(i) = errlist_t(2);
    err_tz(i) = errlist_t(3);
    time(i)  = i*dt; 
end

f2 = figure;
subplot(3,1,1);
plot(time, err_tx)
title('X err')
ylim([-1 1])
grid on
subplot(3,1,2);
plot(time, err_ty)
title('Y err')
ylim([-1 1])
grid on
subplot(3,1,3);
plot(time, err_tz)
title('Z err')
ylim([-1 1])
grid on
sgtitle('Task space error')
