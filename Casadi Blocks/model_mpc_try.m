% Trajectory Tracking + Multiple shooting
clear all
close all
clc

% CasADi v3.4.5
% addpath('C:\Users\mehre\OneDrive\Desktop\CasADi\casadi-windows-matlabR2016a-v3.4.5')
% CasADi v3.5.5
addpath('C:\Users\z041288\Downloads\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

T = 0.1; %[s]
N = 10; % prediction horizon
rob_diam = 0.3;

Ts = 0.1; 
m = 1412; 
Iz = 1536.7; 
kf = -128916; 
kr = -85944;
lf = 1.06; 
lr = 1.85; 


v_max = 0.6; v_min = -v_max; vl_max = 0.1; vl_min = -vl_max;
omega_max = pi/4; omega_min = -omega_max;
a_min = -0.2; a_max = 0.2; 
s_min = -0.2; s_max = 0.2;

x = SX.sym('x'); y = SX.sym('y'); theta = SX.sym('theta');
vlongi = SX.sym('vlongi'); vlat = SX.sym('vlat'); omega = SX.sym('omega');
states = [x;y;theta;vlongi;vlat;omega]; n_states = length(states);

a = SX.sym('a'); steer = SX.sym('steer');
controls = [a;steer]; n_controls = length(controls);
rhs = [vlongi*cos(theta) - vlat*sin(theta);vlat*cos(theta) + vlongi*sin(theta);...
    omega;a;((m*vlongi*vlat + Ts*(lf*kf-lr*kr)*omega - Ts*kf*steer*vlongi - Ts*m*vlongi*vlongi*omega)/(m*vlongi + Ts*(kf+kr)) - vlat)/Ts;...
    ((Iz*vlongi*omega + Ts*vlat*(lf*kf - lr*kr) - Ts*kf*lf*steer*vlongi)/(Iz*vlongi - Ts*(lf*lf*kf + lr*lr*kr)) - omega)/Ts]; % system r.h.s



f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_controls,N); % Decision variables (controls)
%P = SX.sym('P',n_states + n_states);
P = SX.sym('P',n_states + N*(n_states+n_controls));
% parameters (which include the initial state and the reference along the
% predicted trajectory (reference states and reference controls))

X = SX.sym('X',n_states,(N+1));
% A vector that represents the states over the optimization problem.

obj = 0; % Objective function
g = [];  % constraints vector

%Q = zeros(3,3); Q(1,1) = 1;Q(2,2) = 1;Q(3,3) = 0.5; % weighing matrices (states)
%R = zeros(2,2); R(1,1) = 0.5; R(2,2) = 0.05; % weighing matrices (controls)
Q = diag([1 1 0.5 0.01 0.01 0.01]);
R = diag([0.5 0.05]);

st  = X(:,1); % initial state
g = [g;st-P(1:n_states)]; % initial condition constraints
num = n_states + n_controls;
for k = 1:N
    st = X(:,k);  con = U(:,k);    
    %obj = obj+(st-P(4:6))'*Q*(st-P(4:6)) + con'*R*con; % calculate obj
    obj = obj+(st-P(num*k-1:num*k+4))'*Q*(st-P(num*k-1:num*k+4)) + ...
              (con-P(num*k+5:num*k+6))'*R*(con-P(num*k+5:num*k+6)) ; % calculate obj
    % the number 5 is (n_states+n_controlss)
    st_next = X(:,k+1);
    f_value = f(st,con);
    st_next_euler = st+ (T*f_value);
    g = [g;st_next-st_next_euler]; % compute constraints
end
% make the decision variable one column  vector
OPT_variables = [reshape(X,6*(N+1),1);reshape(U,2*N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;

args.lbg(1:n_states*(N+1)) = 0;  % -1e-20  % Equality constraints
args.ubg(1:n_states*(N+1)) = 0;  % 1e-20   % Equality constraints

args.lbx(1:n_states:n_states*(N+1),1) = -20; %state x lower bound % new - adapt the bound
args.ubx(1:n_states:n_states*(N+1),1) = 20; %state x upper bound  % new - adapt the bound
args.lbx(2:n_states:n_states*(N+1),1) = -2; %state y lower bound
args.ubx(2:n_states:n_states*(N+1),1) = 2; %state y upper bound
args.lbx(3:n_states:n_states*(N+1),1) = -inf; %state theta lower bound
args.ubx(3:n_states:n_states*(N+1),1) = inf; %state theta upper bound
args.lbx(4:n_states:n_states*(N+1),1) = v_min; %state vlongi lower bound
args.ubx(4:n_states:n_states*(N+1),1) = v_max; %state vlongi upper bound
args.lbx(5:n_states:n_states*(N+1),1) = vl_min; %state vlat lower bound
args.ubx(5:n_states:n_states*(N+1),1) = vl_max; %state vlat upper bound 
args.lbx(6:n_states:n_states*(N+1),1) = omega_min; %state omega lower bound
args.ubx(6:n_states:n_states*(N+1),1) = omega_max; %state omega upper bound

args.lbx(n_states*(N+1)+1:n_controls:n_states*(N+1)+n_controls*N,1) = a_min; %a lower bound
args.ubx(n_states*(N+1)+1:n_controls:n_states*(N+1)+n_controls*N,1) = a_max; %a upper bound
args.lbx(n_states*(N+1)+2:n_controls:n_states*(N+1)+n_controls*N,1) = s_min; %steer lower bound
args.ubx(n_states*(N+1)+2:n_controls:n_states*(N+1)+n_controls*N,1) = s_max; %steer upper bound
%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
x0 = [0; 0; 0; 0; 0; 0];    % initial condition.
xs = [1.5 ; 1.5; 0; 0; 0; 0]; % Reference posture.

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,2);        % two control inputs for each robot
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables

sim_tim = 50; % Maximum simulation time

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];

% the main simulaton loop... it works as long as the error is greater
% than 10^-6 and the number of mpc steps is less than its maximum
% value.
main_loop = tic;
while(mpciter < sim_tim / T) % new - condition for ending the loop
    current_time = mpciter*T;  %new - get the current time
    % args.p   = [x0;xs]; % set the values of the parameters vector
    %----------------------------------------------------------------------
    args.p(1:n_states) = x0; % initial condition of the robot posture
    for k = 1:N %new - set the reference to track
        t_predict = current_time + (k-1)*T; % predicted time instant
        x_ref = 0.5*t_predict; y_ref = 1; theta_ref = 0; vlongi_ref = 0.5; vlat_ref = 0; omega_ref = 0;
        a_ref = 0.1; s_ref = 0;
        if x_ref >= 12 % the trajectory end is reached
            x_ref = 12; y_ref = 1; theta_ref = 0; omega_ref = 0; vlongi_ref = 0;
            a_ref = 0; s_ref = 0;
        end
        args.p(num*k-1:num*k+4) = [x_ref, y_ref, theta_ref, vlongi_ref, vlat_ref, omega_ref];
        args.p(num*k+5:num*k+6) = [a_ref, s_ref];
    end
    %----------------------------------------------------------------------    
    % initial value of the optimization variables
    args.x0  = [reshape(X0',n_states*(N+1),1);reshape(u0',n_controls*N,1)];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    u = reshape(full(sol.x(n_states*(N+1)+1:end))',n_controls,N)'; % get controls only from the solution
    xx1(:,1:n_states,mpciter+1)= reshape(full(sol.x(1:n_states*(N+1)))',n_states,N+1)'; % get solution TRAJECTORY
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    % Apply the control and shift the solution
    [t0, x0, u0] = shift(T, t0, x0, u,f);
    xx(:,mpciter+2) = x0;
    X0 = reshape(full(sol.x(1:n_states*(N+1)))',n_states,N+1)'; % get solution TRAJECTORY
    % Shift trajectory to initialize the next step
    X0 = [X0(2:end,:);X0(end,:)];
    mpciter
    mpciter = mpciter + 1;
end;
main_loop_time = toc(main_loop);
average_mpc_time = main_loop_time/(mpciter+1)

Draw_MPC_tracking_v1 (t,xx,xx1,u_cl,xs,N,rob_diam)