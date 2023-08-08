classdef casadi_block_MPC_left_turn < matlab.System & matlab.system.mixin.Propagates & matlab.system.mixin.SampleTime
    % untitled Add summary here
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.
    
    properties
        % Public, tunable properties.
        
    end
    
    properties (DiscreteState)
    end
    
    properties (Access = private)
        % Pre-computed constants.
        N
        Nc
        casadi_solver
        x0
        X0
        mv
        lbx
        ubx
        lbg
        ubg
        path
        target
        target_dist
        tocs
        errs
        er_dis
        effs
        effa
        err_v
        c1
        c2
        c3
        c4
        c5
        comb_x
        comb_y
        theta_arr
    end
    
    
    methods (Access = protected)
        function sts = getSampttpleTimeImpl(casobj)
            sts = createSampleTime(casobj,'Type','Discrete','SampleTime',0.01);
        end
        function num = getNumInputsImpl(~)
            num = 7;
        end
        function num = getNumOutputsImpl(~)
            num = 5;
        end
        function [dt1,dt2,dt3,dt4,dt5] = getOutputDataTypeImpl(~)
            dt1 = 'double';
            dt2 = 'double';
            dt3 = 'double';
            dt4 = 'double';
            dt5 = 'double';
        end
        function [dt1,dt2,dt3,dt4,dt5,dt6,dt7] = getInputDataTypeImpl(~)
            dt1 = 'double';
            dt2 = 'double';
            dt3 = 'double';
            dt4 = 'double';
            dt5 = 'double';
            dt6 = 'double';
            dt7 = 'double';
        end
        function [sz1,sz2,sz3,sz4,sz5] = getOutputSizeImpl(~)
            sz1 = [1,1];
            sz2 = [1,1];
            sz3 = [1,1];
            sz4 = [1,1];
            sz5 = [1,1];
        end
        function [sz1,sz2,sz3,sz4,sz5,sz6,sz7] = getInputSizeImpl(~)
            sz1 = 1;
            sz2 = 1;
            sz3 = 1;
            sz4 = 1;
            sz5 = 1;
            sz6 = 1;
            sz7 = 1;
        end
        function [cp1,cp2,cp3,cp4,cp5,cp6,cp7] = isInputComplexImpl(~)
            cp1 = false;
            cp2 = false;
            cp3 = false;
            cp4 = false;
            cp5 = false;
            cp6 = false;
            cp7 = false;
        end
        function [cp1,cp2,cp3,cp4,cp5] = isOutputComplexImpl(~)
            cp1 = false;
            cp2 = false;
            cp3 = false;
            cp4 = false;
            cp5 = false;
        end
        function [fz1,fz2,fz3,fz4,fz5,fz6,fz7] = isInputFixedSizeImpl(~)
            fz1 = true;
            fz2 = true;
            fz3 = true;
            fz4 = true;
            fz5 = true;
            fz6 = true;
            fz7 = true;
        end
        function [fz1,fz2,fz3,fz4,fz5] = isOutputFixedSizeImpl(~)
            fz1 = true;
            fz2 = true;
            fz3 = true;
            fz4 = true;
            fz5 = true;
        end
        function setupImpl(casobj)
            % Implement tasks that need to be performed only once, such as pre-computed constants.
            
            addpath('C:\Users\z041288\Desktop\MPC_implementation')
            
            % CasADi v3.5.5
            addpath('C:\Users\z041288\Desktop\MPC_implementation\casadi-windows-matlabR2016a-v3.5.5')
            %load('casadi_objects.mat');
            load('left_turn_data.mat');
            import casadi.*
            
            %car pars
            m = 1900;
            Iz = 750;
            lf = 1.2;
            lr = 1.6;
            cf = 70000;
            cr = 70000;
            
            %controller pars
            T =0.01; %[s]
            casobj.N = 15; %% 15 cuts now delays a bit
            casobj.Nc = 3;
            
            d_max = 0.5; d_min = -0.5;
            ax_max=1; ax_min=-1;
            
            
            import casadi.*
            % symbols
            xdot=MX.sym('xdot');
            ydot = MX.sym('ydot');
            yaw=MX.sym('yaw');
            yawdot = MX.sym('yawdot');
            X=MX.sym('X');
            Y=MX.sym('Y');
            delta = MX.sym('delta');
            ax=MX.sym('ax');
            
            % calculate lateral forces linear
            alphaf=atan2((yawdot*lf+ydot),(xdot))-delta;
            alphar=atan2((ydot-yawdot*lr),(xdot));
            fyf=-cf*alphaf;
            fyr=-cr*alphar;
            
            % states ,equations and input
            states = [xdot;ydot;yawdot;yaw;X;Y];
            n_states = length(states);
            controls =[delta ax];
            n_controls = length(controls);
            
            rhs = [ ax+ yawdot*ydot
                (1/m)*(2*fyr+2*fyf*cos(delta)-m*xdot*yawdot)
                (1/Iz)*(2*lf*fyf*cos(delta)-2*lr*fyr)
                yawdot
                xdot*cos(yaw)-ydot*sin(yaw)
                xdot*sin(yaw)+ydot*cos(yaw) ]; % system r.h.s
               
            f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
            U = SX.sym('U',n_controls,casobj.Nc); % Decision variables (controls)
            P = SX.sym('P',n_states + casobj.N*(n_states)); % parameters (which include the initial state and the reference state)
            X = SX.sym('X',n_states,(casobj.N+1)); % A vector that represents the states over the optimization problem.
            
            obj = 0; % Objective funcion
            g = [];  % constraints vector
            
            Q = diag([50;0;0;50;240;300]); % weighing matrices (states) % high speed penalty causes spike in yaw rate at fast segment
            R =diag([1510;125]); % weighing matrices (controls rates)
            R1=diag([0;0.01]);
            st  = X(:,1); % initial state
            g = [g;st-P(1:6)]; % initial condition constraints
            for k = 1:casobj.N
                st = X(:,k);
                if k<casobj.Nc
                    con=U(:,k);
                else
                    con=U(:,casobj.Nc);
                end
                obj = obj+((st-P(6*k+1:6*k+6))'*Q*(st-P(6*k+1:6*k+6))); %%%targets function
                st_next = X(:,k+1);
                f_value = f(st,con);
                st_next_euler = st+ (0.01*f_value);
                g = [g;st_next-st_next_euler;]; % compute constraints
            end
            
            con_prev=0;
            for k=1:casobj.Nc
                if k<casobj.Nc
                    con=U(:,k);
                else
                    con=U(:,casobj.Nc);
                end
                con_rate=con-con_prev;
                obj = obj+con'*R1*con+con_rate'*R*con_rate;
                con_prev=con;
            end
            
            OPT_variables = [reshape(X,6*(casobj.N+1),1);reshape(U,2*casobj.Nc,1)];
            nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);
            
            opts = struct;
            opts.ipopt.max_iter = 2000;
            opts.ipopt.print_level =0;%0,3
            opts.print_time = 0;
            opts.ipopt.acceptable_tol =1e-8;
            opts.ipopt.acceptable_obj_change_tol = 1e-6;
            
            solver = nlpsol('solver', 'ipopt', nlp_prob,opts);
            
            
            args = struct;
            args.lbg(1:6*(casobj.N+1)) =0; % Equality constraints
            args.ubg(1:6*(casobj.N+1)) = 0  ; % Equality constraints
            
            args.lbx(1:6:6*(casobj.N+1),1) = 0; %states lower bound - xdot
            args.ubx(1:6:6*(casobj.N+1),1) = 10; %states upper bound - xdot
            args.lbx(2:6:6*(casobj.N+1),1) = -inf; %states lower bound - ydot
            args.ubx(2:6:6*(casobj.N+1),1) = inf; %states upper bound - ydot
            args.lbx(3:6:6*(casobj.N+1),1) = -inf; %states lower bound - yawrate
            args.ubx(3:6:6*(casobj.N+1),1) = inf; %states upper bound - yawrate
            args.lbx(4:6:6*(casobj.N+1),1) = -inf; %states lower bound - yaw
            args.ubx(4:6:6*(casobj.N+1),1) = inf; %states upper bound - yaw
            args.lbx(5:6:6*(casobj.N+1),1) = -inf; %states lower bound - X
            args.ubx(5:6:6*(casobj.N+1),1) = inf; %states upper bound - X
            args.lbx(6:6:6*(casobj.N+1),1) = -inf; %states lower bound - Y
            args.ubx(6:6:6*(casobj.N+1),1) = inf; %states upper bound - Y
            args.lbx(6*(casobj.N+1)+1:2:6*(casobj.N+1)+2*casobj.Nc,1) = d_min; %delta lower bound
            args.ubx(6*(casobj.N+1)+1:2:6*(casobj.N+1)+2*casobj.Nc,1) = d_max; %delta upper bound
            args.lbx(6*(casobj.N+1)+2:2:6*(casobj.N+1)+2*casobj.Nc,1) = ax_min; %acceleration lower bound
            args.ubx(6*(casobj.N+1)+2:2:6*(casobj.N+1)+2*casobj.Nc,1) = ax_max; %acceleration upper bound
            
            casobj.casadi_solver = solver;
            casobj.lbx = args.lbx;
            casobj.ubx = args.ubx;
            casobj.lbg = args.lbg;
            casobj.ubg = args.lbg;
            % intial states, paths
            % state : [xdot;ydot;yawdot;yaw;X;Y]
            % control : [delta;ax]
            casobj.x0=[10*(5/18);0;0;0;1.82;-60];
            casobj.comb_x = comb_x; casobj.comb_y = comb_y;
            casobj.path=[casobj.comb_x',casobj.comb_y'];
            
            casobj.mv = zeros(casobj.Nc,2);       % 2 control inputs for each state
            casobj.X0 = repmat(casobj.x0,1,casobj.N+1)';
            casobj.target =1;
            casobj.target_dist=100;
            casobj.tocs=[];
            casobj.c1 = 0; casobj.c2 = 0; casobj.c3 = 0; casobj.c4 = 0; casobj.c5 = 0;
            casobj.theta_arr = theta_arr; 
            
     end
        function [delta_out,acc_out,break_out,vref,err] = stepImpl(casobj,xd_in,yd_in,yawd_in,yaw_in,x_in,y_in,t)
             err=0;
             err_ve=0;
            
            tic
            if t==0
                x_0=[0.001;yd_in;yawd_in;yaw_in;x_in;y_in];
                
            else
                
                x_0=[xd_in;yd_in;yawd_in;yaw_in;x_in;y_in];
            end
            % x_in = x_in + 2.666*cos(yaw_in); y_in = y_in + 2.666*sin(yaw_in);
            p(1:6)= x_0; % set the values of the parameters vector
            
            %% main loop %%
            
                       
            path_ref = zeros(casobj.N,2);
            if y_in < -5                
                next_N_ind = find(casobj.path(:,2)>y_in, casobj.N);
                if size(next_N_ind,1) < casobj.N 
                    len_N = size(next_N_ind,1);
                    path_ref(1:len_N,1) = casobj.path(next_N_ind,1);
                    path_ref(1:len_N,2) = casobj.path(next_N_ind,2);
                    path_ref(len_N+1:end,1) = casobj.path(end,1);
                    path_ref(len_N+1:end,2) = casobj.path(end,2);
                else
                    path_ref(1:casobj.N,1) = casobj.path(next_N_ind,1);
                    path_ref(1:casobj.N,2) = casobj.path(next_N_ind,2);
                end
            else if y_in >= -5
                    next_N_ind = find(casobj.path(:,1)<x_in, casobj.N);
                    if size(next_N_ind,1) < casobj.N 
                        len_N = size(next_N_ind,1);
                        path_ref(1:len_N,1) = casobj.path(next_N_ind,1);
                        path_ref(1:len_N,2) = casobj.path(next_N_ind,2);
                        path_ref(len_N+1:end,1) = casobj.path(end,1);
                        path_ref(len_N+1:end,2) = casobj.path(end,2);
                    else
                        path_ref(1:casobj.N,1) = casobj.path(next_N_ind,1);
                        path_ref(1:casobj.N,2) = casobj.path(next_N_ind,2);
                    end
                end
            end
                    
            
            
            %plot(path_ref(:,1), path_ref(:,2));
            
            for i=1:casobj.N
                
                yr = 0;
                if i == 1 
                    if path_ref(i,1)-x_in == 0
                        yr = pi/2;
                    else
                        yr = atan(((path_ref(i,2)-y_in)/(path_ref(i,1)-x_in)));
                    end                    
                else
                    if path_ref(i,1)-path_ref(i-1,1) == 0
                        yr = pi/2;
                    else
                        yr = atan(((path_ref(i,2)-path_ref(i-1,2))/(path_ref(i,1)-path_ref(i-1,1))));
                    end
                end
                if yr <= 0
                    yr = yr + pi;
                end
                yaw_ref = yr;
                yawdot_ref = 0;
                xdot_ref = 10*(5/18);
                
                x_ref = path_ref(i,1); 
                y_ref = path_ref(i,2) ; 
                ydot_ref = 0;
                %x_ref = path_ref(9,1);y_ref = path_ref(9,2);

                
                p(6*i+1:6*i+6) = [xdot_ref, ydot_ref, yawdot_ref, yaw_ref, x_ref, y_ref];
                
                %x_ref 
                %y_ref 
                %yaw_ref
                %casobj.exit
                vref = 10*(5/18);
%                 k=j-1;
%                 if k==0
%                     k=size(casobj.path,1);
%                 end
%                 p(6*i+4:6*i+6) = [atan2((casobj.path(j,2)-casobj.path(k,2)),(casobj.path(j,1)-casobj.path(k,1))) , casobj.path(j,1), casobj.path(j,2)];
            end
%            pos_flag
            %casobj.c5
            y_ref
            y_in
            casobj.err_v=[casobj.err_v;err_ve];
            % initial value of the optimization variables
            x0_solve  = [reshape(casobj.X0',6*(casobj.N+1),1);reshape(casobj.mv',2*casobj.Nc,1)];
            sol = casobj.casadi_solver('x0', x0_solve, 'lbx', casobj.lbx, 'ubx', casobj.ubx,...
                'lbg', casobj.lbg, 'ubg', casobj.ubg,'p',p);
            u = reshape(full(sol.x(6*(casobj.N+1)+1:end))',2,casobj.Nc)'; % get controls only from the solution
            casobj.mv = [u(2:size(u,1),:);u(size(u,1),:)];
            casobj.X0 = reshape(full(sol.x(1:6*(casobj.N+1)))',6,casobj.N+1)'; % get solution TRAJECTORY
            casobj.X0 = [casobj.X0(2:end,:);casobj.X0(end,:)];
            %u(1,1)
           
            if t<0.01
                acc_out=0.01;
                delta_out=0;
                break_out=0;
            else
                delta_out=rad2deg(u(1,1))*30;
                if (u(1,2) >=0 )
                    acc_out=u(1,2); % engine acceleration
                    break_out =0; % braking acceleration
                else
                    acc_out= 0;
                    break_out=-u(1,2) ;
                end
            end

             toc  
          casobj.tocs=[casobj.tocs;toc];
          casobj.errs=[casobj.errs;err];  
          casobj.tocs=[casobj.tocs;toc];
          casobj.effs=[casobj.effs;u(1,1)];
          casobj.effa=[casobj.effa;u(1,2)];
          if (t>408.5)
        mean( casobj.errs)
        std(casobj.errs)
        max(casobj.errs)
        mean( casobj.tocs)
        std(casobj.tocs)
        max(casobj.tocs)
        eff_st=sum(casobj.effs.^2)
        eff_ca=sum(casobj.effa.^2)
          end
        end   
    end
end