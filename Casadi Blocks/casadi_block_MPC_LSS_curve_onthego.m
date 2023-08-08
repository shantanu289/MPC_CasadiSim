classdef casadi_block_MPC_LSS_curve_onthego < matlab.System & matlab.system.mixin.Propagates & matlab.system.mixin.SampleTime
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
        path_ref
    end
    
    
    methods (Access = protected)
        function sts = getSampttpleTimeImpl(casobj)
            sts = createSampleTime(casobj,'Type','Discrete','SampleTime',0.01);
        end
        function num = getNumInputsImpl(~)
            num = 10;
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
        function [dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8,dt9,dt10] = getInputDataTypeImpl(~)
            dt1 = 'double';
            dt2 = 'double';
            dt3 = 'double';
            dt4 = 'double';
            dt5 = 'double';
            dt6 = 'double';
            dt7 = 'double';
            dt8 = 'double';
            dt9 = 'double';
            dt10 = 'double';
        end
        function [sz1,sz2,sz3,sz4,sz5] = getOutputSizeImpl(~)
            sz1 = [1,1];
            sz2 = [1,1];
            sz3 = [1,1];
            sz4 = [1,1];
            sz5 = [1,1];
        end
        function [sz1,sz2,sz3,sz4,sz5,sz6,sz7,sz8,sz9,sz10] = getInputSizeImpl(~)
            sz1 = 1;
            sz2 = 1;
            sz3 = 1;
            sz4 = 1;
            sz5 = 1;
            sz6 = 1;
            sz7 = 1;
            sz8 = 1;
            sz9 = 1;
            sz10 = 1;
        end
        function [cp1,cp2,cp3,cp4,cp5,cp6,cp7,cp8,cp9,cp10] = isInputComplexImpl(~)
            cp1 = false;
            cp2 = false;
            cp3 = false;
            cp4 = false;
            cp5 = false;
            cp6 = false;
            cp7 = false;
            cp8 = false; 
            cp9 = false; 
            cp10 = false;
        end
        function [cp1,cp2,cp3,cp4,cp5] = isOutputComplexImpl(~)
            cp1 = false;
            cp2 = false;
            cp3 = false;
            cp4 = false;
            cp5 = false;
        end
        function [fz1,fz2,fz3,fz4,fz5,fz6,fz7,fz8,fz9,fz10] = isInputFixedSizeImpl(~)
            fz1 = true;
            fz2 = true;
            fz3 = true;
            fz4 = true;
            fz5 = true;
            fz6 = true;
            fz7 = true;
            fz8 = true; 
            fz9 = true;
            fz10 = true;
        end
        function [fz1,fz2,fz3,fz4,fz5] = isOutputFixedSizeImpl(~)
            fz1 = true;
            fz2 = true;
            fz3 = true;
            fz4 = true;
            fz5 = true;
        end
        function setupImpl(casobj,xd_in,yd_in,yawd_in,yaw_in,x_in,y_in,~,~,~,vlongi_ref)
            % Implement tasks that need to be performed only once, such as pre-computed constants.
            
            addpath('C:\Users\z041288\Desktop\MPC_implementation')
            
            % CasADi v3.5.5
            addpath('C:\Users\z041288\Downloads\casadi-windows-matlabR2016a-v3.5.5')
            %load('casadi_objects.mat');
            %load('mpc_lss_data_r1200_lat_0.5_v40.mat');
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
            casobj.N = 6; %% 15 cuts now delays a bit
            casobj.Nc = 2;
            
            
            d_max = 0.34; d_min = -0.34;
            ax_max=0.3; ax_min=-1;           
            
            
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
            rhs = rhs + 0.1*rhs.*rand(size(rhs)) ;   
            f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
            U = SX.sym('U',n_controls,casobj.Nc); % Decision variables (controls)
            P = SX.sym('P',n_states + casobj.N*(n_states)); % parameters (which include the initial state and the reference state)
            X = SX.sym('X',n_states,(casobj.N+1)); % A vector that represents the states over the optimization problem.
            
            obj = 0; % Objective funcion
            g = [];  % constraints vector
            
            Q = diag([50*vlongi_ref/20;0;0;50;240;300]); % weighing matrices (states) % high speed penalty causes spike in yaw rate at fast segment
            R =diag([1510;125]); % weighing matrices (controls rates)
            R1=diag([0;0]);
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
            args.ubx(1:6:6*(casobj.N+1),1) = 30; %states upper bound - xdot
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
            casobj.x0=[xd_in;yd_in;yawd_in;yaw_in;x_in;y_in];
            %casobj.comb_x = comb_x; casobj.comb_y = comb_y;
            %casobj.path=[casobj.comb_x',casobj.comb_y'];
            
            casobj.mv = zeros(casobj.Nc,2);       % 2 control inputs for each state
            casobj.X0 = repmat(casobj.x0,1,casobj.N+1)';
            casobj.target =1;
            casobj.target_dist=100;
            casobj.tocs=[];
            casobj.c1 = 0; casobj.c2 = 0; casobj.c3 = 0; casobj.c4 = 0; casobj.c5 = 0;
            
            
     end
        function [delta_out,acc_out,break_out,vref,err] = stepImpl(casobj,xd_in,yd_in,yawd_in,yaw_in,x_in,y_in,t,R,vlat_ref,vlongi_ref)
            %% rear axle to front axle transformation
            x_in = x_in + 2.66*cos(yaw_in);
            y_in = y_in + 2.66*sin(yaw_in);
            %% generate path %%
            if casobj.c1 == 0 
                theta = atan(vlat_ref/vlongi_ref); % radians 
                psi = atan(y_in/x_in);             % radians
                ang_vel = yawd_in; 
                R_new = sqrt(vlongi_ref^2)/ang_vel;
                x_centre = x_in - R_new*cos(psi + theta);
                y_centre = y_in - R_new*sin(psi + theta);
                theta_final = psi + theta + pi/2;
                npts_arc = 1.25*round(R_new*(pi/2)*3/vlongi_ref);                               
                x_start = x_in;
                theta_lin = linspace(psi+theta, theta_final, npts_arc);                
                x_arc = x_centre + R_new*cos(theta_lin);
                y_arc = y_centre + R_new*sin(theta_lin);
                casobj.comb_x = x_arc; casobj.comb_y = y_arc;                 
                casobj.path = [casobj.comb_x', casobj.comb_y'];
                plot(casobj.comb_x, casobj.comb_y, 'og');
                x_in 
                y_in
                casobj.c1 = 1;
            end
            
            
            %% --- %%
            err=0;
            err_ve=0;            
            tic
            if t==0
                x_0=[xd_in;yd_in;yawd_in;yaw_in;x_in;y_in];                
            else                
                x_0=[xd_in;yd_in;yawd_in;yaw_in;x_in;y_in];
            end            
            p(1:6)= x_0; % set the values of the parameters vector
            
            %% main loop %%         
                     
            casobj.path_ref = zeros(casobj.N, 2);
            dist_arr = sqrt((x_in-casobj.comb_x).^2 + (y_in-casobj.comb_y).^2);
            [min_dis, min_ind] = min(dist_arr);
            start_ind = min_ind + 1;
            if start_ind > size(casobj.path, 1) - casobj.N + 1 
                casobj.path_ref(1:size(casobj.path,1)+1-start_ind,:) = casobj.path(start_ind:size(casobj.path,1),:);
                casobj.path_ref(size(casobj.path,1)+1-start_ind+1:casobj.N) = casobj.path(end,:);
            else 
                casobj.path_ref(1:casobj.N, :) = casobj.path(start_ind:start_ind + casobj.N - 1, :);
            end          
            %plot(path_ref(:,1), path_ref(:,2));            
            for i=1:casobj.N                
                yr = 0;
                if i == 1 
                    if casobj.path_ref(i,1)-x_in == 0
                        yr = pi/2;
                    else
                        yr = atan(((casobj.path_ref(i,2)-y_in)/(casobj.path_ref(i,1)-x_in)));
                    end                    
                else
                    if casobj.path_ref(i,1)-casobj.path_ref(i-1,1) == 0
                        yr = pi/2;
                    else
                        yr = atan(((casobj.path_ref(i,2)-casobj.path_ref(i-1,2))/(casobj.path_ref(i,1)-casobj.path_ref(i-1,1))));
                    end
                end
                yaw_ref = yr;
                yawdot_ref = 0;
                xdot_ref = vlongi_ref;                
                x_ref = casobj.path_ref(i,1); 
                y_ref = casobj.path_ref(i,2) ; 
                ydot_ref = 0;
                p(6*i+1:6*i+6) = [xdot_ref, ydot_ref, yawdot_ref, yaw_ref, x_ref, y_ref];               
                vref = vlongi_ref;
            end
            %% perform optimization %%

            casobj.err_v=[casobj.err_v;err_ve];
            % initial value of the optimization variables
            x0_solve  = [reshape(casobj.X0',6*(casobj.N+1),1);reshape(casobj.mv',2*casobj.Nc,1)];
            sol = casobj.casadi_solver('x0', x0_solve, 'lbx', casobj.lbx, 'ubx', casobj.ubx,'lbg', casobj.lbg, 'ubg', casobj.ubg,'p',p);
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