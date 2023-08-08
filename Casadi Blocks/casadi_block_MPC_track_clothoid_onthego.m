classdef casadi_block_MPC_track_clothoid_onthego < matlab.System & matlab.system.mixin.Propagates & matlab.system.mixin.SampleTime
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
            num = 14;
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
        function [dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8,dt9,dt10,dt11,dt12,dt13,dt14] = getInputDataTypeImpl(~)
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
            dt11 = 'double';
            dt12 = 'double';
            dt13 = 'double';
            dt14 = 'double';
        end
        function [sz1,sz2,sz3,sz4,sz5] = getOutputSizeImpl(~)
            sz1 = [1,1];
            sz2 = [1,1];
            sz3 = [1,1];
            sz4 = [1,1];
            sz5 = [1,1];
        end
        function [sz1,sz2,sz3,sz4,sz5,sz6,sz7,sz8,sz9,sz10,sz11,sz12,sz13,sz14] = getInputSizeImpl(~)
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
            sz11 = 1;
            sz12 = 1;
            sz13 = 1;
            sz14 = 1;
        end
        function [cp1,cp2,cp3,cp4,cp5,cp6,cp7,cp8,cp9,cp10,cp11,cp12,cp13,cp14] = isInputComplexImpl(~)
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
            cp11 = false;
            cp12 = false;
            cp13 = false;
            cp14 = false;
        end
        function [cp1,cp2,cp3,cp4,cp5] = isOutputComplexImpl(~)
            cp1 = false;
            cp2 = false;
            cp3 = false;
            cp4 = false;
            cp5 = false;
        end
        function [fz1,fz2,fz3,fz4,fz5,fz6,fz7,fz8,fz9,fz10,fz11,fz12,fz13,fz14] = isInputFixedSizeImpl(~)
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
            fz11 = true;
            fz12 = true;
            fz13 = true;
            fz14 = true;
        end
        function [fz1,fz2,fz3,fz4,fz5] = isOutputFixedSizeImpl(~)
            fz1 = true;
            fz2 = true;
            fz3 = true;
            fz4 = true;
            fz5 = true;
        end
        function setupImpl(casobj,xd_in,yd_in,yawd_in,yaw_in,x_in,y_in,~,~,~,~,~,~,~,vv_ref)
            % Implement tasks that need to be performed only once, such as pre-computed constants.
            
            addpath('C:\Users\z041288\Desktop\MPC_implementation');
            
            % CasADi v3.5.5
            addpath('C:\Users\z041288\Desktop\MPC_implementation\casadi-windows-matlabR2016a-v3.5.5');
            addpath('C:\Users\z041288\Desktop\MPC_implementation\clothoid_fitting_git\G1fitting');
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
            casobj.N = 9; %% 15 cuts now delays a bit
            casobj.Nc = 3;
            
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
            
            Q = diag([50*vv_ref/20;0;0;50;240;300]); % weighing matrices (states) % high speed penalty causes spike in yaw rate at fast segment
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
            args.ubx(1:6:6*(casobj.N+1),1) = 15; %states upper bound - xdot
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
            casobj.mv = zeros(casobj.Nc,2);       % 2 control inputs for each state
            casobj.X0 = repmat(casobj.x0,1,casobj.N+1)';
            casobj.target =1;
            casobj.target_dist=100;
            casobj.tocs=[];
            casobj.c1 = 0; casobj.c2 = 0; casobj.c3 = 0; casobj.c4 = 0; casobj.c5 = 0;
        end
     
        %% function executed at every iteration %%
        function [delta_out,acc_out,break_out,vref,err] = stepImpl(casobj,xd_in,yd_in,yawd_in,yaw_in,x_in,y_in,t,x_s_in, y_s_in, t_s_in, x_s_out, y_s_out, t_s_out, vv_ref)
            %% transformation to front axle %%           
            x_in = x_in + 2.66*cos(yaw_in);
            y_in = y_in + 2.66*sin(yaw_in);            
            %% generate reference path %%
            if casobj.c1 == 0                  
                degtorad = pi/180;                % degree to radian conversion;                 

                t_in_rad = t_s_in*degtorad; 
                t_out_rad = t_s_out*degtorad;               
                
                [k1,dk1,L1] = buildClothoid(x_s_in, y_s_in, t_in_rad, x_s_out, y_s_out, t_out_rad); 
                npts_cl = 3*round(L_clothoid*3*1.33/vv_ref);
                XY1 = pointsOnClothoid(x_s_in, y_s_in, t_in_rad, k1, dk1, L1, npts_cl) ; 
                casobj.comb_x = XY1(1,:);
                casobj.comb_y = XY1(2,:);
                casobj.path = [casobj.comb_x', casobj.comb_y'];
                casobj.c1 = 1;
                plot(casobj.comb_x, casobj.comb_y, 'or');                         
                % one can add initial and final padding based on the known
                % initial and final heading and co-ordinates - linear only
            end
            
            err=0;
            err_ve=0;
            
            tic
            if t==0
                x_0=[0.001;yd_in;yawd_in;yaw_in;x_in;y_in];    
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
            yawflag = 0;
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
                yawdot_ref = -(pi/72)/0.01;
                xdot_ref = vv_ref;
                ydot_ref = 0;
%                 if abs(yr-yaw_in) > 5*pi/180 && yawflag == 0
%                     yawflag = i;
%                     x_ref = casobj.path_ref(yawflag,1);
%                     y_ref = casobj.path_ref(yawflag,2);                     
%                 else if yawflag ~= 0
%                         x_ref = casobj.path_ref(yawflag,1);
%                         y_ref = casobj.path_ref(yawflag,2);
%                     else                                             
%                         x_ref = casobj.path_ref(i,1);
%                         y_ref = casobj.path_ref(i,2); 
%                     end
%                 end
                x_ref = casobj.path_ref(i,1); y_ref = casobj.path_ref(i,2);
                p(6*i+1:6*i+6) = [xdot_ref, ydot_ref, yawdot_ref, yaw_ref, x_ref, y_ref];
                
                vref = vv_ref;
            end
%             figure(3);
%             plot(p(11:6:end), p(12:6:end), 'og'); hold on;
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