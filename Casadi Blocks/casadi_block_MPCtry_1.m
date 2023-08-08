classdef casadi_block_MPCtry_1 < matlab.System & matlab.system.mixin.Propagates & matlab.system.mixin.SampleTime
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
        thetax
        ppx1
        ppy1
        ppx2
        ppy2
        ppx3
        ppy3
        ppx4
        ppy4
        ppx5
        ppy5
        exit
        XY1
        XY2
        XY3
        XY4
        XY5
        XY6
        XY7
        XY8
        XY9
        XY10
        arc1_x
        arc1_y
        arc2_x
        arc2_y
        arc3_x
        arc3_y
        arc4_x
        arc4_y
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
            addpath('C:\Users\z041288\Downloads\casadi-windows-matlabR2016a-v3.5.5')
            %load('casadi_objects.mat');
            load('casadi_objects_ex1.mat');
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
            casobj.Nc = 4;
            
            d_max = 0.34; d_min = -0.34;
            ax_max=0.3; ax_min=-1;
            
            
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
            
            Q = diag([30;0;0;50;165;165]); % weighing matrices (states) % high speed penalty causes spike in yaw rate at fast segment
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
            casobj.x0=[6;0;0;0;1.5;-57];
            casobj.path=[comb_x',comb_y'];
            casobj.mv = zeros(casobj.Nc,2);       % 2 control inputs for each state
            casobj.X0 = repmat(casobj.x0,1,casobj.N+1)';
            casobj.target =1;
            casobj.target_dist=100;
            casobj.tocs=[];
            casobj.c1 = 0; casobj.c2 = 0; casobj.c3 = 0; casobj.c4 = 0; casobj.c5 = 0;
            casobj.ppx1 = ppx1; casobj.ppy1 = ppy1; casobj.thetax = thetax; casobj.ppx2 = ppx2;
            casobj.ppy2 = ppy2; casobj.ppx3 = ppx3; casobj.ppy3 = ppy3; casobj.ppx4 = ppx4; casobj.ppy4 = ppy4; casobj.ppx5 = ppx5; 
            casobj.ppy5 = ppy5;  casobj.XY1 = XY1; casobj.XY2 = XY2; casobj.XY3 = XY3; casobj.XY4 = XY4; casobj.XY5 = XY5; casobj.XY6 = XY6; casobj.XY7 = XY7; casobj.XY8 = XY8;casobj.XY9 = XY9; casobj.XY10 = XY10;
            casobj.arc1_x = arc1_x; casobj.arc1_y = arc1_y; casobj.arc2_x = arc2_x; casobj.arc2_y = arc2_y; casobj.arc3_x = arc3_x; casobj.arc3_y=arc3_y; 
            casobj.arc4_x = arc4_x; casobj.arc4_y = arc4_y; %casobj.exit = exit;
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
            
%             casobj.target_dist=sqrt((x_0(5)-casobj.path( casobj.target,1))^2 + (x_0(6)-casobj.path( casobj.target,2))^2);
            p(1:6)= x_0; % set the values of the parameters vector
            % search for closest point in the next 10 points
%             for goal=casobj.target:casobj.target+1
%                 if goal>size(casobj.path,1)
%                     i=goal-size(casobj.path,1);
%                 else
%                     i=goal;
%                 end
%                 dist=sqrt((x_0(5)-casobj.path(i,1))^2 + (x_0(6)-casobj.path(i,2))^2);
%                 if dist<casobj.target_dist
%                     casobj.target=i;
%                     err=casobj.target_dist;
%                     casobj.target_dist=dist;
%                 end
%                 if  dist<0.07
%                     casobj.target= casobj.target+1;
% %                     if casobj.target> size(casobj.path,1)
% %                         casobj.target=casobj.target-size(casobj.path,1);
% %                     end
%                     err=casobj.target_dist;
%                     break;
%                 end
%             end
            for i=1:casobj.N
%                 if  casobj.target+i>size(casobj.path,1)
%                     j=casobj.target+i-size(casobj.path,1);
%                 else
%                     j= casobj.target+i-1;
%                 end 
                yaw_ref = polyval(casobj.thetax, t);
                yawdot_ref = -(pi/72)/0.01;
                xdot_ref = 5.5;
                ydot_ref = 0;
                x_ref = 0; y_ref = 0; pos_flag = 0;
                
                %% exit = 1 %%
                if 3>2  
                    if x_in <= casobj.XY1(1,end)
                        pos_flag = 1;
                    else if (x_in > casobj.XY1(1,end) && x_in <= casobj.XY2(1,end))  || (y_in > casobj.XY1(2,end) && y_in <= casobj.XY2(2,end))
                            pos_flag = 2;
                        else if x_in > casobj.XY2(1,end) && x_in <= casobj.arc1_x(end,1)
                                pos_flag = 3;
                            else if (x_in > casobj.arc1_x(end,1) && x_in <= casobj.XY3(1,end)) %|| (y_in > casobj.arc1_y(end,1) && y_in <= casobj.XY3(2,end))
                                    pos_flag = 4;
                                else if x_in > casobj.XY3(1,end) && x_in < casobj.XY4(1, end)
                                        pos_flag = 5;
                                    else if x_in >= casobj.XY4(1,end)
                                            pos_flag = 6;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                %% exit = 2 %%
                if 2>3
                    if y_in <= casobj.XY1(2,end)
                        pos_flag = 1;
                    else if y_in > casobj.XY1(2,end) && y_in <= casobj.XY2(2,end)
                            pos_flag = 2;
                        else if (y_in > casobj.XY2(2,end) && y_in <= casobj.arc2_y(end,1)) %|| (x_in >= casobj.XY2(1,end))
                                pos_flag = 3;
                            else if(y_in > casobj.arc2_y(end,1) && y_in <= casobj.XY5(2,end))
                                    pos_flag = 4;
                                else if y_in > casobj.XY5(2,end) && y_in < casobj.XY6(2, end)
                                        pos_flag = 5;
                                    else if y_in >= casobj.XY6(2,end)
                                            pos_flag = 6;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                %% exit = 3 %%
                if 2>3
                    if y_in <= casobj.XY1(2,end)
                        pos_flag = 1;
                    else if y_in > casobj.XY1(2,end) && y_in <= casobj.XY2(2,end)
                            pos_flag = 2;
                        else if (y_in > casobj.XY2(2,end) && x_in > 0) || (y_in > casobj.arc3_y(end,1) && x_in < 0)
                                pos_flag = 3;
                            else if  (y_in < casobj.arc3_y(end,1) && x_in < 0) && y_in > casobj.XY7(2,end)
                                    pos_flag = 4;
                                else if (x_in <= casobj.XY7(1,end) && x_in >= casobj.XY8(1,end)) || ((y_in < casobj.XY7(2,end) && y_in > casobj.XY8(2,end)) && x_in < 0)
                                        pos_flag = 5;
                                    else if (x_in <= casobj.XY8(1,end)) || ((y_in <= casobj.XY8(2,end)) && x_in < 0)
                                            pos_flag = 6;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                
                %% exit = 4 %%
                if 2>3
                    if y_in <= casobj.XY1(2,end) && x_in > 0
                        pos_flag = 1;
                    else if (y_in > casobj.XY1(2,end) && y_in <= casobj.XY2(2,end)) && (x_in > 0)
                            pos_flag = 2;
                        else if y_in > casobj.XY2(2,end)
                                pos_flag = 3;
                            else if (y_in < casobj.arc4_y(end,1) && y_in >= casobj.XY9(2,end)) && (x_in < 0)
                                    pos_flag = 4;
                                else if (y_in >= casobj.XY10(2,end) && y_in < casobj.XY9(2,end)) && (x_in < 0)
                                        pos_flag = 5;
                                    else if y_in < casobj.XY10(2,end) && x_in < 0
                                            pos_flag = 6;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                
                %% exit 1 %% 
                if 3>2
                    if pos_flag == 1  % the trajectory end is reached
                        casobj.c1 = casobj.c1 + 0.01/2;
                        x_ref = polyval(casobj.ppx1, casobj.c1); y_ref = polyval(casobj.ppy1, casobj.c1);
                        x_ref = (x_ref<=casobj.XY1(1,1))*casobj.XY1(1,1) + ((x_ref <= casobj.XY1(1,end)) && x_ref>casobj.XY1(1,1))*x_ref + ((x_ref > casobj.XY1(1,end)) || casobj.c1>15)*casobj.XY1(1,end);
                        y_ref = (y_ref<=casobj.XY1(2,1) && casobj.c1<15)*casobj.XY1(2,end) + ((y_ref <= casobj.XY1(2,end)) && y_ref>casobj.XY1(2,1) && casobj.c1<15)*y_ref + (casobj.c1>15)*casobj.XY1(2,end);
                    else if pos_flag == 2
                            casobj.c2 = casobj.c2 + 0.01/2.2;
                            x_ref = polyval(casobj.ppx2, casobj.c2); y_ref = polyval(casobj.ppy2, casobj.c2);
                            x_ref = (x_ref <= casobj.XY2(1,end))*x_ref + (x_ref > casobj.XY2(1,end))*casobj.XY1(1,end);
                            y_ref = ((y_ref <= casobj.XY2(2,end)) && casobj.c2<=15)*y_ref + ((y_ref > casobj.XY2(2,end))||casobj.c2>15)*casobj.XY1(2,end);
                        else if pos_flag == 3
                                casobj.c3 = casobj.c3 + 0.01/2.2 ;
                                x_ref = polyval(casobj.ppx3, casobj.c3); y_ref = polyval(casobj.ppy3, casobj.c3);
                                x_ref = (x_ref <= casobj.arc1_x(end,1))*x_ref + (x_ref > casobj.arc1_x(end,1))*casobj.arc1_x(end,1);
                                y_ref = (y_ref <= casobj.arc1_y(end,1))*y_ref + (y_ref > casobj.arc1_y(end,1))*casobj.arc1_y(end,1);
                            else if pos_flag == 4
                                    casobj.c4 = casobj.c4 + 0.01/2 ;
                                    x_ref = polyval(casobj.ppx4, casobj.c4); y_ref = polyval(casobj.ppy4, casobj.c4);
                                    x_ref = (x_ref <= casobj.XY3(1,end))*x_ref + (x_ref > casobj.XY3(1,end))*casobj.XY3(1,end);
                                    y_ref = (y_ref <= casobj.XY3(2,end))*y_ref + (y_ref > casobj.XY3(2,end))*casobj.XY3(2,end);
                                else if pos_flag == 5
                                        casobj.c5 = casobj.c5 + 0.01/2.5 ;
                                        x_ref = polyval(casobj.ppx5, casobj.c5); y_ref = polyval(casobj.ppy5, casobj.c5) ;
                                    else if pos_flag == 6
                                            x_ref = casobj.XY4(1,end); y_ref = casobj.XY4(2,end); yaw_ref = 0; xdot_ref = 0; ydot_ref = 0; yawdot_ref = 0;
                                            flag = 1;
                                            break;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                %% exit 2 %%
                if 2>3
                    if pos_flag == 1  % the trajectory end is reached
                        casobj.c1 = casobj.c1 + 0.01/3;
                        x_ref = polyval(casobj.ppx1, casobj.c1); y_ref = polyval(casobj.ppy1, casobj.c1);
                        x_ref = (x_ref <= casobj.XY1(1,end))*x_ref + (x_ref > casobj.XY1(1,end))*casobj.XY1(1,end);
                        y_ref = (y_ref <= casobj.XY1(2,end))*y_ref + (y_ref > casobj.XY1(2,end))*casobj.XY1(2,end);
                    else if pos_flag == 2
                            casobj.c2 = casobj.c2 + 0.01/2.5;
                            x_ref = polyval(casobj.ppx2, casobj.c2); y_ref = polyval(casobj.ppy2, casobj.c2);
                            x_ref = (x_ref <= casobj.XY2(1,end))*x_ref + (x_ref > casobj.XY2(1,end))*casobj.XY2(1,end);
                            y_ref = (y_ref <= casobj.XY2(2,end))*y_ref + (y_ref > casobj.XY2(2,end))*casobj.XY2(2,end);
                        else if pos_flag == 3
                                casobj.c3 = casobj.c3 + 0.01/5 ;
                                x_ref = polyval(casobj.ppx3, casobj.c3); y_ref = polyval(casobj.ppy3, casobj.c3);
                                x_ref = (x_ref <= min(casobj.arc2_x))*min(casobj.arc2_x) + (x_ref > max(casobj.arc2_x))*max(casobj.arc2_x) + ((x_ref > min(casobj.arc2_x))&&(x_ref <= max(casobj.arc2_x)))*x_ref;
                                y_ref = (y_ref <= casobj.arc2_y(end,1))*y_ref + (y_ref > casobj.arc2_y(end,1))*casobj.arc2_y(end,1);
                            else if pos_flag == 4
                                    casobj.c4 = casobj.c4 + 0.01/3 ;
                                    x_ref = polyval(casobj.ppx4, casobj.c4); y_ref = polyval(casobj.ppy4, casobj.c4);
                                    x_ref = (x_ref <= min(casobj.XY5(1,:)))*min(casobj.XY5(1,:)) + (x_ref > max(casobj.XY5(1,:)))*max(casobj.XY5(1,:)) + ((x_ref > min(casobj.XY5(1,:)))&&(x_ref <= max(casobj.XY5(1,:))))*x_ref;
                                    y_ref = (y_ref <= min(casobj.XY5(2,:)))*min(casobj.XY5(2,:)) + (y_ref > max(casobj.XY5(2,:)))*max(casobj.XY5(2,:)) + ((y_ref > min(casobj.XY5(2,:)))&&(y_ref <= max(casobj.XY5(2,:))))*y_ref;
                                else if pos_flag == 5
                                        casobj.c5 = casobj.c5 + 0.01/2 ;
                                        x_ref = polyval(casobj.ppx5, casobj.c5); 
                                        y_ref = polyval(casobj.ppy5, casobj.c5) ;
                                        x_ref = (x_ref <= min(casobj.XY6(1,:)))*min(casobj.XY6(1,:)) + (x_ref > max(casobj.XY6(1,:)))*max(casobj.XY6(1,:)) + ((x_ref > min(casobj.XY6(1,:)))&&(x_ref <= max(casobj.XY6(1,:))))*x_ref;
                                        y_ref = (y_ref <= min(casobj.XY6(2,:)))*min(casobj.XY6(2,:)) + (y_ref > max(casobj.XY6(2,:)))*max(casobj.XY6(2,:)) + ((y_ref > min(casobj.XY6(2,:)))&&(y_ref <= max(casobj.XY6(2,:))))*y_ref;
                                    else if pos_flag == 6
                                            x_ref = casobj.XY6(1,end); y_ref = casobj.XY6(2,end); yaw_ref = 0; xdot_ref = 0; ydot_ref = 0; yawdot_ref = 0;
                                            flag = 1;
                                            break;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                
                %% exit 3 %%
                if 2>3
                    if pos_flag == 1  % the trajectory end is reached
                        casobj.c1 = casobj.c1 + 0.01/0.005;
                        x_ref = polyval(casobj.ppx1, casobj.c1); y_ref = polyval(casobj.ppy1, casobj.c1);
                        x_ref = (x_ref <= casobj.XY1(1,end))*x_ref + (x_ref > casobj.XY1(1,end))*casobj.XY1(1,end);
                        y_ref = (y_ref <= casobj.XY1(2,end))*y_ref + (y_ref > casobj.XY1(2,end))*casobj.XY1(2,end);
                    else if pos_flag == 2
                            casobj.c2 = casobj.c2 + 0.01/1;
                            x_ref = polyval(casobj.ppx2, casobj.c2); y_ref = polyval(casobj.ppy2, casobj.c2);
                            x_ref = (x_ref <= casobj.XY2(1,end))*x_ref + (x_ref > casobj.XY2(1,end))*casobj.XY2(1,end);
                            y_ref = (y_ref <= casobj.XY2(2,end))*y_ref + (y_ref > casobj.XY2(2,end))*casobj.XY2(2,end);
                        else if pos_flag == 3
                                casobj.c3 = casobj.c3 + 0.01/5 ;
                                x_ref = polyval(casobj.ppx3, casobj.c3); y_ref = polyval(casobj.ppy3, casobj.c3);
                                x_ref = (x_ref <= min(casobj.arc3_x))*min(casobj.arc3_x) + (x_ref > max(casobj.arc3_x))*max(casobj.arc3_x) + ((x_ref > min(casobj.arc3_x))&&(x_ref <= max(casobj.arc3_x)))*x_ref;
                                y_ref = (x_ref < 0 && y_ref <= casobj.arc3_y(end,1))*casobj.arc3_y(end,1) + (x_ref < 0 && y_ref > casobj.arc3_y(end,1))*y_ref + (x_ref > 0 && y_ref > min(casobj.arc3_y))*y_ref;
                            else if pos_flag == 4
                                    casobj.c4 = casobj.c4 + 0.01/3 ;
                                    x_ref = polyval(casobj.ppx4, casobj.c4); y_ref = polyval(casobj.ppy4, casobj.c4);
                                    x_ref = (x_ref <= min(casobj.XY7(1,:)))*min(casobj.XY7(1,:)) + (x_ref > max(casobj.XY7(1,:)))*max(casobj.XY7(1,:)) + ((x_ref > min(casobj.XY7(1,:)))&&(x_ref <= max(casobj.XY7(1,:))))*x_ref;
                                    y_ref = (y_ref <= min(casobj.XY7(2,:)))*min(casobj.XY7(2,:)) + (y_ref > max(casobj.XY7(2,:)))*max(casobj.XY7(2,:)) + ((y_ref > min(casobj.XY7(2,:)))&&(y_ref <= max(casobj.XY7(2,:))))*y_ref;
                                else if pos_flag == 5
                                        casobj.c5 = casobj.c5 + 0.01/2 ;
                                        x_ref = polyval(casobj.ppx5, casobj.c5);
                                        y_ref = polyval(casobj.ppy5, casobj.c5) ;
                                        x_ref = (x_ref <= min(casobj.XY8(1,:)))*min(casobj.XY8(1,:)) + (x_ref > max(casobj.XY8(1,:)))*max(casobj.XY8(1,:)) + ((x_ref > min(casobj.XY8(1,:)))&&(x_ref <= max(casobj.XY8(1,:))))*x_ref;
                                        y_ref = ((y_ref <= min(casobj.XY8(2,:)))|| casobj.c5>=15)*min(casobj.XY8(2,:)) + (y_ref > max(casobj.XY8(2,:)))*max(casobj.XY8(2,:)) + (((y_ref > min(casobj.XY8(2,:)))&&(y_ref <= max(casobj.XY8(2,:))))&&(casobj.c5<15))*y_ref;
                                    else if pos_flag == 6
                                            x_ref = casobj.XY8(1,end); y_ref = casobj.XY8(2,end); yaw_ref = 0; xdot_ref = 0; ydot_ref = 0; yawdot_ref = 0;
                                            flag = 1;
                                            break;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                %% exit 4 %%
                if 2>3
                    if pos_flag == 1  % the trajectory end is reached
                        casobj.c1 = casobj.c1 + 0.01/3;
                        x_ref = polyval(casobj.ppx1, casobj.c1); y_ref = polyval(casobj.ppy1, casobj.c1);
                        x_ref = (x_ref <= casobj.XY1(1,end))*x_ref + (x_ref > casobj.XY1(1,end))*casobj.XY1(1,end);
                        y_ref = (y_ref <= casobj.XY1(2,end))*y_ref + (y_ref > casobj.XY1(2,end))*casobj.XY1(2,end);
                    else if pos_flag == 2
                            casobj.c2 = casobj.c2 + 0.01/2.5;
                            x_ref = polyval(casobj.ppx2, casobj.c2); y_ref = polyval(casobj.ppy2, casobj.c2);
                            x_ref = (x_ref <= casobj.XY2(1,end))*x_ref + (x_ref > casobj.XY2(1,end))*casobj.XY2(1,end);
                            y_ref = (y_ref <= casobj.XY2(2,end))*y_ref + (y_ref > casobj.XY2(2,end))*casobj.XY2(2,end);
                        else if pos_flag == 3
                                casobj.c3 = casobj.c3 + 0.01/3.5 ;
                                x_ref = polyval(casobj.ppx3, casobj.c3); y_ref = polyval(casobj.ppy3, casobj.c3);
                                x_ref = (x_ref <= min(casobj.arc4_x))*min(casobj.arc4_x) + (x_ref > max(casobj.arc4_x))*max(casobj.arc4_x) + ((x_ref > min(casobj.arc4_x))&&(x_ref <= max(casobj.arc4_x)))*x_ref;
                                y_ref = (x_ref < 0 && y_ref <= casobj.arc4_y(end,1))*casobj.arc4_y(end,1) + (x_ref < 0 && y_ref > casobj.arc4_y(end,1))*y_ref + (x_ref > 0 && y_ref > min(casobj.arc4_y))*y_ref;
                            else if pos_flag == 4
                                    casobj.c4 = casobj.c4 + 0.01/1.5 ;
                                    x_ref = polyval(casobj.ppx4, casobj.c4); y_ref = polyval(casobj.ppy4, casobj.c4);
                                    x_ref = (x_ref <= min(casobj.XY9(1,:)))*min(casobj.XY9(1,:)) + (x_ref > max(casobj.XY9(1,:)))*max(casobj.XY9(1,:)) + ((x_ref > min(casobj.XY9(1,:)))&&(x_ref <= max(casobj.XY9(1,:))))*x_ref;
                                    y_ref = (y_ref <= min(casobj.XY9(2,:)))*min(casobj.XY9(2,:)) + (y_ref > max(casobj.XY9(2,:)))*max(casobj.XY9(2,:)) + ((y_ref > min(casobj.XY9(2,:)))&&(y_ref <= max(casobj.XY9(2,:))))*y_ref;
                                else if pos_flag == 5
                                        casobj.c5 = casobj.c5 + 0.01/1 ;
                                        x_ref = polyval(casobj.ppx5, casobj.c5);
                                        y_ref = polyval(casobj.ppy5, casobj.c5) ;
                                        x_ref = (x_ref <= min(casobj.XY10(1,:)))*min(casobj.XY10(1,:)) + (x_ref > max(casobj.XY10(1,:)))*max(casobj.XY10(1,:)) + ((x_ref > min(casobj.XY10(1,:)))&&(x_ref <= max(casobj.XY10(1,:))))*x_ref;
                                        y_ref = ((y_ref <= min(casobj.XY10(2,:)))|| casobj.c5>=15)*min(casobj.XY10(2,:)) + (y_ref > max(casobj.XY10(2,:)))*max(casobj.XY10(2,:)) + (((y_ref > min(casobj.XY10(2,:)))&&(y_ref <= max(casobj.XY10(2,:))))&&(casobj.c5<15))*y_ref;
                                    else if pos_flag == 6
                                            x_ref = casobj.XY8(1,end); y_ref = casobj.XY8(2,end); yaw_ref = 0; xdot_ref = 0; ydot_ref = 0; yawdot_ref = 0;
                                            flag = 1;
                                            break;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                p(6*i+1:6*i+6) = [xdot_ref, ydot_ref, yawdot_ref, yaw_ref, x_ref, y_ref];
                
                %x_ref 
                %y_ref 
                %yaw_ref
                %casobj.exit
                vref = 6;
%                 k=j-1;
%                 if k==0
%                     k=size(casobj.path,1);
%                 end
%                 p(6*i+4:6*i+6) = [atan2((casobj.path(j,2)-casobj.path(k,2)),(casobj.path(j,1)-casobj.path(k,1))) , casobj.path(j,1), casobj.path(j,2)];
            end
            pos_flag
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
                break_out=1 ;
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