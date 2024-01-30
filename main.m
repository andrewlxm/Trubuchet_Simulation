%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  o                                    %
%                                 /|                                    %
%                                / |    <--- counterweight length       %
%                short arm ---> /  |                                    %
%                              /   C  <--- counterweight                %
%                             /                                         %
%                            o                                          %
%                           /|                                          %
%    arm gravity center -> . |                                          %
%                         /  |                                          %
%                        /   |                                          %
%      long arm --->    /    |  <----- pole                             %
%                      /     |                                          %
%                     /      |                                          %
%                    /       |                                          %
%                   /________|__________P <--- projectile               %
%                                                                       %
%                  |<-- sling length -->|                               %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters Definition

h     =   258/100;      % Pole height, from ground to axis of rotation

l_al  =   365/100;      % Long arm length
l_acg =    60/100;      % Arm center of gravity, distance measured from the rotation axis 
l_as  =    93/100;      % Short arm length
l_c   =   114/110;      % Length of counterweight
l_s   =   313/100;      % Sling length

m_a   =  66.00;      % Mass of the arm
m_c   = 173.00;      % Mass of the counterweight
m_p   =   3.85;      % Mass of the projectile

I_a   =   300.0;     % Inertia of arm
I_c   =   100.0;     % Ineria of counterweight

r_p   =   10/100;    % Projectile radius, used for display and animations
h = h - r_p;

release_angle = 190+20+65*0;
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  o                                            %
%                                 /|                                            %
%                                / |  <--- alpha angle                          %
%                               /  | (between counterweight line and vertical)  %
%                              /   C                                            %
%                             /\                                                %
%                            o | <- theta angle                                 %
%                           /|_/                                                %
%                          / |                                                  %
%                         /  |                                                  %
%                        /   |                                                  %
%                       /    |                                                  %
%                      /     |                                                  %
%                     /      |                                                  %
%                    /       |                                                  %
%                   /._______|__________P                                       %
%                    ^                                                          %
%                    beta angle between sling and horizontal                    %
%                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All angles are positive along the trigonometric rotation

stage1 = true;
stage2 = true;
stage3 = true;

plotVectors = false;

% Forces
ScaleFactor = 1;  % Scale Factor to display forces acting upon the counterweight VC the projectile
g = 9.81;

F_norm = g * m_p;
F_g_c  = g * m_c;
F_g_a  = F_g_c * l_as / l_al;
F_g_p  = g * m_p;

%% Motion Stage 1
% Set Initial angles and angle velocities
theta_0      = asin(h/l_al)+pi/2;
alpha_0      = pi - theta_0;
beta_0       = theta_0 - 0.5*pi;

theta_dot_0  = 0;
alpha_dot_0  = 0;
beta_dot_0   = 0;

init_cond_1  = [theta_0, theta_dot_0, alpha_0, alpha_dot_0]';

% Solve system with Runge-Kutta 45 
% and find theta, alpha, theta_dot and alpha_dot
t_interval_1 = [0,10];
ode45options = odeset('AbsTol',1e-10,'RelTol',1e-10,'MaxStep',3e-2,'Refine',1);
[t_stage1,Y] = ode45(@(t,Y) odeFctTrebuchetPart1Final( t, Y, m_c, m_p, m_a, l_c, l_as, l_al, l_acg, I_a, I_c) , t_interval_1 , init_cond_1, ode45options);
theta         = Y(:,1);
theta_dot     = Y(:,2);
alpha         = Y(:,3);
alpha_dot     = Y(:,4);

t        = t_stage1;

y_a      = h + l_al * cos(theta);
beta     = zeros(length(theta),1);
beta_dot = zeros(length(theta),1);
for k=1:length(theta)
    if y_a(k) < l_s
        beta(k) = 2*pi - asin(y_a(k)/l_s);
    else
        beta(k) = beta(k-1);
    end
    if k>1
        beta_dot(k) = (beta(k) - beta(k-1) ) / (t(k) - t(k-1));
    end
end    
delta = beta - pi - theta;



% Forces Stage 1

F_a_s  = F_g_a * cos(delta); % Force acting upon point a along sling line (projection of F_g_a along l_s)
F_v_b  = m_c * l_as * theta_dot.^2;

% Find end of Stage 1
projectileTakeOff = false;
takeOffFactor = 2.83;
for k=1: length(t)
    if F_a_s(k)*sin(beta(k))*cos(beta(k)+3*pi/2) > takeOffFactor * F_g_p && projectileTakeOff == false
        disp('End of stage 1 at')
        disp(t(k))
        k_end_stage1 = k;
        projectileTakeOff = true;
    end
end

t          = t(1:k_end_stage1);
theta      = theta(1:k_end_stage1);
theta_dot  = theta_dot(1:k_end_stage1);
alpha      = alpha(1:k_end_stage1);
alpha_dot  = alpha_dot(1:k_end_stage1);
y_a        = y_a(1:k_end_stage1);
beta       = beta(1:k_end_stage1);
beta_dot   = beta_dot(1:k_end_stage1);
delta      = delta(1:k_end_stage1);

%% Motion Stage 2
if stage2
    % Set Initial angles and angle velocities
    theta_1      = theta(k_end_stage1);
    alpha_1      = alpha(k_end_stage1);
    beta_1       = beta(k_end_stage1);

    theta_dot_1  = theta_dot(k_end_stage1);
    alpha_dot_1  = alpha_dot(k_end_stage1);
    beta_dot_1   = beta_dot(k_end_stage1);

    init_cond_2  = [theta_1, theta_dot_1, alpha_1, alpha_dot_1, beta_1, beta_dot_1]';

    % Solve system with Runge-Kutta 45 
    % and find theta, alpha, beta, theta_dot, alpha_dot and beta_dot
    t_interval_2 = [0,10];
    ode45options = odeset('AbsTol',1e-10,'RelTol',1e-10,'MaxStep',3e-2,'Refine',1);
    [t_stage2,Y] = ode45(@(t,Y) odeFctTrebuchetPart2Final( t,Y, m_c, m_p, m_a, l_c, l_as, l_al, l_s, l_acg, I_a, I_c) , t_interval_2 , init_cond_2, ode45options);
    t            = cat(1, t(1:end-1),         t_stage2+t(k_end_stage1));
    theta        = cat(1, theta(1:end-1),     Y(:,1));
    theta_dot    = cat(1, theta_dot(1:end-1), Y(:,2));
    alpha        = cat(1, alpha(1:end-1),     Y(:,3));
    alpha_dot    = cat(1, alpha_dot(1:end-1), Y(:,4));
    beta         = cat(1, beta(1:end-1),      Y(:,5));
    beta_dot     = cat(1, beta_dot(1:end-1),  Y(:,6));

    % Forces during Stage 2
    delta  = beta - pi - theta;
    F_a_s  = F_g_a * cos(delta); % Force acting upon point a along sling line (projection of F_g_a along l_s)
    F_v_b  = m_c * l_as * theta_dot.^2;

    % Find end of Stage 2
    projectileReleased = false;
    for k=1: length(t)
        k_end_stage2 = k;
        if  (beta(k)+pi/2-theta(k))*180/pi < release_angle && projectileReleased == false
            projectileReleased = true;
            disp('End of stage 2 at')
            disp(t(k_end_stage2))
            break
        end
    end
    k_end_stage2 = k_end_stage2;
    
    t          = t(1:k_end_stage2);
    theta      = theta(1:k_end_stage2);
    theta_dot  = theta_dot(1:k_end_stage2);
    alpha      = alpha(1:k_end_stage2);
    alpha_dot  = alpha_dot(1:k_end_stage2);
    beta       = beta(1:k_end_stage2);
    beta_dot   = beta_dot(1:k_end_stage2);
    delta      = delta(1:k_end_stage2);
end

x_a =   - l_al * sin(theta);
y_a = h + l_al * cos(theta);

x_p = x_a + l_s * cos(beta);
y_p = y_a + l_s * sin(beta);

v_p_x = (x_p(end) - x_p(end-1)) /  ( t(end) - t(end-1) ) ;
v_p_y = (y_p(end) - y_p(end-1)) /  ( t(end) - t(end-1) ) ;

gamma_0   = atan(v_p_y/v_p_x);
v_p_0     = sqrt( (v_p_x).^2 + (v_p_y).^2 );
if gamma_0 < 0
    gamma_0   = atan(v_p_y/v_p_x);
    v_p_0     = - sqrt( (v_p_x).^2 + (v_p_y).^2 );
end
x_release = x_p(end);
y_release = y_p(end)+r_p;

%% Motion Stage 3
if stage3
    % Set Initial angles and angle velocities
    theta_2      = theta(k_end_stage2);
    alpha_2      = alpha(k_end_stage2);
    beta_2       = beta(k_end_stage2);

    theta_dot_2  = theta_dot(k_end_stage2);
    alpha_dot_2  = alpha_dot(k_end_stage2);
    beta_dot_2   = beta_dot(k_end_stage2);

    init_cond_3  = [theta_2, theta_dot_2, alpha_2, alpha_dot_2, beta_2, beta_dot_2]';

    % Solve system with Runge-Kutta 45 
    % and find theta, alpha, beta, theta_dot, alpha_dot and beta_dot
    t_interval_3 = [0,4.5];
    ode45options = odeset('AbsTol',1e-10,'RelTol',1e-10,'MaxStep',3e-2,'Refine',1);
    [t_stage3,Y] = ode45(@(t,Y) odeFctTrebuchetPart3Final( t,Y, m_c, m_p*0.001, m_a, l_c, l_as, l_al, l_s, l_acg, I_a, I_c) , t_interval_3 , init_cond_3, ode45options);
    t            = cat(1, t(1:end-1),         t_stage3+t(k_end_stage2));
    theta        = cat(1, theta(1:end-1),     Y(:,1));
    theta_dot    = cat(1, theta_dot(1:end-1), Y(:,2));
    alpha        = cat(1, alpha(1:end-1),     Y(:,3));
    alpha_dot    = cat(1, alpha_dot(1:end-1), Y(:,4));
    beta         = cat(1, beta(1:end-1),      Y(:,5));
    beta_dot     = cat(1, beta_dot(1:end-1),  Y(:,6));
    
    theta_dot_dot = finiteDerivative(t,theta_dot);
    alpha_dot_dot = finiteDerivative(t,alpha_dot);
    beta_dot_dot  = finiteDerivative(t,beta_dot);

    % Forces during Stage 3
    delta  = beta - pi - theta;
    F_a_s  = F_g_a * cos(delta); % Force acting upon point a along sling line (projection of F_g_a along l_s)
    F_v_b  = m_c * l_as * theta_dot.^2;

    k_end_stage3 = length(t_stage3);
end



x_a =   - l_al * sin(theta);
y_a = h + l_al * cos(theta)+r_p;

x_b =     l_as * sin(theta);
y_b = h - l_as * cos(theta)+r_p;

x_c = x_b - l_c * sin(alpha);
y_c = y_b - l_c * cos(alpha);

x_p = x_a + l_s * cos(beta);
y_p = y_a + l_s * sin(beta);

x_p_dot_dot = zeros(k_end_stage2,1);
y_p_dot_dot = zeros(k_end_stage2,1);
a_p         = zeros(k_end_stage2,1);
for k=1:k_end_stage2
    x_p_dot_dot(k) =                - l_al * ( theta_dot_dot(k)*cos(theta(k)) - theta_dot(k).^2*sin(theta(k)) ) ;
    x_p_dot_dot(k) = x_p_dot_dot(k) -  l_s * ( beta_dot_dot(k)*sin(beta(k))   + beta_dot(k).^2*cos(beta(k)) ) ;

    y_p_dot_dot(k) =                - l_al * ( theta_dot_dot(k)*sin(theta(k)) + theta_dot(k).^2*cos(theta(k)) ) ;
    y_p_dot_dot(k) = x_p_dot_dot(k) +  l_s * ( beta_dot_dot(k)*cos(beta(k))   - beta_dot(k).^2*sin(beta(k)) ) ;
    
    a_p(k) = sqrt (x_p_dot_dot(k).^2 + y_p_dot_dot(k).^2)/g;
end

x_c_dot_dot = zeros(k_end_stage3,1);
y_c_dot_dot = zeros(k_end_stage3,1);
a_c         = zeros(k_end_stage3,1);
for k=1:k_end_stage3
    x_c_dot_dot(k) =                  l_as * ( theta_dot_dot(k)*cos(theta(k)) - theta_dot(k).^2*sin(theta(k)) ) ;
    x_c_dot_dot(k) = x_c_dot_dot(k) -  l_c * ( alpha_dot_dot(k)*cos(alpha(k)) - alpha_dot(k).^2*sin(alpha(k)) ) ;

    y_c_dot_dot(k) =                  l_as * ( theta_dot_dot(k)*sin(theta(k)) + theta_dot(k).^2*cos(theta(k)) ) ;
    y_c_dot_dot(k) = x_c_dot_dot(k) +  l_c * ( alpha_dot_dot(k)*sin(alpha(k)) + alpha_dot(k).^2*cos(alpha(k)) ) ;
    
    a_c(k) = sqrt (x_c_dot_dot(k).^2 + y_c_dot_dot(k).^2)/g;
end

%% Motion Stage 4
c_p = 0.1;
[t_p_f, x_p_f, y_p_f] = projectileAirDrag(m_p, c_p, gamma_0, v_p_0, x_release, y_release, -r_p);

% t_p        = cat(1, t(1:k_end_stage2), t_p_f);
% x_p        = cat(1, x_p(1:end-1), x_p_f);
% y_p        = cat(1, y_p(1:end-1), y_p_f);

%% Animation and Plots

f = figure;


NSubPlot = 4;
subPlotId(1,:) = [ 1, 4];
subPlotId(2,:) = [ 2,  6];
subPlotId(3,:) = [ 7, 10];
subPlotId(4,:) = [ 8, 12];

for n=1:NSubPlot
    if n~=2
        sph = subplot(4,3,subPlotId(n,:));
        daspect([1 1 1])
    end
    xlim = 1.05*(l_al+l_s);
    if n~=4
        axis([ -xlim , xlim*0.7, -1.0, h+l_al+l_s+1])
    else
        axis([ min(-xlim, min(x_p_f)-1) , max(xlim, max(x_p_f)+1), -1.0,  max(y_p_f)+1])
    end
    hold on

    if n~=2 && n~=4
        grid on
        % Fixed Lines
        poleLine           = line('XData',[ 0      ,      0], 'YData',[0      , h+0    ], 'linewidth',1.2, 'color',[.7 .7 .7]);
        groundLine         = line('XData',[ -xlim  ,   xlim], 'YData',[0      ,   0    ], 'linewidth',1.2, 'color',[.7 .7 .7]);

        % Moving Lines and Points
        armLine(n)            = line('XData',[ x_a(1) , x_b(1)], 'YData',[y_a(1) , y_b(1)], 'linewidth',1*max((m_a)^.25,0.1));
        counterweightLine(n)  = line('XData',[ x_b(1) , x_c(1)], 'YData',[y_b(1) , y_c(1)], 'linewidth',2);
        slingLine(n)          = line('XData',[ x_a(1) , x_p(1)], 'YData',[y_a(1) , y_p(1)], 'linewidth',0.1);
        counterWeightPoint(n) = plot( x_c(1), y_c(1), '.', 'MarkerSize', 20*max((abs(m_c))^(1/4),0.1));
    end
    if n==4
        % Fixed Lines
        poleLine           = line('XData',[ 0      ,            0], 'YData',[0      , h+r_p     ], 'linewidth',1.2, 'color',[.7 .7 .7]);
        groundLine         = line('XData',[ min(-xlim, min(x_p_f)-1) , max(xlim, max(x_p_f)+1),], 'YData',[0      , 0     ], 'linewidth',1.2, 'color',[.7 .7 .7]);

        % Moving Lines and Points
        armLine(n)            = line('XData',[ x_a(1) , x_b(1)], 'YData',[y_a(1) , y_b(1)], 'linewidth',1*max((m_a)^.25,0.1));
        counterweightLine(n)  = line('XData',[ x_b(1) , x_c(1)], 'YData',[y_b(1) , y_c(1)], 'linewidth',2);
        slingLine(n)          = line('XData',[ x_a(1) , x_p(1)], 'YData',[y_a(1) , y_p(1)], 'linewidth',0.1);
%         projectilePoint(n)    = plot( x_p(1), y_p(1), '.', 'MarkerSize', 20*max((abs(m_p))^(1/4),0.3));
    end
    projectilePoint(n)    = plot( x_p(1), y_p(1), '.', 'MarkerSize', 20*max((abs(m_p))^(1/4),0.3));
end

% Forces Vectors
if plotVectors
    vectXb = quiver(x_b(1), y_b(1), 0, -F_g_c/F_norm/ScaleFactor, 'linewidth',2*sqrt(ScaleFactor), 'color','b');
    vectXb.MaxHeadSize = 2;

    vectXvb = quiver(x_b(1), y_b(1), F_v_b(1)/F_norm * sin(theta(1)), -F_v_b(1)/F_norm * cos(theta(1)), 'linewidth',2, 'color','r');
    vectXvb.MaxHeadSize = 2;

    vectXa = quiver(x_a(1), y_a(1), F_g_a/F_norm*cos(theta(1))*sin(theta(1)), F_g_a/F_norm*sin(theta(1))*sin(theta(1)), 'linewidth',2, 'color','r');
    vectXa.MaxHeadSize = 2;

    vectXc = quiver(x_c(1), y_c(1), 0, -F_g_c/F_norm/ScaleFactor, 'linewidth',2*sqrt(ScaleFactor), 'color','k');
    vectXc.MaxHeadSize = 2;

    vectXa_s = quiver(x_p(1), y_p(1), -F_a_s(1)/F_norm*cos(beta(1))*sin(theta(1)), -F_a_s(1)/F_norm*sin(beta(1))*sin(theta(1)), 'linewidth',2, 'color','c');
    vectXa_s.MaxHeadSize = 2;

    vectXa_p = quiver(x_p(1), y_p(1), 0, F_a_s(1)/F_norm*sin(beta(1))*cos(beta(1)+3*pi/2), 'linewidth',2, 'color','g');
    vectXa_p.MaxHeadSize = 2;

    vectXp = quiver(x_p(1), y_p(1), 0, -F_g_p/F_norm, 'linewidth',2);
    vectXp.MaxHeadSize = 2;
end




step    = 5;
start   = step*2;
t_pause = 1e-8;
stopMotion1 = false;

subplot(4,3,subPlotId(1,:));
[x_p_down, y_p_down, x_p_up, y_p_up] = safeLine(x_p(1:100), y_p(1:100), r_p*3);
for k=1:length(y_p_down)
    y_p_down(k) = max(0, y_p_down(k));
end
dangerDown = plot(x_p_down, y_p_down, 'Color', [.7 .7 .7]);
dangerUp   = plot(x_p_up,   y_p_up,   'Color', [.7 .7 .7]);
uistack(dangerUp,'bottom')


coord_up      = [x_p_down, y_p_down];
coord_low     = [x_p_up,y_p_up];
coord_combine = [coord_up;flipud(coord_low)];
dangerArea = fill(coord_combine(:,1),coord_combine(:,2),[.7 .7 .7]);
uistack(dangerArea,'bottom')


coord_circle_weight_x = linspace(-2.44, 2.44, 1001)';
coord_circle_weight_y_down = zeros(length(coord_circle_weight_x),1);
coord_circle_weight_y_up = zeros(length(coord_circle_weight_x),1);
for k=1 : length(coord_circle_weight_x)
    coord_circle_weight_y_down(k) = h - sqrt(2.44.^2 - coord_circle_weight_x(k).^2);
    coord_circle_weight_y_up(k)   = h + sqrt(2.44.^2 - coord_circle_weight_x(k).^2);
end
dangerCircle_down = plot(coord_circle_weight_x, coord_circle_weight_y_down, 'Color', [.7 .7 .7]); 
dangerCircle_up   = plot(coord_circle_weight_x, coord_circle_weight_y_up, 'Color', [.7 .7 .7]); 

coord_circle_up      = [coord_circle_weight_x, coord_circle_weight_y_down];
coord_circle_low     = [coord_circle_weight_x,coord_circle_weight_y_up];
coord_circle_combine = [coord_circle_up;flipud(coord_circle_low)];
dangerCircleArea = fill(coord_circle_combine(:,1),coord_circle_combine(:,2),[.7 .7 .7]);
uistack(dangerCircleArea,'bottom')

pause(0);
for k=start:step: length(t)

    tic
    
    for n=1:NSubPlot
        if n~=2
            subplot(4,3,subPlotId(n,:));
        end
        if t(k) < t(k_end_stage2)
            plot(x_p(k), y_p(k), 'r.-');
            set(projectilePoint(n), 'XData',[ x_p(k)],          'YData',[y_p(k) ] );
        else
            
            stopMotion1 = true;
            if n==4 && k-k_end_stage2+1 < length(x_p_f)
                plot(x_p_f(k-k_end_stage2+1), y_p_f(k-k_end_stage2+1), 'r.-');
                set(projectilePoint(n), 'XData',[ x_p_f(k-k_end_stage2+1)], 'YData',[y_p_f(k-k_end_stage2+1) ] );
            end
        end
        
        if (n~=1 && n~=2 && n~=4) || ( n==1 && ~stopMotion1)
            if x_a(k) < x_p(k) && ~stopMotion1
                subplot(4,3,subPlotId(1,:));
                armArea = fill([ -4 -4 x_a(1) 0 x_a(1) ], [h*1.5 0 0 h*1.5 h*1.5 ], [.7 .7 .7]);
                armArea.EdgeColor = [.7 .7 .7];
                uistack(armArea, 'bottom')
                slingArea = fill([x_a(k) x_a(k) x_p(k) x_p(k)], [0 y_a(k) y_p(k) 0], [.7 .7 .7]);
                slingArea.EdgeColor = [.7 .7 .7];
                uistack(slingArea, 'bottom')
            end
            % Moving Lines and Points
            set(armLine(n),            'XData',[ x_a(k) , x_b(k)],  'YData',[y_a(k) , y_b(k)] );
            set(counterweightLine(n),  'XData',[ x_b(k) , x_c(k)],  'YData',[y_b(k) , y_c(k)] );
            set(slingLine(n),          'XData',[ x_a(k) , x_p(k)],  'YData',[y_a(k) , y_p(k)] );
            set(counterWeightPoint(n), 'XData',[ x_c(k)],           'YData',[y_c(k) ] );

            % Forces Vectors
            if plotVectors
                vectXb.XData = x_b(k); vectXb.YData = y_b(k); 
                vectXb.UData = -F_g_c/F_norm*sin(alpha(k))*cos(alpha(k))/ScaleFactor; 
                vectXb.VData = -F_g_c/F_norm*cos(alpha(k))*cos(alpha(k))/ScaleFactor;

                vectXvb.XData = x_b(k); vectXvb.YData = y_b(k); 
                vectXvb.UData =  F_v_b(k)/F_norm * sin(theta(k)); 
                vectXvb.VData = -F_v_b(k)/F_norm * cos(theta(k));

                vectXa.XData = x_a(k); vectXa.YData = y_a(k); 
                vectXa.UData = F_g_a/F_norm*cos(theta(k))*sin(theta(k)); 
                vectXa.VData = F_g_a/F_norm*sin(theta(k))*sin(theta(k));

                vectXc.XData = x_c(k); vectXc.YData = y_c(k); vectXc.UData = 0; vectXc.VData = -F_g_c/F_norm/ScaleFactor;

                vectXa_s.XData = x_p(k); vectXa_s.YData = y_p(k); 
                vectXa_s.UData = -F_a_s(k)/F_norm*cos(beta(k))*sin(theta(k));
                vectXa_s.VData = -F_a_s(k)/F_norm*sin(beta(k))*sin(theta(k));

                vectXa_p.XData = x_p(k); vectXa_p.YData = y_p(k); 
                vectXa_p.UData = 0;
                vectXa_p.VData = F_a_s(k)/F_norm*sin(beta(k))*cos(beta(k)+3*pi/2);

                vectXp.XData = x_p(k); vectXp.YData = y_p(k); vectXp.UData = 0; vectXp.VData = -F_g_p/F_norm;
            end
        elseif n==4
            % Moving Lines and Points
            set(armLine(n),            'XData',[ x_a(k) , x_b(k)],  'YData',[y_a(k) , y_b(k)] );
            set(counterweightLine(n),  'XData',[ x_b(k) , x_c(k)],  'YData',[y_b(k) , y_c(k)] );
            set(slingLine(n),          'XData',[ x_a(k) , x_p(k)],  'YData',[y_a(k) , y_p(k)] );  
        end
    end
    drawnow
    tStop = toc;
    pause(max( (t(k) - t(k-step) - tStop)*0.2  , t_pause))
end