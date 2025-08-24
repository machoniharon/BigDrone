# BigDrone
%% Define the artificial potential field parameters
% Points of attraction
xa = 70;
ya = 50;
% Location of obstacles
nobs = 5;
xo = [20 30 35 50 65];
yo = [13 25 35 38 45];
% Size of obstacles (rho_0)
rho_0 = [7 3 6 4 4];
% Constants.
Ka = 1; 
Kv = 50; 
Kr = [500 250 400 300 300];
mass = 5;
%% Integration step.
% initial conditions
x0 = 5; y0 = 5;
% Initialise simulation parameters
x = x0; y = y0;
vx = 0; vy = 0;
terminate = false;
t = 0;
dt = 0.01;
% Store variables for later assessment
traj = [];
Fac = []; 
Frc = []; 
Fdc = [];
Fc = [];
while(~terminate)
    %% Check if reached destination
    rho = sqrt((xa-x)^2+(ya-y)^2);
    if(rho < 0.1)
        terminate = true;
    end
    %
    %% Calculate the conservative forces
    % Attractive force (Fa)
    Fa = zeros(2,1);

    Fa(1,1) = Fa(1,1) + Ka *(xa - x);
    Fa(2,1) = Fa(2,1) + Ka *(ya - y);
        
    Fac = [Fac Fa];
    % Repulsive force (Fr)
    Fr = zeros(2,1);
    for i = 1:nobs
       rho_obs = sqrt((x-xo(i))^2+(y-yo(i))^2);
        if rho_obs < rho_0(i)
           Fr(1,1) = Fr(1,1) + Kr(i) * (1/rho_obs - 1/rho_0(i)) * rho_obs^1.01 * (x-xo(i)); % Equation for repulsive force
           Fr(2,1) = Fr(2,1) + Kr(i) * (1/rho_obs - 1/rho_0(i)) * rho_obs^1.01 * (y-yo(i)); % Equation for repulsive force
        end
    end
    Frc = [Frc Fr];
    % Damping force (Fd)
    Fd = zeros(2,1);
    Fd(1,1) = -Kv * vx;
    Fd(2,1) = -Kv * vy;

    Fdc = [Fdc Fd];
    % Total force (F)
    
    F(1,1) = Fa(1,1) + Fr(1,1) + Fd(1,1);
    F(2,1) = Fa(2,1) + Fr(2,1) + Fd(2,1);
   
    Fc = [Fc F];
    % Velocity and Position (vx, vy, x, y) through integrating the equations

    vx = vx + (F(1,1)/mass) * dt; % Update velocity in x-direction
    vy = vy + (F(2,1)/mass) * dt; % Update velocity in y-direction
    x = x + vx * dt; % Update position in x-direction
    y = y + vy * dt; % Update position in y-direction

    traj = [traj [x;y]]; %plot the trajectory
end
%% Display the results
% Plot the trajectory
figure(1);
% Now the attractive potential contour
xmin = 0; xmax = 80;
ymin = 0; ymax = 60;
[X,Y] = meshgrid(xmin:5:xmax,ymin:5:ymax);
Z = 1.0*(0.5*Ka*sqrt((X-xa).^2+(Y-ya).^2));
Zmax = max(max(Z));
Z = 10*Z/Zmax;
contour(X,Y,Z,50);hold on; 
plot(traj(1,1),traj(2,1),'bo','LineWidth',6);
plot(xa,ya,'bx','LineWidth',10);
xlabel('x(m)');
ylabel('y(m)');
axis('equal');
% The obstacle contours
for kk = 1:nobs
    r = 0:0.01:rho_0(kk);
    tht = 0:0.1:2*pi+0.1;
    [THT,R] = meshgrid(tht,r);
    Z = 0.5*Kr(kk)*((sqrt((R.*cos(THT)).^2+(R.*sin(THT)).^2))-(rho_0(kk))).^2;
    Zmax = max(max(Z));
    Z = 10*Z/Zmax;
    contour(xo(kk)+R.*cos(THT),yo(kk)+R.*sin(THT),Z,50);
end
% Draw the obstacle boundaries
for jj = 1:nobs
    tht = 0:0.01:2*pi;
    xobs = xo(jj)+rho_0(jj)*cos(tht);
    yobs = yo(jj)+rho_0(jj)*sin(tht);
    plot(xobs,yobs,':r','LineWidth',2.0);
end
%% Finally, plot the trajectory
plot(traj(1,:),traj(2,:),'g','LineWidth',2.0);
hold off;
axis equal;
axis tight;
