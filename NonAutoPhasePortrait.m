%------------------------%
% Phase Portrait Example %
%------------------------%

clc
clear all
close all

% Phase Plane Portraits
f = @(t,x) [x(2); -sin(x(1))*t];
% f = @(t,x) [-10*x(1); x(1)*t];
% f = @(t,x) [x(2); -0.2*(x(1)^2-1)*x(2) - x(1)];


% Creates two matrixies one for all the x-values on the grid, and one for
% all the y-values on the grid. Note that x and y are matrixies of the same
% size and shape. This is determined by nx and ny.
nx = 40;
ny = 40;
nt = 20;
tvec = linspace(0,10,nt);
y1 = linspace(-5,5,nx);
y2 = linspace(-5,5,ny);
[x,y,t] = meshgrid(y1,y2,tvec);

% We can use a single loop over each element to comput the derivatives at
% each point (y1,y2)
u = zeros(nx,nt);     % Preallocate the x-quiver direction
v = zeros(ny,nt);     % Preallocate the y-quiver direction
t = zeros(nt); % Derivatives at the starting time t=0

for j = 1:nt
    for i = 1:nx
        Yprime = f(t(j),[x(i) y(i)]);
        u(i,j) = Yprime(1);
        v(i,j) = Yprime(2);
    end
end
% Normalize the quivers
for j = 1:nt
    for i = 1:nx
    Vmod = sqrt(u(i,j)^2 + v(i,j)^2);
    u(i,j) = u(i,j)/Vmod;
    v(i,j) = v(i,j)/Vmod;
    end
end


% Plot the vector field
figure
X = meshgrid(u,v);
[u_mesh,v_mesh] = meshgrid
for j = 1:nt
    hold on
    quiver(x(:,:,j),y(:,:,j),u(:,j),v(:,j),'r')
end
quiver(x,y,u,v,'r')
xlabel('x')
ylabel('dx/dt')
title('Phase Portrait')
axis equal
xlim([min(y1) max(y1)])
ylim([min(y2) max(y2)])

% Let us plot some solutions over the vector field
hold on
t0 = 0; tf = 3; time = [t0,tf];
for x0 = linspace(-4,4,4)
for y0 = linspace(-4,4,4) % initial conditions
    [ts,ys] = ode45(f,time,[x0;y0]);
    for t = time
    plot(ys(:,1),ys(:,2),'k-','LineWidth',2);
    plot(ys(1,1),ys(1,2),'bo'); % Starting Point
    end
end
end
grid on

%% END OF LINE %