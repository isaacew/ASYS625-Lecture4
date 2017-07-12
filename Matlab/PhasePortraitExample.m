%------------------------%
% Phase Portrait Example %
%------------------------%

clc
clear all
close all

% Phase Plane Portraits
f = @(t,x) [x(2); -sin(x(1))];

% Creates two matrixies one for all the x-values on the grid, and one for
% all the y-values on the grid. Note that x and y are matrixies of the same
% size and shape. This is determined by nx and ny.
nx = 50;
ny = 20;
y1 = linspace(-2,8,nx);
y2 = linspace(-2,2,ny);
[x,y] = meshgrid(y1,y2);

% We can use a single loop over each element to comput the derivatives at
% each point (y1,y2)
u = zeros(size(x));     % Preallocate the x-quiver direction
v = zeros(size(y));     % Preallocate the y-quiver direction
t = 0; % Derivatives at the starting time t=0
for i = 1:numel(x)
    Yprime = f(t,[x(i) y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end
% Normalize the quivers
for i = 1:numel(x)
Vmod = sqrt(u(i)^2 + v(i)^2);
u(i) = u(i)/Vmod;
v(i) = v(i)/Vmod;
end
% Plot the vector field
figure
quiver(x,y,u,v,'r')
xlabel('x')
ylabel('dx/dt')
title('Phase Portrait')
axis equal
xlim([min(y1) max(y1)])
ylim([min(y2) max(y2)])

% Let us plot some solutions over the vector field
hold on
t0 = 0; tf = 10; time = [t0,tf];
for x0 = linspace(0,2,10) % initial conditions
    [ts,ys] = ode45(f,time,[0;x0]);
    plot(ys(:,1),ys(:,2),'k-','LineWidth',2);
    plot(ys(1,1),ys(1,2),'bo'); % Starting Point
end
grid on

%% END OF LINE %