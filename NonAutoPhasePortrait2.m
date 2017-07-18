%------------------------%
% Phase Portrait Example %
%------------------------%

clc
clear all
close all

% Phase Plane Portraits
%f = @(t,x) [x(2); -sin(x(1))*t];
% f = @(t,x) [-10*x(1); x(1)*t];
% f = @(t,x) [x(2); -0.2*(x(1)^2-1)*x(2) - x(1)];
f = @(t,x) [-x(1) - exp(-0.2*t)*x(2) ; -x(1) - x(2)];
%f = @(t,x) [cos(2*t)-cos(5/9*t) ; sin(2*t)-sin(5/9*t)];

% Creates two matrixies one for all the x-values on the grid, and one for
% all the y-values on the grid. Note that x and y are matrixies of the same
% size and shape. This is determined by nx and ny.
nt = 100;
t0 = 0.1;
tf = 25;
dt = tf/nt;
x0 = 5; xt = x0;
y0 = 2; yt = y0;
for t = linspace(t0,tf,nt);
    nx = 30;
    ny = 30;
    y1 = linspace(-5,5,nx);
    y2 = linspace(-5,5,ny);
    [x,y] = meshgrid(y1,y2);

    % We can use a single loop over each element to comput the derivatives at
    % each point (y1,y2)
    u = zeros(size(x));     % Preallocate the x-quiver direction
    v = zeros(size(y));     % Preallocate the y-quiver direction
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
    figure(100)
    clf
    quiver(x,y,u,v,'r')
    xlabel('x')
    ylabel('dx/dt')
    title('Phase Portrait')
    axis equal
    xlim([min(y1) max(y1)])
    ylim([min(y2) max(y2)])


% Let us plot some solutions over the vector field
hold on
    [ts,ys] = ode45(f,[t,t+dt],[xt;yt]);
    [ts2,ys2] = ode45(f,[t0,t+dt],[x0;y0]);
    plot(ys(:,1),ys(:,2),'k-','LineWidth',2);
    plot(ys(1,1),ys(1,2),'bo'); % Starting Point
    plot(ys2(:,1),ys2(:,2),'k-','LineWidth',2);
    xt = ys(end,1);
    yt = ys(end,2);
drawnow
end

disp('Done')
%% END OF LINE %