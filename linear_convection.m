% linear convection equation, u_t + au_x = 0
clc; clear; close all

nx = 3000;
dx = 2 / (nx - 1);
CFL = 0.7;
a = 1; % wavespeed
dt = CFL*dx/a;
nt = round(0.1/dt);

%% upwind scheme

u = initial(nx,dx);
for n = 1:nt
    un = u;
    for j = 2:nx
        u(j) = un(j) - a*dt/dx*(un(j) - un(j-1)); 
    end
end

hold on
plot(linspace(0,2,nx),u,':','LineWidth',3);
xlabel('x'); ylabel('u'); ylim([0 3])
hold off
%% Lax-Friedrichs scheme

u = initial(nx,dx);
for n = 1:nt
    un = u;
    for j = 2:nx-1
        u(j) = 1/2*(un(j+1)+un(j-1)) - a*dt/dx/2*(un(j+1) - un(j-1)); 
    end
end

hold on
plot(linspace(0,2,nx),u,':','LineWidth',3);
xlabel('x'); ylabel('u'); ylim([0 3])
hold off

%% MacCormack

u = initial(nx,dx);
u_star = ones(1,nx);
for n = 1:nt
    un = u;
    for j = 1:nx-1
        u_star(j) = un(j) - dt/dx*a*(un(j+1)-un(j));
    end
    for j = 2:nx-1
        u(j) = (un(j)+u_star(j))/2 - dt/dx/2*a*(u_star(j)-u_star(j-1));
    end
end

hold on
plot(linspace(0,2,nx),u,':','LineWidth',3);
xlabel('x'); ylabel('u'); ylim([0 3])
hold off


function u = initial(nx,dx)
    % initial data
    u = ones(1, nx);
    u(round(0.5/dx):round(1/dx+1)) = 2;
    figure
    plot(linspace(0,2,nx),u,'LineWidth',3);
    xlabel('x'); ylabel('u'); ylim([0 3])
end