%EULER 1D Euler equation
close all; clc; clear;
% U: conservative variables: rho, rho u, E
% result: primitive variables: rho, u, p

%% LF
% nx = 5001;
% CFL = 0.9;
% LF('A',nx,CFL)
% nx = 40001;
% LF('B',nx,CFL)

%% MacCormack
% nx = 5001;
% CFL = 0.95;
% MacCormack('A',nx,CFL)
% nx = 40001;
% MacCormack('B',nx,CFL)

%% ROE
% nx = 5001;
% CFL = 0.4;
% Roe('A',nx,CFL)
nx = 40001;
CFL = 0.25;
Roe('B',nx,CFL)

function [U,dx,endtime,interval] = init_A(scheme,nx)
%INIT_A initialize problem A
          
    interval = [0, 1];
    discontinuity = 0.3;
    UL = [1 0 2.5];
    UR = [0.125 0 0.25];
    endtime = 0.2;
    
    dx = (interval(2)-interval(1)) / (nx - 1);  % calculation interval is [0,1]
    U = ones(3,nx);     % 3-by-nx: rho, rho*u, E
    
    for i = 1:3
        U(i,1:round((discontinuity-interval(1))/dx)) = UL(i); 
        U(i,round((discontinuity-interval(1))/dx)+1:end) = UR(i);       
    end
    plot_line(U,interval);
    title('Initial U=(\rho,\rho u,E)^T')
    xlabel('x')
    legend('\rho','\rho u','E')
    set(gca,'FontSize',16,'FontName','Times New Roman'); axis tight
    
    saveas(gcf,['tex/img/init_',scheme,'_A'],'epsc')
end

function [U,dx,endtime,interval] = init_B(scheme,nx)
%INIT_B initialize problem B

    interval = [-5,5];
    discontinuity = -4;
    UL = [3.857143 2.629369 10.33333];
    UR = [0 0 1];
    endtime = 1.8;
    
    dx = (interval(2)-interval(1)) / (nx - 1);  % calculation interval is [0,1]
    result = ones(3,nx);
    
    for i = 1:3
        result(i,1:round((discontinuity-interval(1))/dx)) = UL(i);
        if i == 1
            len = length(result(i,round((discontinuity-interval(1))/dx)+1:end));
            result(i,round((discontinuity-interval(1))/dx)+1:end) = 1 + 0.2*sin(5*linspace(discontinuity,interval(2),len));
        else
            result(i,round((discontinuity-interval(1))/dx)+1:end) = UR(i);
        end        
    end
    plot_line(result,interval);
    title('(\rho,u,p)(x,0)')
    xlabel('x')
    legend('\rho','u','p')
    set(gca,'FontSize',16,'FontName','Times New Roman'); axis tight  
    saveas(gcf,['tex/img/init_',scheme,'_B'],'epsc')
    
    U = primitive2conservative(result);
end

function F = flux(U)
%FLUX flux F

    gamma = 1.4;
    F = ones(size(U));
    F(1,:) = U(2,:);
    F(2,:) = 1/2*(3-gamma)*U(2,:).^2./U(1,:) + (gamma-1)*U(3,:);
    F(3,:) = gamma*U(2,:)./U(1,:).*U(3,:) - 1/2*(gamma-1)*U(2,:).^3./U(1,:).^2;
end

function result = conservative2primitive(U)
%CONSERVATIVE2PRIMITIVE
    
    gamma = 1.4;
    result = ones(size(U));
    result(1,:) = U(1,:);
    result(2,:) = U(2,:)./U(1,:);
    result(3,:) = (U(3,:)-1/2*result(1,:).*result(2,:).^2)*(gamma-1);
end

function U = primitive2conservative(result)
%PRIMITIVE2CONSERVATIVE

    gamma = 1.4;
    U = ones(size(result));
    U(1,:) = result(1,:);
    U(2,:) = result(1,:).*result(2,:);
    U(3,:) = 1/2*result(1,:).*result(2,:).^2 + result(3,:)/(gamma-1);
end

function result = solution(U,interval,scheme,problem,time)
%PRIMITIVE non-conservative results: rho, u, p

    nx = size(U,2);
    result = conservative2primitive(U);        
    plot_line(result,interval);
    title(['t = ',num2str(time),' result=(\rho,u,p)'])
    xlabel('x')
    legend('\rho','u','p')
    set(gca,'FontSize',16,'FontName','Times New Roman'); axis tight
    
    saveas(gcf,['tex/img/',scheme,'_',problem],'epsc')
end

function plot_line(result,interval)
%PLOT_INTERVAL lines
    
    figure
    hold on
    nx = size(result,2);
    x = linspace(interval(1),interval(2),nx);
    plot(x,result(1,:),'-','LineWidth',3)
    plot(x,result(2,:),'-.','LineWidth',3)
    plot(x,result(3,:),':','LineWidth',3)
    hold off
end

function dt = time_step(U,dx,CFL)
%TIME_STEP calculate time step, dt = CFL*dx/Smax, Smax = max{|u|+a}

    gamma = 1.4;
    result = conservative2primitive(U);
    sound = sqrt(gamma*abs(result(3,:)./result(1,:)));
    Smax = max(abs(result(2,:))+sound);
    dt = CFL*dx/Smax;
end

function LF(problem,nx,CFL)
%LF 
    if problem == 'A'
        [U,dx,endtime,interval] = init_A('LF',nx);
    elseif problem == 'B'
        [U,dx,endtime,interval] = init_B('LF',nx);
    end
    
    t = 0;
    n = 0;
    while t < endtime + 1e-8 % float correction
        Un = U;
        Fn = flux(Un);
        dt = time_step(Un,dx,CFL);
        for j = 2:nx-1
            U(:,j) = (Un(:,j+1)+Un(:,j-1))/2 - dt/dx/2*(Fn(:,j+1)-Fn(:,j-1));
        end
        t = t + dt;
        n = n + 1;
    end
    fprintf('Lax-Friedrichs:\n')
    fprintf('dx = %f\tCFL=%f\ttimestep num=%d\n',dx,CFL,n)
    solution(U,interval,'LF',problem,t);
end

function MacCormack(problem,nx,CFL)
%MACCORMACK
    
    if problem == 'A'
        [U,dx,endtime,interval] = init_A('LF',nx);
    elseif problem == 'B'
        [U,dx,endtime,interval] = init_B('LF',nx);
    end
    
    t = 0;
    n = 0;
    while t < endtime + 1e-8 % float correction
        Un = U;
        Fn = flux(Un);
        dt = time_step(Un,dx,CFL);
        U_star = ones(3,nx-1);
        % predictor step, forward difference, 1:nx-1
        for j = 1:nx-1
            U_star(:,j) = Un(:,j) - dt/dx*(Fn(:,j+1)-Fn(:,j));
        end
        Fn_star = flux(U_star);

        % corrector step, backward difference, 2:nx
        for j = 2:nx-1
            U(:,j) = (Un(:,j)+U_star(:,j))/2 - dt/dx/2*(Fn_star(:,j)-Fn_star(:,j-1));
        end
        t = t + dt;
        n = n + 1;
    end
    fprintf('MacCormack:\n')
    fprintf('dx = %f\tCFL=%f\ttimestep num=%d\n',dx,CFL,n)
    solution(U,interval,'Mac',problem,t);
end

function Roe(problem,nx,CFL)
%Roe
    
    if problem == 'A'
        [U,dx,endtime,interval] = init_A('Roe',nx);
    elseif problem == 'B'
        [U,dx,endtime,interval] = init_B('Roe',nx);
    end
    
    t = 0;
    n = 0;
    while t < endtime + 1e-8 % float correction
        Un = U;        
        dt = time_step(Un,dx,CFL);
        F = Roe_flux(Un);
        

        for j = 2:nx-1
            U(:,j) = Un(:,j) - dt/dx*(F(:,j)-F(:,j-1));
        end      
        
        t = t + dt;
        n = n + 1;
    end
    
    fprintf('Roe:\n')
    fprintf('dx = %f\tCFL=%f\ttimestep num=%d\n',dx,CFL,n)
    solution(U,interval,'Roe',problem,t);
end

function F = Roe_flux(Un)
%ROE_FLUX numerical flux of Roe scheme

    gamma = 1.4;
    % the Roe averaged variables: u, H, a;
    resultn = conservative2primitive(Un); % rho, u, p
    rho_L = abs(resultn(1,1:end-1));
    rho_R = abs(resultn(1,2:end));
    u_L = resultn(2,1:end-1);
    u_R = resultn(2,2:end);
    E_L = Un(3,1:end-1);
    E_R = Un(3,2:end);
    velocity = (sqrt(rho_L).*u_L+sqrt(rho_R).*u_R)./(sqrt(rho_L)+sqrt(rho_R));
    H_L = gamma*E_L./rho_L - (gamma-1)/2*u_L.^2;
    H_R = gamma*E_R./rho_R - (gamma-1)/2*u_R.^2;
    enthalpy = (sqrt(rho_L).*H_L+sqrt(rho_R).*H_R)./(sqrt(rho_L)+sqrt(rho_R));
    sound = sqrt((gamma-1)*abs(enthalpy-1/2*velocity.^2));
    
    nx = size(Un,2);
    % the averaged eigenvalues lambda
    lambda = ones(3,nx-1);
    lambda(1,:) = velocity - sound;
    lambda(2,:) = velocity;
    lambda(3,:) = velocity + sound;
    
    % the averaged eigenvectors K
    K = ones(3,3,nx-1); % vecter index, component index, spatial index
    K(1,2,:) = velocity - sound;
    K(1,3,:) = enthalpy - velocity.*sound;
    K(2,2,:) = velocity;
    K(2,3,:) = 1/2*velocity.^2;
    K(3,2,:) = velocity + sound;
    K(3,3,:) = enthalpy + velocity.*sound;
     
    % jumps in the conserved quantity
    Delta_u = ones(3,nx-1);
    Delta_u(1,:) = rho_R - rho_L;
    Delta_u(2,:) = rho_R.*u_R - rho_L.*u_L;
    Delta_u(3,:) = E_R - E_L;
    % the wave strengths alpha
    alpha = ones(3,nx-1);
    alpha(2,:) = (gamma-1)./sound.^2.*(Delta_u(1,:).*(enthalpy-velocity.^2)+velocity.*Delta_u(2,:)-Delta_u(3,:));
    alpha(1,:) = 1/2./sound.*(Delta_u(1,:).*(velocity+sound)-Delta_u(2,:)-sound.*alpha(2,:));
    alpha(3,:) = Delta_u(1,:) - (alpha(1,:)+alpha(2,:));
    
    % numerical flux
    sum = zeros(3,nx-1);
    multi = sum;
    for i = 1:3
        multi(i,:) = 1/2*alpha(i,:).*abs(lambda(i,:));
    end
    for i = 1:3
        sum = sum + multi.*squeeze(K(i,:,:)); 
    end
    Fn = flux(Un);
    F_L = Fn(:,1:end-1);
    F_R = Fn(:,2:end);
    F = 1/2*(F_L+F_R) - sum;
end

