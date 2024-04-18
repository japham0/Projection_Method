clc; clear all; close all

% This script tests the averaging operators Ex and Ey 
% used to compute the advection term

Nx = 10;
L = 1;
h = L/Nx;

% Initialize u, v, P, S
u_s = zeros(Nx+1, Nx); 
for i = 1:Nx+1
    for j = 1:Nx
        u_s(i,j) = u_star((i-1)*h,j*h-h/2);
    end
end
v_s = zeros(Nx, Nx+1); 
for i = 1:Nx
    for j = 1:Nx+1
        v_s(i,j) = v_star(i*h-h/2, (j-1)*h);
    end
end

P_t = zeros(Nx,Nx);
for i = 1:Nx
    for j = 1:Nx
        P_t(i,j) = P(i*h-h/2, j*h-h/2);
    end
end

S_t = zeros(Nx+1, Nx+1);
for i = 1:Nx+1
    for j = 1:Nx+1
        S_t(i,j) = S((i-1)*h, (j-1)*h);
    end
end

% Test Ey function
% Grid options: v = 1, u = 2, P = 3, S = 4
fprintf("Testing Ey: \nNorm(Ey(v)-Ey_true(v)): ")
max(max(abs(Ey(1,v_s,Nx)-Ey_true(1,v_s,Nx, h))))

fprintf("Norm(Ey(u)-Ey_true(u)): ")
max(max(abs(Ey(2,u_s,Nx)-Ey_true(2,u_s,Nx, h))))

fprintf("Norm(Ey(P)-Ey_true(P)): ")
max(max(abs(Ey(3,P_t,Nx)-Ey_true(3,P_t,Nx,h))))

fprintf("Norm(Ey(S)-Ey_true(S)): ")
max(max(abs(Ey(4,S_t,Nx)-Ey_true(4,S_t,Nx,h))))

% Test Ex function
fprintf("Testing Ex: \nNorm(Ex(v)-Ex_true(v)")
max(max(abs(Ex(1, v_s, Nx)-Ex_true(1,v_s,Nx,h))))

fprintf("Norm(Ex(u)-Ex_true(u)): ")
max(max(abs(Ex(2, u_s, Nx)-Ex_true(2,u_s,Nx,h))))

fprintf("Norm(Ex(P)-Ex_true(P)): ")
max(max(abs(Ex(3, P_t, Nx)-Ex_true(3,P_t,Nx,h))))

fprintf("Norm(Ex(S)-Ex_true(S)): ")
max(max(abs(Ex(4, S_t, Nx)-Ex_true(4,S_t,Nx,h))))



% Compute Ey for different grids
function Ey = Ey(grid, w, Nx)

if grid == 1
% v-velocity
% stored on P-grid
    v = w;
    Ey_v = zeros(Nx,Nx);
    for i = 1:Nx
        for j = 1:Nx
            Ey_v(i,j) = (v(i,j+1) + v(i,j))/2;
            
        end
    end
    Ey = Ey_v;
end

if grid == 2
% u-velocity
% stored on S-grid
    u = w;
    Ey_u = zeros(Nx+1,Nx);
    for i = 1:Nx+1
        for j = 1:Nx-1
            Ey_u(i,j) = (u(i,j+1) + u(i,j))/2;
            % Ey_u(i,j) = 1;
        end
    end
    Ey = Ey_u;
end 

if grid == 3
% Pressure field
% stored on v-grid
    Pt = w;
    Ey_P = zeros(Nx,Nx+1);
    for i = 1:Nx
        for j = 1:Nx-1
            Ey_P(i,j+1) = (Pt(i,j+1) + Pt(i,j))/2;
            %Ey_P(i,j+1) = 1;
        end
    end
    Ey = Ey_P;
end

if grid == 4
% Stress field
% stored on u-grid
    St = w;
    Ey_S = zeros(Nx+1,Nx);
    for i = 1:Nx+1
        for j = 1:Nx
            Ey_S(i,j) = (St(i,j+1) + St(i,j))/2;
        end
    end
    Ey = Ey_S;
end
end


% Compute Ex for different grids
function Ex = Ex(grid, w, Nx)
if grid == 1
% v-velocity
% stored on Stress-grid
    v = w;
    Ex_v = zeros(Nx+1,Nx+1);
    for i = 1:Nx-1
        for j = 1:Nx+1
            Ex_v(i+1,j) = (v(i+1,j) + v(i,j))/2;
            %Ex_v(i+1,j)=1;
        end
    end
    Ex = Ex_v;
end
if grid == 2
% u-velocity
% stored on P-grid
    u = w;
    Ex_u = zeros(Nx,Nx);
    for i = 1:Nx
        for j = 1:Nx
            Ex_u(i,j) = (u(i+1,j) + u(i,j))/2;
            %Ex_u(i+1,j)=1;
        end
    end
    Ex = Ex_u;
end
if grid == 3
% Pressure-field
% stored on u-grid
    P_t = w;
    Ex_P = zeros(Nx+1,Nx);
    for i = 1:Nx-1
        for j = 1:Nx
            Ex_P(i+1,j) = (P_t(i+1,j) + P_t(i,j))/2;
            %Ex_P(i+1,j)=1;
        end
    end
    Ex = Ex_P;
end
if grid == 4
% Stress-field
% stored on v-grid
    S_t = w;
    Ex_S = zeros(Nx,Nx+1);
    for i = 1:Nx
        for j = 1:Nx+1
            Ex_S(i,j) = (S_t(i+1,j) + S_t(i,j))/2;
            %Ex_P(i+1,j)=1;
        end
    end
    Ex = Ex_S;
end
end

function Ex_true = Ex_true(grid, w, Nx, h)
if grid == 1
% v-velocity
% stored on S-grid
    Ex_v_true = zeros(Nx+1,Nx+1);
    for i = 1:Nx-1
        for j = 1:Nx+1
            Ex_v_true(i+1,j) = ...
                (v_star((i+1)*h-h/2, (j-1)*h) + v_star(i*h-h/2, (j-1)*h))/2;
        end
    end
    Ex_true = Ex_v_true;
end
if grid == 2
% u-velocity
% stored on P-grid
    Ex_v_true = zeros(Nx,Nx);
    for i = 1:Nx
        for j = 1:Nx
            Ex_v_true(i,j) = ...
                (u_star(i*h,j*h-h/2) + u_star((i-1)*h,j*h-h/2))/2;
        end
    end
    Ex_true = Ex_v_true;
end 
if grid == 3
% Pressure field
% stored on u-grid
    Ex_P_true = zeros(Nx+1,Nx);
    for i = 1:Nx-1
        for j = 1:Nx
            Ex_P_true(i+1,j) = ...
                (P((i+1)*h-h/2, j*h-h/2)+P(i*h-h/2,j*h-h/2))/2;
        end
    end
    Ex_true = Ex_P_true;
end 
if grid == 4
% Stress field
% stored on v-grid
    Ex_S_true = zeros(Nx,Nx+1);
    for i = 1:Nx
        for j = 1:Nx+1
            Ex_S_true(i,j) = ...
                (S(i*h, (j-1)*h)+S((i-1)*h, (j-1)*h))/2;
        end
    end
    Ex_true = Ex_S_true;
end 
end



function Ey_true = Ey_true(grid, w, Nx, h)
if grid == 1
% v-velocity
% stored on P-grid
    Ey_v_true = zeros(Nx,Nx);
    for i = 1:Nx
        for j = 1:Nx
            Ey_v_true(i,j) = ...
                (v_star(i*h-h/2, j*h) + v_star(i*h-h/2, (j-1)*h))/2;
        end
    end
    Ey_true = Ey_v_true;
end
if grid == 2
% u-velocity
% stored on S-grid
    Ey_u_true = zeros(Nx+1,Nx);
    for i = 1:Nx+1
        for j = 1:Nx-1
            Ey_u_true(i,j) = ...
                (u_star((i-1)*h,(j+1)*h-h/2) + u_star((i-1)*h,j*h-h/2))/2;
        end
    end
    Ey_true = Ey_u_true;
end
if grid == 3
% Pressure field
% stored on v-grid
    Ey_P_true = zeros(Nx, Nx+1);
    for i = 1:Nx
        for j = 1:Nx-1
            Ey_P_true(i,j+1) =...
            (P(i*h-h/2, (j+1)*h-h/2)+P(i*h-h/2,j*h-h/2))/2;
        end
    end
    Ey_true = Ey_P_true;
end

if grid == 4
% Stress field
% stored on u-grid
    Ey_S_true = zeros(Nx+1, Nx);
    for i = 1:Nx+1
        for j = 1:Nx
            Ey_S_true(i,j) =...
            (S((i-1)*h, j*h)+S((i-1)*h, (j-1)*h))/2;
        end
    end
    Ey_true = Ey_S_true;
end
end





function u = u_star(x,y)
u = sin(pi*x)*sin(pi*y);
end

function v = v_star(x,y)
v = sin(pi*x)*sin(pi*y);
end
function P = P(x,y)
P = sin(pi*x)*sin(pi*y);
end
function S = S(x,y)
S = sin(pi*x)*sin(pi*y);
end