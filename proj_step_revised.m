clc; close all; clear;

proj_step(4)

function out = proj_step(m)

% This function executes one projection step on 
% the domain [0,1]x[0,1] with a grid size m x m. 
% It follows the calculation:

% u* = u_df + (nabla phi)
% div u* = div u_df + nabla^2 phi
% div u* = nabla^2 phi
% solves poisson equation

% And returns the error between u_divfree and u_true

L = 1;
h = L / m;

u_s = zeros(m+1, m);
u_t = zeros(m+1, m);

v_s = zeros(m, m+1);
v_t = zeros(m, m+1);

dudx = zeros(m,m); % Can only have derivatives on m x m grid using FD
dvdy = zeros(m,m);


%% Initialize u*
for i = 1:m+1
    for j = 1:m
        u_s(i,j) = u_star((i-1)*h,j*h-h/2);
        u_t(i,j) = u_true((i-1)*h,j*h-h/2);
    end
end

for i = 1:m
    for j = 1:m+1
        v_s(i,j) = v_star(i*h-h/2, (j-1)*h);
        v_t(i,j) = v_true(i*h-h/2, (j-1)*h);
    end
end
u_s
v_s

%% Solve Poisson Equation
A = A_matrix(m-1); % Construct Poisson matrix
dudx = div_u(u_s,m,h);
dvdy = div_v(v_s,m,h);

% POISSON SOLVER
% nabla^2 u = b
% b = dudx + dvdy
sol = poisson_solver(A,m,dudx,dvdy);

% u_df = u* - grad(phi)
phi = pressure(m,sol);

%% Compute divergence-free part of u
for j = 1:m
    for i = 2:m
        % rows 2 through m
        u_df(i,j) = u_s(i,j) - (phi(i,j)-phi(i-1,j))/h;
    end
    % row 1
    u_df(1,j) = u_s(1,j) - phi(1,j)/h;
    % row m + 1
    u_df(m+1,j) = u_s(m+1,j) - phi(m,j)/h;
end

for j = 2:m
    % columns 2 through m
    for i = 1:m
        v_df(i,j) = v_s(i,j) - (phi(i,j)-phi(i,j-1))/h;
        % column 1
        v_df(i,1) = v_s(i,1) - (phi(i,1))/h;
        % column m+1
        v_df(i,m+1) = v_s(i,m+1) - phi(i,m)/h;
    end
end

% figure(2)
% mesh(v_t-v_df)
% xlabel("x")
% ylabel("y")
% figure(3)
% mesh(u_t-u_df)
% xlabel("x")
% ylabel("y")
% figure(4)
% plot(b)

% Difference between u_div_free and u_true
u_e = max(max(u_df-u_t));
v_e = max(max(v_df-v_t));

out = [u_e, v_e];
end

function out = pressure(m, sol)
% This function returns the pressure field on
% the mxm grid given the vector from the solution 
% to the poisson equation

phi = zeros(m*m);
% update phi (matrix)
for j = 1:m
    for i = 1:m
        phi(i,j) = sol(i+(j-1)*m);
    end
end

out = phi;
end

% div_u calculation
function du = div_u(u_s,m,h)
for j = 1:m
    for i = 1:m
        du(i,j) = (u_s(i+1,j)-u_s(i,j))/h;
    end
end
end

% div_v calculation
function dv = div_v(v_s,m,h)
for j = 1:m
    for i = 1:m
        dv(i,j) = (v_s(i,j+1)-v_s(i,j))/h;
    end
end
end

function u = u_star(x,y)
u = sin(x)*cos(y);
%u = (x^3/3) * y^2;
end

function u = u_true(x,y)
u = sin(x)*cos(y);
%u = (x^3/3) * y^2;
end
function v = v_star(x,y)
v = -cos(x)*sin(y);
%v = (-y^3/3) * x^2;
end
function v = v_true(x,y)
v = -cos(x)*sin(y);
%v = (-y^3/3) * x^2;
end


function A = A_matrix(n)
% This function returns the matrix for solving the 
% 2D Poisson equation with Neumann boundary conditions
h = 1/n;
x = linspace(0,1, n+1);
y = linspace(0,1, n+1);

I = eye(n+1);
I(1,1) = 1/2;
I(n+1,n+1) = 1/2;
D = kron(I,I);

e = ones((n+1)^2,1);
v = D*e;
L = spdiags([-e 4*e -e],[-1 0 1],(n+1),(n+1));
L(1,2)  = -2;
L(n+1,n)= -2;
A = L;
for j = 1:n
    A = [A,spalloc(j*(n+1),n+1,0)
        spalloc(n+1,j*(n+1),0),L];
end
A = A + spdiags([-e -e],[-n-1 n+1],(n+1)^2,(n+1)^2);
for j=1:n+1
    A(j,j+n+1) = -2;
    A((n+1)^2+1-j,(n+1)*n+1-j) = -2;
end
A = n^2*A;

%% SET CORNERS & FIX 1 EQN
A = full(A);
A(1,1) = 1;
A(1,2) = 0;
A(1,n+2) = 0;
end

function sol = poisson_solver(A,m,dudx,dvdy)

% u* = u_df + (nabla phi)
% div u* = div u_df + nabla^2 phi
% div u* = nabla^2 phi
% need to solve poisson equation

% POISSON SOLVER
% nabla^2 u = b
% b = dudx + dvdy

b = zeros(m*m,1);
sol = zeros(m*m,1);

% initialize b
for j = 1:m
    for i = 1:m
        b(i+(j-1)*m) = dudx(i,j) + dvdy(i,j);
    end
end

sol = A \ b;
end
