clc; close all; clear;

% This script solves the heat equation u_t = epsilon * u_xx
% on a MAC grid using the Crank-Nicholson scheme. It produces
% the plot showing O(h^2)convergence. The grid follows the 
% storage structure for the x-component of the velocity field

h_list = [];
Nx_list = [5,10,20,40];
iter = 1;
mac_err = [];
while iter < 5
Nx = Nx_list(iter);
L = 1;       % Domain [0,1] x [0,1]
h = L / Nx;
h_list(iter) = h;
Ny = Nx;

Nt = 400;     % Number of timesteps
dt = 0.001;  
epsilon = 1;

r = epsilon * dt / (2*h^2);
u_e = [];
v_e = [];

u_s = zeros(Nx+1, Nx); % Initialize u_star

for i = 1:Nx+1
    for j = 1:Nx
        u_s(i,j) = u_star((i-1)*h,j*h-h/2);
    end
end
u_s_init = u_s;

% True solution
x = linspace(0, L, Nx+1);
y = linspace(h/2, L-h/2, Nx);
[X, Y] = meshgrid(x,y);
u_true = zeros(size(u_s));
u_true_time = [];

for t = 1:Nt
    u_true = exp(-2*t*dt*pi^2)*sin(pi*X).*sin(pi*Y);
    u_true_time = [u_true_time, u_true(:)];
end

surf(x,y,reshape(u_true, Nx, Nx+1))

bu = zeros(Nx*(Nx+1), 1);
Au = zeros(Nx*(Nx+1));

for t = 1:Nt

    if t == 1
        u_s = u_s(:);
    end

    bu = build_Ru(u_s,r, Nx);

    Au = build_Au(Nx,dt,h, epsilon);
    u_step = Au \ bu;
    u_s = u_step;

end

u_true_end = u_true_time(:,end);
mac_err(iter) = max(abs(u_true_end-u_step));
iter = iter + 1;

end
figure()
loglog(h_list, 10*h_list.^2);
hold on
loglog(h_list, mac_err)
grid on
xlabel("$h$", Interpreter="latex", FontSize=20)
ylabel("$\Vert u_{true}-u_{approx} \Vert$", Interpreter="latex", FontSize=20)
title("Convergence of Solution to Heat Equation on MAC Grid", Interpreter="latex", FontSize=20)
legend("$O(h^2)$", "error", Interpreter="latex", FontSize=20)

figure()
surf(x,y,reshape(u_true_end-u_step, Nx, Nx+1))



function out = build_Ru(u_s, r, Nx)
bu = zeros(Nx*(Nx+1),1);
for i = 2:Nx
    bu(i) = u_s(i-1) - 5*u_s(i) + u_s(i+1) + u_s(i+(Nx+1));
    n = i + (Nx+1) * (Nx-1);
    bu(n) = u_s(n-1) - 5*u_s(n) + u_s(n+1) + u_s(n-(Nx+1));
end

for j = 1:Nx
    n = 1 + (j-1)*(Nx+1);
    bu(n) = 0;
    
    m = (Nx+1)+ (j-1)*(Nx+1);
    bu(m) = 0;
end

for i = 2:Nx
    for j = 2:Nx-1
        n = i + (j-1)*(Nx+1);
        bu(n) = u_s(n) + r*(u_s(n-1) - 2*(u_s(n)) + u_s(n+1)) + ...
                         r*(u_s(n-(Nx+1)) - 2*(u_s(n)) + u_s(n-(Nx+1)));
    end
end

out = bu;
end


function out = build_Rv(v_s, r, Nx)
bv = zeros(Nx*(Nx+1),1);
for j = 2:Nx
    m = 1 + (j-1)*Nx;
    bv(m) = v_s(m-Nx) - 5*v_s(m) + v_s(m+1) + v_s(m+Nx);
    
    bv(j*Nx) = v_s(j*Nx-Nx) + v_s(j*Nx-1) - 5*v_s(j*Nx) + v_s(j*Nx+Nx);
end

for i = 1:Nx
    bv(i) = 0;
    m = i + Nx*Nx;
    bv(m) = 0;
end

for i = 2:Nx-1
    for j = 2:Nx
        n = i + (j-1) * Nx;
        bv(n) = v_s(n) + r*(v_s(n-1) - 2*v_s(n) + v_s(n+1))...
                       + r*(v_s(n-Nx) - 2*v_s(n) + v_s(n+Nx));
    end
end

out = bv;
end



function out = build_Au(Nx, dt, h, epsilon)
% This function constructs the Crank-Nicholson
% LHS matrix with homogenous Dirichlet BC's

r = epsilon * dt / (2*h^2);
Ny = Nx;

% Initialize matrices for the implicit scheme
A = zeros(Nx*(Nx+1), Nx*(Nx+1));

% Set up the coefficients for the matrix A
for i = 1:Nx+1
    for j = 1:Nx
        n = i + (j - 1) * (Nx+1); 

        if i == 1 
            % Homogeneous Dirichlet boundary condition
            A(n, n) = 1;
            b(n) = 0; % Boundary value is set to zero
        
        end
        if i == Nx+1
            A(n,n) = 1;
            
        end 
        if j == 1 && i ~= 1 && i ~= Nx+1
            A(n, n) = -5;
            A(n, n-1) = 1;
            A(n, n+1) = 1;
            A(n, n+(Nx+1))= 1;
        
        end

        if j ~= 1 && i ~= 1 && i ~= Nx+1 && j~=Nx
            A(n, n) = 1 + 4*r;
            A(n, n-1) = -r;
            A(n, n+1) = -r;
            A(n, n-(Nx+1))= -r;
            A(n, n+(Nx+1))= -r;
        end

        if j == Nx && i ~= 1 && i ~= Nx+1 
            A(n, n) = -5;
            A(n, n-1) = 1;
            A(n, n+1) = 1;
            A(n, n-(Nx+1))= 1;
        end
    end
end

out = A;
end


function out = build_Av(Nx, dt, h, epsilon)
% This function constructs the Crank-Nicholson
% LHS matrix with homogenous Dirichlet BC's

r = epsilon * dt / (2*h^2);
Ny = Nx;

% Initialize matrices for the implicit scheme
A = zeros(Nx*(Nx+1), Nx*Nx);


% Set up the coefficients for the matrix A
for i = 1:Nx
    for j = 1:Nx+1
        n = i + (j - 1) * Nx; % Mapping grid to Matrix A in LHS
        if  j == 1 || j == Nx+1

            A(n, n) = 1;

        else
            A(n, n) = 1 + 4*r;
            A(n, n-1) = -r;
            A(n, n+1) = -r;
            A(n, n-Nx) = -r;
            A(n, n+Nx) = -r;
        end
        if (i == 1 || i == Nx) && (j~=1) && (j~=Nx+1)
            A(n,n) = -5;
            if i == 1
                A(n, n-1)=0;
                A(n, n+1)=1;
                A(n, n-Nx)=1;
                A(n, n+Nx)=1;
            end
            if i == Nx
                A(n,n+1)=0;
                A(n, n+Nx)=1;
                A(n, n-1)=1;
                A(n, n-Nx)=1;
            end
        end
    end
end
out = A;
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
%u = sin(pi*x)*cos(pi*y);
%u = ((x^3)/3)*y^2;
%u = -pi*sin(pi*x)*cos(pi*y);
u = sin(pi*x)*sin(pi*y);
end

function v = v_star(x,y)
%v = -cos(pi*x)*sin(pi*y);
%v = -((y^3)/3)*x^2;
%v = -pi*cos(pi*x)*sin(pi*y);
v = sin(pi*x)*sin(pi*y);
end



function P = P_matrix(n)
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
P = L;
for j = 1:n
    P = [P,spalloc(j*(n+1),n+1,0)
        spalloc(n+1,j*(n+1),0),L];
end
P = P + spdiags([-e -e],[-n-1 n+1],(n+1)^2,(n+1)^2);
for j=1:n+1
    P(j,j+n+1) = -2;
    P((n+1)^2+1-j,(n+1)*n+1-j) = -2;
end
% for neumann only test use:
P = -full(P);

%P = -n^2*full(P);


%% SET CORNERS & FIX 1 EQN
P = full(P);
P(1,1) = 1;
P(1,2) = 0;
P(1,n+2) = 0;
%full(P)
%rank(P)
end

function out = A_n(n)

% This function returns the n+1 x n+1
% matrix A to solve the homogenous 
% Poisson equation with Neumann boundary conditions

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
A = -full(A);
A
rank(A)
A(1,1) = 1;
A(1,2) = 0;
A(1,n+2) = 0;
A
rank(A)
out = A;
end

function P_m = P_m(P_m,m)

% interior nodes
for j = 2:m-1
    for i = 2:m-1
        P_m(i+(j-1)*m,i+(j-1)*m-m) = 1;
        P_m(i+(j-1)*m,i+(j-1)*m-1) = 1;
        P_m(i+(j-1)*m,i+(j-1)*m) = -4;
        P_m(i+(j-1)*m,i+(j-1)*m+1) = 1;
        P_m(i+(j-1)*m,i+(j-1)*m+m) = 1;
    end
end

% west wall
for j = 2:m-1
    i = 1;
    P_m(i+(j-1)*m,i+(j-1)*m-m) = 1;
    P_m(i+(j-1)*m,i+(j-1)*m) = -3;
    P_m(i+(j-1)*m,i+(j-1)*m+1) = 1;
    P_m(i+(j-1)*m,i+(j-1)*m+m) = 1;
end

% east wall
for j = 2:m-1
    i = m;
    P_m(i+(j-1)*m,i+(j-1)*m-m) = 1;
    P_m(i+(j-1)*m,i+(j-1)*m-1) = 1;
    P_m(i+(j-1)*m,i+(j-1)*m) = -3;
    P_m(i+(j-1)*m,i+(j-1)*m+m) = 1;
end

% south wall
for i = 2:m-1
    j = 1;
    P_m(i+(j-1)*m,i+(j-1)*m-1) = 1;
    P_m(i+(j-1)*m,i+(j-1)*m) = -3;
    P_m(i+(j-1)*m,i+(j-1)*m+1) = 1;
    P_m(i+(j-1)*m,i+(j-1)*m+m) = 1;
end

% north wall
for i = 2:m-1
    j = m;
    P_m(i+(j-1)*m,i+(j-1)*m-m) = 1;
    P_m(i+(j-1)*m,i+(j-1)*m-1) = 1;
    P_m(i+(j-1)*m,i+(j-1)*m) = -3;
    P_m(i+(j-1)*m,i+(j-1)*m+1) = 1;
end

% west-south corner
P_m(1,1) = 1;

% south-east corner
i = m;
j = 1;

P_m(i+(j-1)*m,i+(j-1)*m-1) = 1;
P_m(i+(j-1)*m,i+(j-1)*m) = -2;
P_m(i+(j-1)*m,i+(j-1)*m+m) = 1;

% east-north corner
i = m;
j = m;

P_m(i+(j-1)*m,i+(j-1)*m-m) = 1;
P_m(i+(j-1)*m,i+(j-1)*m-1) = 1;
P_m(i+(j-1)*m,i+(j-1)*m) = -2;

% north-west corner
i = 1;
j = m;

P_m(i+(j-1)*m,i+(j-1)*m-m) = 1;
P_m(i+(j-1)*m,i+(j-1)*m) = -2;
P_m(i+(j-1)*m,i+(j-1)*m+1) = 1;

end