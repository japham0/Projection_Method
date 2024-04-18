clc; close all; clear;

% This script solves the heat equation u_t = epsilon*u_xx
% on a MAC grid, using the Crank-Nicholson method. The 
% grid follows the structure for the y-component of the
% velocity field. It shows the error decreasing with 
% O(h^2) convergence on a loglog plot.

test_divergence = [];
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
    v_e = [];

    v_s = zeros(Nx, Nx+1); % Initialize v_star
    for i = 1:Nx
        for j = 1:Nx+1
            v_s(i,j) = v_star(i*h-h/2, (j-1)*h);
        end
    end

    % True solution on MAC Grid
    x = linspace(h/2, L-h/2, Nx);
    y = linspace(0, L, Nx+1);
    [X, Y] = meshgrid(x,y);
    v_true = zeros(size(v_s));
    v_true_time = [];

    for t = 1:Nt
        v_true = exp(-2*t*dt*pi^2)*sin(pi*X).*sin(pi*Y);
        v_true_time = [v_true_time, v_true(:)];
    end
    surf(x,y,reshape(v_true, Nx+1, Nx))

    bv = zeros(Nx*(Nx+1), 1);
    Av = zeros(Nx*(Nx+1));   

    for t = 1:Nt

        if t == 1
            v_s = v_s(:);
        end
        bv = build_Rv(v_s,r, Nx);
        Av = build_Av(Nx,dt,h, epsilon);
        v_step = Av \ bv;
        v_s = v_step;

    end

    v_true_end = v_true_time(:,end);
    mac_err(iter) = max(abs(v_true_end-v_step))
    iter = iter + 1;

end

figure()
loglog(h_list, 10*h_list.^2);
hold on
loglog(h_list, mac_err)
grid on
xlabel("$h$", Interpreter="latex", FontSize=20)
ylabel("$\Vert v_{true}-v_{approx} \Vert$", Interpreter="latex", FontSize=20)
title("Convergence of Solution to Heat Equation on MAC Grid", Interpreter="latex", FontSize=20)
legend("$O(h^2)$", "error", Interpreter="latex", FontSize=20)

figure()
surf(x,y,reshape(v_step, Nx+1, Nx))



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
