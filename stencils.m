clc; clear all; close all;

n = 10;

test_neumann2(n)

m = 20;

test_neumann2(m)

p = 40;

test_neumann2(p)

% test_neumann2(n)

% test_neumann1(n)

% u_star = CN(n);

% u_star_reshape = reshape(u_star, n-1, n-1);


function out = CN(n)

% This function returns u* from the Crank-Nicholson
% scheme (2I - 0.5dt*mu P)u* = (2I + 0.5dt*mu P)u
n = 5;
dt = 0.01;
mu = 0.01; 
P = A_d(n);

A_cn = eye((n-1)^2) - 0.5*dt*mu*P;
R    = eye((n-1)^2) + 0.5*dt*mu*P;

x = 1/n * (1:n-1);
y = 1/n * (1:n-1);
u = [];
for j = 1:n-1
    u = [u, (x.^2 + 5*x - 2)*cos(y(j))];
end
u = u';

u_star = A_cn \ R*u ;
out = u_star
end 

function out = test_neumann1(n)

% This function solves the poisson equation 
% u_xx + u_yy = -13sin(3x + 2y) on the square
% [0,1] x [0,1] with Neumann boundary conditions 
% and a step size dx = dy = h = 1/n.

h = 1/n; 

%% NEUMANN BOUNDARY
% sin(3x+2y)
% du / dy = 2cos(3x+2y)
% du / dx = 3cos(3x+2y)

x = linspace(0,1, n+1);
y = linspace(0,1, n+1);

% lower boundary: y=0  -du/dy
low   = -2*cos(3*x);

% left boundary:  x=0  -du/dx
left  = -3*cos(2*y);

% right boundary: x=1   du/dx
right = 3*cos(3+2*y);

% upper boundary: y=1   du/dy
up    = 2*cos(3*x+2);

%% RHS

F=[];
for j=1:n+1
   F =[F, 13*sin(3*x+2*y(j))];
end

b = F;
%b'
% lower border
b(1:n+1) = b(1:n+1)+2*(1/h)*low;
% left border
b(1:n+1:(n+1)*n+1) = b(1:n+1:(n+1)*n+1) + 2*(1/h)*left;
% right border
b(n+1:n+1:(n+1)^2) = b(n+1:n+1:(n+1)^2) + 2*(1/h)*right;
% upper border
b((n+1)*n+1:(n+1)^2) = b((n+1)*n+1:(n+1)^2) + 2*(1/h)*up;
b = h^2 * b';

%rank(A)
b(1) = sin(3*x(1)+2*y(1)); % b(1) = 0

%% SOLVE SYSTEM
% u = A \ b;
u = A_n(n) \ b;


%% COMPARE TO TRUE SOLUTION

true = [];
for j = 1:n+1
   true = [true, sin(3*x+2*y(j))];
end

true = reshape(true, n+1, n+1);

% PLOTTING
figure;

% Plotting the true solution
subplot(1, 2, 1);
surf(x, y, true);
title('True Solution');
xlabel('x');
ylabel('y');
zlabel('u');

% Plotting the approximate solution
subplot(1, 2, 2);
surf(x, y, reshape(u, n+1, n+1));
title('Approximate Solution');
xlabel('x');
ylabel('y');
zlabel('u');

% return abs error from true solution
% out = norm(u' - true(:));
end

function out = test_neumann2(n)

% This function solves the poisson equation 
% u_xx + u_yy = (2-5x-x^2)cos(y) on the square
% [0,1] x [0,1] with Neumann boundary conditions 
% and a step size dx = dy = h = 1/n.

h = 1/n; 

%% NEUMANN BOUNDARY

x = linspace(0,1, n+1);
y = linspace(0,1, n+1);

% lower boundary: y=0  -du/dy
low   = 0*x;

% left boundary:  x=0  -du/dx
left  = -5*cos(y);

% right boundary: x=1   du/dx
right = 7*cos(y);

% upper boundary: y=1   du/dy
up    = -(x.^2+5*x).*sin(1);

%% RHS
F=[];
for j=1:n+1
   F =[F, (-x.^2-5*x+2)*(cos(y(j)))];
end

b = F;

% lower border
b(1:n+1) = b(1:n+1)+2*(1/h)*low;
% left border
b(1:n+1:(n+1)*n+1) = b(1:n+1:(n+1)*n+1) + 2*(1/h)*left;
% right border
b(n+1:n+1:(n+1)^2) = b(n+1:n+1:(n+1)^2) + 2*(1/h)*right;
% upper border
b((n+1)*n+1:(n+1)^2) = b((n+1)*n+1:(n+1)^2) + 2*(1/h)*up;

b = h^2 * b';

b(1) = 0; % (x^2 + 5x)cos(y)=0 at x=y=0




u = A_n(n)\b;
true = [];
for j = 1:n+1
   true = [true, (x.^2+5*x)*(cos(y(j)))];
end
true = reshape(true, n+1, n+1);

% PLOTTING
figure;

% Plotting the true solution
subplot(1, 2, 1);
surf(x, y, true);
title('True Solution');
xlabel('x');
ylabel('y');
zlabel('u');

% Plotting the approximate solution
subplot(1, 2, 2);
surf(x, y, reshape(u, n+1, n+1));
title('Approximate Solution');
xlabel('x');
ylabel('y');
zlabel('u');

% return abs error from true solution
% out = norm(u' - true(:));

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
A = full(A);

A(1,1) = 1;
A(1,2) = 0;
A(1,n+2) = 0;

out = A;
end



function out = test_dirichlet(n)

% This function solves the poisson equation 
% u_xx + u_yy = (2-5x-x^2)cos(y) on the square
% [0,1] x [0,1] with Dirichlet boundary conditions 
% and a step size dx = dy = h = 1/n.
% n is input, number of points in one spacial dim

h = 1/n; 
x = 1/n*(0:n);
y = 1/n*(0:n);

% lower boundary
g_low = x.^2 + 5*x;

% upper boundary
g_up  = (x.^2 + 5*x) * cos(1);

% left boundary
g_left = 0 .* y;

% right boundary
g_right = 6 * cos(y);


F = [];
x = 1/n * (1:n-1);
y = 1/n * (1:n-1);
for j = 1:n-1
    F = [F, (x.^2 + 5*x - 2)*cos(y(j))];
end

b = F;

% low
b(1:n-1) = b(1:n-1) + (1/h^2) * g_low(2:n);

%left
b(1:n-1:(n-1)*(n-2)+1) = b(1:n-1:(n-1)*(n-2)+1)+(1/h^2) * g_left(2:n);

% right
b((n-1):n-1:(n-1)^2) = b(n-1:n-1:(n-1)^2) + (1/h^2) * g_right(2:n);

% up
b((n-1)*(n-2)+1:(n-1)^2) = b((n-1)*(n-2)+1:(n-1)^2) + (1/h^2) * g_up(2:n);

b = h^2 * b;
b = b';

% Approximate solution
% approx_solution = A_d(n,b);
approx_solution = A_d(n) \ b;
 
% True solution
x = linspace(1/n, 1-1/n, n-1);
y = linspace(1/n, 1-1/n, n-1);

true_solution = [];
for j = 1:n-1
    true_solution = [true_solution, (5*x+x.^2)*cos(y(j))];
end

% norm((approx_solution'-true_solution), "inf")
% dA = h^2
out = norm((approx_solution'-true_solution)*h, 2);

% Reshape solutions for surface plotting
approx_solution_reshaped = reshape(approx_solution, [n-1, n-1]);
true_solution_reshaped = reshape(true_solution, [n-1, n-1]);
% 
% max(max(approx_solution_reshaped-true_solution_reshaped))
% Plot approximate solution
% figure;
% subplot(1, 2, 1);
% surf(x, y, approx_solution_reshaped - true_solution_reshaped);
% title('Approximate Solution');
% xlabel('x');
% ylabel('y');
% zlabel('Solution');

% Plot true solution
% subplot(1, 2, 2);
% surf(x, y, true_solution_reshaped);
% title('True Solution');
% xlabel('x');
% ylabel('y');
% zlabel('Solution');
end


function out = A_d(n, b)

% This function returns the (n-1) x (n-1) 
% matrix A which solves the homogeneous
% Poisson equation with Dirichlet boundary conditions

off_diag = ones((n-1)^2,1);
A = spdiags([-off_diag -off_diag 4*off_diag -off_diag -off_diag],...
    [-(n-1) -1 0 1 n-1],(n-1)^2,(n-1)^2);

for j=1:n-2
    A(j*(n-1),j*(n-1)+1) = 0;
    A(j*(n-1)+1,j*(n-1)) = 0;
end

A = full(A);

out = A;

end

