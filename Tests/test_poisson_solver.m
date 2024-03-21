clc; close all; clear;

% This scripts solves the poisson equation 
% nabla^2 phi = -2pi^2cos(pi*x)*sin(pi*x)
% with homogeneous Neumann boundary conditions.
% It produces the plot showing O(h^2) convergence 

h_list = [];
Nx_list = [4,10,20,40];
iter = 1;
test_poisson = [];


while iter < 5
    Nx = Nx_list(iter);
    test_poisson(iter) = test_neumann1(Nx);
    h_list(iter) = 1/Nx;
    iter = iter + 1;
end
loglog(h_list,h_list.^2);
hold on
loglog(h_list, test_poisson)
grid on
xlabel("$h$", Interpreter="latex", FontSize=20)
ylabel("$\Vert \phi_{true}-\phi_{approx} \Vert$", Interpreter="latex", FontSize=20)
title("$O(h^2)$ Poisson Solver", Interpreter="latex", FontSize=20)
legend("$O(h^2)$", "Error", Interpreter="latex", FontSize=14)


% Functions
function out = test_neumann1(n)
h = 1/n; 

%% NEUMANN BOUNDARY

% True: phi = cos(pi*x)cos(pi*y)
% du / dx = -pi*sin(pi*x)cos(pi*y)
% du / dy = -pi*cos(pi*x)sin(pi*y)

x = linspace(0,1, n+1);
y = linspace(0,1, n+1);

% lower boundary: y=0  -du/dy
low   = 0*x;

% left boundary:  x=0  -du/dx
left  = 0*x;

% right boundary: x=1   du/dx
right = 0*x;

% upper boundary: y=1   du/dy
up    = 0*x;

%% RHS

F=[];
for j=1:n+1
   F =[F, -2*pi^2*cos(pi*x)*cos(pi*y(j))];
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

%rank(A)
b(1) = cos(pi*x(1))*cos(pi*y(1)); % b(1) = 0

%% SOLVE SYSTEM
% u = A \ b;
% u = A_n(n) \ b;
 phi = P_matrix(n) \ b;
% u = P_m(zeros(n+1), n+1) \ b;
% u = pinv(P_matrix(n))*b;

%% COMPARE TO TRUE SOLUTION

true = [];
for j = 1:n+1
   true = [true, cos(pi*x)*cos(pi*y(j))];
end

true = reshape(true, n+1, n+1);
out = max(max(abs(phi-true(:))));

% PLOTTING
% Plotting the true solution
figure(1)
surf(x, y, true);
title('True Solution');
xlabel('x');
ylabel('y');
zlabel('u');

% Plotting the approximate solution
figure(2)
surf(x, y, reshape(phi, n+1, n+1));
title('Approximate Solution');
xlabel('x');
ylabel('y');
zlabel('u');

% Plotting the difference
% figure(3)
% surf(x, y, true-reshape(u, n+1, n+1));
% title('True - Approximate Solution');
% xlabel('x');
% ylabel('y');
% zlabel('u');

% return abs error from true solution
%out = norm(u' - true(:));
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