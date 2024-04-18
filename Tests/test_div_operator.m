clc; close all; clear;
% This file tests the divergence operator
% computed by 

% 3/22/24 shows order 2 

test_divergence = [];
h_list = [];
Nx_list = [5,10,20,40];
iter = 1;
test_poisson = [];


 while iter < 5

    Nx = Nx_list(iter);
    L = 1;       % Domain [0,1] x [0,1]
    h = L / Nx;
    h_list(iter) = h;
    Ny = Nx;

    Nt = 10;     % Number of timesteps
    dt = 0.001;
    epsilon = 1;

    r = epsilon * dt / (2*h^2);
    u_e = [];
    v_e = [];

    u_s = zeros(Nx+1, Nx); % Initialize u_star
    v_s = zeros(Nx, Nx+1); % Initialize v_star

    for i = 1:Nx+1
        for j = 1:Nx
            u_s(i,j) = u_star((i-1)*h,j*h-h/2);
        end
    end
    for i = 1:Nx
        for j = 1:Nx+1
            v_s(i,j) = v_star(i*h-h/2, (j-1)*h);
        end
    end
    
    % u_s
    % v_s
    dudx = div_u(u_s,Nx,h);
    dvdy = div_v(v_s,Nx,h);
    
    % True Divergence
    x = linspace(h/2, L-h/2, Nx);
    y = linspace(h/2, L-h/2, Nx);
    [X, Y] = meshgrid(x,y);
    true = -2*pi^2*cos(pi*X).*cos(pi*Y);

    numerical_divergence = dudx(:) + dvdy(:);
    %true_divergence = zeros(size(numerical_divergence));
    true_divergence = true(:);
    test_divergence(iter) = max(abs(numerical_divergence-true_divergence));
    iter = iter + 1;

 end
 test_divergence
 figure(1)
 subplot(1,2,1);
 surf(X, Y, reshape(numerical_divergence, [Nx, Nx]), 'EdgeColor', 'none');
 title('Approximate Divergence');
 xlabel('x');
 ylabel('y');
 zlabel('Divergence');

 subplot(1,2,2);
 surf(X, Y, true, 'EdgeColor', 'none');
 title('True Divergence');
 xlabel('x');
 ylabel('y');
 zlabel('Divergence');
 figure()

 figure(2)
 loglog(h_list,h_list.^2);
 hold on
 loglog(h_list, test_divergence, Marker="*")
 grid on
 xlabel("$h$", Interpreter="latex", FontSize=20)
 ylabel("$\Vert (\nabla \cdot u)_{approx}-(\nabla \cdot u)_{true} \Vert$", Interpreter="latex", FontSize=20)
 title("$O(h^2)$ Convergence of Divergence Operator", Interpreter="latex", FontSize=20)
 legend("$O(h^2)$", "error", Interpreter="latex", FontSize=20)


% du/dx calculation
function du = div_u(u_s,m,h)
for j = 1:m
    for i = 1:m
        du(i,j) = (u_s(i+1,j)-u_s(i,j))/(h);
    end
end
end

% dv/dy calculation
function dv = div_v(v_s,m,h)
for j = 1:m
    for i = 1:m
        dv(i,j) = (v_s(i,j+1)-v_s(i,j))/(h);
    end
end
end

function u = u_star(x,y)
%u = ((x^3)/3)*y^2;
%u = 2*cos(pi*x)*y^2;
u = -pi*sin(pi*x)*cos(pi*y);
end

function v = v_star(x,y)
%v = -((y^3)/3)*x^2;
%v = 2*pi*sin(pi*x)*(y^3/3);
v = -pi*cos(pi*x)*sin(pi*y);
end

