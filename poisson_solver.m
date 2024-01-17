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
