function outmat = hs_matc4(x,lambda,omega,rho,psi)
    % This is the 4x4 derivative matrix for the complex Ginzburg-Landau
    % equation, linearized about the Hocking-Stewartson Pulse.  The matrix
    % is given in a block form, with 2x2 sub blocks.  omega, rho, and psi
    % are parameters for the system, while lambda is the spectral parameter
    % and x is the non-autonomous coupling of the derivative to the pulse.
    outmat = zeros(4,4);

    % define the main components of the matrix, as per Afendikov and
    % Bridges 2001
    B    = [1, -omega; omega, 1];
    R    = [cos(psi), -sin(psi); sin(psi), cos(psi)];
    temp = omega*log(cosh(x));
    Q = [cos(temp)/cosh(x); -sin(temp)/cosh(x)];
    P = Q'*Q*eye(2) + 2*(Q*Q');
    P = B*(eye(2)+B)*P;
    P = B*B - P;
    
    % fill the off diagonal blocks
    outmat(3:4,1:2) = P + lambda*rho*R;
    outmat(1:2,3:4) = eye(2);
end