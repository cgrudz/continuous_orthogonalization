function dYdx = stiefel_derivative(x,Y,lambda,omega,rho,psi)

    % This is the stiefel manifold derivative of the CGL equation.
    A = zeros(4,4);

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
    A(3:4,1:2) = P + lambda*rho*R;
    A(1:2,3:4) = eye(2);
    
    % reshape the frame from col vector into a two frame
    Y = reshape(Y,4,2);
    YH = Y';
    
    % define dYdx
    M = skewc( A - 2*Y*YH*symc(A) );
    dYdx = M*Y;

    dYdx = dYdx(:);
end