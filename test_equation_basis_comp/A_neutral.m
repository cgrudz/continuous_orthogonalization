function Y_prime = A_neutral(xi,Y,lambda)

    A = [0,                        1;
        lambda + 1 - 6/cosh(xi)^2, 0];
    
    sigma = sqrt(lambda+1);    
    Y_prime = (A - sigma*eye(2))*Y;