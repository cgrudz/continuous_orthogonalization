function A = A_neutral(xi,lambda)

    A = [0,                        1;
        lambda + 1 - 6/cosh(xi)^2, 0];
    
    sigma = sqrt(lambda+1);    
    A = (A - sigma*eye(2));