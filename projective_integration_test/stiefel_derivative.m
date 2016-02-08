function dQdx = stiefel_derivative(x,Q,derivative_handle)

A = derivative_handle(x);
Q = reshape(Q,4,2);
QH = conj(Q)';

M = skewc( A - 2*Q*QH*symc(A) );
dQdx = M*Q;

dQdx = dQdx(:);
