function out = A_bistable(x,lambda)

out = [0,                           1;
       lambda + 1 - 6/cosh(x)^2, 0];
end