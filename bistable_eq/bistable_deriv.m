function output = bistable_deriv(t,Y,A,lambda)
    output = (A(t,lambda)-eye(2)*lambda)*Y;
end