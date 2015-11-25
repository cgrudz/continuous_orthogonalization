function deriv=drury_edit(t,y,lambda,A_mat,n,k)
% ydot=drury_edit(t,y,lambda,A,s,p,n,k,mu,damping)
%
% Returns the ODE output for the polar method using the method of Drury
%
% Input "t" and "y" are provided by ode45, "A" is a function handle to the
% desired Evans matrix, p is a given parameter,
%"n" is the dimension of the system and "k" is the
% dimension of the manifold.

W = reshape(y,n,k);
A_temp = A_mat(t,lambda);
deriv = [reshape((eye(n)-W*W')*A_temp*W,n*k,1)];
end
