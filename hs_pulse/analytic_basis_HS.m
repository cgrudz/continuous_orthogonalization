function [out, projects] = analytic_basis_HS(projection,x,preimage, ...
                                             A,posneg,eps,par)
% [out, projects] = analytic_basis(projection,x,preimage,s,p,A, ...
%                                   posneg,eps,Q,p_old_in)
%
% Returns an analytic basis at specified infinity using the method of Kato.
%
% Input "projection" is a function handle to the projection function to be
% used, "x" is the numerical value of infinity, "preimage" is the contour
% on which the Evans function is computed, "A" is a function handle to the
% Evans matrix, "posneg" is 1 or -1 determining which space the
% projection function should return, and "eps" is the tolerance in the
% projection function. 

iterations = size(preimage,2);
[p_old, Q1] = projection(A(x,preimage(1),par.o,par.r,par.p),posneg,eps);
[n,k]=size(Q1);

[U,T] = schur(A(x,preimage(1),par.o,par.r,par.p),'complex');
E = ordeig(T);
k = length(find(posneg*real(E)>eps));
US = ordschur(U,T,posneg*real(E)>eps);

out = zeros(n,k,iterations);
projects = zeros(size(p_old,1),size(p_old,2),iterations);

% if imag(preimage(1)) == 0
%     out(:,:,1) = real(US(:,1:k));
% else
    out(:,:,1) = US(:,1:k);
% end

for j=2:iterations
    [U,T] = schur(A(x,preimage(j),par.o,par.r,par.p),'complex');
    E = ordeig(T);
    k = length(find(posneg*real(E)>eps));
    US = ordschur(U,T,posneg*real(E)>eps);
    out(:,:,j) =  US(:,1:k);
end

