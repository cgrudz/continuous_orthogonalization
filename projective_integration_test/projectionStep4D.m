function [x,Q] = projectionStep4D(x,Q)

Q = reshape(Q,4,2);

for i = 1:2
  P = Q'*Q-eye(2);
  Q = Q - 0.5*Q*P;
end

Q = Q(:);
