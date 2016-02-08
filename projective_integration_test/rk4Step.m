function [tNew,yNew] = rk4Step(fHandle,t,y,h);

% Runge-Kutta 4 butcher matrix
aButch = diag([1/2 1/2 1],-1);
bButch = [1/6 1/3 1/3 1/6];
cButch = [0 1/2 1/2 1]';

% First Runge Kutta stages
k = feval(fHandle,t,y);
sumKutta = h*bButch(1)*k;

% Other Runge Kutta stages
for stage = 2:4
  k = feval(fHandle, t+h*cButch(stage) , y+h*aButch(stage,stage-1)*k);
  sumKutta = sumKutta + h*bButch(stage)*k;
end

% One step forward with the solution
tNew = t+h;
yNew = y+sumKutta;
