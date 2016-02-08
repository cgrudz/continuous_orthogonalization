%%%%%%%%%   REAL 2D *****************
%clear all; close all; clc;
%
%%% Eigenvalue
%n = 1;
%lambda = (n+0.5)^2;
%
%%% ODE Matrix
%A = [ 0 1; -lambda 0 ];
%
%%% Initial condition
%Q0 = [ 0; 1 ];
%
%%% ODE Right-hand side
%rhs = @(x,Q) qODE2D(x,Q,A);
%
%% Integrate
%numSteps = 1000;
%options.postProcess = 'projectionStep2D';
%sol = rk4(rhs,0,pi/numSteps,numSteps,Q0(:),options);
%Q = sol.y;
%X = sol.t;
%figure, plot(X,Q(:,1),'k'), hold on, plot(X,Q(:,2),'r'),hold off; 
%figure, plot(Q(:,1),Q(:,2)), axis square;
%
%%% Verify
%C  = [ 0 1 ];
%Q1 = reshape(Q(end,:),2,1);
%Q1
%det( C*Q1 )
%
%P = zeros(size(X));
%for n = 1:length(X)
%  q = reshape(Q(n,:),2,1);
%  P(n) = norm( conj(q')*q - 1 );
%end
%figure, plot(X,P,'-o')

%%%%%%%%%%    COMPLEX 2D *****************
%clear all; close all; clc;
%
%%% Eigenvalue
%%n = 1;
%%ambda = (n+0.5)^2, pause
%lambda = 0.5 + 1i*0.5, pause
%
%disp(' Results in C2 ');
%
%%% ODE Matrix
%A = [ 0 1; -lambda 0 ];
%
%%% Initial condition
%q0 = [ 0; 1 ];
%
%%% ODE Right-hand side
%rhs = @(x,q) qODE2D(x,q,A);
%
%% Integrate
%numSteps = 10000;
%options.postProcess = 'projectionStep2D';
%sol = rk4(rhs,0,pi/numSteps,numSteps,q0,options);
%q = sol.y;
%x = sol.t;
%%figure, plot(x,q(:,1),'k'), hold on, plot(x,q(:,2),'r'),hold off; 
%%figure, plot(q(:,1),q(:,2)), axis square;
%
%%% Verify
%q(end,:)
%
%P = zeros(size(x));
%for n = 1:length(x)
%  qi = reshape(q(n,:),2,1);
%  P(n) = conj(qi')*qi - 1;
%end
%figure, plot(x,real(P))
%figure, plot(x,imag(P))
%
%disp(' Results in R4 ');
%
%%% ODE Matrix
%A = [ 0 1; -real(lambda) 0 ];
%B = [ 0 0; -imag(lambda) 0 ];
%
%%% Initial condition
%q0 = [ 0; 1; 0; 0 ];
%
%%% ODE Right-hand side
%rhs = @(x,q) qODE2DComplex(x,q,A,B);
%
%% Integrate
%numSteps = 1000;
%options.postProcess = 'projectionStep2DComplex';
%sol = rk4(rhs,0,pi/numSteps,numSteps,q0,options);
%q = sol.y;
%x = sol.t;
%%figure, plot(x,q(:,1),'k'), hold on, plot(x,q(:,2),'r'),hold off; 
%%figure, plot(q(:,1),q(:,2)), axis square;
%
%%% Verify
%q(end,:)
%
%P = zeros(size(x));
%for n = 1:length(x)
%  qi = reshape(q(n,:),4,1);
%  P(n) = sum(qi.^2) - 1;
%end
%figure, plot(x,real(P))

%%%%%%%%%%   REAL 4D *****************
%%% Eigenvalue
%lambda = 4.7300407;
%
%%% ODE Matrix
%A = [ 0 1 0 0; 0 0 1 0; 0 0 0 1; lambda^4 0 0 0];
%
%%% Initial condition
%Q0 = [ 0 0; 0 0 ; 1 0; 0 1];
%
%%% ODE Right-hand side
%rhs = @(x,Q) qODE4D(x,Q,A);
%
%%% Integrate
%numSteps = 1000;
%options.postProcess = 'projectionStep4D';
%sol = rk4(rhs,0,1/numSteps,numSteps,Q0(:),options);
%Q = sol.y;
%X = sol.t;
%
%%% Verify
%C  = [ 1 0 0 0; 0 1 0 0 ];
%Q1 = reshape(Q(end,:),4,2);
%Q1
%det( C*Q1 )
%
%P = zeros(size(X));
%for n = 1:length(X)
%  QN = reshape(Q(n,:),4,2);
%  P(n) = norm( conj(QN')*QN - eye(2) , 'fro' );
%end
%plot(X,P)
%
%%%%%%%%%%   COMPLEX 4D *****************
%%% Eigenvalue
%lambda = 4.7300407;
%
%%% ODE Matrix
%A = [ 0 1 0 0; 0 0 1 0; 0 0 0 1; lambda^4 0 0 0];
%
%%% Initial condition
%Q0 = [ 0 0; 0 0 ; 1 0; 0 1];
%
%%% ODE Right-hand side
%rhs = @(x,Q) qODE4DComplex(x,Q,A);
%
%%% Integrate
%numSteps = 1000;
%options.postProcess = 'projectionStep4DComplex';
%sol = rk4(rhs,0,1/numSteps,numSteps,split(Q0(:)),options);
%Q = sol.y;
%X = sol.t;
%
%%% Verify
%C  = [ 1 0 0 0; 0 1 0 0 ];
%Q1 = reshape(gather(Q(end,:)'),4,2);
%Q1
%det( C*Q1 )
%
%P = zeros(size(X));
%for n = 1:length(X)
%  QN = reshape(gather(Q(n,:)'),4,2);
%  P(n) = norm( conj(QN')*QN - eye(2) , 'fro' );
%end
%figure; plot(X,P)
%
%%%%%%%%    COMPLEX 2D *****************
%clear all; close all; clc;
%
%% Eigenvalue
%n = 1;
%lambda = (n+0.5)^2;
%lambda = lambda*(1 + 1i);
%
%disp(' Results in C2 ');
%
%%% ODE Matrix
%A = [ 0 1; -lambda 0 ];
%
%%% Initial condition
%q0 = [ 0; 1 ];
%
%%% ODE Right-hand side
%rhs = @(x,q) qODE2D(x,q,A);
%
%% Integrate
%numSteps = 1000;
%options.postProcess = 'projectionStep2D';
%sol = rk4(rhs,0,pi/numSteps,numSteps,q0,options);
%q = sol.y;
%x = sol.t;
%%figure, plot(x,q(:,1),'k'), hold on, plot(x,q(:,2),'r'),hold off; 
%%figure, plot(q(:,1),q(:,2)), axis square;
%
%%% Verify
%q(end,:)
%
%P = zeros(size(x));
%for n = 1:length(x)
%  qi = reshape(q(n,:),2,1);
%  P(n) = conj(qi')*qi - 1;
%end
%figure, plot(x,real(P))
%figure, plot(x,imag(P))
%
%disp(' Results in R4 ');
%
%%% ODE Matrix
%A = [ 0 1; -real(lambda) 0 ];
%B = [ 0 0; -imag(lambda) 0 ];
%
%%% Initial condition
%q0 = [ 0; 1; 0; 0 ];
%
%%% ODE Right-hand side
%rhs = @(x,q) qODE2DComplex(x,q,A,B);
%
%% Integrate
%numSteps = 1000;
%options.postProcess = 'projectionStep2DComplex';
%sol = rk4(rhs,0,pi/numSteps,numSteps,q0,options);
%q = sol.y;
%x = sol.t;
%%figure, plot(x,q(:,1),'k'), hold on, plot(x,q(:,2),'r'),hold off; 
%%figure, plot(q(:,1),q(:,2)), axis square;
%
%%% Verify
%q(end,:)
%
%P = zeros(size(x));
%for n = 1:length(x)
%  qi = reshape(q(n,:),4,1);
%  P(n) = sum(qi.^2) - 1;
%end
%figure, plot(x,real(P))
%
%disp(' Results in C2 With ODE45');
%
%%% ODE Matrix
%A = [ 0 1; -lambda 0 ];
%
%%% Initial condition
%q0 = [ 0; 1 ];
%
%%% ODE Right-hand side
%rhs = @(x,q) qODE2D(x,q,A);
%
%% Integrate
%sol = ode45(rhs,[0 pi],q0);
%q = (sol.y)';
%x = (sol.x)';
%
%%% Verify
%q(end,:)
%
%P = zeros(size(x));
%for n = 1:length(x)
%  qi = reshape(q(n,:),2,1);
%  P(n) = conj(qi')*qi - 1;
%end
%figure, plot(x,real(P))
%figure, plot(x,imag(P))
%
%disp(' Results in R4 with ODE45');
%
%%% ODE Matrix
%A = [ 0 1; -real(lambda) 0 ];
%B = [ 0 0; -imag(lambda) 0 ];
%
%%% Initial condition
%q0 = [ 0; 1; 0; 0 ];
%
%%% ODE Right-hand side
%rhs = @(x,q) qODE2DComplex(x,q,A,B);
%
%% Integrate
%sol = ode45(rhs,[0,pi],q0);
%q = (sol.y)';
%x = (sol.x)';
%
%%% Verify
%q(end,:)
%
%P = zeros(size(x));
%for n = 1:length(x)
%  qi = reshape(q(n,:),4,1);
%  P(n) = sum(qi.^2) - 1;
%end
%figure, plot(x,real(P))



%%%%%%%%    COMPLEX 4D ALT *****************
%clear all; close all; clc;
%disp('Results COMPLEX 4D ALT, real lambda');
%lambda = 4.7300407;
%r = 0.0;
%theta = 0;
%a = realPartLambda4(real(lambda),imag(lambda),r,theta);
%b = imagPartLambda4(real(lambda),imag(lambda),r,theta);
%
%%% ODE Matrix
%A = [ 0 1 0 0; 0 0 1 0; 0 0 0 1; lambda^4 0 0 0];
%AR = [ 0 1 0 0; 0 0 1 0; 0 0 0 1; a 0 0 0];
%AI = [ 0 0 0 0; 0 0 0 0; 0 0 0 0; b 0 0 0];
%
%% Initial condition
%Q0 = [ 0 0; 0 0 ; 1 0; 0 1];
%
%%% ODE Right-hand side
%rhs = @(x,Q) qODE4DComplexAlt(x,Q,AR,AI);
%
%%% Integrate
%numSteps = 10000;
%options.postProcess = 'projectionStep4DComplexAlt';
%sol = rk4(rhs,0,1/numSteps,numSteps,split(Q0(:)),options);
%Q = sol.y;
%X = sol.t;
%
%%% Verify
%C  = [ 1 0 0 0; 0 1 0 0 ];
%Q1 = reshape(gather(Q(end,:)'),4,2);
%Q1
%det( C*Q1 )
%
%%P = zeros(size(X));
%%for n = 1:length(X)
%%  QN = reshape(gather(Q(n,:)'),4,2);
%%  P(n) = norm( conj(QN')*QN - eye(2) , 'fro' );
%%end
%%figure; plot(X,P); drawnow;
%
%%% Verify
%Q = sol.y;
%P = zeros(size(X));
%R = zeros(size(X));
%for n = 1:length(X)
%  U  = reshape(Q(n,1:8),4,2);
%  V  = reshape(Q(n,9:16),4,2);
%  %QN = U + i*V;
%  %P(n) = norm( conj(QN')*QN - eye(2) , 'fro' );
%  P(n) = norm( U'*U + V'*V - eye(2) , 'fro' );
%  R(n) = norm( U'*V - V'*U , 'fro' );
%end
%figure; plot(X,P); drawnow;
%figure; plot(X,R); drawnow;
%
%disp('Results COMPLEX 4D ALT, complex lambda');
%lambda = 4.7300407;
%r = 0.2;
%theta = 1/16*pi;
%a = realPartLambda4(real(lambda),imag(lambda),r,theta);
%b = imagPartLambda4(real(lambda),imag(lambda),r,theta);
%
%%% ODE Matrix
%A = [ 0 1 0 0; 0 0 1 0; 0 0 0 1; lambda^4 0 0 0];
%AR = [ 0 1 0 0; 0 0 1 0; 0 0 0 1; a 0 0 0];
%AI = [ 0 0 0 0; 0 0 0 0; 0 0 0 0; b 0 0 0];
%
%% Initial condition
%Q0 = [ 0 0; 0 0 ; 1 0; 0 1];
%
%%% ODE Right-hand side
%rhs = @(x,Q) qODE4DComplexAlt(x,Q,AR,AI);
%
%%% Integrate
%numSteps = 10000;
%options.postProcess = 'projectionStep4DComplexAlt';
%sol = rk4(rhs,0,1/numSteps,numSteps,split(Q0(:)),options);
%Q = sol.y;
%X = sol.t;
%
%%% Verify
%P = zeros(size(X));
%R = zeros(size(X));
%for n = 1:length(X)
%  U  = reshape(Q(n,1:8),4,2);
%  V  = reshape(Q(n,9:16),4,2);
%  %QN = U + i*V;
%  %P(n) = norm( conj(QN')*QN - eye(2) , 'fro' );
%  P(n) = norm( U'*U + V'*V - eye(2) , 'fro' );
%  R(n) = norm( U'*V - V'*U , 'fro' );
%end
%figure; plot(X,P); drawnow;
%figure; plot(X,R); drawnow;

%%%%%%%%%  Orr Sommerfeld *****************
%% Eigenvalue

alpha  = 1.020547;
R      = 5772.2218;
c0     = 0.264003;
lambda = -i*alpha*c0;

%% Initial condition
Q0 = [ 0 0; 0 0 ; 1 0; 0 1];

%% ODE Right-hand side
rhs = @(x,Q) qOrrSommerfeld(x,Q,alpha,R,lambda);

%% Integrate
numSteps = 2000000;
options.postProcess  = 'projectionStep4D';
options.saveSolution = false;
sol = rk4(rhs,-1,2/numSteps,numSteps,Q0(:),options);
Q = sol.y;
X = sol.t;

%% Verify
C  = [ 1 0 0 0; 0 1 0 0 ];
Q1 = reshape(Q,4,2);
Q1
det( C*Q1 )

%%% Verify
%C  = [ 1 0 0 0; 0 1 0 0 ];
%Q1 = reshape(Q(end,:),4,2);
%Q1
%det( C*Q1 )

%P = zeros(size(X));
%for n = 1:length(X)
%  QN = reshape(Q(n,:),4,2);
%  P(n) = norm( conj(QN')*QN - eye(2) , 'fro' );
%end
%plot(X,P)
