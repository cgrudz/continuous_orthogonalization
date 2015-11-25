%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the test script for running the Drury continuous
% orthogonalization method on the Hocking-Stewartson Pulse of the 
% linear, complex Ginzburg-Landau equation.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Opening lines
clc; clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the problem parameters

% parameters for the equation
omega = 3;
rho = 1/sqrt(5);
psi = atan(2);

par.o = omega;
par.r = rho;
par.p = psi;

% wave parameter integration parameters

% endpoints
L = -10.0;
R = 10.0;
xi_steps = 1000;
xi_range = linspace(L,R,xi_steps);

%define the spectral contour
lam_steps = 2000;
preimage = 0 + exp(1i*linspace(0,2*pi,lam_steps+1));
preimage = preimage(1:end-1);
delta_s = (2*pi)/lam_steps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the k dimensional frame of initial conditions

% % choose the basis un/stable subspace
pm = 1;
 
% % choose the precision for the eigenvalue threshold
eps = .0001;

% define the function handle for the matrix equation
deriv = @hs_matc4;
% %define the function handle for the subspace-projection routine
proj = @projection_edit;
 
% define the inital frame for the given the projection funciton handle,
% the numerical infinity, the contour, matrix function handle, choice
% of un/stable subspace, and epsilon tolerance of the projection funciton
[init,projects]= analytic_basis_HS(proj,L,preimage,deriv,pm,eps,par);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Integrate the frame by two methods for comparison
%  we integrate via the Drury method and via standard ODE45

% define the storage matrix for the path of frames
[state_dim,frame_dim,lam_steps] = size(init);
frame = zeros(xi_steps,state_dim,frame_dim,lam_steps);

for i = 1:lam_steps
    i/lam_steps
    [steps, path] = ode15s(@(t,y) drury_edit(t,y,...
                                     preimage(i),deriv,state_dim, ...
                                     frame_dim,par), xi_range,init(:,:,i));
    
    frame(:,:,:,i) =reshape(path,xi_steps,state_dim,frame_dim);   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate geometric phase of the un/stable manifold at plotting steps

% Gram-Schmidt
nrm = sqrt(sum(conj(frame).*frame,2));
nrm = repmat(nrm,[1,state_dim,1,1]);
Y = frame./nrm;
Y_temp = sum(conj(Y(:,:,1,:)).*Y(:,:,2,:),2);
Y_temp = repmat(Y_temp,[1,state_dim,1,1]);
Y_temp = Y(:,:,2,:) - Y(:,:,1,:).*Y_temp;
Y(:,:,2,:) = Y_temp;

% Y = frame;

% Assuming a closed contour, reindex one step forward and one step backward
% in order to approximate the lambda derivative via the difference equation
% for each vector in the frame
shift_f = cat(4,Y(:,:,:,2:end),Y(:,:,:,1));
shift_b = cat(4,Y(:,:,:,end),Y(:,:,:,1:end-1));

% The lambda derivatives are approximated via
% \frac{X(\lambda(s_2),t) - X(\lambda(s_0),t)}{2\delta s}
Y_lam = (shift_f - shift_b)/(2*delta_s);

% We create storage for the connection at each point of the contour and
% each time step we wish to plot
con = connection_multi_dim(Y,Y_lam);
%con = raw_connection_form(Y,Y_lam);

geo_phase = imag(squeeze(sum(con*delta_s,2)/(2*pi)));

figure(1)
plot(xi_range,geo_phase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%