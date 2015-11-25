%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the test script for running the Drury continuous
% orthogonalization method on the Boussinesq equation.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Opening lines
clc; clear all; close all; beep off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the problem parameters

% wave parameter integration parameters

% endpoints
L = -10.0;
R = 10.0;
xi_steps = 2000;
xi_range = linspace(L,R,xi_steps);

%define the spectral contour
lam_steps = 500;
preimage = .1*exp(1i*linspace(0,2*pi,lam_steps+1));
preimage = preimage(1:end-1);
delta_s = (2*pi)/lam_steps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the k dimensional frame of initial conditions

% choose the basis un/stable subspace
pm = 1;

% choose the precision for the eigenvalue threshold
eps = .0001;

%define the function handle for the matrix equation
deriv = @A_bistable;
%define the function handle for the subspace-projection routine
proj = @projection_edit;

% define the inital frame for the given the projection funciton handle,
% the numerical infinity, the contour, matrix function handle, choice
% of un/stable subspace, and epsilon tolerance of the projection funciton
[init,projects]= analytic_basis_bistab(proj,L,preimage,deriv,pm,eps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Integrate the contour
[state_dim,frame_dim,lam_steps] = size(init);

Y = zeros([xi_steps,state_dim,lam_steps]);
for i = 1:lam_steps
    i/lam_steps
    [steps,path] = ode45(@(t,y) drury_edit(t,y,preimage(i),deriv, ...
                                           state_dim,frame_dim), ... 
                                           xi_range,init(:,i));

    Y(:,:,i) = path;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the connection and plot phase

step_f = cat(3,Y(:,:,2:end),Y(:,:,1));
step_b = cat(3,Y(:,:,end),Y(:,:,1:end-1));

Y_lam = (step_f - step_b)/(2*delta_s);

connection_M = connection_multi_dim(Y,Y_lam);
connection_S = sum(Y_lam.*conj(Y),2);

geo_phase_M = sum(connection_M,2)*delta_s;
geo_phase_M = imag(squeeze(geo_phase_M/(2*pi)));

geo_phase_S = sum(connection_S,3)*delta_s;
geo_phase_S = imag(squeeze(geo_phase_S/(2*pi)));

figure()
subplot(1,2,1)
plot(xi_range,geo_phase_M)
subplot(1,2,2)
plot(xi_range,geo_phase_S)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%