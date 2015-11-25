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
% %% Define the contour of initial conditions
% 
% contour = zeros([2,lam_steps]);
% for i = 1:lam_steps
%     contour(:,i) = [1,sqrt(preimage(i)+1)];
% end

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
init = squeeze(init);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Integrate the contour

traj = zeros([xi_steps,2,lam_steps]);
for i = 1:lam_steps
    i/lam_steps
    [steps,path] = ode45(@(t,y) A_neutral(t,y,preimage(i)),...
                           xi_range,init(:,i));
              
    traj(:,:,i) = path;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the connection and plot phase

step_f = cat(3,traj(:,:,2:end),traj(:,:,1));
step_b = cat(3,traj(:,:,end),traj(:,:,1:end-1));

Y_lam = (step_f - step_b)/(2*delta_s);

connection = imag(sum(Y_lam.*conj(traj),2));
connection = connection./sum(traj.*conj(traj),2);

geo_phase = sum(connection,3)*delta_s;
geo_phase = squeeze(geo_phase/(2*pi));

plot(xi_range,geo_phase)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%