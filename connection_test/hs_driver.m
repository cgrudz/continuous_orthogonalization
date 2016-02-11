%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the test script for the projective unitary integration
% method on the Hocking-Stewartson Pulse of the 
% linear, complex Ginzburg-Landau equation.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Opening lines
clc; clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the problem parameters

% experiment title, will concat the center of the contour to the name
exp_name = 'PERK_test_';

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
R = -9.9;
xi_steps = 500;
xi_range = linspace(L,R,xi_steps);
delta_xi = (R-L)/xi_steps;

%define the spectral contour
lam = 0;
lam_steps = 1000;
preimage = lam + .5*exp(1i*linspace(0,2*pi,lam_steps+1));
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
%% Integrate the frame by PERK

% define the storage matrix for the path of frames
[state_dim,frame_dim,lam_steps] = size(init);
frame = zeros(xi_steps,state_dim,frame_dim,lam_steps);

for i = 1:lam_steps
    % define the initial condition for the lambda step
    temp = init(:,:,i);
    % define the stiefel derivative for the parameter value
    dQdx = @(x,Q)stiefel_derivative(x,Q,lam,omega,rho,psi);
    % integrate the initial frame and store
    frame(:,:,:,i) = rk4(dQdx,L,delta_xi,xi_steps,temp(:));  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate geometric phase of the unstable manifold

% Assuming a closed contour, reindex one step forward and one step backward
% in order to approximate the lambda derivative via the difference equation
% for each vector in the frame
shift_f = cat(4,frame(:,:,:,2:end),frame(:,:,:,1));
shift_b = cat(4,frame(:,:,:,end),frame(:,:,:,1:end-1));

% The lambda derivatives are approximated via
% \frac{X(\lambda(s_2),t) - X(\lambda(s_0),t)}{2\delta s}
Y_lam = (shift_f - shift_b)/(2*delta_s);

% We create storage for the connection at each point of the contour and
% each time step we wish to plot
con = connection_multi_dim(frame,Y_lam);

% compute the geometric phase and relative phase 
geo_phase = imag(squeeze(sum(con*delta_s,2)/(2*pi)));
rel_phase = geo_phase - geo_phase(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save data for post processing

data = struct('GP',geo_phase,'RP',rel_phase,'par',par,'K',preimage, ...
              'con',con, 'frame',frame,'plt_points',xi_range);

filename = strcat(exp_name,num2str(lam));
save(filename,'data')

gp_data = struct('GP',geo_phase,'RP',rel_phase);
filename_1 = strcat(filename,'_GP_RP');
save(filename_1,'gp_data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%