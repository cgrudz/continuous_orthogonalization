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
%% Integrate the frame by two methods for comparison
%  we integrate via the Drury method and via standard ODE45

% define the storage matrix for the path of frames
[state_dim,frame_dim,lam_steps] = size(init);
drury_frame = zeros(state_dim,frame_dim,xi_steps,lam_steps);

rk_frame = zeros(state_dim,frame_dim,xi_steps,lam_steps);

for i = 1:lam_steps
    i
%      sigma = sqrt(preimage(i)+1);
%      vec = [1;sigma];
%      vec_o = vec/sqrt(vec'*vec);
%      [rk_steps, rk_path] = ode45(@(t,y) bistable_deriv(t,y,...
%                                    deriv,preimage(i)),xi_range,vec_o);
    [rk_steps, rk_path] = ode45(@(t,y) bistable_deriv(t,y,deriv,...
                                     preimage(i)),xi_range,init(:,:,i));
    rk_frame(:,:,:,i) = rk_path';
    
%     [drury_steps, drury_path] = ode45(@(t,y) drury_edit(t,y,...
%                                     preimage(i),deriv,state_dim, ...
%                                      frame_dim), xi_range,vec_o);
    [drury_steps, drury_path] = ode45(@(t,y) drury_edit(t,y,...
                                       preimage(i),deriv,state_dim, ...
                                       frame_dim), xi_range,init(:,:,i));
    drury_frame(:,:,:,i) = drury_path';
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate geometric phase of the un/stable manifold at plotting steps

%All safe and debugged from here

% Assuming a closed contour, reindex one step forward and one step backward
% in order to approximate the lambda derivative via the difference equation
% for each vector in the frame
shift_f_1 = cat(4,drury_frame(:,:,:,2:end),drury_frame(:,:,:,1));
shift_b_1 = cat(4,drury_frame(:,:,:,end),drury_frame(:,:,:,1:end-1));


shift_f_2 = cat(4,rk_frame(:,:,:,2:end),rk_frame(:,:,:,1));
shift_b_2 = cat(4,rk_frame(:,:,:,end),rk_frame(:,:,:,1:end-1));

% The lambda derivatives are approximated via
% \frac{X(\lambda(s_2),t) - X(\lambda(s_0),t)}{2\delta s}
D_paths_lam = (shift_f_1 - shift_b_1)/(2*delta_s);
rk_paths_lam = (shift_f_2 - shift_b_2)/(2*delta_s);

% We create storage for the connection at each point of the contour and
% each time step we wish to plot
connection_1 = connection_multi_dim(drury_frame,D_paths_lam);
connection_2 = connection(drury_frame,D_paths_lam);
%connection_2 = connection(rk_frame,rk_paths_lam);

phase_1 = squeeze(sum(connection_1*delta_s,2)/(2*pi));
phase_2 = squeeze(sum(connection_2*delta_s,2)/(2*pi));

figure(1)
plot(xi_range,phase_1)
figure(2)
plot(xi_range,phase_2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%