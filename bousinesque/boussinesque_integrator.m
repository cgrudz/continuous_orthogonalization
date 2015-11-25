%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the test script for running the Drury continuous
% orthogonalization method on the Boussinesq equation.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Opening lines
clc; clear all; close all; beep off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the problem parameters

% Boussinesq parameters
p = 0.4;

% wave parameter integration parameters

% endpoints
L = -10.0;
R = 10.0;
% timestep of phase plot
h = 2;
% timesteps to return the geometric phase
tsteps = linspace(L,R,((R-L)/h));

%define the spectral contour
lam_steps = 3200;
preimage = 10+0.5*exp(1i*linspace(0,2*pi,lam_steps));
%delta_s = (2*pi)/lam_steps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the k dimensional frame of initial conditions

% choose the basis un/stable subspace
pm = 1;

% choose the precision for the eigenvalue threshold
eps = .0001;

% define the function handle for the matrix equation
deriv = @A;
%define the function handle for the subspace-projection routine
proj = @projection_edit;

% define the inital frame for the given the projection funciton handle,
% the numerical infinity, the contour, matrix function handle, choice
% of un/stable subspace, and epsilon tolerance of the projection funciton
[init,projects]= analytic_basis_edit(proj,L,preimage,deriv,pm,eps,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Implement the Drury Continuous Orthogonalization integral of the frame

% define the storage matrix for the path of frames
[state_dim,frame_dim,lam_steps] = size(init);
frame_paths = zeros(state_dim,frame_dim,lam_steps,length(tsteps));
frame_paths(:,:,:,1) = init;
for i = 1:lam_steps
    i
    for j = 1:length(tsteps)-1
        %integrate the frame and store the values for each point tstep
        %where the phase will be ploted at
        temp_struc =ode45(@(t,y) ... 
                          drury_edit(t,y,preimage(i),deriv,state_dim,frame_dim,p), ...
                          [tsteps(j),tsteps(j+1)],frame_paths(:,:,i,j));
        frame_plot_point = reshape(temp_struc.y(:,end),state_dim,frame_dim);
        frame_paths(:,:,i,j+1) = frame_plot_point;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate geometric phase of the un/stable manifold at plotting steps

% Assuming a closed contour, reindex one step forward and one step backward
% in order to approximate the lambda derivative via the difference equation
% for each vector in the frame
shift_f = cat(3,frame_paths(:,:,2:end,:),frame_paths(:,:,1,:));
shift_b = cat(3,frame_paths(:,:,end,:),frame_paths(:,:,1:end-1,:));

% The lambda derivatives are approximated via
% \frac{X(\lambda(s_2),t) - X(\lambda(s_0),t)}{2\delta s}
frame_paths_lam = (shift_f - shift_b)*(lam_steps/4*pi);

% We create storage for the connection at each point of the contour and
% each time step we wish to plot
connection = zeros(lam_steps,length(tsteps));

for i=1:lam_steps
    % For each point in the contour
    for j = 1:length(tsteps)
        % and each plot point, we will compute the connection via the
        % matrix equation
        temp = frame_paths(:,:,i,j)'*frame_paths_lam(:,:,i,j);
        for l = 1:frame_dim
            % by the alternating sum of the elements of the matrix
            % at contour point i and plot point j, computed for each
            % row l of the frame_dim X frame_dim matrix
            connection(i,j) = (connection(i,j) + ...
                              (-1)^(l+1)*sum(temp(l,1:2:end)) + ...
                              (-1)^(l)*sum(temp(l,2:2:end)));
        end
    end
end
phase = sum(connection,1)/(1i*lam_steps);
plot(tsteps,phase)
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%