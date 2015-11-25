function output = connection_multi_dim(frame,frame_lam)
    % This calculates the connection for the orthonormal frame of spanning
    % vectors for the un/stable bundle, as would be given by their wedge in
    % the exterior algebra.  The frame storage dimensions should be in the
    % form described below:
    [xi_steps, state_dim, frame_dim, lam_steps] = size(frame);
    output = zeros(xi_steps,lam_steps);
    M = (-1).^(1:frame_dim);
    M = M'*M;
    
    for i=1:lam_steps
        % For each point in the contour
        for j = 1:xi_steps
            % and each plot point, we will compute the connection via the
            % matrix equation
            % by the alternating sum of the elements of the matrix
            % at contour point i and plot point j, computed for each
            % row l of the frame_dim X frame_dim matrix         
            omega     = squeeze(frame(j,:,:,i));
            omegaH    = omega';
            omega_lam = squeeze(frame_lam(j,:,:,i));
            temp      = (omegaH*omega_lam).*M;
            output(j,i) = sum(temp(:));
        end
    end
end