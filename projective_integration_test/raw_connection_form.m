function output = raw_connection_form(frame,frame_lam)
    % This calculates the connection for the orthonormal frame of spanning
    % vectors for the un/stable bundle, as would be given by their wedge in
    % the exterior algebra.  The frame storage dimensions should be in the
    % form described below:
    [xi_steps, state_dim, frame_dim, lam_steps] = size(frame);
    output = zeros(xi_steps,lam_steps);
    
    for i=1:lam_steps
        % For each point in the contour
        for j = 1:xi_steps
            % and each plot point, we will compute the connection via the
            % matrix equation
            % by the alternating sum of the elements of the matrix
            % at contour point i and plot point j, computed for each
            % row l of the frame_dim X frame_dim matrix         
            o1     = squeeze(frame(j,:,1,i));
            o2     = squeeze(frame(j,:,2,i));
            o1_l = squeeze(frame_lam(j,:,1,i));
            o2_l = squeeze(frame_lam(j,:,2,i));
            temp = (o1_l*o1')*(o2*o2') + (o2_l*o2')*(o1*o1');
            temp = temp - (o1_l*o2')*(o2*o1') -(o1*o2')*(o2_l*o1');
            output(j,i) = temp;
        end
    end
end