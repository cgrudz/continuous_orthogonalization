function yOut = rk4(fHandle,t0,h,numSteps,y0)
    % This is an edit of Daniele Avitabile's code for the test run of the
    % projective integration step for the numerical connection on the frame

    %% Initialising
    t = t0;
    y = y0;

    %% Allocating
    yOut = zeros(numSteps,4,2); 
    yOut(1,:,:) = reshape(y0,4,2);
    %% Time loop
    for it = 2:numSteps
        % RK step
        [t,y] = rk4Step(fHandle,t,y,h);
        % Euler projection step
        [t,y] = projectionStep4D(t,y);
        
        yOut(it,:,:) = reshape(y,4,2);
    end
    
end
