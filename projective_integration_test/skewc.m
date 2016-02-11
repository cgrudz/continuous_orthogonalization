function M = skewc(F)
    % compute the skew symmetric part of the matrix F
    M = 0.5*( F - F' );
end