function output = connection(Y,Y_prime)
    % state dimension is the first index so we collapse along this
    % dimension, following the form of the connection as given by
    % \frac{Imag(<Y_\lambda,Y>}{<Y,Y>}
    output = squeeze(imag(sum(conj(Y).*Y_prime,1)));
    output = output./squeeze(sum(conj(Y).*Y,1));
end