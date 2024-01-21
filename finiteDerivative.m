

function [y_dot] = finiteDerivative(t,y)
    N = length(y);
    y_dot = zeros(N,1);
    for k=1:N-1
        y_dot(k+1) = ( y(k+1)-y(k) ) / ( t(k+1) - t(k) );
    end
end