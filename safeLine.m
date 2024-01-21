

function [x_down, y_down, x_up, y_up] = safeLine(x,y, safeD)
N = length(x);
x_down = zeros(N,1);
y_down = zeros(N,1);
x_up   = zeros(N,1);
y_up   = zeros(N,1);

for k=1:N
    if k < N
        slope     = (y(k+1) - y(k)) / (x(k+1) - x(k));
    end
    if abs(slope) < 1e-4
        x_down(k) =  x(k) ;
        y_down(k) =  y(k) -safeD;

        x_up(k)   =  x(k) ;
        y_up(k)   =  y(k) +safeD;
    else if abs(slope) > 1e8
        x_down(k) =  x(k) -safeD;
        y_down(k) =  y(k) ;

        x_up(k)   =  x(k) +safeD;
        y_up(k)   =  y(k) ;
    else
        a         = -1/slope;
        x_down(k) = ( -sqrt(a.^2+1) * safeD + a.^2 * x(k) + x(k) ) / ( a.^2 + 1 );
        y_down(k) = ( - a * sqrt(a.^2+1) * safeD + a.^2 * y(k) + y(k) ) / ( a.^2 + 1 );

        x_up(k)   = (  sqrt(a.^2+1) * safeD + a.^2 * x(k) + x(k) ) / ( a.^2 + 1 );
        y_up(k)   = (   a * sqrt(a.^2+1) * safeD + a.^2 * y(k) + y(k) ) / ( a.^2 + 1 );
    end
end

end