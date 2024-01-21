
function [t,x,y] = projectileAirDrag(m_p, c_p, gamma_0, v_0, x_release, y_release, y_gnd)

t = linspace(0,20,2001);
g = 9.81;
% c = 1;
% m = 1;
% gamma = 45 * pi /180;
% v_0 = 30;
v_t = m_p * g / c_p;           % terminal speed
x = x_release + v_0 * v_t/g * cos(gamma_0) * (1- exp(-g*t/v_t) );
y = y_release + v_t / g * (v_0*sin(gamma_0) + v_t) * (1- exp(-g*t/v_t) ) - v_t*t;

for k=1:length(y)
    if y(k) < y_gnd
        k_stop = k-1;
        break
    end
end
x = x(1:k_stop);
y = y(1:k_stop);
t = t(1:k_stop);
end