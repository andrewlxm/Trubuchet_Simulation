function dYdt = odeFctTrebuchetPart3Final(t,Y, m_c, m_p, m_a, l_c, l_as, l_al, l_s, l_acg, I_a, I_c)
% Y(1) : theta
% Y(2) : theta_dot = d theta / dt
% Y(3) : alpha
% Y(4) : alpha_dot = d alpha / dt
% Y(5) : beta
% Y(6) : beta_dot = d alpha / dt

g = 9.81;

J = m_c * l_as * l_c;
K = m_p * l_al * l_s;
L = g * (m_c*l_as - m_p*l_al - m_a * l_acg);
M = g * m_c * l_c;
N = g * m_p * l_s;

A =  m_p * l_al.^2 + m_c * l_as.^2  + I_a;
B = -J * cos( Y(1) + Y(3) );
C =  K * sin( Y(5) - Y(1) );
D =  J * Y(4).^2 * sin(Y(1)+Y(3)) + K * Y(6).^2 * cos(Y(5)-Y(1)) + L*sin(Y(1));

E =  m_c * l_c.^2 + I_c;
F =  J * Y(4).^2 * sin(Y(1)+Y(3)) + M*sin(Y(3));

G =  m_p * l_s.^2;
H = -K * Y(2).^2 * cos(Y(5)-Y(1)) + N * cos(Y(5));


DEN1 = -A*G*E + B.^2*G + C.^2*E;
DEN2 = E*(C.^2 - A*G) + B.^2*G;

dYdt = [ Y(2);
        ( -B*F*G - C*H*E + D*G*E ) / DEN1 - 0*10*Y(1) ;
        Y(4);
        ( A*F*G + B*C*H - B*D*G - C.^2*F ) / DEN2 ;
        Y(6);
        ( A*H*E - B.^2*H + B*C*F - C*D*E ) / DEN1 ;
];
end