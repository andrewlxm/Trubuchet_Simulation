function dYdt = odeFctTrebuchetPart1Final(t,Y, m_c, m_p, m_a, l_c, l_as, l_al, l_acg, I_a, I_c)
% Y(1) : theta
% Y(2) : theta_dot = d theta / dt
% Y(3) : alpha
% Y(4) : alpha_dot = d alpha / dt

g = 9.81;

F = m_c * l_as * l_c;
G = g * (m_c*l_as - m_p*l_al - m_a * l_acg);
H = g * m_c * l_c;


A =  m_p * l_al^2 + m_c * l_as^2  + I_a;
B = -F * cos( Y(1) + Y(3) );
C =  F * Y(4).^2 * sin(Y(1)+Y(3)) + G*sin(Y(1));
D =  m_c * l_c.^2 + I_c;
E =  F * Y(2).^2 * sin(Y(1)+Y(3)) + H*sin(Y(3));

DEN = D*A - B^2;

dYdt = [ Y(2);
        ( B*E - D*C ) /DEN ;
        Y(4);
        ( B*C - A*E ) /DEN ;
];
end