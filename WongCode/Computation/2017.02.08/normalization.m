    m_pi = 0.13957;
    eta_jet = 0 ; 
    s_nn = 200;    
    m_N = 0.93957;
    y_N = acosh(s_nn/(2*m_N));
%-------------------------        
% a set of parameters
    pT_trig = 1.5;    
    N_jet = 0.15 + 0.1 * pT_trig;     
    T_jet = 0.19 + 0.06 * pT_trig;
    sigma0 = 0.5;
    m_a = 1.1;
    q = 0.8;   
    frNk = 3;
    f_J = 0.82;
    a = 0.5;
    T = 0.5;
    m_d = 1;    
    
    A_ridge = 1 / normalisation(T, m_pi, m_d, y_N, a);
    disp(A_ridge);
    
function [norm] = normalisation(T, m_pi, m_d, y_N, a)
    
%   fx = @(d_pit, d_yi) (1 - sqrt(m_pi^2 + d_pit.^2) .* exp(abs(d_yi) - y_N) /m_pi).^a ;
    fA = @(d_pit, d_yi) d_pit.* (1 - sqrt(m_pi^2 + d_pit.^2) .* exp(abs(d_yi) - y_N) /m_pi).^a .* exp(-sqrt(m_pi^2 + d_pit.^2)./T)./(sqrt(m_d^2 + d_pit.^2)) ;
    
    pit_min = 0;
    pit_max = 10;
    yi_min = @(d_pit) -(log(m_pi./sqrt(m_pi^2 + d_pit.^2)) + y_N);
    yi_max = @(d_pit) (log(m_pi./sqrt(m_pi^2 + d_pit.^2)) + y_N);
    
    norm = 2*pi * integral2(fA, pit_min, pit_max, yi_min, yi_max);    
    
    end