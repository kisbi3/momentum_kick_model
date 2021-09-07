% Drawing the ridge and jet structure
% constants
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
%-------------------------    
    A_ridge = 1 / normalisation(T, m_pi, m_d, y_N, a);
    
    phi = -2:0.01:2;
    deta = -1:0.01:1;
    pft = 2:0.01:3;
    
    [PHI, DETA, PFT] = meshgrid(phi, deta, pft);
    
    eta = DETA + eta_jet;

    pf1 = PFT .* cos(PHI);
    pf2 = PFT .* sin(PHI);
    pf3 = PFT .* sinh(eta);
    
    pi1 = pf1 - q/cosh(eta_jet);
    pi2 = pf2;
    pi3 = pf3 - q*sinh(eta_jet)/cosh(eta_jet);
    pit = sqrt(pi1.^2 + pi2.^2);
    
    Ef = sqrt(pf1.^2 + pf2.^2 + pf3.^2 + m_pi^2);
    Ei = sqrt(pi1.^2 + pi2.^2 + pi3.^2 + m_pi^2);
    
    yf = log((Ef + pf3)./(Ef - pf3))./2;
    yi = log((Ei + pi3)./(Ei - pi3))./2;
    
    mtf = sqrt(m_pi^2 + pf1.^2 + pf2.^2);
    mti = sqrt(m_pi^2 + pi1.^2 + pi2.^2);
    
    x = sqrt(m_pi^2 + pit.^2) .* exp(abs(yi) - y_N) ./m_pi ;
    sigma_pi = sigma0 * m_a ./(sqrt(m_a^2 + PFT.^2));
    
%   etaf = log((sqrt(pf1.^2 + pf2.^2 + pf3.^2) + pf3)./((sqrt(pf1.^2 + pf2.^2 + pf3.^2) - pf3)));
%   etai = log((sqrt(pi1.^2 + pi2.^2 + pi3.^2) + pi3)./((sqrt(pi1.^2 + pi2.^2 + pi3.^2) - pi3)));
   
    Nridge =  frNk * 2/3 * A_ridge * (1 - x).^a .* exp(-mti./T)./(sqrt(m_d^2 + pit.^2)) .*  Ef./Ei .* sqrt(1 - m_pi^2./(mtf.^2 .* (cosh(yf)).^2));
    
    Njet = f_J * N_jet .* exp((m_pi - mtf) ./ T_jet) ./(T_jet .* (m_pi + T_jet)) .* exp(-(PHI.^2 + DETA.^2) ./ (2 * sigma_pi.^2)) ./ (2 * pi * sigma_pi.^2);
    
    N = Nridge + Njet;
    
    Q = trapz(pft, N, 3);
    [gPHI, gDETA] = meshgrid(phi, deta);
    mesh(gPHI, gDETA, Q);
   
function [norm] = normalisation(T, m_pi, m_d, y_N, a)
    
%   fx = @(d_pit, d_yi) (1 - sqrt(m_pi^2 + d_pit.^2) .* exp(abs(d_yi) - y_N) /m_pi).^a ;
    fA = @(d_pit, d_yi) d_pit.* (1 - sqrt(m_pi^2 + d_pit.^2) .* exp(abs(d_yi) - y_N) /m_pi).^a .* exp(-sqrt(m_pi^2 + d_pit.^2)./T)./(sqrt(m_d^2 + d_pit.^2)) ;
    
    pit_min = 0;
    pit_max = 10;
    yi_min = @(d_pit) -(log(m_pi./sqrt(m_pi^2 + d_pit.^2)) + y_N);
    yi_max = @(d_pit) (log(m_pi./sqrt(m_pi^2 + d_pit.^2)) + y_N);
    
    norm = 2*pi * integral2(fA, pit_min, pit_max, yi_min, yi_max);    
    
    end