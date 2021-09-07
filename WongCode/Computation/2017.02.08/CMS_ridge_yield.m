% Calculating ridge yield
% constants
    m_pi = 0.13957;
    eta_jet = 0 ; 
    s_nn = 200;    
    m_N = 0.93957;
    y_N = acosh(s_nn/(2*m_N));
%-------------------------        
% a set of parameters
    sigma0 = 0.5;
    m_a = 1.1;
    q = 0.7;   
    f_J = 0.82;
    a = 0.5;
    m_d = 1;    
%-------------------------   
for i = 1 : 5
    xx = [-1.1781, -0.9817, -0.7854, -0.589, -0.3927, -0.1963, 0, 0.1963, 0.3927, 0.589, 0.7854, 0.9817, 1.1781];

    switch i
        case 1
            C_ZYAM = 0;
            pT_trig = 0.65;
            ry = [0.0069, 0.0216, 0.0405, 0.0612, 0.0802, 0.0931, 0.0983, 0.0931, 0.0802, 0.0612, 0.0405, 0.0216, 0.0069];
        case 2
            C_ZYAM = 0;
            pT_trig = 1.5;
            ry = [0.0134, 0.0432, 0.0812, 0.1253, 0.1686, 0.199, 0.2091, 0.199, 0.1686, 0.1253, 0.0812, 0.0432, 0.0134];
        case 3
            C_ZYAM = 0;
            pT_trig = 3;
            ry = [0.013, 0.0462, 0.0943, 0.1492, 0.2059, 0.2428, 0.2542, 0.2428, 0.2059, 0.1492, 0.0943, 0.0462, 0.013];
        case 4
            C_ZYAM = 0;
            pT_trig = 5;
            ry = [0.0091, 0.0258, 0.0808, 0.1221, 0.1592, 0.177, 0.204, 0.177, 0.1592, 0.1221, 0.0808, 0.0258, 0.0091];
        case 5
            C_ZYAM = 0;
            pT_trig = 9;
            ry = [-0.0101, 0.0073, 0.0362, 0.0451, 0.0612, 0.1078, 0.113, 0.1078, 0.0612, 0.0451, 0.0362, 0.0073, -0.0101];
    end
    subplot(1,5,i);
    
% a set of variables   
    N_jet = 0.6017 + 0.1704 * pT_trig;     
    T_jet = 0.228 + 0.072 * pT_trig;
    fr = exp(-1.137/pT_trig);
    Nk = 10.69652482 * exp(-0.1692 * pT_trig);
    %7.937
    frNk = fr * Nk;
    T = 0.6;
    A_ridge = 1 / normalisation(T, m_pi, m_d, y_N, a);
%-------------------------   
    phi = -pi:0.01:pi;
    deta = -1:0.01:1;
    pft = 1:0.01:2;
    
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
    
%   Ptotal = trapz(pft, N, 3);
%   Qtotal = trapz(deta, Ptotal, 1) - C_ZYAM; %Total
    
    Pridge = trapz(pft, Nridge, 3);
    Qridge = trapz(deta, Pridge, 1) - C_ZYAM; %Ridge
    
    hold on    
    plot(phi, Qridge, 'r -', 'linewidth', 2)
    plot(xx, ry,'*','marker','s','markerfacecolor','b', 'markeredgecolor','b');    
    xlim([-pi,pi]);
    ylim([0,1]);
    
end   

function [norm] = normalisation(T, m_pi, m_d, y_N, a)
    
%   fx = @(d_pit, d_yi) (1 - sqrt(m_pi^2 + d_pit.^2) .* exp(abs(d_yi) - y_N) /m_pi).^a ;
    fA = @(d_pit, d_yi) d_pit.* (1 - sqrt(m_pi^2 + d_pit.^2) .* exp(abs(d_yi) - y_N) /m_pi).^a .* exp(-sqrt(m_pi^2 + d_pit.^2)./T)./(sqrt(m_d^2 + d_pit.^2)) ;
    
    pit_min = 0;
    pit_max = 10;
    yi_min = @(d_pit) -(log(m_pi./sqrt(m_pi^2 + d_pit.^2)) + y_N);
    yi_max = @(d_pit) (log(m_pi./sqrt(m_pi^2 + d_pit.^2)) + y_N);
    
    norm = 2*pi * integral2(fA, pit_min, pit_max, yi_min, yi_max);    
    
end