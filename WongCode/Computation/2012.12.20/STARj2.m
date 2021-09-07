
function STARj2
    

    q = 1;
    mpi = 0.13957;
    TMP = 0.5;
    etajet = 0 ;
    Aridge = 1 ;
    fRNK = 3;
    fJ = 0.632;
    ma = 1.1;
    md = 1;
    a = 0.5;
    sigmapizero = 0.5;
    snn = 200;
    mN = 0.93957;
    
    
    yb = acosh(snn/(2*mN));
    
    
    for phi = -1:0.01:1.5
        
    inte2 = 0;    
        
    for pt = 2:0.005:4
    
    deta = -1:0.005:1;
    
    eta = deta + etajet;

    pf1 = pt .* cos(phi);
    pf2 = pt .* sin(phi);
    pf3 = pt .* sinh(eta);
    
    pi1 = pf1 - q/cosh(etajet);
    pi2 = pf2;
    pi3 = pf3 - q*sinh(etajet)/cosh(etajet);
    
    Ef = sqrt(pf1.^2 + pf2.^2 + pf3.^2 + mpi^2);
    Ei = sqrt(pi1.^2 + pi2.^2 + pi3.^2 + mpi^2);
    
    yf = log((Ef + pf3)./(Ef - pf3))./2;
    yi = log((Ei + pi3)./(Ei - pi3))./2;
    
    mtf = sqrt(mpi^2 + pf1.^2 + pf2.^2);
    mti = sqrt(mpi^2 + pi1.^2 + pi2.^2);
  
    pit = sqrt(pi1.^2 + pi2.^2);
    
    x = sqrt(mpi^2 + pit.^2) .* exp(abs(yi) - yb) ./mpi ;
    
    etaf = log((sqrt(pf1.^2 + pf2.^2 + pf3.^2) + pf3)./((sqrt(pf1.^2 + pf2.^2 + pf3.^2) - pf3)));
    etai = log((sqrt(pi1.^2 + pi2.^2 + pi3.^2) + pi3)./((sqrt(pi1.^2 + pi2.^2 + pi3.^2) - pi3)));

    NJF = 0.75;
    TJF = 0.55;
    
    sigmapi = sigmapizero * ma ./(sqrt(ma^2 + pt.^2));
    
    
    Nridge =  fRNK * 2/3 * Aridge * (1 - x).^a .* exp(-sqrt(mpi^2 + pit.^2)/TMP)./(sqrt(md^2 + pit.^2)) .*  Ef./Ei .* sqrt(1 - mpi^2./(mtf.^2 .* (cosh(yf)).^2));
    
    Njet = fJ * NJF .* (exp((mpi - sqrt(mpi^2 + pt.^2)) ./ TJF) ./(TJF .* (mpi + TJF)))  .* exp(-(phi^2 + deta.^2) ./ (2 * sigmapi.^2)) ./ (2 * pi * sigmapi.^2);

    N =  0.18 * Njet / 0.115;
    
    inte1 = sum(N(:) * 0.005);

    inte2 = inte2 + inte1 * pt * 0.005;
   
    end
    
    plot(phi, inte2)
    
    hold all
    
    end
    
end