function Ridge

    Ni = 1;
    m = 1;
    T = 5;
    etajet = 0 ;
    sigmay = 5.5;
    q = 0.8;
    
    [phi,deta,pt] = meshgrid(-2:0.1:2,-1:0.1:1,-10:0.1:10);
    
    eta = deta + etajet;

    pf1 = pt.*cos(phi);
    pf2 = pt.*sin(phi);
    pf3 = pt.*sinh(eta);
    
    pi1 = pf1 - q/cosh(etajet);
    pi2 = pf2;
    pi3 = pf3 - q*sinh(etajet)./cosh(etajet);
    
    Ef = sqrt(pf1.^2 + pf2.^2 + pf3.^2 + m^2);
    Ei = sqrt(pi1.^2 + pi2.^2 + pi3.^2 + m^2);
    
    yf = log((Ef + pf3)./(Ef - pf3))./2;
    yi = log((Ei + pi3)./(Ei - pi3))./2;
    
    mtf = sqrt(m^2 + pf1.^2 + pf2.^2);
    mti = sqrt(m^2 + pi1.^2 + pi2.^2);
  
    pft = sqrt(pf1.^2 + pf2.^2);
    pit = sqrt(pi1.^2 + pi2.^2);
    
    etaf = log((sqrt(pf1.^2 + pf2.^2 + pf3.^2) + pf3)./((sqrt(pf1.^2 + pf2.^2 + pf3.^2) - pf3)));
    etai = log((sqrt(pi1.^2 + pi2.^2 + pi3.^2) + pi3)./((sqrt(pi1.^2 + pi2.^2 + pi3.^2) - pi3)));
    
    detaf = etaf - etajet;
    dphi = atan(pf2./pf1);
    
    Ai = Ni*exp(m/T)/(((2*pi)^(3/2))*(sigmay)*T); 
    
    N = Ai * exp((-yi.^2)/(2*sigmay^2)) .* exp(-mti/T)./mti .*  Ef./Ei .* sqrt(1 - m^2/(mtf.^2 .* (cosh(yf)).^2));
    
    disp(N);
    
    mesh(phi,deta,N)
    
    
end