
cc    pkickden.f, to calcaulte the distribtuion of partons 
cc    the denominator is modified to have a different mass in transverse mass
cc       dN/ptdpt ~ dexp(-sqrt(am**2+pt**2)/T)/ sqrt(amden**2+pt**2)
cc    we also plot dN/ptdpt (for pickden.07 for the case adding the jet compnt
cc      after a momentum kick  P(q)=delta (qvec - q1)
cc      q1 is along jet direction, input, allow a distribution in jet 
cc      sigyi, Ti are input parameters  
cc    f77 -o pkickden.x pkickden.f

      IMPLICIT REAL*8 (A-H,O-Z)

      dimension Py(100,100,100),Peta(100,100,100)
      dimension fiv  (100),etav(100),ptv(100)
      dimension wetaj(100),weta(100),wfi(100)
      dimension dNptdpt(100),dNdpt(100)
      dimension dNdfi(100),dNdfij(100),dNdfitot(100),dNdeta(100)
      dimension dNptdptj(100)

	open ( 5,file='pkickden.in',status='unknown')
        open ( 6,file='pkickden.xn',status='unknown') 
	open ( 7,file='pkickden.07',status='unknown')  !plot yield vs pt
                                                          !gives ptdisn.agr
	open ( 8,file='pkickden.08',status='unknown')  !plot dndfi @y=0 & 1.7
                                                          !gives part of fitcombinen.agr
	open (10,file='pkickden.10',status='unknown')  !plot dn/deta(eta),@phi=0,dphi 
                                                          !gives part of fitcombinen.agr
                                                          !      all  @ pt=2 GeV 
	open ( 9,file='pkick.09',status='unknown')  !prepr yld(fi,eta) 4 mathca 


 1     read (5,*,end=999)  ietatype, itri2rec
                ! if ietatype=1, out eta=etaLabi-etajet=Deltaeta(STAR)
                ! if ietatype=2, out eta=etaLabi
                ! if itri2rec=0, no triangle->rec corecntn made
                ! if itri2rec=1, tri->rec extrapolted enhanced area ~2x
       write (6,*) ' Input: ietatype =',ietatype
       write (6,*) ' if ietatype=1,eta=etaLabi-etajet=Deltaeta(STAR)'
       write (6,*) ' if ietatype=2,eta=etaLabi'   
       write (6,*) ' Weshould mostlyuse ietatype=2,for detctr etalimits'
       write (6,*) itri2rec,' [if =0,no tri->rec|if=1,tri->rec,~2x]'

       read (5,*)  etajmn,etajmx,netaj,nsetaj ! etajet trigger,mn,mx,ietaj,#sides
       if (netaj .gt. 1)   go to  15
       detaj=0.d0
       wetaj(1)=1.d0*nsetaj
       sumnetaj=1.d0
       go to 20
 15    detaj=(etajmx-etajmn)/dfloat(netaj-1)
          wetaj    (1)=0.5d0*nsetaj           ! weights for etaj integration
       do ietaj=2,netaj-1
          wetaj(ietaj)= 1.d0*nsetaj
       enddo
          wetaj(netaj)=0.5d0*nsetaj
       sumnetaj=dfloat(netaj-1)

 20    write(6,*)  '(Jet) etajmn,etajmx,netaj=',etajmn,etajmx,netaj
       write(6,*)  '(Jet)               detaj=',detaj
       write(6,*)  '# of above sides of (Jet),nsetaj=',nsetaj

       read (5,*)  am,q1,sigyi,Ti
cc     initial pt distribution=exp(-amT/T)/amT; amT=sqrt(ampi^2+pt^2)
       write(6,*)  ' am(GeV)=',am
       write(6,*)  ' kick in jet directn(GeV)q1=',q1
       write(6,*)  ' initial sigyi=',sigyi
       write(6,*)  ' initial T_i(GeV)=',Ti    
       write(6,*)  ' initial pt dis=exp(-amT/T)/amT'
       write(6,*)  '                amT=sqrt(ampi^2+pt^2) '
       write(6,*)  ' transfer q-distribution 3--beam, 2 -- beam x jet'
       read (5,*)  amden
       write(6,*)  'mass-denomintr is changed to amden=',amden
       read (5,*)  anjetpt, Tjet, anjetfi, sigfi  
            ! anjetpt,Tjet =# and T of jet partons,(integrated ovr eta, fi)
            ! anjetfi,sigfi=# and wdth of jet partons,(intg ovr eta, pt sectn)
       write(6,*)  ' parton # in jet=anjetpt(intg eta fi)=',anjetpt
       write(6,*)  ' Tempertaure of jet=Tjet(GeV)  =',Tjet 
       write(6,*)  ' parton # in jet=anjetfi(intg eta, part-pt=',anjetfi
       write(6,*)  ' parton fi-width for pt secctions  =',sigfi 
       anormjetpt=1.d0/( Tjet * (am+Tjet) * dexp(am/Tjet) )  ! normalizn const
       anormjetfi=1.d0/(dsqrt(6.2831853d0)*sigfi)            ! normalizn const

       read (5,*)  etamn,etamx,neta,nseta ! fnl eta output depndg on ietatype:
         ! if ietatype=1, eta=Deta=etaLabi-etajet or etaLabi=eta+etajet
         ! if ietatype=2, eta=etalab              or etaLabi=eta
         ! nseta = 1 only 1 eta(positive) side included, =2 both side incld

       if (neta .gt. 1)   go to  25
       deta=0.d0
       weta(1)=1.d0*nseta
       go to 30

 25    deta=(etamx-etamn)/dfloat(neta-1)
       weta   (1)=0.5d0*nseta           ! weights for etaj integration
          do ieta=2,neta-1
          weta(ieta)=1.0d0*nseta
          enddo
       weta(neta)=0.5d0*nseta
   
 30    write (6,*)  '(Output eta) etamn,etamx,neta=', etamn,etamx,neta 
       write (6,*)  '(Output eta)             deta=', deta
       write(6,*)   '# of sides of above eta nseta=', nseta

       read (5,*)  fimin, dfi, nfi, nsfi  ! final particle phi
       write (6,*)  'fimin, dfi, nfi =', fimin, dfi, nfi    
       write(6,*)   '# of above eta sides  nsfi=', nsfi

       if (nfi .gt. 1)   go to  35
       wfi(1)=1.d0*nsfi
       go to 40
 35    wfi  (1)=0.5d0*nsfi           ! weights for fij integration
       do ifi=2,nfi-1
          wfi(ifi)=1.0d0*nsfi
       enddo
          wfi(nfi)=0.5d0*nsfi

 40    read (5,*)  ptmin, dpt, npt ! final particle pt 
       write (6,*)  'ptmin, dpt, npt =', ptmin, dpt, npt


       read (5,*)  iptplt,ifiplt,ietaplt,ifietaplt
                  ! if above=0, no plot
                  ! if 1 for iptplt,ifiplt,ietaplt, plot integrated plot
                  ! if 2 for        ifiplt, single plot fixed eta&pt
                  ! if 2 for        ietaplt,single plot fixed fi&pt
                  ! if ifietaplt=0 no plot of yield(phi,eta)
                  ! if ifietaplt=1 plot dn/dphideta(phi,eta), fi=x, eta=y
                  ! if ifietaplt=2 plot dn/dphidets(phi,eta), eta=x, fi=y
       write (6,*)  'iptplt=   ',iptplt    ! if=1 plot yield vs pt, intg eta&fi
       write (6,*)  'ifiplt=   ',ifiplt    ! if=1 plot yield vs phi,@ e.g.eta=0?
       write (6,*)  'ietaplt=  ',ietaplt   ! if=1 plot yield vs eta,@ e.g.fi=0?
       write (6,*)  'ifietaplt=',ifietaplt ! if=1 get yield vs fi-eta,@e.g.pt=2GeV?
       write (6,*) 'iptplt,ifiplt,ietaplt,ifietaplt=',
     >              iptplt,  ifiplt, ietaplt, ifietaplt
       read (5,*)          fiplteta,etapltfi,fietapltpt 
       write(6,*) 'locations of plot is fixed for (xx=)'
       write(6,*) '        fiplt(eta)=,ietaplt(fi),fietaplt(pt)',
     >                     fiplteta,etapltfi,fietapltpt

       read (5,*)  Aptplt,Afiplt,Aetaplt,Afietaplt   ! constants for plotting 
       write (6,*)  'Aptplt,constant for plotng=',Aptplt
       write (6,*)  'Afiplt,constant for plotng=',Afiplt
       write (6,*)  'Aetaplt,constant for plotng=',Aetaplt
       write (6,*)  'Afietaplt,const for plt fi-eta=',Afietaplt

       write (7,701)
       write (8,701)
       write (10,701)
 701   format('@type xy')
       sum=0.d0
       do ipt=1,400
          pt=dfloat(ipt-1)*0.01
          amti=dsqrt(am**2+pt**2)
          term=dexp(-amti/Ti)/dsqrt(amden**2+pt**2)
          sum=sum+term*pt
       enddo
       sum=6.2831853*sum*0.01d0
       anormpt=1.d0/sum
 
c       aN=dexp(am/Ti)/((2.d0*3.14159d0)**1.5*sigyi*Ti)
       aN=anormpt/(dsqrt(2.d0*3.14159d0)*sigyi)

         do 100 ieta=1,neta
            eta=etamn+(ieta-1)*deta
            etav(ieta)=eta

            do 90  ifi=1,nfi
               fi=fimin+(ifi-1)*dfi
               fiv(ifi)=fi

               do 80  ipt=1,npt
                  pt=ptmin+(ipt-1)*dpt
                  ptv(ipt)=pt
                                    ! these are final p
                                    ! to beam3-z,1&2 are Trnsv-perp
                  sumdndy=0.d0
                  sumdndeta=0.d0

               do 70 ietaj=1,netaj   ! average over all etajet distributn
                  etaj=etajmn+(ietaj-1)*detaj    ! eta of jet trigger
                  if (ietatype .eq. 1)  etaLab=eta+etaj !eta=output=etaLab-etaj
                  if (ietatype .eq. 2)  etaLab=eta      !eta=output=etaLab

                  p3=pt*dsinh(etaLab)  ! p in 3 (z-direction, beam) uses etaLab
                  p2=pt*dsin(fi)       ! p in 2 (y-direction, beamXjet)
                  p1=pt*dcos(fi)       ! p in 1 (x-direction, jet)
                  amt=dsqrt(am**2+p1**2+p2**2)
                  Ef =dsqrt(am**2+p1**2+p2**2+p3**2)
                  y  =0.5d0*dlog((Ef+p3)/(Ef-p3))
                              ! p(final)=p=p(init)+q, p1i i for initial
                p1i=p1-q1/dcosh(etaj)   ! etaj is eta of jet
                p2i=p2
                p3i=p3-q1*dsinh(etaj)/dcosh(etaj)
                amti=dsqrt(am**2+p1i**2+p2i**2)  ! to beam3,1&2 are Trnsv-perp
                amtiden=dsqrt(amden**2+p1i**2+p2i**2)     ! denominator amden 
                Ei  =dsqrt(am**2+p1i**2+p2i**2+p3i**2)
                yi  =0.5d0*dlog((Ei+p3i)/(Ei-p3i)) 
c                termy=aN*dexp(-yi**2/(2.d0*sigyi**2))
c     >               *(dexp(-amti/Ti)/amti)
c     >               *Ef/Ei
                termy=aN*dexp(-yi**2/(2.d0*sigyi**2))
     >               *(dexp(-amti/Ti)/amtiden)
     >               *Ef/Ei
                termeta=termy*dsqrt(1.d0-am**2/(amt**2*dcosh(y)**2))
                sumdndy  =sumdndy  +wetaj(ietaj)*termy     ! multiplied by weight
                sumdndeta=sumdndeta+wetaj(ietaj)*termeta   ! multiplied by weight
 70             continue

                sumdndy  =sumdndy  /sumnetaj  ! assume dN/detajet=uniform
                sumdndeta=sumdndeta/sumnetaj  ! assume dN/detajet=uniform

         Py  (ieta,ifi,ipt)=sumdndy          ! dN/dy dphi pt dpt
         Peta(ieta,ifi,ipt)=sumdndeta        ! dN/deta dphi pt dpt
c         write (6,*) '    y,fi,pt',  y,fi,pt,Py  (ieta,ifi,ipt)
         write (6,*) ' eta,fi,pt', eta,fi,pt,Peta(ieta,ifi,ipt)
 80      continue
!         write (6,*) 
 90      continue
         write (6,*) 
 100     continue

       if (iptplt  .eq. 0)   go to 105  
       do ipt=1,npt
          pt=ptmin+(ipt-1)*dpt
          sumdNptdpt=0.d0
          do ieta=1,neta
            do ifi=1,nfi
            sumdNptdpt=sumdNptdpt+weta(ieta)*wfi(ifi)*Peta(ieta,ifi,ipt)
            enddo
          enddo
          sumtm  =sumdNptdpt*deta*dfi ! un-normalized      dN/pt dpt
          dNptdpt(ipt)=Aptplt   *sumtm     ! normalized(byhand) dN/pt dpt 
          dNdpt  (ipt)=Aptplt*pt*sumtm     ! normalized(byhand) dN/   dpt 
          if (ipt .eq. 1)  Amatch=dNdpt(ipt)    ! Amatch is dNdpt at ipt=1 (pt=ptmin)
       enddo

       write (7,710)
 710   format('#  pt, dN/ptdpt(pt) integrated over fi&eta')
       do ipt=1,npt
          write (7,*) ptv(ipt),dNptdpt(ipt)      ! plot out integrated yields dN/ptdpt
       enddo
       write (7,703)
       write (7,711)
 711   format('#  pt, dN/dpt(pt) integrated over fi&eta')
       do ipt=1,npt
          write (7,*) ptv(ipt),dNdpt(ipt)      ! plot out integrated yields dN/ptdpt
       enddo
       write (7,703)
       write (7,712)
 712   format('#  pt, dN/ptdpt(pt) of jet')
       do ipt=1,npt
          amti=dsqrt(am**2+ptv(ipt)**2)
          dnptdptj(ipt)=anjetpt*anormjetpt*dexp(-amti/Tjet)
          write (7,*) ptv(ipt),dnptdptj(ipt)      ! plot out jet dN/ptdpt
       enddo
       write (7,703)
       write (7,713)
 713   format('#  pt,dN/ptdpt(pt), jet+ridge integ  fi&eta+jetdN/ptdpt')
       do ipt=1,npt
          write (7,*) ptv(ipt),dnptdpt(ipt)+dnptdptj(ipt)      
                              ! plot out integrated yields dN/ptdpt
       enddo
       write (7,703)

       write (7,715)
 715   format('#  pt, dNdpt(pt) init formla noramled to first pt pt')
       amti=dsqrt(am**2+ptmin**2)  ! to beam3,1&2 are Trnsv-perp
       amtiden=dsqrt(amden**2+ptmin**2)     ! denominator amden 
c       dNdpttm=ptmin*dexp(-amti/Ti)/amti
       dNdpttm=ptmin*dexp(-amti/Ti)/amtiden
       Anorm=Amatch/dNdpttm        ! normalize the curves
       do ipt=1,npt
          pt=ptmin+(ipt-1)*dpt
           amti=dsqrt(am**2+pt**2)  ! to beam3,1&2 are Trnsv-perp
           amtiden=dsqrt(amden**2+pt**2)     ! denominator amden 
           dNdptinit=Anorm*pt*dexp(-amti/Ti)/amtiden
          write (7,*) pt,dNdptinit    ! plot initial pt dis
       enddo
       write (7,703)
 703   format('&')


 105   if (ifiplt  .eq. 0)   go to 130  
       if (ifiplt  .eq. 1)   go to 110     ! plot out integrtd eta, integrtd pt
       if (ifiplt  .eq. 2)   go to 120     ! plot out a single (eta, pt)

 110  do ifi=1,nfi
          fi=fimin+(ifi-1)*dfi
          sumdNdfi=0.d0
          do ieta=1,neta
           do ipt=1,npt
            pt=ptmin+(ipt-1)*dpt
            term=weta(ieta)*pt*Peta(ieta,ifi,ipt)
            if (ipt .eq. 1)  term=0.5d0*term     ! the weight for ipt=1 is 0.5
            if (itri2rec .eq. 0)  go to 112
            etaetaj=2-etav(ieta)        ! sum etamx+etaj(deltaeta)=2*etamx-delta_eta
            if (dabs(etaetj) .le. 1.d-5)  etaetaj=0.5d0*deta  
            term=term*2.d0/etaetaj      ! the tri2rec scaling is 2/(2-delta_eta)  
 112        sumdNdfi=sumdNdfi+term
           enddo
          enddo
          sumtm     =sumdNdfi*deta*dpt ! dN/dfi
          dNdfi(ifi)=Afiplt*sumtm      ! put in norm Aptplt to normlze 
          dNdfij(ifi)=Anjetfi*anormjetfi*dexp(-fi**2/(2.d0*sigfi**2))
          dNdfitot(ifi)=dNdfi(ifi)+dNdfij(ifi)
       enddo
       write (8,720) 
 720   format('#  fi, dN/dfi integrated over eta&pt')
       write (8,721) etamn,etamx,neta,nseta
 721   format('# etamn,etamx,neta,nseta=',2f10.4,I5,f10.4)
       write (8,722) fimin, dfi, nfi, nsfi
 722   format('# fimin, dfi, nfi, nsfi=',2f10.4,I5,f10.4)

       if (nsfi.eq.1)  go to 115
       do ifi=nfi,2,-1
          write (8,*) -fiv(ifi),dNdfi(ifi)  !plot negative part of phi
       enddo
 115   do ifi=1,nfi
          write (8,*)  fiv(ifi),dNdfi(ifi)  !plot positive part of phi
       enddo
       write (8,703)
       if (nsfi.eq.1)  go to 117
       do ifi=nfi,2,-1
          write (8,*) -fiv(ifi),dNdfij(ifi)  !plot negative part of phi
       enddo
 117   do ifi=1,nfi
          write (8,*)  fiv(ifi),dNdfij(ifi)  !plot positive part of phi
       enddo
       write (8,703)
       if (nsfi.eq.1)  go to 119
       do ifi=nfi,2,-1
          write (8,*) -fiv(ifi),dNdfitot(ifi)  !plot negative part of phi
       enddo
 119   do ifi=1,nfi
          write (8,*)  fiv(ifi),dNdfitot(ifi)  !plot positive part of phi
       enddo
       write (8,703)



       go to 130

 120   write (8,810) fiplteta,ptv(1)
 810   format('#  fi, single(eta,pt)dN/ptdpt@eta=',f10.5,',pt=',f10.5)
       ieta=(fiplteta-etamn)/deta+1
       do ifi=nfi,2,-1
          fi=fimin-(ifi-1)*dfi
          write (8,*) fi,Afiplt*Peta(ieta,ifi,1) !plot eta=fiplteta,pt=2GeV yld
       enddo
       do ifi=1,nfi
          fi=fimin+(ifi-1)*dfi
          write (8,*) fi,Afiplt*Peta(ieta,ifi,1) !plot eta=fiplteta,pt=2GeV yld
       enddo
       write (8,703)

 130   if (ietaplt .eq. 0)  go to 160
       if (ietaplt .eq. 1)  go to 140  !plot dn/deta, integrated over fi, pt 
       if (ietaplt .eq. 2)  go to 150  !plot dn/deta at a single (fi,pt)

cc     plot dn/deta here, integrated over fi, pt 
 140   do ieta=1,neta
          sumdNdeta=0.d0
          do ifi=1,nfi
           do ipt=1,npt
            pt=ptmin+(ipt-1)*dpt
            term=wfi(ifi)*pt*Peta(ieta,ifi,ipt)
            if (ipt .eq. 1)  term=0.5d0*term     ! the weight for ipt=1 is 0.5
            sumdNdfi=sumdNdfi+term
           enddo
          enddo
          sumtm     =sumdNdeta*dfi*dpt ! dN/dfi
          dNdeta(ieta)=Aetaplt*sumtm      ! put in norm Aptplt to normlze 
       enddo
       write (10,740) 
 740   format('#  dN/deta, integrated over fi&pt')
       if (nseta.eq.1)  go to 145
       do ieta=neta,2,-1
          write (10,*) -etav(ieta),dNdeta(ieta)  !plot negative part of eta
       enddo
 145   do ieta=1,neta
          write (10,*)  etav(ieta),dNdeta(ieta)  !plot positive part of eta
       enddo
       write (10,703)
       go to 160


 150   write (10,730)  fiv(1),ptv(1)
 730   format('#  eta varying,single(pt,fi)@pt=2GeV,fi=',f10.4)
       do ieta=1,neta
          write (10,*)  etav(ieta),Aetaplt*Peta(ieta,1,1)
       enddo
       write (10,703)


 160   if (ifietaplt  .eq. 0)   go to 170  

       if (ifietaplt  .eq. 2)   go to 165
cc     plot  fi=x, eta=y, yield=z

         do ieta=1,neta
            write (9,*) (Afietaplt*Peta(ieta,ifi,1),ifi=nfi,2,-1),
     >                  (Afietaplt*Peta(ieta,ifi,1),ifi=1,nfi)    ! at pt=2 GeV
         enddo                ! fi is x, and eta is y.

         go to 170

cc     plot  eta=x, fi=y, yield=z
 165        do ifi=-nfi+2,nfi
            if (ifi .le. 0)  ifiabs=iabs(ifi)+2
            if (ifi .gt. 0)  ifiabs=ifi
           write (9,*) (Afietaplt*Peta(ieta,ifiabs,1),ieta=1,neta)    
                              ! at pt=2 GeV
         enddo                ! eta is x, and fi is y.

 170   continue
 999   stop
       end

