cc    program pardis.f, to calcualte parton distribution function
cc       and dN/dy=aNi*dF/dy  
cc    dF/dy=(normalized dN/dy) after integrating pt dpt
cc    the momentum distribution dF/ptdptdy is normalized to one
cc    initial parton, inclusive, and ppjet


      implicit real*8 (a-h,o-z)
      dimension yyv(500),ptv(500),dNdyptdpt(500,500)
      dimension dFdy(500)
 
        open ( 5,file='pardis.in',status='unknown') ! input data pt,J pt,R 
	open ( 6,file='pardis.xn',status='unknown') ! standard output
	open ( 7,file='pardis.07',status='unknown') ! aNi*dFdy(y),pt-intgtd
	open ( 8,file='pardis.08',status='unknown') ! plt dFdy(pt,y) 
	open ( 9,file='pardis.09',status='unknown') ! 3dplotting


 5      read (5,*,end=999)  Ti,aa,amden,amparent  !amparent=parentparton mass
                                   ! amparent=0.14 for 200GeV, ~0.4 for 62GeV
        write(6,*) 'Ti,aa,amden,amparent', Ti,aa,amden,amparent
        read (5,*)  ppsqrts,am,aNi  ! sqrts(pp),mass detect part
        write(6,*) 'ppsqrts,am,aNi',ppsqrts,am,aNi  ! dndy=aNi*dFdy    
        read (5,*)  ktrldndy,ktrlpt,ktrlyy,ktrl3d !plt dndy?F(pt)_y?F(y)_pt?3D?
        write(6,*) 'ktrldndy,ktrlpt,ktrlyy,ktrl3d',
     >              ktrldndy,ktrlpt,ktrlyy,ktrl3d 
        read (5,*)  iptmx,iptplt,ipt3dplt
        write(6,*) 'iptmx,iptplt,ipt3dplt',iptmx,iptplt,ipt3dplt
        read (5,*)  iyymx,iyyplt,iyy3dplt
        write(6,*) 'iyymx,iyyplt,iyy3dplt',iyymx,iyyplt,iyy3dplt

       emrat=(ppsqrts/2.d0)/0.93890595d0   ! ppsqrts=bpothpp sqrts 
                                           ! E(p) in CMS=ppsqrts/2
       yyb=dlog(emrat+dsqrt(emrat**2-1.d0))  ! yy of beam in CMS=mN*cosh(yb)
       write(6,*) 'yyb=',yyb
cc  examples:        Ti=0.5d0 (GeV);  aa=0.5d0; amden=1.d0(GeV); aNi=1.d0
cc                   am=0.14;  sqrts=200 (in GeV)

       bottom=1.d-7     ! for pt plotting, what is the bottom of scales
       ptmx=4.4d0
       dyy=5.4d0/dfloat(iyymx)
       dpt=ptmx/dfloat(iptmx) 

       SUM=0.D0

       do 40 iyy=1,iyymx
          yy=0.d0+(dfloat(iyy)-1)*dyy
          yyv(iyy)=yy
          sumdFdy=0.d0
       do 30 ipt=1,iptmx
          pt=0.d0+(dfloat(ipt)-1.d0)*dpt
          ptv(ipt)=pt
          amti =dsqrt(am   **2+pt**2)
          denom=dsqrt(amden**2+pt**2)
          xx=(amti/amparent)*dexp(yy-yyb)
          if (xx .le. 1.d0)  go to 10
          term=0.d0
          go to 20
 10       term=(1.d0-xx)**aa*dexp(-amti/Ti)/denom
 20       dndyptdpt(iyy,ipt)=term
          sum=sum+pt*term
          sumdFdy=sumdFdy+pt*term
 30    enddo
       dFdy(iyy)=sumdFdy*dpt                ! nromalized dN/dy
 40   enddo

          sum=2.d0*sum*dyy*dpt

       do iyy=1,iyymx
          dFdy(iyy)=aNi*dFdy(iyy)/sum
       do ipt=1,iptmx
             dndyptdpt(iyy,ipt)=dndyptdpt(iyy,ipt)/sum   ! multiply by Ani
       enddo
       enddo

        if (ktrldndy .ne. 1)  go to 100 
        write(7,701)
        write(7,705) Ti,aa,amden,aNi
 705    format('# Ti=',f10.4,' aa=',F10.4,' amden=',F10.4,' Ani=',F10.4)
        write(7,706) ppsqrts,am   
 706    format('# sqrts=',f10.4,' GeV,  am=',F10.4)
        write(7,709) 
 709    format('# yy | dF/dy =',f10.3)
       do iyy=1,iyymx
       write (7,*) yyv(iyy),dFdy(iyy)
       enddo
       write (7,703)

 100   if (ktrlpt .ne. 1)  go to 200
        write(8,701)
        write(8,705) Ti,aa,amden,aNi
        write(8,706) ppsqrts,am
        do iyy=1,iyymx,iyyplt
        write(8,710) yyv(iyy)
 710    format('# pt | dN/dyptdpt for yy=',f10.3)

          do ipt=1,iptmx
             term=dndyptdpt(iyy,ipt)
             if (term.ne.0.d0)  write (8,*) ptv(ipt),term
          enddo
       amtmx=am*dexp(-(yyv(iyy)-yyb))
       ptmx=dsqrt(amtmx**2-am**2)
       write (8,*) ptmx, bottom

       write (8,703)
       enddo

 200   if (ktrlyy .ne. 1)  go to 300
        write(8,701)
        write(8,705) Ti,aa,amden,aNi
        write(8,706) ppsqrts,am
        do ipt=1,iptmx,iptplt
        write(8,715) ptv(ipt)
 715    format('# y | dN/dyptdpt for pt=',f10.3,' GeV')

        yymx=yyb-dlog(dsqrt(am**2+ptv(ipt)**2)/am)
        write (8,*) -yymx, bottom
          do iyy=iyymx,2,-1
             term=dndyptdpt(iyy,ipt)
             if (term.ne.0.d0)  write (8,*) -yyv(iyy),term
          enddo
          do iyy=1,iyymx
             term=dndyptdpt(iyy,ipt)
             if (term.ne.0.d0)  write (8,*)  yyv(iyy),term
          enddo
        write (8,*)  yymx, bottom
       write (8,703)
       enddo

 300   if (ktrl3d .ne. 1)   go to 400
       do ipt=1,iptmx,ipt3dplt
       write (9,*) (dndyptdpt(iyy,ipt),iyy=iyymx,2,-iyy3dplt),
     >             (dndyptdpt(iyy,ipt),iyy=1,iyymx, iyy3dplt)    
       enddo

 400   go to 5
 999   stop

 701   format('@type xy')
 703   format('&')
       end
