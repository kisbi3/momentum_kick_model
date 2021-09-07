cc    program calxydy.f
      implicit real*8 (a-h,o-z)
      dimension x(100),y(100),xup(100),yup(100),err(100)

        open ( 6,file='calxydy.dat',status='unknown') 
	open ( 7,file='calxy.dat',   status='unknown')  !data
	open ( 8,file='caldy.dat',status='unknown')  !upper data

        write(6,701)
 701    format('@type xydy')

        do i=1,100
           read (7,*,end=20) x(i), y(i)
c           write(6,*) x(i), y(i)
       enddo

 20     imx=i-1

        do i=1,imx
           read (8,*,end=30) xup(i), yup(i)
           err(i)=yup(i)-y(i)
           write(6,*) x(i), y(i), err(i)
        enddo

 30     write(6,703)
 703    format('&')
        stop

         end
