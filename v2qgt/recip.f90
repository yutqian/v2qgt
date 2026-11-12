!!$*  subroutine for computing reciprocal lattice vector
      subroutine recilatt(b1,b2,b3,dSkxky, a1,a2,a3,nkx,nky)
      implicit real*8 (a-h,o-z)
      dimension a1(3),a2(3),a3(3), b1(3),b2(3),b3(3),a2xa3(3)
      dimension b1xb2(3)
      pi=4.*atan(1.)
      call vcross(a2xa3,a2,a3)
      Vcell=dot_product(a1,a2xa3)
      a3mag=dsqrt(dot_product(a3,a3))
      call vcross(b1,a2,a3);call vcross(b2,a3,a1);call vcross(b3,a1,a2)
      b1=2.*pi*b1/Vcell ; b2=2.*pi*b2/Vcell ; b3=2.*pi*b3/Vcell

      call vcross(b1xb2,b1,b2)
      dSkxky=dsqrt(dot_product(b1xb2,b1xb2))/real(nkx)/real(nky)
      return
      end subroutine recilatt

!!$*  subroutine for computing vector cross-product
      subroutine vcross(a,b,c)
      implicit real*8(a-h,o-z)
      dimension a(3),b(3),c(3)
      a(1)=b(2)*c(3)-b(3)*c(2)
      a(2)=b(3)*c(1)-b(1)*c(3)
      a(3)=b(1)*c(2)-b(2)*c(1)
      return
      end subroutine vcross

!!$*  subroutine for computing reciprocal properties
      subroutine reciproperty(nbmax,npmax, b1,b2,b3,ecut,ispinor)
      implicit real*8(a-h,o-z)
      dimension b1(3),b2(3),b3(3),vtmp(3),nbmax(3)
      data c/0.262465831d0/
      pi=4.*atan(1.)

      b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
      b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
      b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)

      phi12=acos((b1(1)*b2(1)+b1(2)*b2(2)+b1(3)*b2(3))/(b1mag*b2mag))
      call vcross(vtmp,b1,b2)
      vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
      sinphi123=(b3(1)*vtmp(1)+b3(2)*vtmp(2)+b3(3)*vtmp(3))/(vmag*b3mag)
      nb1maxA=(dsqrt(ecut*c)/(b1mag*abs(sin(phi12))))+1
      nb2maxA=(dsqrt(ecut*c)/(b2mag*abs(sin(phi12))))+1
      nb3maxA=(dsqrt(ecut*c)/(b3mag*abs(sinphi123)))+1
      npmaxA=nint(4.*pi*nb1maxA*nb2maxA*nb3maxA/3.)

      phi13=acos((b1(1)*b3(1)+b1(2)*b3(2)+b1(3)*b3(3))/(b1mag*b3mag))
      call vcross(vtmp,b1,b3)
      vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
      sinphi123=(b2(1)*vtmp(1)+b2(2)*vtmp(2)+b2(3)*vtmp(3))/(vmag*b2mag)
      phi123=abs(asin(sinphi123))
      nb1maxB=(dsqrt(ecut*c)/(b1mag*abs(sin(phi13))))+1
      nb2maxB=(dsqrt(ecut*c)/(b2mag*abs(sinphi123)))+1
      nb3maxB=(dsqrt(ecut*c)/(b3mag*abs(sin(phi13))))+1
      npmaxB=nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)

      phi23=acos((b2(1)*b3(1)+b2(2)*b3(2)+b2(3)*b3(3))/(b2mag*b3mag))
      call vcross(vtmp,b2,b3)
      vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
      sinphi123=(b1(1)*vtmp(1)+b1(2)*vtmp(2)+b1(3)*vtmp(3))/(vmag*b1mag)
      phi123=abs(asin(sinphi123))
      nb1maxC=(dsqrt(ecut*c)/(b1mag*abs(sinphi123)))+1
      nb2maxC=(dsqrt(ecut*c)/(b2mag*abs(sin(phi23))))+1
      nb3maxC=(dsqrt(ecut*c)/(b3mag*abs(sin(phi23))))+1
      npmaxC=nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)

      nbmax(1)=max0(nb1maxA,nb1maxB,nb1maxC)       ! maximum 
      nbmax(2)=max0(nb2maxA,nb2maxB,nb2maxC)
      nbmax(3)=max0(nb3maxA,nb3maxB,nb3maxC)

      !! multiply 'ispinor' to handle two component spinors
      npmax=ispinor*min0(npmaxA,npmaxB,npmaxC)
      return
      end subroutine reciproperty
