!!$*  subroutine for computing planewave G index
      subroutine plindx(ig,ncnt, &
         ispinor,wk,b1,b2,b3,nbmax,np,ecut,npmax)
      implicit real*8(a-h,o-z)
      dimension wk(3),sumkg(3),b1(3),b2(3),b3(3),nbmax(3)
      integer :: ig(3,npmax)
      data c/0.262465831d0/
      ncnt=0
      do ig3=0,2*nbmax(3)
       ig3p=ig3
       if (ig3.gt.nbmax(3)) ig3p=ig3-2*nbmax(3)-1
       do ig2=0,2*nbmax(2)
       ig2p=ig2
        if (ig2.gt.nbmax(2)) ig2p=ig2-2*nbmax(2)-1
        do ig1=0,2*nbmax(1)
        ig1p=ig1
         if (ig1.gt.nbmax(1)) ig1p=ig1-2*nbmax(1)-1
         do j=1,3
          sumkg(j)=(wk(1)+ig1p)*b1(j)+ &
                   (wk(2)+ig2p)*b2(j)+ &
                   (wk(3)+ig3p)*b3(j)
         enddo
         gtot=sqrt(dot_product(sumkg,sumkg))
         etot=gtot**2/c
         if (etot.lt.ecut) then
          ncnt=ncnt+1
          ig(1,ncnt)=ig1p
          ig(2,ncnt)=ig2p
          ig(3,ncnt)=ig3p
         end if
        enddo
       enddo
      enddo
      if (ispinor*ncnt.ne.np) then
       write(0,*) '*** error - computed ispinor*ncnt=',ispinor*ncnt, &
                 ' != input nplane=',np;stop
      endif
      return
      end subroutine plindx
