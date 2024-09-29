!!$*  subroutine for extending data_set distribution over extended BZ: 2D array and its related 1D variable (rvari)
      subroutine get_ext_variable(xrecivec,xrecilat,xrvari,kext, &
                       recivec,recilat,rvari,nnk,kperiod,nk,iz2, &
                       b1,b2,b3) !iz2=4 for z2, other 1
      implicit real*8 (a-h,o-z)
      real*8  recivec(3,nk*iz2),recilat(3,nk*iz2)
      real*8  xrecivec(3,kperiod*2*kperiod*2*nk*iz2)
      real*8  xrecilat(3,kperiod*2*kperiod*2*nk*iz2)
      real*8  rvari(nk*iz2),xrvari(kperiod*2*kperiod*2*nk*iz2)
      dimension b1(3),b2(3),b3(3)

      kk=0   ! extend variable distribution over extended BZ 
      do ib2=-1*(kperiod-1)+1,kperiod
         do ib1=-1*(kperiod-1)+1,kperiod
           do ik=1,nnk
              kk=kk+1
              xrecivec(1,kk)=recivec(1,ik) +(ib1-1)*b1(1)+(ib2-1)*b2(1)
              xrecivec(2,kk)=recivec(2,ik) +(ib1-1)*b1(2)+(ib2-1)*b2(2)
              xrecivec(3,kk)=recivec(3,ik) +(ib1-1)*b1(3)+(ib2-1)*b2(3)
              xrecilat(1,kk)=recilat(1,ik) +(ib1-1)
              xrecilat(2,kk)=recilat(2,ik) +(ib2-1)
              xrecilat(3,kk)=recilat(3,ik)
              xrvari(kk)=rvari(ik)
              kext=kk
           enddo
         enddo
      enddo

      return
      end subroutine 
