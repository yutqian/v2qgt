!!$*  subroutine for sorting
      subroutine get_sorted_xrvari(xrecivec,xrecilat,xrvari,&
                             kext,kperiod,nk,iz2)
      implicit real*8 (a-h,o-z)
      real*8  xrecivec(3,kperiod*2*kperiod*2*nk*iz2)
      real*8  xrecilat(3,kperiod*2*kperiod*2*nk*iz2)
      real*8  xrvari(kperiod*2*kperiod*2*nk*iz2)
      real*8  xb(3),xtemp

      do k=kext-1,1,-1  ! sorting kx
          do j=1,k
            if(xrecivec(1,j+1) .gt. xrecivec(1,j))then
               xb(:)=xrecivec(:,j)
               xrecivec(:,j)=xrecivec(:,j+1)
               xrecivec(:,j+1)=xb(:)
               xb(:)=xrecilat(:,j)
               xrecilat(:,j)=xrecilat(:,j+1)
               xrecilat(:,j+1)=xb(:)
               xtemp=xrvari(j)
               xrvari(j)=xrvari(j+1)
               xrvari(j+1)=xtemp
            endif
          enddo
      enddo
      do k=kext-1,1,-1  ! sorting ky
          do j=1,k
            if(xrecivec(1,j+1) .eq. xrecivec(1,j))then
                if(xrecivec(2,j+1) .gt. xrecivec(2,j))then
                xb(:)=xrecivec(:,j)
                xrecivec(:,j)=xrecivec(:,j+1)
                xrecivec(:,j+1)=xb(:)
                xb(:)=xrecilat(:,j)
                xrecilat(:,j)=xrecilat(:,j+1)
                xrecilat(:,j+1)=xb(:)
                xtemp=xrvari(j)
                xrvari(j)=xrvari(j+1)
                xrvari(j+1)=xtemp
                endif
            endif
          enddo
      enddo
      return
      end subroutine get_sorted_xrvari
