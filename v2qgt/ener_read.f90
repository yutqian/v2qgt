!!$*  subroutine for finding energy(nn) in the given k index
      subroutine ener_read(ener,isp,k,nk,nband)
         implicit real*8(a-h,o-z)
         dimension wk(3)
         complex*16 cener(nband)
         real*8 occ(nband)
         real*8 ener(nband)
         irec=3+(k-1)*(nband+1)+nk*(nband+1)*(isp-1)  !record addres for "k"-point
         read(10,rec=irec) xnplane,(wk(i),i=1,3), &
         (cener(nn),occ(nn),nn=1,nband)
         ener=real(cener)
!      write(6,*) 'check: nn=1 energies',real(cener(1)),ener(1),cener(1)
         return
      end subroutine ener_read 
