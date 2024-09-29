!!$*  subroutine for writing results
      subroutine write_result(isp,ispin,ispinor,fonameo,foname,filename,&
                       irecl,ecut,nk,nkx,nky,nband,b1,b2,b3,kperiod,   &
                       dSkxky,nini,nmax,xrecivec,xrecilat,kext,        &
                       xrvari,rvari,rvari2,rvari3,rvari4,              &
                       iz2,icd,iz,ivel,ikubo, nprocs)
      implicit real*8 (a-h,o-z)
      dimension b1(3),b2(3),b3(3)
      real*8  xrecivec(3,kperiod*2*kperiod*2*nk*iz2)
      real*8  xrecilat(3,kperiod*2*kperiod*2*nk*iz2)
      real*8  xrvari(kperiod*2*kperiod*2*nk*iz2,3)
      real*8  rvari,rvari2
      real*8  xb(3),xtemp,rvari3(4),rvari4(4)
      character*256 filename,foname,fonameo
      character*256 foname_
      character*20,external :: int2str
       
      if(ikubo .ge. 1 .and. nini .eq. nmax) then
         write(foname_, '(3A )')TRIM(foname),'.EIG-', &
                               trim(ADJUSTL(int2str(nini)))
      else
         foname_ = foname
      endif

      if(isp .eq. 1 .and. ispinor .eq. 1 .and. ispin.eq.2) then
         write(fonameo,'(A,A)')TRIM(foname_),'.UP.dat'
      else if (isp .eq. 2 .and. ispinor .eq. 1 .and. ispin.eq.2)then
          write(fonameo,'(A,A)')TRIM(foname_),'.DN.dat'
      else if (isp .eq. 1 .and. ispinor .eq. 2) then
          write(fonameo,'(A,A)')TRIM(foname_),'.dat'
      else if (isp .eq. 1 .and. ispinor .eq. 1 .and. ispin.eq.1)then
          write(fonameo,'(A,A)')TRIM(foname_),'.dat'
      endif

      open(32,file=fonameo,status='unknown')
      write(32,'(A,I4,A)')"# Job running on ",nprocs," total cores"
      write(32,'(A,A)')   "# File reading...  : ",filename
      write(32,'(A,I9)')"# TOTAL RECORD LENGTH = ",irecl
      if (ispinor .eq. 2)then
        write(32,'(A,I6,A)')"# ISPIN            : ",ispin, &
                           " (LSORBIT = .TRUE.)"
      else
         write(32,'(A,I6,A)')"# ISPIN            : ",ispin,&
                            " (LSORBIT = .FALSE.)"
      endif
      write(32,'(A,F11.4)')  "# ENCUT (eV)       : ",ecut
      write(32,'(A,I6)')     "# NKPOINT          : ",nk
      write(32,'(A,I6,A,I4)')"#  K-GRID          : ",nkx,"   X",nky
      write(32,'(A,I6)')     "# NBANDS           : ",nband
      write(32,'(A,3F13.6)') "# RECIVEC B1 (A^-1): ",(b1(i),i=1,3)
      write(32,'(A,3F13.6)') "# RECIVEC B2       : ",(b2(i),i=1,3)
      write(32,'(A,3F13.6)') "# RECIVEC B3       : ",(b3(i),i=1,3)
      write(32,'(A,F13.6)')   "#  dk^2 = |dk1xk2| = ",dSkxky
      write(32,*)" "

      if(nini .eq. nmax) then
         if(iz == 1)then
          write(32,'(A,I4)')"# Z2 invariant for the BAND : ",nini
         elseif(iz+ivel+icd .eq. 0)then
           write(32,'(A,I4)')"# Chern Number for the BAND : ",nmax
         endif
      else
        if(iz == 1)then
         write(32,'(A,I4,A,I4)')"# Z2 invariant for the BANDS : ",nini,&
                              " - ",nmax
        elseif(iz+ivel+icd .eq. 0)then
         write(32,'(A,I4,A,I4)')"# Chern Number for the BANDS : ",nini,&
                              "    -  ",nmax
        endif
      endif
 
      if(iz == 1)then !Z2 INVARIANT
         write(32,'(A,I2)')"# Z2 Invariant (top) =    ",& 
                       mod((nint(rvari)),2)
         write(32,'(A,I2)')"# Z2 Invariant (bottom) =    ",&
                        mod((nint(rvari2)),2)
         write(32,'(A)')"# Berry Curvature F (A^2) : &
          -Im[logPI_S(det(S(K_s,K_s+1)))]/dk^2"
         write(32,'(A)')"# N-field strength       : &
          Sum_s{Im[log(det(S(K_s,k_s+1)))]}/2pi - F/2pi "
         write(32,'(A)')"# (cart) kx        ky        kz(A^-1) &
            n-field strength      ,   (recip)kx        ky        kz"

      elseif(icd == 1)then !OPTICAL SELECTIVITY
       write(32,'(A,I4,A,I4)')"# OPTICAL SELECTIVITY BETWEEN BANDS: ",&
                             nini,"    -  ",nmax
       write(32,'(A)')"# n(k,w_cv)= |P(k,s,cv,+)|^2 - |P(k,s,cv,-)|^2"
       write(32,'(A)')"#            ---------------------------------"
       write(32,'(A)')"#            |P(k,s,cv,+)|^2 + |P(k,s,cv,-)|^2"
       write(32,'(A)')"#  The TRANSITION MATRIX ELEMENT P ="
       write(32,'(A)')"#   P(k,s,cv,+ or -) = 1/sqrt(2)[p_x(k,cv,s) +&
          (or -) i*p_y(k,cv,s)]"                                        
       write(32,'(A)')"#  THE INTERBAND TRANSITION MATRIX p_x,y ="      
       write(32,'(A)')"#   p_x,y(k,cv,s)=<psi(k,c,s)|-i*hbar*1/dx(y)|  &
          psi(k,v,s>"                                                   
       write(32,'(A,4F16.6)')"# MAXVAL of SELECTIVITY at kx,ky,kz      &
          (in reci)= ",(rvari3(i),i=1,4)                                
       write(32,'(A,4F16.6)')"# MINVAL of SELECTIVITY at kx,ky,kz      &
          (in reci)= ",(rvari4(i),i=1,4)                                
       write(32,'(A)')"# (cart) kx        ky        kz(A^-1)           &
          selectivity(n(k)),        (recip)kx        ky        kz"

      elseif(ikubo .ge. 1)then !Berry curvature using Kubo formula
         write(32,'(A)')"# Chern Number is sum of &
             Berry Curvature over 1BZ"
         write(32,'(A,F16.4)')"# Chern Number =   ",rvari
         write(32,'(A,I4)')"# Berry curvature using kubo for BAND &
             index n ", nini
         write(32,'(A)')"# Omega_n=     |P(k,s,nm,+)|^2-|P(k,s,nm,-)|^2"
         write(32,'(A)')"#        sum_n -------------------------------"
         write(32,'(A)')"#       (n/=m) |energy(k,s,n)-energy(k,s,m)|^2"
         write(32,'(A)')"#  The TRANSITION MATRIX ELEMENT P ="
         write(32,'(A)')"#  P(k,s,nm,+ or -) = 1/sqrt(2)[p_x(k,nm,s) + &
             (or -) i*p_y(k,nm,s)]"
         write(32,'(A)')"#  THE INTERBAND TRANSITION MATRIX p_x,y ="
         write(32,'(A)')"#  p_x,y(k,nm,s)=<psi(k,n,s)|-i*hbar*1/dx(y)| &
            psi(k,m,s)>"
!       write(32,'(A,4F16.6)')"# MAXVAL of SELECTIVITY at kx,ky,kz 
!    &(in reci)= ",(rvari3(i),i=1,4)
!       write(32,'(A,4F16.6)')"# MINVAL of SELECTIVITY at kx,ky,kz 
!    &(in reci)= ",(rvari4(i),i=1,4)
         write(32,'(A)')"# (cart) kx        ky        kz(A^-1) &
          BERRYKUBO (A^2),        Qxx       Qxy        Qyy    &
                    (recip)kx        ky        kz"

      else !BERRYCURVATURE
         write(32,'(A)')"# Chern Number is sum of                 &
             Berry Curvature over 1BZ"                                 
         write(32,'(A,F16.4)')"# Chern Number =   ",rvari           
         write(32,'(A,4F16.4)')"# MAXVAL of BERRYCURV at kx,ky,kz  &
              (in reci)= ",(rvari3(i),i=1,4)                           
         write(32,'(A,4F16.4)')"# MINVAL of BERRYCURV at kx,ky,kz  &
              (in reci)= ",(rvari4(i),i=1,4)                           
         write(32,'(A)')"# Berry Curvature (A^2) :                 &
              -Im[logPI_S(det(S(K_s,K_s+1)))]/dk^2"                    
         write(32,'(A)')"# (cart) kx        ky        kz(A^-1)     &
              Berry Curvature (A^2), Qxx    Qxy     Qyy    &
              (recip)kx        ky        kz"
      endif

      do ik=1,kext
         write(32,'(3F11.6,A,3F20.4,A,3F11.6)')(xrecivec(i,ik),i=1,3),&
             "     ",(xrvari(ik,i),i=1,3),"            ",                     &
             (xrecilat(i,ik),i=1,3)
      enddo

      close(32)
      return
      end subroutine write_result

      subroutine write_special_kpoint(b1,b2,b3)
         implicit real*8(a-h,o-z)
         dimension SKP(3),SP(3),b1(3),b2(3),b3(3)
         open(41,file='SKP.dat',status='unknown')
         do i=1,3;SKP(i)=( 0.                     * b1(i))+     &
                         ( 0.                     * b2(i))+     &
                         ( 0.                     * b3(i));enddo 
         write(41,'(3F12.6,A)')(SKP(i),i=1,3),"  # G"            
                                                                 
         do i=1,3;SKP(i)=( 2. / 3.                * b1(i))+     &
                         ( 1. / 3.                * b2(i))+     &
                         ( 0.                     * b3(i));enddo 
         write(41,'(3F12.6,A)')(SKP(i),i=1,3),"  # K1"           
         do i=1,3;SKP(i)=(-1. / 3.                * b1(i))+     &
                         ( 1. / 3.                * b2(i))+     &
                         ( 0.                     * b3(i));enddo 
         write(41,'(3F12.6,A)')(SKP(i),i=1,3),"  # K1"           
                                                                 
         do i=1,3;SKP(i)=( 1. / 3.                * b1(i)) +    &
                         ( 2. / 3.                * b2(i)) +    &
                         ( 0.                     * b3(i));enddo 
         write(41,'(3F12.6,A)')(SKP(i),i=1,3),"  # K2"           
         do i=1,3;SKP(i)=( 1. / 3.                * b1(i))+     &
                         (-1. / 3.                * b2(i))+     &
                         ( 0.                     * b3(i));enddo 
         write(41,'(3F12.6,A)')(SKP(i),i=1,3),"  # K2"           
                                                                 
         do i=1,3;SKP(i)=( 1. / 2.                * b1(i))+     &
                         ( 0.                     * b2(i))+     &
                         ( 0.                     * b3(i));enddo 
         write(41,'(3F12.6,A)')(SKP(i),i=1,3),"  # M1"           
                                                                 
         do i=1,3;SKP(i)=( 1. / 2.                * b1(i))+     &
                         ( 1. / 2.                * b2(i))+     &
                         ( 0.                     * b3(i));enddo 
         write(41,'(3F12.6,A)')(SKP(i),i=1,3),"  # M2"           
                                                                 
         do i=1,3;SKP(i)=( 0.                     * b1(i))+     &
                        ( 1. / 2.                * b2(i))+      &
                        ( 0.                     * b3(i));enddo
         write(41,'(3F12.6,A)')(SKP(i),i=1,3),"  # M3"

         close(41)
      end subroutine write_special_kpoint 
