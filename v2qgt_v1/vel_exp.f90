!!$*  subroutine for computing velocity expectation value for state psi(n,k)
   subroutine vel_expectation(a1,a2,a3,b1,b2,b3, &
        kperiod,nbmax,npmax,ecut,ispinor,ispin,nband,ne,nk,wklist, &
        nplist,ni,nj,irecl,filename,foname)
      implicit real*8 (a-h, o-z)
      complex*8  coeff(npmax)
      complex*16 coeffi(npmax),coeffj(npmax)
      complex*16 coeffiu(npmax),coeffid(npmax)
      complex*16 coeffju(npmax),coeffjd(npmax)
      complex*16 vel_x,vel_y
      real*8  vel_expt(2,nk,ispin),wklist(3,nk) 
      real*8  recivec(3,nk),recilat(3,nk)
      real*8  xvel_expt(2,kperiod*2*kperiod*2*nk,ispin)
      real*8  xrecivec(3,kperiod*2*kperiod*2*nk)
      real*8  xrecilat(3,kperiod*2*kperiod*2*nk)
      real*8  vel_x_expt_max(4),vel_y_expt_max(4)
      real*8  vel_x_expt_min(4),vel_y_expt_min(4)
      dimension a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),a2xa3(3),sumkg(3)
      dimension wk(3),nbmax(3),xb(3)
      integer ig(3,npmax),nplist(nk)
      integer ni,nj,ne,nk,nband,np,npmax,kperiod,ispin,irecl
      character*256 filename,foname,fonameo
      data c/0.262465831d0/ ! constant c = 2m/hbar**2 [1/eV Ang^2]
      data hbar/6.58211928E-16/ !h/2pi [eV * s]
      data xm/0.510998910E+6/ ! electron mass (eV/c^2)
      data anginv/1.0E+10/     ! inverse angstrom (1/Ang)
      pi=4.*atan(1.)
      do isp=1,ispin
         do ik=1, nk
            coeff=(0.,0.)
            coeffi=(0.,0.);coeffj=(0.,0.)
            coeffiu=(0.,0.);coeffid=(0.,0.)
            coeffju=(0.,0.);coeffjd=(0.,0.)
            wk(:)=wklist(:,ik)
            np=nplist(ik)
            if(ni .ne. nj)then
             call plindx(ig,ncnt, ispinor,wk,b1,b2,b3,nbmax,np,ecut,npmax)
             read(10,rec=(3+(ik-1)*(nband+1)+ &
                         nk*(nband+1)*(isp-1)+ni))(coeff(i),i=1,np)
             coeffi=coeff;coeff=(0.,0.)
             read(10,rec=(3+(ik-1)*(nband+1)+ &
                         nk*(nband+1)*(isp-1)+nj))(coeff(i),i=1,np)
             coeffj=coeff;coeff=(0.,0.)
             else if(ni .eq. nj)then
              call plindx(ig,ncnt, ispinor,wk,b1,b2,b3,nbmax,np,ecut,npmax)
              read(10,rec=(3+(ik-1)*(nband+1)+ &
                         nk*(nband+1)*(isp-1)+ni))(coeff(i),i=1,np)
              coeffi=coeff;coeff=(0.,0.)
            endif
            vel_x=(0.,0.);vel_y=(0.,0.)
            do iplane=1,ncnt
             xkgx=(wk(1)+ig(1,iplane))*b1(1)+ &
                  (wk(2)+ig(2,iplane))*b2(1)+ &
                  (wk(3)+ig(3,iplane))*b3(1)
             xkgy=(wk(1)+ig(1,iplane))*b1(2)+ &
                  (wk(2)+ig(2,iplane))*b2(2)+ &
                  (wk(3)+ig(3,iplane))*b3(2)
             if(ispinor .eq. 2) then
              coeffiu(iplane)=coeffi(iplane)
              coeffid(iplane)=coeffi(iplane+ncnt)
              if(ni .eq. nj)then
               coeffju(iplane)=coeffi(iplane)
               coeffjd(iplane)=coeffi(iplane+ncnt)
               else
               coeffju(iplane)=coeffj(iplane)
               coeffjd(iplane)=coeffj(iplane+ncnt)
              endif
              ! -i*hbar/m <psi|d/dx|psi>
              vel_x=vel_x+     &
              anginv*hbar/xm*(conjg(coeffiu(iplane))*xkgx*coeffju(iplane)+ &
                             conjg(coeffid(iplane))*xkgx*coeffjd(iplane))
              ! -i*hbar/m <psi|d/dy|psi> 
              vel_y=vel_y+     &
              anginv*hbar/xm*(conjg(coeffiu(iplane))*xkgy*coeffju(iplane)+ &
                              conjg(coeffid(iplane))*xkgy*coeffjd(iplane))
              else if (ispinor .eq. 1) then
               coeffiu(iplane)=coeffi(iplane)
               if(ni .eq. nj)then
                coeffju(iplane)=coeffi(iplane)
                else
                 coeffju(iplane)=coeffj(iplane)
               endif
               vel_x=vel_x+ &
              anginv*hbar/xm*conjg(coeffiu(iplane))*xkgx*coeffju(iplane)
               vel_y=vel_y+ &
              anginv*hbar/xm*conjg(coeffiu(iplane))*xkgy*coeffju(iplane)
             endif ! ispinor
            enddo  ! iplane
            vel_expt(1,ik,isp)=real(vel_x,8)
            vel_expt(2,ik,isp)=real(vel_y,8)
            write(6,'(A,I4,5F11.6)')"# IK, K(reci), VEL_EXPT(x,y) : ", &
                                  ik,wk,(vel_expt(i,ik,isp),i=1,2)
         enddo  !ik

         if (ik .eq. 1)then
          vel_x_expt_max(4)=vel_expt(1,ik,isp)
          vel_x_expt_max(1)=wklist(1,ik)
          vel_x_expt_max(2)=wklist(2,ik)
          vel_x_expt_max(3)=wklist(3,ik)
          vel_x_expt_min(4)=vel_expt(1,ik,isp)
          vel_x_expt_min(1)=wklist(1,ik)
          vel_x_expt_min(2)=wklist(2,ik)
          vel_x_expt_min(3)=wklist(3,ik)
          else if(ik.ge.2.and.vel_expt(1,ik,isp).ge.vel_x_expt_max(4))then
           vel_x_expt_max(4)=vel_expt(1,ik,isp)
           vel_x_expt_max(1)=wklist(1,ik)
           vel_x_expt_max(2)=wklist(2,ik)
           vel_x_expt_max(3)=wklist(3,ik)
          else if(ik.ge.2.and.vel_expt(1,ik,isp).le.vel_x_expt_min(4))then
           vel_x_expt_min(4)=vel_expt(1,ik,isp)
           vel_x_expt_min(1)=wklist(1,ik)
           vel_x_expt_min(2)=wklist(2,ik)
           vel_x_expt_min(3)=wklist(3,ik)
         endif

         if (ik .eq. 1)then
          vel_y_expt_max(4)=vel_expt(2,ik,isp)
          vel_y_expt_max(1)=wklist(1,ik)
          vel_y_expt_max(2)=wklist(2,ik)
          vel_y_expt_max(3)=wklist(3,ik)
          vel_y_expt_min(4)=vel_expt(2,ik,isp)
          vel_y_expt_min(1)=wklist(1,ik)
          vel_y_expt_min(2)=wklist(2,ik)
          vel_y_expt_min(3)=wklist(3,ik)
          else if(ik.ge.2.and.vel_expt(2,ik,isp).ge.vel_y_expt_max(4))then
           vel_y_expt_max(4)=vel_expt(2,ik,isp)
           vel_y_expt_max(1)=wklist(1,ik)
           vel_y_expt_max(2)=wklist(2,ik)
           vel_y_expt_max(3)=wklist(3,ik)
          else if(ik.ge.2.and.vel_expt(2,ik,isp).le.vel_y_expt_min(4))then
           vel_y_expt_min(4)=vel_expt(2,ik,isp)
           vel_y_expt_min(1)=wklist(1,ik)
           vel_y_expt_min(2)=wklist(2,ik)
           vel_y_expt_min(3)=wklist(3,ik)
         endif
         do ik=1,nk
          do j=1,3
           recivec(j,ik)=wklist(1,ik)*b1(j)+ &
                         wklist(2,ik)*b2(j)+ & 
                         wklist(3,ik)*b3(j)
           recilat(j,ik)=wklist(j,ik)
          enddo
         enddo

         kk=0   ! extend berry curvature distribution over extended BZ 
         do ib2=-1*(kperiod-1)+1,kperiod
          do ib1=-1*(kperiod-1)+1,kperiod  ! you may adjust these values as you wish..
           do ik=1,nk
            kk=kk+1
            xrecivec(1,kk)=recivec(1,ik) +(ib1-1)*b1(1)+(ib2-1)*b2(1)
            xrecivec(2,kk)=recivec(2,ik) +(ib1-1)*b1(2)+(ib2-1)*b2(2)
            xrecivec(3,kk)=recivec(3,ik) +(ib1-1)*b1(3)+(ib2-1)*b2(3)
            xrecilat(1,kk)=recilat(1,ik) +(ib1-1)
            xrecilat(2,kk)=recilat(2,ik) +(ib2-1)
            xrecilat(3,kk)=recilat(3,ik)
            xvel_expt(1,kk,isp)=vel_expt(1,ik,isp)
            xvel_expt(2,kk,isp)=vel_expt(2,ik,isp)
            kext=kk
           enddo
          enddo
         enddo

!$$*  sorting k-points and the corresponding optical selectivity on the
!     periodically repeated data grid
         write(6,*)" "
         write(6,'(A)')"# SORTING K-grids..."
         do k=kext-1,1,-1  ! sorting kx
          do j=1,k
           if(xrecivec(1,j+1) .gt. xrecivec(1,j))then
            xb(:)=xrecivec(:,j)
            xrecivec(:,j)=xrecivec(:,j+1)
            xrecivec(:,j+1)=xb(:)
            xb(:)=xrecilat(:,j)
            xrecilat(:,j)=xrecilat(:,j+1)
            xrecilat(:,j+1)=xb(:)
            xtemp=xvel_expt(1,j,isp)
            xvel_expt(1,j,isp)=xvel_expt(1,j+1,isp)
            xvel_expt(1,j+1,isp)=xtemp
            xtemp=xvel_expt(2,j,isp)
            xvel_expt(2,j,isp)=xvel_expt(2,j+1,isp)
            xvel_expt(2,j+1,isp)=xtemp
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
             xtemp=xvel_expt(1,j,isp)
             xvel_expt(1,j,isp)=xvel_expt(1,j+1,isp)
             xvel_expt(1,j+1,isp)=xtemp
             xtemp=xvel_expt(2,j,isp)
             xvel_expt(2,j,isp)=xvel_expt(2,j+1,isp)
             xvel_expt(2,j+1,isp)=xtemp
             endif
            endif
           enddo
          enddo
 
         if(isp .eq. 1 .and. ispinor .eq. 1 .and. ispin.eq.2) then
          write(fonameo,'(A,A)')TRIM(foname),'.UP.dat'
          else if (isp .eq. 2 .and. ispinor .eq. 1 .and. ispin.eq.2)then
           write(fonameo,'(A,A)')TRIM(foname),'.DN.dat'
          else if (isp .eq. 1 .and. ispinor .eq. 2) then
           write(fonameo,'(A,A)')TRIM(foname),'.dat'
          else if (isp .eq. 1 .and. ispinor .eq. 1 .and. ispin.eq.1)then
           write(fonameo,'(A,A)')TRIM(foname),'.dat'
         endif
         open(61,file=fonameo,status='unknown')
         write(61,'(A,A)')"# File reading... : ",filename
         write(61,'(A,I9)')"# TOTAL RECORD LENGTH = ",irecl
         if (ispinor .eq. 2)then
          write(61,'(A,I6,A)')"# ISPIN            : ",ispin, &
                             " (LSORBIT = .TRUE.)"
          else
           write(61,'(A,I6,A)')"# ISPIN            : ",ispin, &
                              " (LSORBIT = .FALSE.)"
         endif
         write(61,'(A,F11.4)')  "# ENCUT (eV)       : ",ecut
         write(61,'(A,I6)')     "# NKPOINT          : ",nk
         write(61,'(A,I6)')     "# NBANDS           : ",nband
         write(61,'(A,3F13.6)') "# RECIVEC B1 (A^-1): ",(b1(i),i=1,3)
         write(61,'(A,3F13.6)') "# RECIVEC B2       : ",(b2(i),i=1,3)
         write(61,'(A,3F13.6)') "# RECIVEC B3       : ",(b3(i),i=1,3)
         write(61,*)" "

         write(61,'(A,I4)')"# VEOLOCITY EXPECTATION VALUE of BAND:",ni
         write(61,'(A)')"# <v(n,k)>= 1/hbar dE(k)/dk = &
         1/m_e<psi(n,k)|p|psi(n,k)>, p=-i*hbar*d/dx,m_e=elect_rest_mass"
         write(61,'(A,4F16.6)')"# MAXVAL of VEL_EXPT <v_x>  &
         (in reci)= ",(vel_x_expt_max(i),i=1,4)
         write(61,'(A,4F16.6)')"# MINVAL of VEL_EXPT <v_x>  &
         (in reci)= ",(vel_x_expt_min(i),i=1,4)
         write(61,'(A,4F16.6)')"# MAXVAL of VEL_EXPT <v_y> &
         (in reci)= ",(vel_y_expt_max(i),i=1,4)
         write(61,'(A,4F16.6)')"# MINVAL of VEL_EXPT <v_y> &
         (in reci)= ",(vel_y_expt_min(i),i=1,4)
         write(61,'(A)')"# (cart) kx     ky     kz(A^-1) &
         vel_expt(vx(n,k), vy(n,k))(m/s)  (recip)kx      ky      kz"
         do ik=1,kext
           write(61,'(3F11.6,A,2F10.6,A,3F11.6)')(xrecivec(i,ik),i=1,3), &
            "   ",(xvel_expt(i,ik,isp),i=1,2),"       ", &
             (xrecilat(i,ik),i=1,3)
         enddo
         close(61)

      enddo   !isp

!     write(6,'(A)')"# DONE! "
!     do isp=1, ispin
!      if(isp .eq. 1 .and. ispinor .eq. 1 .and. ispin .eq. 2) then
!       write(6,'(A,A,A)')"#  Results are summarized in ",TRIM(foname),
!    &                    ".UP.dat for spin-1"
!       else if (isp.eq.2 .and. ispinor.eq.1 .and. ispin.eq.2)then
!        write(6,'(A,A,A)')"#  Results are summarized in ",
!    &                     TRIM(foname),".DN.dat for spin-2"
!       else if (isp .eq. 1 .and. ispinor .eq. 2) then
!        write(6,'(A,A,A)')"#  Results are summarized in ",
!    &                     TRIM(foname),".dat"
!       else if (isp.eq.1 .and. ispinor.eq.1 .and. ispin.eq.1) then
!        write(6,'(A,A,A)')"#  Results are summarized in ",
!    &                     TRIM(foname),".dat"
!      endif
!     enddo
      close(10)
      return
      end subroutine vel_expectation
