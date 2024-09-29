!!$*  subroutine for computing berry curvature for the nn band in the given k-point using kubo formula
      !!subroutine kubo_quantum_metric(berrycurv_kubo, &
      subroutine kubo_quantum_metric(gxy, &
                b1,b2,b3,wklist,isp,                  &
                nband,ecut,ispinor,nplist,nbmax,npmax,&
                nk,nn,                                &
                nprocs,myrank,mpi_comm_earth)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      dimension nbmax(3),nplist(nk),wk(3),ener(nband)
      real*8    berrycurv_kubo(nk),b1(3),b2(3),b3(3),wklist(3,nk)
      real*8    berrycurv_kubo_(nk)
      ! quantum metric components
      real*8    gxy(nk),gxy_(nk)
      real*8    gxx(nk),gxx_(nk)
      real*8    gyy(nk),gyy_(nk)
      real*8    metertoang 
      ! descarded
      complex*16 ctrans_mtrx_left,ctrans_mtrx_right
      complex*16 ctrans_mtrx_left2,ctrans_mtrx_right2
      ! g_xy
      complex*16 mtrxy_plus_left,mtrxy_plus_right 
      complex*16 mtrxy_minus_left,mtrxy_minus_right
      complex*16 Pplus_Pminus_xy
      ! g_xx
      complex*16 mtrxx_plus_left, mtrxx_plus_right 
      complex*16 mtrxx_minus_left,mtrxx_minus_right
      complex*16 Pplus_Pminus_xx
      ! g_yy
      complex*16 mtryy_plus_left, mtryy_plus_right 
      complex*16 mtryy_minus_left,mtryy_minus_right
      complex*16 Pplus_Pminus_yy

      complex*16 cinter_mtrx_x,cinter_mtrx_y,cinter_mtrx_z
      complex*16 cinter_mtrx_x2,cinter_mtrx_y2,cinter_mtrx_z2
      complex*16 coeffn(npmax),coeffm(npmax)
      complex*8  coeff(npmax)
      complex*16 coeffnu(npmax),coeffnd(npmax)
      complex*16 coeffmu(npmax),coeffmd(npmax)
      integer :: ig(3,npmax)
!     data hbarkg/1.054571596E-34/ !h/2pi [kg m^2/s]
!     data xm/9.1093837015E-31/ ! electron mass [kg]
!     data eV/1.60217663E-19/ ! electron volt [kg c^2=kg m^2/s^2]
!     data ang/1.E-10/ ! Angstrom  [A = 10^-10 m]
!     data anginv/1.E+10/ ! Angstrom^-1  [A^-1 = 10^10 m^-1]
      data hbar/6.582119569e-16/ 
      data c/2.99792458E+08/ ! constant c [m/s]
      data em/0.510998950E+6/ ! electron mass (eV/c^2)
      data metertoang/1.E+10/ ! [1 m = 10^10 A]
      real*8     berrycurv_unit ! [A^4]  without b1, b1: A^-1
      integer :: nprocs, myrank, mpi_comm_earth, mpierr
      integer*4  ourjob(nprocs), ourjob_disp(0:nprocs-1)
   
      call mpi_job_distribution_chain(nk, ourjob, ourjob_disp, nprocs)

!     berrycurv_unit=(hbarkg**4)/((xm*eV)**2)  ! [m^4] without b1, b1: A^-1
!     berrycurv_unit=(hbarkg**4)/((xm*eV)**2)*(anginv**4)  ! [A^4] without b1, b1: A^-1
      berrycurv_unit=(hbar**4)/((em)**2)*(metertoang**4)*c**4  ! [A^4] without b1, b1: A^-1
!     data hbar/1./ !Here, I will set hbar = 1. for the simplicity
!     write(6,*) "start kubo routine..";stop
      berrycurv_kubo = 0.0d0 
      berrycurv_kubo_ = 0.0d0

      g_xy = 0.0d0
      g_xy_ = 0.0d0
      g_xx = 0.0d0
      g_xx_ = 0.0d0
      g_yy = 0.0d0
      g_yy_ = 0.0d0

!     do ik=1,nk
      do ik=sum(ourjob(1:myrank))+1, sum(ourjob(1:myrank+1))
          call ener_read(ener,isp,ik,nk,nband)      
          !write(6,*) "start ener routine..";stop
          !write(6,*) 'checking nn=1 energy',ener(1);stop
          do mm=1,nband      ! band index running from 1 to N except for m=n (the # of total bands = N)
               !if(mm.ne.nn) then 
               if(mm.ne.nn .and. abs(ener(nn)-ener(mm)).gt.epsilon(1d0)) then 
               cinter_mtrx_x=(0.,0.)
               cinter_mtrx_y=(0.,0.)
               cinter_mtrx_z=(0.,0.)

               cinter_mtrx_x2=(0.,0.)
               cinter_mtrx_y2=(0.,0.)
               cinter_mtrx_z2=(0.,0.)


               coeff=(0.,0.)
               coeffn=(0.,0.);coeffm=(0.,0.) ! between band indices n and m
               coeffnu=(0.,0.);coeffmu=(0.,0.)
               coeffnd=(0.,0.);coeffmd=(0.,0.)
               wk(:)=wklist(:,ik)
               np=nplist(ik)
               call plindx(ig,ncnt, ispinor,wk,b1,b2,b3,nbmax,np,ecut,npmax)
               !write(6,*)ik,nband,isp,nn,np;stop
               read(10,rec=(3+(ik-1)*(nband+1)+                      &
                          nk*(nband+1)*(isp-1)+nn))(coeff(i),i=1,np)  
               coeffn=coeff;coeff=(0.,0.)                             
               read(10,rec=(3+(ik-1)*(nband+1)+                      &
                          nk*(nband+1)*(isp-1)+mm))(coeff(i),i=1,np)  
               coeffm=coeff;coeff=(0.,0.)                             
!              write(6,*) coeffn,coeffm;stop                          
               do iplane=1,ncnt                                       
                  xkgx=(wk(1)+ig(1,iplane))*b1(1)+                     &
                       (wk(2)+ig(2,iplane))*b2(1)+                     &
                       (wk(3)+ig(3,iplane))*b3(1)                       
                  xkgy=(wk(1)+ig(1,iplane))*b1(2)+                     &
                       (wk(2)+ig(2,iplane))*b2(2)+                     &
                       (wk(3)+ig(3,iplane))*b3(2)
                  xkgz=(wk(1)+ig(1,iplane))*b1(3)+                     &
                       (wk(2)+ig(2,iplane))*b2(3)+                     &
                       (wk(3)+ig(3,iplane))*b3(3)


! computing P_x^{nm} and P_y^{nm}
! soc  case:
                  if(ispinor .eq. 2) then
                     coeffnu(iplane)=coeffn(iplane)
                     coeffnd(iplane)=coeffn(iplane+ncnt)
                     coeffmu(iplane)=coeffm(iplane)
                     coeffmd(iplane)=coeffm(iplane+ncnt)

                     cinter_mtrx_x=cinter_mtrx_x+                     &
                         conjg(coeffnu(iplane))*xkgx*coeffmu(iplane)+ & 
                         conjg(coeffnd(iplane))*xkgx*coeffmd(iplane)   
                     cinter_mtrx_y=cinter_mtrx_y+                     &
                         conjg(coeffnu(iplane))*xkgy*coeffmu(iplane)+ &
                         conjg(coeffnd(iplane))*xkgy*coeffmd(iplane)
                     cinter_mtrx_z=cinter_mtrx_z+                     &
                         conjg(coeffnu(iplane))*xkgz*coeffmu(iplane)+ & 
                         conjg(coeffnd(iplane))*xkgz*coeffmd(iplane)   

!                  cinter_mtrx_x2=cinter_mtrx_x2+                     &
!                         conjg(coeffmu(iplane))*xkgx*coeffnu(iplane)+& 
!                         conjg(coeffmd(iplane))*xkgx*coeffnd(iplane)  
!                  cinter_mtrx_y2=cinter_mtrx_y2+                      &
!                          conjg(coeffmu(iplane))*xkgy*coeffnu(iplane)+&
!                          conjg(coeffmd(iplane))*xkgy*coeffnd(iplane)
!                  cinter_mtrx_z2=cinter_mtrx_z2+                      &
!                          conjg(coeffmu(iplane))*xkgz*coeffnu(iplane)+&
!                          conjg(coeffmd(iplane))*xkgz*coeffnd(iplane)

! nsoc  case:
                  else if(ispinor .eq. 1) then
                      !write(6,*)coeffn(iplane), xkgx
                     cinter_mtrx_x=cinter_mtrx_x+                   &
                           conjg(coeffn(iplane))*xkgx*coeffm(iplane)  ! P_x^{nm}=<u_{n,k}|p_x|u_{m,k}>
                     cinter_mtrx_y=cinter_mtrx_y+                     &
                           conjg(coeffn(iplane))*xkgy*coeffm(iplane) ! P_y^{nm}=<u_{n,k}|p_y|u_{m,k}>
                     cinter_mtrx_z=cinter_mtrx_z+                     &
                           conjg(coeffn(iplane))*xkgz*coeffm(iplane) ! P_z^{nm}=<u_{n,k}|p_z|u_{m,k}>


!                    cinter_mtrx_x2=cinter_mtrx_x2+                 &
!                           conjg(coeffm(iplane))*xkgx*coeffn(iplane)  ! P_x^{mn}=<u_{m,k}|p_x|u_{n,k}>
!                    cinter_mtrx_y2=cinter_mtrx_y2+                  &
!                           conjg(coeffm(iplane))*xkgy*coeffn(iplane) ! P_y^{mn}=<u_{m,k}|p_y|u_{n,k}>
!                    cinter_mtrx_z2=cinter_mtrx_z2+                  &
!                           conjg(coeffm(iplane))*xkgz*coeffn(iplane) ! P_z^{mn}=<u_{m,k}|p_z|u_{n,k}>
                  endif
               enddo ! iplane loop end
       !stop
               ! g_xy
               mtrxy_plus_left  = cinter_mtrx_x + (0.,1.)*cinter_mtrx_y
               mtrxy_plus_right = conjg(cinter_mtrx_x) + &
                                         (0.,1.)*conjg(cinter_mtrx_y)

               mtrxy_minus_left = cinter_mtrx_x - (0.,1.)*cinter_mtrx_y
               mtrxy_minus_right  = conjg(cinter_mtrx_x) -&
                                         (0.,1.)*conjg(cinter_mtrx_y)

               Pplus_Pminus_xy = (mtrxy_plus_left*mtrxy_plus_right - &
                                  mtrxy_minus_left*mtrxy_minus_right)

               ! g_xx
               mtrxx_plus_left  = cinter_mtrx_x + (0.,1.)*cinter_mtrx_x
               mtrxx_plus_right = conjg(cinter_mtrx_x) + &
                                         (0.,1.)*conjg(cinter_mtrx_x)

               mtrxx_minus_left = cinter_mtrx_x - (0.,1.)*cinter_mtrx_x
               mtrxx_minus_right  = conjg(cinter_mtrx_x) -&
                                         (0.,1.)*conjg(cinter_mtrx_x)

               Pplus_Pminus_xx = (mtrxx_plus_left*mtrxx_plus_right - &
                                  mtrxx_minus_left*mtrxx_minus_right)

               ! g_yy
               mtryy_plus_left  = cinter_mtrx_y + (0.,1.)*cinter_mtrx_y
               mtryy_plus_right = conjg(cinter_mtrx_y) + &
                                         (0.,1.)*conjg(cinter_mtrx_y)

               mtryy_minus_left = cinter_mtrx_y - (0.,1.)*cinter_mtrx_y
               mtryy_minus_right  = conjg(cinter_mtrx_y) -&
                                         (0.,1.)*conjg(cinter_mtrx_y)

               Pplus_Pminus_yy = (mtryy_plus_left*mtryy_plus_right - &
                                  mtryy_minus_left*mtryy_minus_right)


! quantum metric method1
             ! gxy
             g_xy(ik) = g_xy(ik) + (0.,1.)* &
                (Pplus_Pminus_xy /(ener(nn)-ener(mm))**2)           &
            *berrycurv_unit

             ! gxx
             g_xx(ik) = g_xx(ik) + (0.,1.)* &
                (Pplus_Pminus_xx /(ener(nn)-ener(mm))**2)           &
            *berrycurv_unit

             ! gyy
             g_yy(ik) = g_yy(ik) + (0.,1.)* &
                (Pplus_Pminus_yy /(ener(nn)-ener(mm))**2)           &
            *berrycurv_unit

! quantum metric method2
!              berrycurv_kubo(ik) = berrycurv_kubo(ik) + real(  &
!             (cinter_mtrx_x*cinter_mtrx_y2)                    &
!              /(ener(nn)-ener(mm))**2 )*berrycurv_unit 

               endif ! if mm is not equal to nn
           enddo ! band index mm loop end 

      write(6,'(2A,I4,4F16.6)')"# IK, K(reci),", &
             "gxx (A^2, Kubo) : ", &
                               ik,wk,g_xx(ik)
      write(6,'(2A,I4,4F16.6)')"# IK, K(reci),", &
             "gxy (A^2, Kubo) : ", &
                               ik,wk,g_xy(ik)
      write(6,'(2A,I4,4F16.6)')"# IK, K(reci),", &
             "gyy (A^2, Kubo) : ", &
                               ik,wk,g_yy(ik)

      enddo !ik loop end

#ifdef MPI_USE
      call MPI_ALLREDUCE(g_xy,g_xy_,nk, MPI_REAL8,&
                         MPI_SUM, mpi_comm_earth, mpierr)
      call MPI_ALLREDUCE(g_xx,g_xx_,nk, MPI_REAL8,&
                         MPI_SUM, mpi_comm_earth, mpierr)
      call MPI_ALLREDUCE(g_yy,g_yy_,nk, MPI_REAL8,&
                         MPI_SUM, mpi_comm_earth, mpierr)
      g_xx  = g_xx_
      g_xy  = g_xy_
      g_yy  = g_yy_
#endif

      return
      end subroutine kubo_quantum_metric



