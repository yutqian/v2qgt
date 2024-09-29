!!$*  subroutine for computing berry curvature for the nn band in the given k-point using kubo formula
      subroutine kubo_berry_curvature1(berrycurv_kubo, &
                b1,b2,b3,wklist,isp,                  &
                nband,ecut,ispinor,nplist,nbmax,npmax,&
                nk,nn,                                &
                nprocs,myrank,mpi_comm_earth)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      dimension nbmax(3),nplist(nk),wk(3),ener(nband)
      real*8    berrycurv_kubo(nk),b1(3),b2(3),b3(3),wklist(3,nk)
      real*8    berrycurv_kubo_(nk)
      real*8    metertoang 
      complex*16 ctrans_mtrx_left,ctrans_mtrx_right
      complex*16 cinter_mtrx_x,cinter_mtrx_y
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

! computing P_x^{nm} and P_y^{nm}
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
                  else if(ispinor .eq. 1) then
                      !write(6,*)coeffn(iplane), xkgx
                     cinter_mtrx_x=cinter_mtrx_x+                   &
                           conjg(coeffn(iplane))*xkgx*coeffm(iplane)  ! P_x^{nm}=<u_{n,k}|p_x|u_{m,k}>
                     cinter_mtrx_y=cinter_mtrx_y+                     &
                           conjg(coeffn(iplane))*xkgy*coeffm(iplane) ! P_y^{nm}=<u_{n,k}|p_y|u_{m,k}>
                  endif
               enddo ! iplane loop end
       !stop
               ctrans_mtrx_left  = cinter_mtrx_x + (0.,1.)*cinter_mtrx_y
               ctrans_mtrx_right = cinter_mtrx_x - (0.,1.)*cinter_mtrx_y
               !write(6,*) ctrans_mtrx_left
               !write(6,*) ener(nn)-ener(mm),1./(ener(nn)-ener(mm))**2

               !! method 1
              berrycurv_kubo(ik) = berrycurv_kubo(ik)                 &
             -((abs(ctrans_mtrx_left))**2-(abs(ctrans_mtrx_right))**2)&
              /(ener(nn)-ener(mm))**2                                 &
              *berrycurv_unit


               endif ! if mm is not equal to nn
           enddo ! band index mm loop end 
      write(6,'(2A,I4,4F16.6)')"# IK, K(reci),", &
             "Berry Curvature (A^2, Kubo) : ", &
                               ik,wk,berrycurv_kubo(ik)
      enddo !ik loop end

#ifdef MPI_USE
      call MPI_ALLREDUCE(berrycurv_kubo,berrycurv_kubo_,nk, MPI_REAL8,&
                         MPI_SUM, mpi_comm_earth, mpierr)
      berrycurv_kubo  = berrycurv_kubo_
#endif

      return
      end subroutine kubo_berry_curvature1



!!$*  subroutine for computing berry curvature for the nn band in the given k-point using kubo formula
      subroutine kubo_berry_curvature2(berrycurv_kubo, &
                b1,b2,b3,wklist,isp,                  &
                nband,ecut,ispinor,nplist,nbmax,npmax,&
                nk,nn,                                &
                nprocs,myrank,mpi_comm_earth)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      dimension nbmax(3),nplist(nk),wk(3),ener(nband)
      real*8    berrycurv_kubo(nk),b1(3),b2(3),b3(3),wklist(3,nk)
      real*8    berrycurv_kubo_(nk)
      real*8    metertoang 
      complex*16 ctrans_mtrx_left,ctrans_mtrx_right
      complex*16 cinter_mtrx_x,cinter_mtrx_y
      complex*16 cinter_mtrx_x2,cinter_mtrx_y2
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

               cinter_mtrx_x2=(0.,0.)
               cinter_mtrx_y2=(0.,0.)


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
 
                   cinter_mtrx_x2=cinter_mtrx_x2+                     &
                          conjg(coeffmu(iplane))*xkgx*coeffnu(iplane)+& 
                          conjg(coeffmd(iplane))*xkgx*coeffnd(iplane)  
                   cinter_mtrx_y2=cinter_mtrx_y2+                      &
                           conjg(coeffmu(iplane))*xkgy*coeffnu(iplane)+&
                           conjg(coeffmd(iplane))*xkgy*coeffnd(iplane)

! nsoc  case:
                  else if(ispinor .eq. 1) then
                      !write(6,*)coeffn(iplane), xkgx
                     cinter_mtrx_x=cinter_mtrx_x+                   &
                           conjg(coeffn(iplane))*xkgx*coeffm(iplane)  ! P_x^{nm}=<u_{n,k}|p_x|u_{m,k}>
                     cinter_mtrx_y=cinter_mtrx_y+                     &
                           conjg(coeffn(iplane))*xkgy*coeffm(iplane) ! P_y^{nm}=<u_{n,k}|p_y|u_{m,k}>


                     cinter_mtrx_x2=cinter_mtrx_x2+                 &
                            conjg(coeffm(iplane))*xkgx*coeffn(iplane)  ! P_x^{mn}=<u_{n,k}|p_x|u_{m,k}>
                     cinter_mtrx_y2=cinter_mtrx_y2+                  &
                            conjg(coeffm(iplane))*xkgy*coeffn(iplane) ! P_y^{mn}=<u_{n,k}|p_y|u_{m,k}>
                  endif
               enddo ! iplane loop end
       !stop
!               ctrans_mtrx_left  = cinter_mtrx_x + (0.,1.)*cinter_mtrx_y
!               ctrans_mtrx_right = cinter_mtrx_x - (0.,1.)*cinter_mtrx_y
               !write(6,*) ctrans_mtrx_left
               !write(6,*) ener(nn)-ener(mm),1./(ener(nn)-ener(mm))**2


              !! method 2
!              berrycurv_kubo(ik) = berrycurv_kubo(ik) -2*aimag(  &
!             (cinter_mtrx_x*cinter_mtrx_y2)                      &
!              /(ener(nn)-ener(mm))**2 )                          &
!             *berrycurv_unit

! quantum metric
              berrycurv_kubo(ik) = berrycurv_kubo(ik) + real(  &
             (cinter_mtrx_x*cinter_mtrx_y2)                    &
              /(ener(nn)-ener(mm))**2 )*berrycurv_unit 

               endif ! if mm is not equal to nn
           enddo ! band index mm loop end 
      write(6,'(2A,I4,4F16.6)')"# IK, K(reci),", &
             "Berry Curvature (A^2, Kubo) : ", & 
                               ik,wk,berrycurv_kubo(ik)
      enddo !ik loop end

#ifdef MPI_USE
      call MPI_ALLREDUCE(berrycurv_kubo,berrycurv_kubo_,nk, MPI_REAL8,&
                         MPI_SUM, mpi_comm_earth, mpierr)
      berrycurv_kubo  = berrycurv_kubo_
#endif

      return
      end subroutine kubo_berry_curvature2




!!$*  subroutine for computing berry curvature for the nn band in the given k-point using kubo formula
      subroutine kubo_berry_curvature_tot(berrycurv_kubo, &
                b1,b2,b3,wklist,isp,                  &
                nband,ecut,ispinor,nplist,nbmax,npmax,&
                nk,nn,nmax,iqm,                       &
                nprocs,myrank,mpi_comm_earth)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      dimension nbmax(3),nplist(nk),wk(3),ener(nband)
      real*8    berrycurv_kubo(nk,3),b1(3),b2(3),b3(3),wklist(3,nk)
      real*8    berrycurv_kubo_(nk,3)
      real*8    berrycurv_kubo2(nk)

!      real*8    quan_metric(nk), quan_metric_(nk)
      real*8    metertoang 
!     complex*16 ctrans_mtrx_left,ctrans_mtrx_right
!     complex*16 ctrans_mtrx_left2,ctrans_mtrx_right2

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
      complex*16 quan_geom_xy,quan_geom_xx,quan_geom_yy
   
      call mpi_job_distribution_chain(nk, ourjob, ourjob_disp, nprocs)

!     berrycurv_unit=(hbarkg**4)/((xm*eV)**2)  ! [m^4] without b1, b1: A^-1
!     berrycurv_unit=(hbarkg**4)/((xm*eV)**2)*(anginv**4)  ! [A^4] without b1, b1: A^-1
      berrycurv_unit=(hbar**4)/((em)**2)*(metertoang**4)*c**4  ! [A^4] without b1, b1: A^-1
!     data hbar/1./ !Here, I will set hbar = 1. for the simplicity
!     write(6,*) "start kubo routine..";stop
      berrycurv_kubo = 0.0d0 
      berrycurv_kubo_ = 0.0d0
      quan_geom_xy = 0.0d0
      quan_geom_xx = 0.0d0
      quan_geom_yy = 0.0d0

!     do ik=1,nk
      do ik=sum(ourjob(1:myrank))+1, sum(ourjob(1:myrank+1))
          call ener_read(ener,isp,ik,nk,nband)      
          !write(6,*) "start ener routine..";stop
          !write(6,*) 'checking nn=1 energy',ener(1);stop
          !do mm=1,nband      ! band index running from 1 to N except for m=n (the # of total bands = N)
          do mm=nmax+1,nband      ! band index running from 1 to N except for m=n (the # of total bands = N)
!               if(mm.ne.nn) then 
                if(mm.ne.nn .and. abs(ener(nn)-ener(mm)).gt.epsilon(1d0)) then 
!               if(mm.ne.nn .and. &
!                     abs(ener(nn)-ener(mm)).gt. (0.00001)) then 
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
               call plindx(ig,ncnt, ispinor,wk,b1,b2,b3,nbmax,np,& 
                             ecut,npmax)
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
!                 xkgz=(wk(1)+ig(1,iplane))*b1(3)+                     &
!                      (wk(2)+ig(2,iplane))*b2(3)+                     &
!                      (wk(3)+ig(3,iplane))*b3(3)

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
!                    cinter_mtrx_z=cinter_mtrx_z+                     &
!                        conjg(coeffnu(iplane))*xkgz*coeffmu(iplane)+ & 
!                        conjg(coeffnd(iplane))*xkgz*coeffmd(iplane)   

                   cinter_mtrx_x2=cinter_mtrx_x2+                     &
                          conjg(coeffmu(iplane))*xkgx*coeffnu(iplane)+& 
                          conjg(coeffmd(iplane))*xkgx*coeffnd(iplane)  
                   cinter_mtrx_y2=cinter_mtrx_y2+                      &
                           conjg(coeffmu(iplane))*xkgy*coeffnu(iplane)+&
                           conjg(coeffmd(iplane))*xkgy*coeffnd(iplane)
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
!                    cinter_mtrx_z=cinter_mtrx_z+                     &
!                          conjg(coeffn(iplane))*xkgz*coeffm(iplane) ! P_z^{nm}=<u_{n,k}|p_z|u_{m,k}>


                     cinter_mtrx_x2=cinter_mtrx_x2+                 &
                            conjg(coeffm(iplane))*xkgx*coeffn(iplane)  ! P_x^{mn}=<u_{m,k}|p_x|u_{n,k}>
                     cinter_mtrx_y2=cinter_mtrx_y2+                  &
                            conjg(coeffm(iplane))*xkgy*coeffn(iplane) ! P_y^{mn}=<u_{m,k}|p_y|u_{n,k}>
!                    cinter_mtrx_z2=cinter_mtrx_z2+                  &
!                           conjg(coeffm(iplane))*xkgz*coeffn(iplane) ! P_z^{mn}=<u_{m,k}|p_z|u_{n,k}>
                  endif
               enddo ! iplane loop end
       !stop


               quan_geom_xx = ((cinter_mtrx_x*cinter_mtrx_x2)/ & 
                            (ener(nn)-ener(mm))**2 )*berrycurv_unit
               quan_geom_xy = ((cinter_mtrx_x*cinter_mtrx_y2)/ & 
                            (ener(nn)-ener(mm))**2 )*berrycurv_unit
               quan_geom_yy = ((cinter_mtrx_y*cinter_mtrx_y2)/ & 
                            (ener(nn)-ener(mm))**2 )*berrycurv_unit

               if(iqm .eq. 1) then
                 berrycurv_kubo(ik,1)=berrycurv_kubo(ik,1)+ & 
                                             real(quan_geom_xx)
                 berrycurv_kubo(ik,2)=berrycurv_kubo(ik,2)+ &
                                             real(quan_geom_xy)
                 berrycurv_kubo(ik,3)=berrycurv_kubo(ik,3)+ &
                                             real(quan_geom_yy)
               else 
                berrycurv_kubo(ik,1)=berrycurv_kubo(ik,1) & 
                                                -2*aimag(quan_geom_xx)
                berrycurv_kubo(ik,2)=berrycurv_kubo(ik,2) & 
                                                -2*aimag(quan_geom_xy)
                berrycurv_kubo(ik,3)=berrycurv_kubo(ik,3) & 
                                                -2*aimag(quan_geom_yy)
               endif           


!               berrycurv_kubo2(ik) = berrycurv_kubo(ik,2)
               endif ! if mm is not equal to nn
           enddo ! band index mm loop end 
      write(6,'(2A,I4,6F16.6)')"# IK, K(reci),", &
             "Berry Curvature (A^2, Kubo) : Qxx      Qxy      Qyy   ", &
                             ik,wk,berrycurv_kubo(ik,:)
      enddo !ik loop end
!
#ifdef MPI_USE
      call MPI_ALLREDUCE(berrycurv_kubo,berrycurv_kubo_,nk, MPI_REAL8,&
                         MPI_SUM, mpi_comm_earth, mpierr)
      berrycurv_kubo  = berrycurv_kubo_

#endif

      return
      end subroutine kubo_berry_curvature_tot




!!$*  subroutine for computing berry curvature for the nn band in the given k-point using kubo formula
      subroutine kubo_berry_curvature(berrycurv_kubo, &
                b1,b2,b3,wklist,isp,                  &
                nband,ecut,ispinor,nplist,nbmax,npmax,&
                nk,nn,nmax,iqm,                       &
                nprocs,myrank,mpi_comm_earth)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      dimension nbmax(3),nplist(nk),wk(3),ener(nband)
      real*8    berrycurv_kubo(nk,3),b1(3),b2(3),b3(3),wklist(3,nk)
      real*8    berrycurv_kubo_(nk,3)

      real*8    berrycurv_kubo2(nk)
      real*8    berrycurv_kubo2_(nk)
!     real*8    quan_metric(nk), quan_metric_(nk)
      real*8    metertoang 
!     complex*16 ctrans_mtrx_left,ctrans_mtrx_right
!     complex*16 ctrans_mtrx_left2,ctrans_mtrx_right2

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
      complex*16 quan_geom_xy,quan_geom_xx,quan_geom_yy
   
      call mpi_job_distribution_chain(nk, ourjob, ourjob_disp, nprocs)

!     berrycurv_unit=(hbarkg**4)/((xm*eV)**2)  ! [m^4] without b1, b1: A^-1
!     berrycurv_unit=(hbarkg**4)/((xm*eV)**2)*(anginv**4)  ! [A^4] without b1, b1: A^-1
      berrycurv_unit=(hbar**4)/((em)**2)*(metertoang**4)*c**4  ! [A^4] without b1, b1: A^-1
!     data hbar/1./ !Here, I will set hbar = 1. for the simplicity
!     write(6,*) "start kubo routine..";stop
      berrycurv_kubo = 0.0d0 
      berrycurv_kubo_ = 0.0d0
      quan_geom_xy = 0.0d0
      quan_geom_xx = 0.0d0
      quan_geom_yy = 0.0d0

!     do ik=1,nk
      do ik=sum(ourjob(1:myrank))+1, sum(ourjob(1:myrank+1))
          call ener_read(ener,isp,ik,nk,nband)      
          !write(6,*) "start ener routine..";stop
          !write(6,*) 'checking nn=1 energy',ener(1);stop
          !do mm=1,nband      ! band index running from 1 to N except for m=n (the # of total bands = N)
          do mm=11,11      ! band index running from 1 to N except for m=n (the # of total bands = N)
!          do mm=nmax+1,nband      ! band index running from 1 to N except for m=n (the # of total bands = N)
!               if(mm.ne.nn) then 
                if(mm.ne.nn .and. abs(ener(nn)-ener(mm)).gt.epsilon(1d0)) then 
!               if(mm.ne.nn .and. &
!                     abs(ener(nn)-ener(mm)).gt. (0.00001)) then 

               cinter_mtrx_x=(0.,0.)
               cinter_mtrx_y=(0.,0.)
!              cinter_mtrx_z=(0.,0.)

               cinter_mtrx_x2=(0.,0.)
               cinter_mtrx_y2=(0.,0.)
!              cinter_mtrx_z2=(0.,0.)

               coeff=(0.,0.)
               coeffn=(0.,0.);coeffm=(0.,0.) ! between band indices n and m
               coeffnu=(0.,0.);coeffmu=(0.,0.)
               coeffnd=(0.,0.);coeffmd=(0.,0.)
               wk(:)=wklist(:,ik)
               np=nplist(ik)
               call plindx(ig,ncnt, ispinor,wk,b1,b2,b3,nbmax,np,& 
                             ecut,npmax)
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
!                 xkgz=(wk(1)+ig(1,iplane))*b1(3)+                     &
!                      (wk(2)+ig(2,iplane))*b2(3)+                     &
!                      (wk(3)+ig(3,iplane))*b3(3)


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
!                    cinter_mtrx_z=cinter_mtrx_z+                     &
!                        conjg(coeffnu(iplane))*xkgz*coeffmu(iplane)+ & 
!                        conjg(coeffnd(iplane))*xkgz*coeffmd(iplane)   

                   cinter_mtrx_x2=cinter_mtrx_x2+                     &
                          conjg(coeffmu(iplane))*xkgx*coeffnu(iplane)+& 
                          conjg(coeffmd(iplane))*xkgx*coeffnd(iplane)  
                   cinter_mtrx_y2=cinter_mtrx_y2+                      &
                           conjg(coeffmu(iplane))*xkgy*coeffnu(iplane)+&
                           conjg(coeffmd(iplane))*xkgy*coeffnd(iplane)
!                  cinter_mtrx_z2=cinter_mtrx_z2+                      &
!                          conjg(coeffmu(iplane))*xkgz*coeffnu(iplane)+&
!                          conjg(coeffmd(iplane))*xkgz*coeffnd(iplane)
! nsoc  case:
                  else if(ispinor .eq. 1) then
                      !!!write(6,*)coeffn(iplane), xkgx
                     cinter_mtrx_x=cinter_mtrx_x+                   &
                           conjg(coeffn(iplane))*xkgx*coeffm(iplane)  ! P_x^{nm}=<u_{n,k}|p_x|u_{m,k}>
                     cinter_mtrx_y=cinter_mtrx_y+                     &
                           conjg(coeffn(iplane))*xkgy*coeffm(iplane) ! P_y^{nm}=<u_{n,k}|p_y|u_{m,k}>


                     cinter_mtrx_x2=cinter_mtrx_x2+                 &
                            conjg(coeffm(iplane))*xkgx*coeffn(iplane)  ! P_x^{mn}=<u_{n,k}|p_x|u_{m,k}>
                     cinter_mtrx_y2=cinter_mtrx_y2+                  &
                            conjg(coeffm(iplane))*xkgy*coeffn(iplane) ! P_y^{mn}=<u_{n,k}|p_y|u_{m,k}>
                  endif
               enddo ! iplane loop end
       !stop


               quan_geom_xy = ((cinter_mtrx_x*cinter_mtrx_y2)/ & 
                            (ener(nn)-ener(mm))**2 )*berrycurv_unit
               quan_geom_xx = ((cinter_mtrx_x*cinter_mtrx_x2)/ & 
                            (ener(nn)-ener(mm))**2 )*berrycurv_unit
               quan_geom_yy = ((cinter_mtrx_y*cinter_mtrx_y2)/ & 
                            (ener(nn)-ener(mm))**2 )*berrycurv_unit

               if(iqm .eq. 1) then
                 berrycurv_kubo(ik,1)=berrycurv_kubo(ik,1)+ & 
                                             real(quan_geom_xx)
                 berrycurv_kubo(ik,2)=berrycurv_kubo(ik,2)+ &
                                             real(quan_geom_xy)
                 berrycurv_kubo(ik,3)=berrycurv_kubo(ik,3)+ &
                                             real(quan_geom_yy)
               else 
                berrycurv_kubo(ik,1)=berrycurv_kubo(ik,1) & 
                                                -2*aimag(quan_geom_xx)
                berrycurv_kubo(ik,2)=berrycurv_kubo(ik,2) & 
                                                -2*aimag(quan_geom_xy)
                berrycurv_kubo(ik,3)=berrycurv_kubo(ik,3) & 
                                                -2*aimag(quan_geom_yy)
               endif           



               endif ! if mm is not equal to nn
           enddo ! band index mm loop end 

!      berrycurv_kubo2 = berrycurv_kubo(ik,3)
      write(6,'(2A,I4,6F16.6)')"# IK, K(reci),", &
             "Berry Curvature (A^2, Kubo) : ", &
                             ik,wk,berrycurv_kubo(ik,:)
      enddo !ik loop end

#ifdef MPI_USE
      call MPI_ALLREDUCE(berrycurv_kubo,berrycurv_kubo_,nk, MPI_REAL8,&
                         MPI_SUM, mpi_comm_earth, mpierr)
      berrycurv_kubo  = berrycurv_kubo_

#endif

      return
      end subroutine kubo_berry_curvature



