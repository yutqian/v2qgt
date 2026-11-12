

!!$*  subroutine for computing berry curvature for the nn band in the given k-point using kubo formula
      subroutine kubo_berry_curvature_non_abelian1(berrycurv_kubo, &
                b1,b2,b3,wklist,isp,                  &
                nband,ecut,ispinor,nplist,nbmax,npmax,&
                nk,nn,mm,nmax,iqm,                       &
                nprocs,myrank,mpi_comm_earth)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
       integer kk, mm
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
      integer :: nskip_gap
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
      !mm = 8
      do ik=sum(ourjob(1:myrank))+1, sum(ourjob(1:myrank+1))
          call ener_read(ener,isp,ik,nk,nband)      
          !write(6,*) "start ener routine..";stop
          !write(6,*) 'checking nn=1 energy',ener(1);stop
          do kk=1,nband      ! band index running from 1 to N except for m=n (the # of total bands = N)
!          do kk=nmax+1,nband      ! band index running from 1 to N except for m=n (the # of total bands = N)
!               if(kk.ne.nn) then 
                if(kk.ne.nn .and. kk .ne. mm .and.                &
                       abs(ener(nn)-ener(kk)).gt.epsilon(1d0)     &
                 .and. abs(ener(mm)-ener(kk)).gt.epsilon(1d0)) then 
!               if(kk.ne.nn .and. &
!                     abs(ener(nn)-ener(kk)).gt. (0.00001)) then 

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
                          nk*(nband+1)*(isp-1)+kk))(coeff(i),i=1,np)  
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
                            (ener(nn)-ener(kk))**2 )*berrycurv_unit
               quan_geom_xx = ((cinter_mtrx_x*cinter_mtrx_x2)/ & 
                            (ener(nn)-ener(kk))**2 )*berrycurv_unit
               quan_geom_yy = ((cinter_mtrx_y*cinter_mtrx_y2)/ & 
                            (ener(nn)-ener(kk))**2 )*berrycurv_unit

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



               endif ! if kk is not equal to nn
           enddo ! band index kk loop end 

!      berrycurv_kubo2 = berrycurv_kubo(ik,3)
      write(6,'(2A,I4,6F16.6)')"# IK, K(reci),", &
             "Berry Curvature (A^2, Kubo) : ", &
                             ik,wk,berrycurv_kubo(ik,:)
      enddo !ik loop end

#ifdef MPI_USE
      call MPI_ALLREDUCE(berrycurv_kubo,berrycurv_kubo_,3*nk, MPI_REAL8,&
                         MPI_SUM, mpi_comm_earth, mpierr)
      berrycurv_kubo  = berrycurv_kubo_

#endif

      return
      end subroutine kubo_berry_curvature_non_abelian1





!!$*  subroutine for computing berry curvature for the nn band in the given k-point using kubo formula
      subroutine kubo_berry_curvature_non_abelian2(berrycurv_kubo, &
                Q_nm, & 
                b1,b2,b3,wklist,isp,                  &
                nband,ecut,ispinor,nplist,nbmax,npmax,&
                nk,nn,mm,abelian,nmax,iqm,                       &
                nprocs,myrank,mpi_comm_earth)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
       integer kk, mm, abelian
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
      complex*16 cinter_mtrx_x_nk,cinter_mtrx_y_nk
      complex*16 cinter_mtrx_x_km,cinter_mtrx_y_km
      complex*16 cinter_mtrx_x_kn,cinter_mtrx_y_kn

      complex*16 coeffn(npmax),coeffm(npmax),coeffk(npmax)
      complex*8  coeff(npmax)
      complex*16 coeffnu(npmax),coeffnd(npmax)
      complex*16 coeffmu(npmax),coeffmd(npmax)
      complex*16 coeffku(npmax),coeffkd(npmax)
      integer :: ig(3,npmax)
      integer :: nskip_gap
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
      complex*16 quan_geom_xy,quan_geom_xx,quan_geom_yy,quan_geom_yx
      complex*16,parameter :: zi=(0.0d0,1.0d0)
   
      call mpi_job_distribution_chain(nk, ourjob, ourjob_disp, nprocs)

!     berrycurv_unit=(hbarkg**4)/((xm*eV)**2)  ! [m^4] without b1, b1: A^-1
!     berrycurv_unit=(hbarkg**4)/((xm*eV)**2)*(anginv**4)  ! [A^4] without b1, b1: A^-1
      berrycurv_unit=(hbar**4)/((em)**2)*(metertoang**4)*c**4  ! [A^4] without b1, b1: A^-1
!     data hbar/1./ !Here, I will set hbar = 1. for the simplicity
!     write(6,*) "start kubo routine..";stop
      berrycurv_kubo = 0.0d0 
      berrycurv_kubo_ = 0.0d0
      quan_geom_xy = 0.0d0
      quan_geom_yx = 0.0d0
      quan_geom_xx = 0.0d0
      quan_geom_yy = 0.0d0

!     do ik=1,nk
      !mm = 8
      write(6,'(2A,3I4,6F16.6)')"# ,", &
                                 "abelian, mm :", &
                                abelian,mm ,Q_nm

      do ik=sum(ourjob(1:myrank))+1, sum(ourjob(1:myrank+1))
          call ener_read(ener,isp,ik,nk,nband)      
          !write(6,*) "start ener routine..";stop
          !write(6,*) 'checking nn=1 energy',ener(1);stop
          do kk=1,nband      ! band index running from 1 to N except for m=n (the # of total bands = N)
!          do kk=nmax+1,nband      ! band index running from 1 to N except for m=n (the # of total bands = N)
!               if(kk.ne.nn) then 
! abelian case
! non-abelian
         
                if((abelian==2 .and. kk.ne.nn .and. kk .ne. mm .and.                &
                       abs(ener(nn)-ener(kk)).gt.epsilon(1d0)     &
                 .and. abs(ener(mm)-ener(kk)).gt.epsilon(1d0)) .or. &
                  (abelian==1 .and. kk.ne.nn .and. &
                     abs(ener(nn)-ener(kk)).gt. (0.00001)) ) then 

               cinter_mtrx_x=(0.,0.)
               cinter_mtrx_y=(0.,0.)
!              cinter_mtrx_z=(0.,0.)

               cinter_mtrx_x2=(0.,0.)
               cinter_mtrx_y2=(0.,0.)
!              cinter_mtrx_z2=(0.,0.)

               cinter_mtrx_x_nk=(0.,0.)
               cinter_mtrx_y_nk=(0.,0.)
               cinter_mtrx_x_km=(0.,0.)
               cinter_mtrx_y_km=(0.,0.)

               cinter_mtrx_x_kn=(0.,0.)
               cinter_mtrx_y_kn=(0.,0.)

               coeff=(0.,0.)
               coeffn=(0.,0.);coeffm=(0.,0.);coeffk=(0.,0.) ! between band indices n and m
               coeffnu=(0.,0.);coeffmu=(0.,0.);coeffku=(0.,0.)
               coeffnd=(0.,0.);coeffmd=(0.,0.);coeffkd=(0.,0.)
               wk(:)=wklist(:,ik)
               np=nplist(ik)
               call plindx(ig,ncnt, ispinor,wk,b1,b2,b3,nbmax,np,& 
                             ecut,npmax)
               !write(6,*)ik,nband,isp,nn,np;stop
               read(10,rec=(3+(ik-1)*(nband+1)+                      &
                          nk*(nband+1)*(isp-1)+nn))(coeff(i),i=1,np)  
               coeffn=coeff;coeff=(0.,0.)                             

               read(10,rec=(3+(ik-1)*(nband+1)+                      &
                          nk*(nband+1)*(isp-1)+kk))(coeff(i),i=1,np)  
               coeffk=coeff;coeff=(0.,0.)                             

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

                  if(abelian ==1 ) then
!abelian case
! computing P_x^{nk}, P_y^{nk}, P_x^{kn}, P_y^{kn}
! Abelian Kubo: matrix elements between target band n and summation band k
! soc  case:
                      if(ispinor .eq. 2) then
                         coeffnu(iplane)=coeffn(iplane)
                         coeffnd(iplane)=coeffn(iplane+ncnt)
                         coeffku(iplane)=coeffk(iplane)
                         coeffkd(iplane)=coeffk(iplane+ncnt)
    
                         cinter_mtrx_x=cinter_mtrx_x+                     &
                             conjg(coeffnu(iplane))*xkgx*coeffku(iplane)+ & 
                             conjg(coeffnd(iplane))*xkgx*coeffkd(iplane)   
                         cinter_mtrx_y=cinter_mtrx_y+                     &
                             conjg(coeffnu(iplane))*xkgy*coeffku(iplane)+ &
                             conjg(coeffnd(iplane))*xkgy*coeffkd(iplane)
!                        cinter_mtrx_z=cinter_mtrx_z+                     &
!                            conjg(coeffnu(iplane))*xkgz*coeffku(iplane)+ & 
!                            conjg(coeffnd(iplane))*xkgz*coeffkd(iplane)   
    
                       cinter_mtrx_x2=cinter_mtrx_x2+                     &
                              conjg(coeffku(iplane))*xkgx*coeffnu(iplane)+& 
                              conjg(coeffkd(iplane))*xkgx*coeffnd(iplane)  
                       cinter_mtrx_y2=cinter_mtrx_y2+                      &
                               conjg(coeffku(iplane))*xkgy*coeffnu(iplane)+&
                               conjg(coeffkd(iplane))*xkgy*coeffnd(iplane)
!                      cinter_mtrx_z2=cinter_mtrx_z2+                      &
!                              conjg(coeffku(iplane))*xkgz*coeffnu(iplane)+&
!                              conjg(coeffkd(iplane))*xkgz*coeffnd(iplane)
! nsoc  case:
                      else if(ispinor .eq. 1) then
                          !!!write(6,*)coeffn(iplane), xkgx
                         cinter_mtrx_x=cinter_mtrx_x+                   &
                               conjg(coeffn(iplane))*xkgx*coeffk(iplane)  ! P_x^{nm}=<u_{n,k}|p_x|u_{m,k}>
                         cinter_mtrx_y=cinter_mtrx_y+                     &
                               conjg(coeffn(iplane))*xkgy*coeffk(iplane) ! P_y^{nm}=<u_{n,k}|p_y|u_{m,k}>
    
    
                         cinter_mtrx_x2=cinter_mtrx_x2+                 &
                                conjg(coeffk(iplane))*xkgx*coeffn(iplane)  ! P_x^{mn}=<u_{n,k}|p_x|u_{m,k}>
                         cinter_mtrx_y2=cinter_mtrx_y2+                  &
                                conjg(coeffk(iplane))*xkgy*coeffn(iplane) ! P_y^{mn}=<u_{n,k}|p_y|u_{m,k}>
                      endif
                  elseif(abelian == 2)then ! non-abelian
                      if(ispinor .eq. 2) then
                       ! <n|Hx|k> = UP + DW
                       cinter_mtrx_x_nk = cinter_mtrx_x_nk + &
                            conjg(coeffn(iplane))*xkgx*coeffk(iplane) + &
                            conjg(coeffn(iplane+ncnt))*xkgx*coeffk(iplane+ncnt)
                    
                       ! <n|Hy|k> = UP + DW
                       cinter_mtrx_y_nk = cinter_mtrx_y_nk + &
                            conjg(coeffn(iplane))*xkgy*coeffk(iplane) + &
                            conjg(coeffn(iplane+ncnt))*xkgy*coeffk(iplane+ncnt)
                    
                       ! <k|Hx|m> = UP + DW
                       cinter_mtrx_x_km = cinter_mtrx_x_km + &
                            conjg(coeffk(iplane))*xkgx*coeffm(iplane) + &
                            conjg(coeffk(iplane+ncnt))*xkgx*coeffm(iplane+ncnt)
                    
                       ! <k|Hy|m> = UP + DW
                       cinter_mtrx_y_km = cinter_mtrx_y_km + &
                            conjg(coeffk(iplane))*xkgy*coeffm(iplane) + &
                            conjg(coeffk(iplane+ncnt))*xkgy*coeffm(iplane+ncnt)
                    
                       ! <k|Hx|n> = UP + DW
                       cinter_mtrx_x_kn = cinter_mtrx_x_kn + &
                            conjg(coeffk(iplane))*xkgx*coeffn(iplane) + &
                            conjg(coeffk(iplane+ncnt))*xkgx*coeffn(iplane+ncnt)
                    
                       ! <k|Hy|n> = UP + DW
                       cinter_mtrx_y_kn = cinter_mtrx_y_kn + &
                            conjg(coeffk(iplane))*xkgy*coeffn(iplane) + &
                            conjg(coeffk(iplane+ncnt))*xkgy*coeffn(iplane+ncnt)
!   ! nsoc  case:
                      else if(ispinor .eq. 1) then

!   
                         cinter_mtrx_x_nk=cinter_mtrx_x_nk+                   &
                               conjg(coeffn(iplane))*xkgx*coeffk(iplane)  ! P_x^{nk}=<u_{n}|p_x|u_{k}>
                         cinter_mtrx_y_nk=cinter_mtrx_y_nk+                     &
                               conjg(coeffn(iplane))*xkgy*coeffk(iplane) ! P_y^{nk}=<u_{n}|p_y|u_{k}>
    
                         cinter_mtrx_x_km=cinter_mtrx_x_km+                 &
                                conjg(coeffk(iplane))*xkgx*coeffm(iplane)! P_x^{km}=<u_{k}|p_x|u_{m}>
                         cinter_mtrx_y_km=cinter_mtrx_y_km+                  &
                                conjg(coeffk(iplane))*xkgy*coeffm(iplane)! P_y^{km}=<u_{k}|p_y|u_{m}>

                        cinter_mtrx_x_kn=cinter_mtrx_x_kn+                 &
                               conjg(coeffk(iplane))*xkgx*coeffn(iplane)  ! P_x^{mn}=<u_{n,k}|p_x|u_{m,k}>
    
                        cinter_mtrx_y_kn=cinter_mtrx_y_kn+                  &
                               conjg(coeffk(iplane))*xkgy*coeffn(iplane) ! P_y^{mn}=<u_{n,k}|p_y|u_{m,k}>

!                        cinter_mtrx_x_mk=cinter_mtrx_x_mk+                   &
!                              conjg(coeffm(iplane))*xkgx*coeffk(iplane)  ! P_x^{nm}=<u_{n,k}|p_x|u_{m,k}>
!                        cinter_mtrx_y_mk=cinter_mtrx_y_mk+                     &
!                              conjg(coeffm(iplane))*xkgy*coeffk(iplane) ! P_y^{nm}=<u_{n,k}|p_y|u_{m,k}>

                      endif ! ispinor

                  endif ! Abelian or non-Abelian
               enddo ! iplane loop end
       !stop

!
               if(abelian==1) then !abelian
                    quan_geom_xy = ((cinter_mtrx_x*cinter_mtrx_y2)/ & 
                                 (ener(nn)-ener(kk))**2 )*berrycurv_unit
                    quan_geom_xx = ((cinter_mtrx_x*cinter_mtrx_x2)/ & 
                                 (ener(nn)-ener(kk))**2 )*berrycurv_unit
                    quan_geom_yy = ((cinter_mtrx_y*cinter_mtrx_y2)/ & 
                                 (ener(nn)-ener(kk))**2 )*berrycurv_unit
     
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
               elseif(abelian==2) then ! non-abelian
!                   aa=  cinter_mtrx_x*cinter_mtrx_x2
!                   bb=  cinter_mtrx_x_nk*cinter_mtrx_x_kn

!                   write(6,'(2A,6F16.6)')"# IK, K(reci),", &
!                       "aa, bb : ", &
!                           aa, bb
                   !> calculate Q_nn or Q_nm?
                    if (Q_nm == 1) then  !calculate Q_nm
                       quan_geom_xy = ((cinter_mtrx_x_nk*cinter_mtrx_y_km)/ & 
                                    ((ener(nn)-ener(kk))*(ener(mm)-ener(kk))))*berrycurv_unit
                       quan_geom_yx = ((cinter_mtrx_y_nk*cinter_mtrx_x_km)/ & 
                                    ((ener(nn)-ener(kk))*(ener(mm)-ener(kk))))*berrycurv_unit
   
                       quan_geom_xx = ((cinter_mtrx_x_nk*cinter_mtrx_x_km)/ & 
                                    ((ener(nn)-ener(kk))*(ener(mm)-ener(kk))))*berrycurv_unit
                       quan_geom_yy = ((cinter_mtrx_y_nk*cinter_mtrx_y_km)/ & 
                                    ((ener(nn)-ener(kk))*(ener(mm)-ener(kk))))*berrycurv_unit
                    elseif(Q_nm ==0)then !calculate Q_nn
                       quan_geom_xy = ((cinter_mtrx_x_nk*cinter_mtrx_y_kn)/ & 
                                    ((ener(nn)-ener(kk))*(ener(nn)-ener(kk))))*berrycurv_unit
                       quan_geom_yx =((cinter_mtrx_y_nk*cinter_mtrx_x_kn)/ & 
                                    ((ener(nn)-ener(kk))*(ener(nn)-ener(kk))))*berrycurv_unit
   
                       quan_geom_xx = ((cinter_mtrx_x_nk*cinter_mtrx_x_kn)/ & 
                                    ((ener(nn)-ener(kk))*(ener(nn)-ener(kk))))*berrycurv_unit
                       quan_geom_yy = ((cinter_mtrx_y_nk*cinter_mtrx_y_kn)/ & 
                                    ((ener(nn)-ener(kk))*(ener(nn)-ener(kk))))*berrycurv_unit

                    endif
     

                   !> calculate quantum metric or berry curvature
                    if(iqm .eq. 1) then
                      berrycurv_kubo(ik,1)=berrycurv_kubo(ik,1)+ & 
                                                  real(quan_geom_xx)
                      berrycurv_kubo(ik,2)=berrycurv_kubo(ik,2)+ &
                                                  real((quan_geom_xy+quan_geom_yx)/2)
                      berrycurv_kubo(ik,3)=berrycurv_kubo(ik,3)+ &
                                                  real(quan_geom_yy)
                    else 
                     berrycurv_kubo(ik,1)=0 !berrycurv_kubo(ik,1) & 
                                            !         abs(aimag(quan_geom_xx-quan_geom_xx))
                     ! using abs for magnitude
                     berrycurv_kubo(ik,2)=berrycurv_kubo(ik,2) + & 
                                                     real(zi*(quan_geom_xy-quan_geom_yx))
                     berrycurv_kubo(ik,3)=0 !berrycurv_kubo(ik,3) &
                                          !           abs(aimag(quan_geom_yy-quan_geom_yy))
                    endif ! iqm
               endif ! abelian or non-abelian

               endif ! if kk is not equal to nn
           enddo ! band index kk loop end 

      write(6,'(2A,I4,6F16.6)')"# IK, K(reci),", &
            "Berry Curvature (A^2, Kubo) : ", &
                            ik,wk,berrycurv_kubo(ik,:)
      enddo !ik loop end

#ifdef MPI_USE
      call MPI_ALLREDUCE(berrycurv_kubo,berrycurv_kubo_,3*nk, MPI_REAL8,&
                         MPI_SUM, mpi_comm_earth, mpierr)
      berrycurv_kubo  = berrycurv_kubo_

#endif

      return
      end subroutine kubo_berry_curvature_non_abelian2






!!$*  subroutine for checking the 2D Kramers-pair QGT inequality
!!$*  Eq. (5): sqrt(det g_2D) >= 1/2 ||F_12||.
!!$*
!!$*  Output ineq2d(nk,8):
!!$*    1 g11
!!$*    2 g22
!!$*    3 g12
!!$*    4 detg
!!$*    5 lhs = sqrt(detg)
!!$*    6 omega_norm2 = 1/2 Tr[Omega_SU2^2]
!!$*    7 rhs = 1/2 sqrt(omega_norm2)
!!$*    8 margin = lhs - rhs
      subroutine kubo_qgt_ineq_2d_kramers_pair(ineq2d, &
                b1,b2,b3,wklist,isp,                  &
                nband,ecut,ispinor,nplist,nbmax,npmax,&
                nk,nn,mm,nmax,                        &
                nprocs,myrank,mpi_comm_earth)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      integer kk, mm, nn
      integer nprocs, myrank, mpi_comm_earth, mpierr
      dimension nbmax(3),nplist(nk),wk(3),ener(nband)
      real*8 ineq2d(nk,18),ineq2d_(nk,18)
      real*8 b1(3),b2(3),b3(3),wklist(3,nk)
      real*8 metertoang
      real*8 berrycurv_unit
      real*8 g11,g22,g12,detg,lhs,rhs,margin
      real*8 omega_norm2,e0,denom,xkgx,xkgy
      real*8 gap,gap_tol,mingap
      integer*4 ourjob(nprocs),ourjob_disp(0:nprocs-1)
      integer :: ig(3,npmax)
      integer :: nskip_gap
      integer :: ext_mode, ext_bmin, ext_bmax, use_kk
      integer :: envlen, envstat
      character*32 ext_env

      complex*8 coeff(npmax)
      complex*16 coeffn(npmax),coeffm(npmax),coeffk(npmax)
      complex*16 vx_nk,vy_nk,vx_mk,vy_mk
      complex*16 vx_kn,vy_kn,vx_km,vy_km
      complex*16 Qxx(2,2),Qyy(2,2),Qxy(2,2),Qyx(2,2)
      complex*16 Omega(2,2),Omega_su2(2,2),trOmega
      complex*16,parameter :: zi=(0.0d0,1.0d0)

      data hbar/6.582119569e-16/
      data c/2.99792458E+08/
      data em/0.510998950E+6/
      data metertoang/1.E+10/

      call mpi_job_distribution_chain(nk,ourjob,ourjob_disp,nprocs)

      berrycurv_unit=(hbar**4)/((em)**2)*(metertoang**4)*c**4
      ineq2d=0.0d0
      ineq2d_=0.0d0
      gap_tol=1.0d-6

!     External-band selection for Cartesian kx,ky mode.
!     QGT_EXT_MODE=all    : all kk except nn,mm
!     QGT_EXT_MODE=near   : kk inside [QGT_EXT_BMIN,QGT_EXT_BMAX]
!     QGT_EXT_MODE=remote : kk outside [QGT_EXT_BMIN,QGT_EXT_BMAX]
      ext_mode=0
      ext_bmin=1
      ext_bmax=nband
      ext_env=' '
      call get_environment_variable('QGT_EXT_MODE',ext_env, &
           length=envlen,status=envstat)
      if(envstat.eq.0 .and. envlen.gt.0) then
         if(ext_env(1:4).eq.'near' .or. ext_env(1:4).eq.'NEAR') &
              ext_mode=1
         if(ext_env(1:6).eq.'remote' .or. ext_env(1:6).eq.'REMOTE') &
              ext_mode=2
      endif
      ext_env=' '
      call get_environment_variable('QGT_EXT_BMIN',ext_env, &
           length=envlen,status=envstat)
      if(envstat.eq.0 .and. envlen.gt.0) read(ext_env,*) ext_bmin
      ext_env=' '
      call get_environment_variable('QGT_EXT_BMAX',ext_env, &
           length=envlen,status=envstat)
      if(envstat.eq.0 .and. envlen.gt.0) read(ext_env,*) ext_bmax
      if(ext_bmin.lt.1) ext_bmin=1
      if(ext_bmax.gt.nband) ext_bmax=nband
      if(ext_bmax.lt.ext_bmin) then
         if(myrank.eq.0) then
            write(6,'(A,2I8)') '# ERROR: bad QGT_EXT_BMIN/BMAX = ', &
                 ext_bmin, ext_bmax
         endif
         stop
      endif
      if(myrank.eq.0) then
         write(6,'(A,ES12.4,A)') '# QGT_CART gap_tol = ', gap_tol, &
              ' eV'
         if(ext_mode.eq.0) then
            write(6,'(A)') '# QGT_EXT_MODE = all, external kk: all bands except nn,mm'
         else if(ext_mode.eq.1) then
            write(6,'(A,2I8)') '# QGT_EXT_MODE = near, active window = ', &
                 ext_bmin, ext_bmax
         else if(ext_mode.eq.2) then
            write(6,'(A,2I8)') '# QGT_EXT_MODE = remote, excluded active window = ', &
                 ext_bmin, ext_bmax
         endif
      endif

      do ik=sum(ourjob(1:myrank))+1,sum(ourjob(1:myrank+1))
         call ener_read(ener,isp,ik,nk,nband)

         Qxx=(0.0d0,0.0d0)
         Qyy=(0.0d0,0.0d0)
         Qxy=(0.0d0,0.0d0)
         Qyx=(0.0d0,0.0d0)

         e0=0.5d0*(ener(nn)+ener(mm))
         mingap=huge(1.0d0)
         nskip_gap=0

         wk(:)=wklist(:,ik)
         np=nplist(ik)
         call plindx(ig,ncnt,ispinor,wk,b1,b2,b3,nbmax,np, &
                      ecut,npmax)

         read(10,rec=(3+(ik-1)*(nband+1)+ &
              nk*(nband+1)*(isp-1)+nn))(coeff(i),i=1,np)
         coeffn=coeff
         coeff=(0.0,0.0)

         read(10,rec=(3+(ik-1)*(nband+1)+ &
              nk*(nband+1)*(isp-1)+mm))(coeff(i),i=1,np)
         coeffm=coeff
         coeff=(0.0,0.0)

         do kk=1,nband
            use_kk=0
            if(kk.ne.nn .and. kk.ne.mm) then
               if(ext_mode.eq.0) then
                  use_kk=1
               else if(ext_mode.eq.1) then
                  if(kk.ge.ext_bmin .and. kk.le.ext_bmax) use_kk=1
               else if(ext_mode.eq.2) then
                  if(kk.lt.ext_bmin .or. kk.gt.ext_bmax) use_kk=1
               endif
            endif
            if(use_kk.eq.1) then
               gap=abs(e0-ener(kk))
               if(gap.lt.mingap) mingap=gap
               if(gap.le.gap_tol) then
                  nskip_gap=nskip_gap+1
               else

               read(10,rec=(3+(ik-1)*(nband+1)+ &
                    nk*(nband+1)*(isp-1)+kk))(coeff(i),i=1,np)
               coeffk=coeff
               coeff=(0.0,0.0)

               vx_nk=(0.0d0,0.0d0)
               vy_nk=(0.0d0,0.0d0)
               vx_mk=(0.0d0,0.0d0)
               vy_mk=(0.0d0,0.0d0)
               vx_kn=(0.0d0,0.0d0)
               vy_kn=(0.0d0,0.0d0)
               vx_km=(0.0d0,0.0d0)
               vy_km=(0.0d0,0.0d0)

               do iplane=1,ncnt
                  xkgx=(wk(1)+ig(1,iplane))*b1(1)+ &
                       (wk(2)+ig(2,iplane))*b2(1)+ &
                       (wk(3)+ig(3,iplane))*b3(1)
                  xkgy=(wk(1)+ig(1,iplane))*b1(2)+ &
                       (wk(2)+ig(2,iplane))*b2(2)+ &
                       (wk(3)+ig(3,iplane))*b3(2)

                  if(ispinor.eq.2) then
                     vx_nk=vx_nk+ &
                        conjg(coeffn(iplane))*xkgx*coeffk(iplane)+ &
                        conjg(coeffn(iplane+ncnt))*xkgx* &
                        coeffk(iplane+ncnt)
                     vy_nk=vy_nk+ &
                        conjg(coeffn(iplane))*xkgy*coeffk(iplane)+ &
                        conjg(coeffn(iplane+ncnt))*xkgy* &
                        coeffk(iplane+ncnt)

                     vx_mk=vx_mk+ &
                        conjg(coeffm(iplane))*xkgx*coeffk(iplane)+ &
                        conjg(coeffm(iplane+ncnt))*xkgx* &
                        coeffk(iplane+ncnt)
                     vy_mk=vy_mk+ &
                        conjg(coeffm(iplane))*xkgy*coeffk(iplane)+ &
                        conjg(coeffm(iplane+ncnt))*xkgy* &
                        coeffk(iplane+ncnt)

                     vx_kn=vx_kn+ &
                        conjg(coeffk(iplane))*xkgx*coeffn(iplane)+ &
                        conjg(coeffk(iplane+ncnt))*xkgx* &
                        coeffn(iplane+ncnt)
                     vy_kn=vy_kn+ &
                        conjg(coeffk(iplane))*xkgy*coeffn(iplane)+ &
                        conjg(coeffk(iplane+ncnt))*xkgy* &
                        coeffn(iplane+ncnt)

                     vx_km=vx_km+ &
                        conjg(coeffk(iplane))*xkgx*coeffm(iplane)+ &
                        conjg(coeffk(iplane+ncnt))*xkgx* &
                        coeffm(iplane+ncnt)
                     vy_km=vy_km+ &
                        conjg(coeffk(iplane))*xkgy*coeffm(iplane)+ &
                        conjg(coeffk(iplane+ncnt))*xkgy* &
                        coeffm(iplane+ncnt)
                  else if(ispinor.eq.1) then
                     vx_nk=vx_nk+conjg(coeffn(iplane))*xkgx* &
                          coeffk(iplane)
                     vy_nk=vy_nk+conjg(coeffn(iplane))*xkgy* &
                          coeffk(iplane)

                     vx_mk=vx_mk+conjg(coeffm(iplane))*xkgx* &
                          coeffk(iplane)
                     vy_mk=vy_mk+conjg(coeffm(iplane))*xkgy* &
                          coeffk(iplane)

                     vx_kn=vx_kn+conjg(coeffk(iplane))*xkgx* &
                          coeffn(iplane)
                     vy_kn=vy_kn+conjg(coeffk(iplane))*xkgy* &
                          coeffn(iplane)

                     vx_km=vx_km+conjg(coeffk(iplane))*xkgx* &
                          coeffm(iplane)
                     vy_km=vy_km+conjg(coeffk(iplane))*xkgy* &
                          coeffm(iplane)
                  endif
               enddo

               denom=(e0-ener(kk))**2

               Qxx(1,1)=Qxx(1,1)+vx_nk*vx_kn/denom*berrycurv_unit
               Qxx(1,2)=Qxx(1,2)+vx_nk*vx_km/denom*berrycurv_unit
               Qxx(2,1)=Qxx(2,1)+vx_mk*vx_kn/denom*berrycurv_unit
               Qxx(2,2)=Qxx(2,2)+vx_mk*vx_km/denom*berrycurv_unit

               Qyy(1,1)=Qyy(1,1)+vy_nk*vy_kn/denom*berrycurv_unit
               Qyy(1,2)=Qyy(1,2)+vy_nk*vy_km/denom*berrycurv_unit
               Qyy(2,1)=Qyy(2,1)+vy_mk*vy_kn/denom*berrycurv_unit
               Qyy(2,2)=Qyy(2,2)+vy_mk*vy_km/denom*berrycurv_unit

               Qxy(1,1)=Qxy(1,1)+vx_nk*vy_kn/denom*berrycurv_unit
               Qxy(1,2)=Qxy(1,2)+vx_nk*vy_km/denom*berrycurv_unit
               Qxy(2,1)=Qxy(2,1)+vx_mk*vy_kn/denom*berrycurv_unit
               Qxy(2,2)=Qxy(2,2)+vx_mk*vy_km/denom*berrycurv_unit

               Qyx(1,1)=Qyx(1,1)+vy_nk*vx_kn/denom*berrycurv_unit
               Qyx(1,2)=Qyx(1,2)+vy_nk*vx_km/denom*berrycurv_unit
               Qyx(2,1)=Qyx(2,1)+vy_mk*vx_kn/denom*berrycurv_unit
               Qyx(2,2)=Qyx(2,2)+vy_mk*vx_km/denom*berrycurv_unit
               endif
            endif
         enddo

         if(mingap.eq.huge(1.0d0)) mingap=0.0d0

         g11=0.5d0*real(Qxx(1,1)+Qxx(2,2))
         g22=0.5d0*real(Qyy(1,1)+Qyy(2,2))
         g12=0.25d0*real(Qxy(1,1)+Qyx(1,1)+ &
                         Qxy(2,2)+Qyx(2,2))

         Omega=zi*(Qxy-Qyx)
         trOmega=Omega(1,1)+Omega(2,2)
         Omega_su2=Omega
         Omega_su2(1,1)=Omega(1,1)-0.5d0*trOmega
         Omega_su2(2,2)=Omega(2,2)-0.5d0*trOmega

         omega_norm2=0.5d0*real( &
              Omega_su2(1,1)*Omega_su2(1,1) + &
              Omega_su2(1,2)*Omega_su2(2,1) + &
              Omega_su2(2,1)*Omega_su2(1,2) + &
              Omega_su2(2,2)*Omega_su2(2,2))

         detg=g11*g22-g12*g12
         if(detg.lt.0.0d0 .and. abs(detg).lt.1.0d-12) detg=0.0d0
         if(omega_norm2.lt.0.0d0 .and. abs(omega_norm2).lt.1.0d-12) &
              omega_norm2=0.0d0

         lhs=sqrt(max(detg,0.0d0))
         rhs=0.5d0*sqrt(max(omega_norm2,0.0d0))
         margin=lhs-rhs

         ineq2d(ik,1)=g11
         ineq2d(ik,2)=g22
         ineq2d(ik,3)=g12
         ineq2d(ik,4)=detg
         ineq2d(ik,5)=lhs
         ineq2d(ik,6)=omega_norm2
         ineq2d(ik,7)=rhs
         ineq2d(ik,8)=margin
         ineq2d(ik,9)=real(Omega_su2(1,1))
         ineq2d(ik,10)=aimag(Omega_su2(1,1))
         ineq2d(ik,11)=real(Omega_su2(1,2))
         ineq2d(ik,12)=aimag(Omega_su2(1,2))
         ineq2d(ik,13)=real(Omega_su2(2,1))
         ineq2d(ik,14)=aimag(Omega_su2(2,1))
         ineq2d(ik,15)=real(Omega_su2(2,2))
         ineq2d(ik,16)=aimag(Omega_su2(2,2))
         ineq2d(ik,17)=mingap
         ineq2d(ik,18)=dble(nskip_gap)

         write(6,'(A,I6,3F12.6,18ES16.6E3)') &
              '# QGT_INEQ_2D_CART ik k frac, values: ', &
              ik,wk,ineq2d(ik,:)
      enddo

#ifdef MPI_USE
      call MPI_ALLREDUCE(ineq2d,ineq2d_,18*nk,MPI_REAL8, &
                         MPI_SUM,mpi_comm_earth,mpierr)
      ineq2d=ineq2d_
#endif

      return
      end subroutine kubo_qgt_ineq_2d_kramers_pair


!!$*  subroutine for checking the 2D Kramers-pair QGT inequality
!!$*  on a fractional reciprocal-coordinate plane (kdir1,kdir2).
!!$*
!!$*  This computes F_{k_i k_j}, not F_{cartesian x, cartesian y}.
!!$*  If kdir1=1 and kdir2=2, the two directions are the fractional
!!$*  reciprocal coordinates k1 and k2, with
!!$*       K = k1*b1 + k2*b2 + k3*b3.
!!$*
!!$*  The directional kinetic matrix element uses
!!$*       dH/dk_i ~ (hbar^2/m) (K+G) dot b_i.
!!$*
!!$*  Output ineq2d(nk,18):
!!$*    1  g11
!!$*    2  g22
!!$*    3  g12
!!$*    4  detg
!!$*    5  lhs = sqrt(detg)
!!$*    6  omega_norm2 = 1/2 Tr[Omega_SU2^2]
!!$*    7  rhs = 1/2 sqrt(omega_norm2)
!!$*    8  margin = lhs - rhs
!!$*    9  Re Omega_SU2(1,1)
!!$*   10  Im Omega_SU2(1,1)
!!$*   11  Re Omega_SU2(1,2)
!!$*   12  Im Omega_SU2(1,2)
!!$*   13  Re Omega_SU2(2,1)
!!$*   14  Im Omega_SU2(2,1)
!!$*   15  Re Omega_SU2(2,2)
!!$*   16  Im Omega_SU2(2,2)
      subroutine kubo_qgt_ineq_2d_kdir_kramers_pair(ineq2d, &
                kdir1,kdir2,                              &
                b1,b2,b3,wklist,isp,                      &
                nband,ecut,ispinor,nplist,nbmax,npmax,    &
                nk,nn,mm,nmax,                            &
                nprocs,myrank,mpi_comm_earth)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      integer kk, mm, nn, kdir1, kdir2
      integer nprocs, myrank, mpi_comm_earth, mpierr
      dimension nbmax(3),nplist(nk),wk(3),ener(nband)
      real*8 ineq2d(nk,18),ineq2d_(nk,18)
      real*8 b1(3),b2(3),b3(3),wklist(3,nk)
      real*8 bdir1(3),bdir2(3),kgcart(3)
      real*8 g11,g22,g12,detg,lhs,rhs,margin
      real*8 omega_norm2,e0,denom,gap,gap_tol,mingap,bscale
      real*8 xkg1,xkg2
      real*8 metertoang
      real*8 berrycurv_unit
      integer*4 ourjob(nprocs),ourjob_disp(0:nprocs-1)
      integer :: ig(3,npmax)
      integer :: nskip_gap
      integer :: ext_mode, ext_bmin, ext_bmax
      integer :: use_kk
      character*64 ext_env

      complex*8 coeff(npmax)
      complex*16 coeffn(npmax),coeffm(npmax),coeffk(npmax)
      complex*16 v1_nk,v2_nk,v1_mk,v2_mk
      complex*16 v1_kn,v2_kn,v1_km,v2_km
      complex*16 Q11(2,2),Q22(2,2),Q12(2,2),Q21(2,2)
      complex*16 Omega(2,2),Omega_su2(2,2),trOmega
      complex*16,parameter :: zi=(0.0d0,1.0d0)

      data hbar/6.582119569e-16/
      data c/2.99792458E+08/
      data em/0.510998950E+6/
      data metertoang/1.E+10/

      if(kdir1.lt.1 .or. kdir1.gt.3 .or. &
         kdir2.lt.1 .or. kdir2.gt.3 .or. kdir1.eq.kdir2) then
         if(myrank.eq.0) then
            write(6,'(A,2I8)') '# ERROR: invalid kdir1,kdir2 = ', &
                 kdir1,kdir2
            write(6,'(A)') '# Use different integer directions among 1,2,3.'
         endif
         stop
      endif

      if(kdir1.eq.1) then
         bdir1(:)=b1(:)
      else if(kdir1.eq.2) then
         bdir1(:)=b2(:)
      else
         bdir1(:)=b3(:)
      endif

      if(kdir2.eq.1) then
         bdir2(:)=b1(:)
      else if(kdir2.eq.2) then
         bdir2(:)=b2(:)
      else
         bdir2(:)=b3(:)
      endif

!     Numerical gap cutoff for the Kubo denominator, in eV.
!     If an outside band is closer than this to the Kramers pair,
!     the two-band subspace is not isolated and the QGT of this pair
!     is not reliable at that k point.
      gap_tol=1.0d-6

!     Optional restriction of the Kubo external-band sum.
!     Default: QGT_EXT_MODE=all, sum over all kk except nn,mm.
!     QGT_EXT_MODE=near   sums only kk in [QGT_EXT_BMIN,QGT_EXT_BMAX]
!                         excluding nn,mm.  For active four-band
!                         nn,mm inside [11,14], this gives kk=13,14.
!     QGT_EXT_MODE=remote sums kk outside [QGT_EXT_BMIN,QGT_EXT_BMAX].
!                         For active four-band [11,14], this gives
!                         remote bands 1..10,15..nband.
      ext_mode=0
      ext_bmin=1
      ext_bmax=nband
      ext_env=' '
      call get_environment_variable('QGT_EXT_MODE',ext_env)
      if(ext_env(1:4).eq.'near' .or. ext_env(1:4).eq.'NEAR') ext_mode=1
      if(ext_env(1:6).eq.'remote' .or. ext_env(1:6).eq.'REMOTE') ext_mode=2
      ext_env=' '
      call get_environment_variable('QGT_EXT_BMIN',ext_env)
      if(len_trim(ext_env).gt.0) read(ext_env,*) ext_bmin
      ext_env=' '
      call get_environment_variable('QGT_EXT_BMAX',ext_env)
      if(len_trim(ext_env).gt.0) read(ext_env,*) ext_bmax
      if(ext_bmin.lt.1) ext_bmin=1
      if(ext_bmax.gt.nband) ext_bmax=nband
      if(ext_bmax.lt.ext_bmin) then
         if(myrank.eq.0) then
            write(6,'(A,2I8)') '# ERROR: bad QGT_EXT_BMIN/BMAX = ', ext_bmin, ext_bmax
         endif
         stop
      endif

!     The original code assumes reciprocal vectors in A^{-1} for the
!     unit factor below.  Some WAVECAR readers store b vectors in m^{-1}.
!     For fractional directions d/dk_i, use (K+G).b_i in A^{-2}.
!     If |b| is much larger than a typical A^{-1} vector, convert m^{-1}
!     to A^{-1} by dividing both factors by 1e10.
      bscale=1.0d0
      if(max(sqrt(sum(b1(:)*b1(:))),sqrt(sum(b2(:)*b2(:))), &
             sqrt(sum(b3(:)*b3(:)))).gt.100.0d0) then
         bscale=metertoang
      endif
      if(myrank.eq.0) then
         write(6,'(A,ES12.4)') '# QGT_KDIR reciprocal-vector scale bscale = ', bscale
         write(6,'(A,ES12.4,A)') '# QGT_KDIR gap_tol = ', gap_tol, ' eV'
         if(ext_mode.eq.0) then
            write(6,'(A)') '# QGT_EXT_MODE = all, external kk: all bands except nn,mm'
         else if(ext_mode.eq.1) then
            write(6,'(A,2I8)') '# QGT_EXT_MODE = near, active window = ', ext_bmin, ext_bmax
         else if(ext_mode.eq.2) then
            write(6,'(A,2I8)') '# QGT_EXT_MODE = remote, excluded active window = ', ext_bmin, ext_bmax
         endif
      endif

      call mpi_job_distribution_chain(nk,ourjob,ourjob_disp,nprocs)

      berrycurv_unit=(hbar**4)/((em)**2)*(metertoang**4)*c**4
      ineq2d=0.0d0
      ineq2d_=0.0d0

      do ik=sum(ourjob(1:myrank))+1,sum(ourjob(1:myrank+1))
         call ener_read(ener,isp,ik,nk,nband)

         Q11=(0.0d0,0.0d0)
         Q22=(0.0d0,0.0d0)
         Q12=(0.0d0,0.0d0)
         Q21=(0.0d0,0.0d0)

         e0=0.5d0*(ener(nn)+ener(mm))
         mingap=huge(1.0d0)
         nskip_gap=0

         wk(:)=wklist(:,ik)
         np=nplist(ik)
         call plindx(ig,ncnt,ispinor,wk,b1,b2,b3,nbmax,np, &
                      ecut,npmax)

         read(10,rec=(3+(ik-1)*(nband+1)+ &
              nk*(nband+1)*(isp-1)+nn))(coeff(i),i=1,np)
         coeffn=coeff
         coeff=(0.0,0.0)

         read(10,rec=(3+(ik-1)*(nband+1)+ &
              nk*(nband+1)*(isp-1)+mm))(coeff(i),i=1,np)
         coeffm=coeff
         coeff=(0.0,0.0)

         do kk=1,nband
            use_kk=0
            if(kk.ne.nn .and. kk.ne.mm) then
               if(ext_mode.eq.0) then
                  use_kk=1
               else if(ext_mode.eq.1) then
                  if(kk.ge.ext_bmin .and. kk.le.ext_bmax) use_kk=1
               else if(ext_mode.eq.2) then
                  if(kk.lt.ext_bmin .or. kk.gt.ext_bmax) use_kk=1
               endif
            endif
            if(use_kk.eq.1) then
               gap=abs(e0-ener(kk))
               if(gap.lt.mingap) mingap=gap
               if(gap.le.gap_tol) then
                  nskip_gap=nskip_gap+1
               else

               read(10,rec=(3+(ik-1)*(nband+1)+ &
                    nk*(nband+1)*(isp-1)+kk))(coeff(i),i=1,np)
               coeffk=coeff
               coeff=(0.0,0.0)

               v1_nk=(0.0d0,0.0d0)
               v2_nk=(0.0d0,0.0d0)
               v1_mk=(0.0d0,0.0d0)
               v2_mk=(0.0d0,0.0d0)
               v1_kn=(0.0d0,0.0d0)
               v2_kn=(0.0d0,0.0d0)
               v1_km=(0.0d0,0.0d0)
               v2_km=(0.0d0,0.0d0)

               do iplane=1,ncnt
                  kgcart(1)=(wk(1)+ig(1,iplane))*b1(1)+ &
                            (wk(2)+ig(2,iplane))*b2(1)+ &
                            (wk(3)+ig(3,iplane))*b3(1)
                  kgcart(2)=(wk(1)+ig(1,iplane))*b1(2)+ &
                            (wk(2)+ig(2,iplane))*b2(2)+ &
                            (wk(3)+ig(3,iplane))*b3(2)
                  kgcart(3)=(wk(1)+ig(1,iplane))*b1(3)+ &
                            (wk(2)+ig(2,iplane))*b2(3)+ &
                            (wk(3)+ig(3,iplane))*b3(3)

                  xkg1=(kgcart(1)*bdir1(1)+kgcart(2)*bdir1(2)+ &
                       kgcart(3)*bdir1(3))/(bscale*bscale)
                  xkg2=(kgcart(1)*bdir2(1)+kgcart(2)*bdir2(2)+ &
                       kgcart(3)*bdir2(3))/(bscale*bscale)

                  if(ispinor.eq.2) then
                     v1_nk=v1_nk+ &
                        conjg(coeffn(iplane))*xkg1*coeffk(iplane)+ &
                        conjg(coeffn(iplane+ncnt))*xkg1* &
                        coeffk(iplane+ncnt)
                     v2_nk=v2_nk+ &
                        conjg(coeffn(iplane))*xkg2*coeffk(iplane)+ &
                        conjg(coeffn(iplane+ncnt))*xkg2* &
                        coeffk(iplane+ncnt)

                     v1_mk=v1_mk+ &
                        conjg(coeffm(iplane))*xkg1*coeffk(iplane)+ &
                        conjg(coeffm(iplane+ncnt))*xkg1* &
                        coeffk(iplane+ncnt)
                     v2_mk=v2_mk+ &
                        conjg(coeffm(iplane))*xkg2*coeffk(iplane)+ &
                        conjg(coeffm(iplane+ncnt))*xkg2* &
                        coeffk(iplane+ncnt)

                     v1_kn=v1_kn+ &
                        conjg(coeffk(iplane))*xkg1*coeffn(iplane)+ &
                        conjg(coeffk(iplane+ncnt))*xkg1* &
                        coeffn(iplane+ncnt)
                     v2_kn=v2_kn+ &
                        conjg(coeffk(iplane))*xkg2*coeffn(iplane)+ &
                        conjg(coeffk(iplane+ncnt))*xkg2* &
                        coeffn(iplane+ncnt)

                     v1_km=v1_km+ &
                        conjg(coeffk(iplane))*xkg1*coeffm(iplane)+ &
                        conjg(coeffk(iplane+ncnt))*xkg1* &
                        coeffm(iplane+ncnt)
                     v2_km=v2_km+ &
                        conjg(coeffk(iplane))*xkg2*coeffm(iplane)+ &
                        conjg(coeffk(iplane+ncnt))*xkg2* &
                        coeffm(iplane+ncnt)
                  else if(ispinor.eq.1) then
                     v1_nk=v1_nk+conjg(coeffn(iplane))*xkg1* &
                          coeffk(iplane)
                     v2_nk=v2_nk+conjg(coeffn(iplane))*xkg2* &
                          coeffk(iplane)

                     v1_mk=v1_mk+conjg(coeffm(iplane))*xkg1* &
                          coeffk(iplane)
                     v2_mk=v2_mk+conjg(coeffm(iplane))*xkg2* &
                          coeffk(iplane)

                     v1_kn=v1_kn+conjg(coeffk(iplane))*xkg1* &
                          coeffn(iplane)
                     v2_kn=v2_kn+conjg(coeffk(iplane))*xkg2* &
                          coeffn(iplane)

                     v1_km=v1_km+conjg(coeffk(iplane))*xkg1* &
                          coeffm(iplane)
                     v2_km=v2_km+conjg(coeffk(iplane))*xkg2* &
                          coeffm(iplane)
                  endif
               enddo

               denom=(e0-ener(kk))**2

               Q11(1,1)=Q11(1,1)+v1_nk*v1_kn/denom*berrycurv_unit
               Q11(1,2)=Q11(1,2)+v1_nk*v1_km/denom*berrycurv_unit
               Q11(2,1)=Q11(2,1)+v1_mk*v1_kn/denom*berrycurv_unit
               Q11(2,2)=Q11(2,2)+v1_mk*v1_km/denom*berrycurv_unit

               Q22(1,1)=Q22(1,1)+v2_nk*v2_kn/denom*berrycurv_unit
               Q22(1,2)=Q22(1,2)+v2_nk*v2_km/denom*berrycurv_unit
               Q22(2,1)=Q22(2,1)+v2_mk*v2_kn/denom*berrycurv_unit
               Q22(2,2)=Q22(2,2)+v2_mk*v2_km/denom*berrycurv_unit

               Q12(1,1)=Q12(1,1)+v1_nk*v2_kn/denom*berrycurv_unit
               Q12(1,2)=Q12(1,2)+v1_nk*v2_km/denom*berrycurv_unit
               Q12(2,1)=Q12(2,1)+v1_mk*v2_kn/denom*berrycurv_unit
               Q12(2,2)=Q12(2,2)+v1_mk*v2_km/denom*berrycurv_unit

               Q21(1,1)=Q21(1,1)+v2_nk*v1_kn/denom*berrycurv_unit
               Q21(1,2)=Q21(1,2)+v2_nk*v1_km/denom*berrycurv_unit
               Q21(2,1)=Q21(2,1)+v2_mk*v1_kn/denom*berrycurv_unit
               Q21(2,2)=Q21(2,2)+v2_mk*v1_km/denom*berrycurv_unit
               endif
            endif
         enddo

         if(mingap.eq.huge(1.0d0)) mingap=0.0d0

         g11=0.5d0*real(Q11(1,1)+Q11(2,2))
         g22=0.5d0*real(Q22(1,1)+Q22(2,2))
         g12=0.25d0*real(Q12(1,1)+Q21(1,1)+ &
                         Q12(2,2)+Q21(2,2))

         Omega=zi*(Q12-Q21)
         trOmega=Omega(1,1)+Omega(2,2)
         Omega_su2=Omega
         Omega_su2(1,1)=Omega(1,1)-0.5d0*trOmega
         Omega_su2(2,2)=Omega(2,2)-0.5d0*trOmega

         omega_norm2=0.5d0*real( &
              Omega_su2(1,1)*Omega_su2(1,1) + &
              Omega_su2(1,2)*Omega_su2(2,1) + &
              Omega_su2(2,1)*Omega_su2(1,2) + &
              Omega_su2(2,2)*Omega_su2(2,2))

         detg=g11*g22-g12*g12
         if(detg.lt.0.0d0 .and. abs(detg).lt.1.0d-12) detg=0.0d0
         if(omega_norm2.lt.0.0d0 .and. abs(omega_norm2).lt.1.0d-12) &
              omega_norm2=0.0d0

         lhs=sqrt(max(detg,0.0d0))
         rhs=0.5d0*sqrt(max(omega_norm2,0.0d0))
         margin=lhs-rhs

         ineq2d(ik,1)=g11
         ineq2d(ik,2)=g22
         ineq2d(ik,3)=g12
         ineq2d(ik,4)=detg
         ineq2d(ik,5)=lhs
         ineq2d(ik,6)=omega_norm2
         ineq2d(ik,7)=rhs
         ineq2d(ik,8)=margin
         ineq2d(ik,9)=real(Omega_su2(1,1))
         ineq2d(ik,10)=aimag(Omega_su2(1,1))
         ineq2d(ik,11)=real(Omega_su2(1,2))
         ineq2d(ik,12)=aimag(Omega_su2(1,2))
         ineq2d(ik,13)=real(Omega_su2(2,1))
         ineq2d(ik,14)=aimag(Omega_su2(2,1))
         ineq2d(ik,15)=real(Omega_su2(2,2))
         ineq2d(ik,16)=aimag(Omega_su2(2,2))
         ineq2d(ik,17)=mingap
         ineq2d(ik,18)=dble(nskip_gap)

         write(6,'(A,I6,2I4,3F12.6,18ES16.6E3)') &
              '# QGT_INEQ_2D_KDIR ik kdir k frac, values: ', &
              ik,kdir1,kdir2,wk,ineq2d(ik,:)
      enddo

#ifdef MPI_USE
      call MPI_ALLREDUCE(ineq2d,ineq2d_,18*nk,MPI_REAL8, &
                         MPI_SUM,mpi_comm_earth,mpierr)
      ineq2d=ineq2d_
#endif

      return
      end subroutine kubo_qgt_ineq_2d_kdir_kramers_pair


      subroutine kubo_qgt_ineq_2d_dvec_kramers_pair(ineq2d, &
                dvec1,dvec2,                              &
                b1,b2,b3,wklist,isp,                      &
                nband,ecut,ispinor,nplist,nbmax,npmax,    &
                nk,nn,mm,nmax,                            &
                nprocs,myrank,mpi_comm_earth)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      integer kk, mm, nn
      integer nprocs, myrank, mpi_comm_earth, mpierr
      dimension nbmax(3),nplist(nk),wk(3),ener(nband)
      real*8 ineq2d(nk,18),ineq2d_(nk,18)
      real*8 b1(3),b2(3),b3(3),wklist(3,nk)
      real*8 bdir1(3),bdir2(3),kgcart(3),dvec1(3),dvec2(3)
      real*8 cross(3),norm1,norm2,crossnorm
      real*8 g11,g22,g12,detg,lhs,rhs,margin
      real*8 omega_norm2,e0,denom,gap,gap_tol,mingap,bscale
      real*8 xkg1,xkg2
      real*8 metertoang
      real*8 berrycurv_unit
      integer*4 ourjob(nprocs),ourjob_disp(0:nprocs-1)
      integer :: ig(3,npmax)
      integer :: nskip_gap
      integer :: ext_mode, ext_bmin, ext_bmax
      integer :: use_kk
      character*64 ext_env

      complex*8 coeff(npmax)
      complex*16 coeffn(npmax),coeffm(npmax),coeffk(npmax)
      complex*16 v1_nk,v2_nk,v1_mk,v2_mk
      complex*16 v1_kn,v2_kn,v1_km,v2_km
      complex*16 Q11(2,2),Q22(2,2),Q12(2,2),Q21(2,2)
      complex*16 Omega(2,2),Omega_su2(2,2),trOmega
      complex*16,parameter :: zi=(0.0d0,1.0d0)

      data hbar/6.582119569e-16/
      data c/2.99792458E+08/
      data em/0.510998950E+6/
      data metertoang/1.E+10/

!     Arbitrary derivative directions in fractional reciprocal coordinates.
!     If dvec=(a,b,c), the derivative direction is
!        bdir = a*b1 + b*b2 + c*b3
!     and the code computes dH/dq for K(q)=K0+q*bdir.
      bdir1(:)=dvec1(1)*b1(:)+dvec1(2)*b2(:)+dvec1(3)*b3(:)
      bdir2(:)=dvec2(1)*b1(:)+dvec2(2)*b2(:)+dvec2(3)*b3(:)

      norm1=sqrt(sum(bdir1(:)*bdir1(:)))
      norm2=sqrt(sum(bdir2(:)*bdir2(:)))
      cross(1)=bdir1(2)*bdir2(3)-bdir1(3)*bdir2(2)
      cross(2)=bdir1(3)*bdir2(1)-bdir1(1)*bdir2(3)
      cross(3)=bdir1(1)*bdir2(2)-bdir1(2)*bdir2(1)
      crossnorm=sqrt(sum(cross(:)*cross(:)))
      if(norm1.lt.1.0d-14 .or. norm2.lt.1.0d-14 .or. &
         crossnorm.lt.1.0d-14*norm1*norm2) then
         if(myrank.eq.0) then
            write(6,'(A)') '# ERROR: invalid QGT_DVEC1/QGT_DVEC2.'
            write(6,'(A,3ES18.8E3)') '# dvec1 = ', dvec1
            write(6,'(A,3ES18.8E3)') '# dvec2 = ', dvec2
            write(6,'(A)') '# Directions must be nonzero and nonparallel.'
         endif
         stop
      endif

!     Numerical gap cutoff for the Kubo denominator, in eV.
!     If an outside band is closer than this to the Kramers pair,
!     the two-band subspace is not isolated and the QGT of this pair
!     is not reliable at that k point.
      gap_tol=1.0d-6

!     Optional restriction of the Kubo external-band sum.
!     Default: QGT_EXT_MODE=all, sum over all kk except nn,mm.
!     QGT_EXT_MODE=near   sums only kk in [QGT_EXT_BMIN,QGT_EXT_BMAX]
!                         excluding nn,mm.  For active four-band
!                         nn,mm inside [11,14], this gives kk=13,14.
!     QGT_EXT_MODE=remote sums kk outside [QGT_EXT_BMIN,QGT_EXT_BMAX].
!                         For active four-band [11,14], this gives
!                         remote bands 1..10,15..nband.
      ext_mode=0
      ext_bmin=1
      ext_bmax=nband
      ext_env=' '
      call get_environment_variable('QGT_EXT_MODE',ext_env)
      if(ext_env(1:4).eq.'near' .or. ext_env(1:4).eq.'NEAR') ext_mode=1
      if(ext_env(1:6).eq.'remote' .or. ext_env(1:6).eq.'REMOTE') ext_mode=2
      ext_env=' '
      call get_environment_variable('QGT_EXT_BMIN',ext_env)
      if(len_trim(ext_env).gt.0) read(ext_env,*) ext_bmin
      ext_env=' '
      call get_environment_variable('QGT_EXT_BMAX',ext_env)
      if(len_trim(ext_env).gt.0) read(ext_env,*) ext_bmax
      if(ext_bmin.lt.1) ext_bmin=1
      if(ext_bmax.gt.nband) ext_bmax=nband
      if(ext_bmax.lt.ext_bmin) then
         if(myrank.eq.0) then
            write(6,'(A,2I8)') '# ERROR: bad QGT_EXT_BMIN/BMAX = ', ext_bmin, ext_bmax
         endif
         stop
      endif

!     The original code assumes reciprocal vectors in A^{-1} for the
!     unit factor below.  Some WAVECAR readers store b vectors in m^{-1}.
!     For fractional directions d/dk_i, use (K+G).b_i in A^{-2}.
!     If |b| is much larger than a typical A^{-1} vector, convert m^{-1}
!     to A^{-1} by dividing both factors by 1e10.
      bscale=1.0d0
      if(max(sqrt(sum(b1(:)*b1(:))),sqrt(sum(b2(:)*b2(:))), &
             sqrt(sum(b3(:)*b3(:)))).gt.100.0d0) then
         bscale=metertoang
      endif
      if(myrank.eq.0) then
         write(6,'(A,ES12.4)') '# QGT_DVEC reciprocal-vector scale bscale = ', bscale
         write(6,'(A,ES12.4,A)') '# QGT_DVEC gap_tol = ', gap_tol, ' eV'
         if(ext_mode.eq.0) then
            write(6,'(A)') '# QGT_EXT_MODE = all, external kk: all bands except nn,mm'
         else if(ext_mode.eq.1) then
            write(6,'(A,2I8)') '# QGT_EXT_MODE = near, active window = ', ext_bmin, ext_bmax
         else if(ext_mode.eq.2) then
            write(6,'(A,2I8)') '# QGT_EXT_MODE = remote, excluded active window = ', ext_bmin, ext_bmax
         endif
      endif

      call mpi_job_distribution_chain(nk,ourjob,ourjob_disp,nprocs)

      berrycurv_unit=(hbar**4)/((em)**2)*(metertoang**4)*c**4
      ineq2d=0.0d0
      ineq2d_=0.0d0

      do ik=sum(ourjob(1:myrank))+1,sum(ourjob(1:myrank+1))
         call ener_read(ener,isp,ik,nk,nband)

         Q11=(0.0d0,0.0d0)
         Q22=(0.0d0,0.0d0)
         Q12=(0.0d0,0.0d0)
         Q21=(0.0d0,0.0d0)

         e0=0.5d0*(ener(nn)+ener(mm))
         mingap=huge(1.0d0)
         nskip_gap=0

         wk(:)=wklist(:,ik)
         np=nplist(ik)
         call plindx(ig,ncnt,ispinor,wk,b1,b2,b3,nbmax,np, &
                      ecut,npmax)

         read(10,rec=(3+(ik-1)*(nband+1)+ &
              nk*(nband+1)*(isp-1)+nn))(coeff(i),i=1,np)
         coeffn=coeff
         coeff=(0.0,0.0)

         read(10,rec=(3+(ik-1)*(nband+1)+ &
              nk*(nband+1)*(isp-1)+mm))(coeff(i),i=1,np)
         coeffm=coeff
         coeff=(0.0,0.0)

         do kk=1,nband
            use_kk=0
            if(kk.ne.nn .and. kk.ne.mm) then
               if(ext_mode.eq.0) then
                  use_kk=1
               else if(ext_mode.eq.1) then
                  if(kk.ge.ext_bmin .and. kk.le.ext_bmax) use_kk=1
               else if(ext_mode.eq.2) then
                  if(kk.lt.ext_bmin .or. kk.gt.ext_bmax) use_kk=1
               endif
            endif
            if(use_kk.eq.1) then
               gap=abs(e0-ener(kk))
               if(gap.lt.mingap) mingap=gap
               if(gap.le.gap_tol) then
                  nskip_gap=nskip_gap+1
               else

               read(10,rec=(3+(ik-1)*(nband+1)+ &
                    nk*(nband+1)*(isp-1)+kk))(coeff(i),i=1,np)
               coeffk=coeff
               coeff=(0.0,0.0)

               v1_nk=(0.0d0,0.0d0)
               v2_nk=(0.0d0,0.0d0)
               v1_mk=(0.0d0,0.0d0)
               v2_mk=(0.0d0,0.0d0)
               v1_kn=(0.0d0,0.0d0)
               v2_kn=(0.0d0,0.0d0)
               v1_km=(0.0d0,0.0d0)
               v2_km=(0.0d0,0.0d0)

               do iplane=1,ncnt
                  kgcart(1)=(wk(1)+ig(1,iplane))*b1(1)+ &
                            (wk(2)+ig(2,iplane))*b2(1)+ &
                            (wk(3)+ig(3,iplane))*b3(1)
                  kgcart(2)=(wk(1)+ig(1,iplane))*b1(2)+ &
                            (wk(2)+ig(2,iplane))*b2(2)+ &
                            (wk(3)+ig(3,iplane))*b3(2)
                  kgcart(3)=(wk(1)+ig(1,iplane))*b1(3)+ &
                            (wk(2)+ig(2,iplane))*b2(3)+ &
                            (wk(3)+ig(3,iplane))*b3(3)

                  xkg1=(kgcart(1)*bdir1(1)+kgcart(2)*bdir1(2)+ &
                       kgcart(3)*bdir1(3))/(bscale*bscale)
                  xkg2=(kgcart(1)*bdir2(1)+kgcart(2)*bdir2(2)+ &
                       kgcart(3)*bdir2(3))/(bscale*bscale)

                  if(ispinor.eq.2) then
                     v1_nk=v1_nk+ &
                        conjg(coeffn(iplane))*xkg1*coeffk(iplane)+ &
                        conjg(coeffn(iplane+ncnt))*xkg1* &
                        coeffk(iplane+ncnt)
                     v2_nk=v2_nk+ &
                        conjg(coeffn(iplane))*xkg2*coeffk(iplane)+ &
                        conjg(coeffn(iplane+ncnt))*xkg2* &
                        coeffk(iplane+ncnt)

                     v1_mk=v1_mk+ &
                        conjg(coeffm(iplane))*xkg1*coeffk(iplane)+ &
                        conjg(coeffm(iplane+ncnt))*xkg1* &
                        coeffk(iplane+ncnt)
                     v2_mk=v2_mk+ &
                        conjg(coeffm(iplane))*xkg2*coeffk(iplane)+ &
                        conjg(coeffm(iplane+ncnt))*xkg2* &
                        coeffk(iplane+ncnt)

                     v1_kn=v1_kn+ &
                        conjg(coeffk(iplane))*xkg1*coeffn(iplane)+ &
                        conjg(coeffk(iplane+ncnt))*xkg1* &
                        coeffn(iplane+ncnt)
                     v2_kn=v2_kn+ &
                        conjg(coeffk(iplane))*xkg2*coeffn(iplane)+ &
                        conjg(coeffk(iplane+ncnt))*xkg2* &
                        coeffn(iplane+ncnt)

                     v1_km=v1_km+ &
                        conjg(coeffk(iplane))*xkg1*coeffm(iplane)+ &
                        conjg(coeffk(iplane+ncnt))*xkg1* &
                        coeffm(iplane+ncnt)
                     v2_km=v2_km+ &
                        conjg(coeffk(iplane))*xkg2*coeffm(iplane)+ &
                        conjg(coeffk(iplane+ncnt))*xkg2* &
                        coeffm(iplane+ncnt)
                  else if(ispinor.eq.1) then
                     v1_nk=v1_nk+conjg(coeffn(iplane))*xkg1* &
                          coeffk(iplane)
                     v2_nk=v2_nk+conjg(coeffn(iplane))*xkg2* &
                          coeffk(iplane)

                     v1_mk=v1_mk+conjg(coeffm(iplane))*xkg1* &
                          coeffk(iplane)
                     v2_mk=v2_mk+conjg(coeffm(iplane))*xkg2* &
                          coeffk(iplane)

                     v1_kn=v1_kn+conjg(coeffk(iplane))*xkg1* &
                          coeffn(iplane)
                     v2_kn=v2_kn+conjg(coeffk(iplane))*xkg2* &
                          coeffn(iplane)

                     v1_km=v1_km+conjg(coeffk(iplane))*xkg1* &
                          coeffm(iplane)
                     v2_km=v2_km+conjg(coeffk(iplane))*xkg2* &
                          coeffm(iplane)
                  endif
               enddo

               denom=(e0-ener(kk))**2

               Q11(1,1)=Q11(1,1)+v1_nk*v1_kn/denom*berrycurv_unit
               Q11(1,2)=Q11(1,2)+v1_nk*v1_km/denom*berrycurv_unit
               Q11(2,1)=Q11(2,1)+v1_mk*v1_kn/denom*berrycurv_unit
               Q11(2,2)=Q11(2,2)+v1_mk*v1_km/denom*berrycurv_unit

               Q22(1,1)=Q22(1,1)+v2_nk*v2_kn/denom*berrycurv_unit
               Q22(1,2)=Q22(1,2)+v2_nk*v2_km/denom*berrycurv_unit
               Q22(2,1)=Q22(2,1)+v2_mk*v2_kn/denom*berrycurv_unit
               Q22(2,2)=Q22(2,2)+v2_mk*v2_km/denom*berrycurv_unit

               Q12(1,1)=Q12(1,1)+v1_nk*v2_kn/denom*berrycurv_unit
               Q12(1,2)=Q12(1,2)+v1_nk*v2_km/denom*berrycurv_unit
               Q12(2,1)=Q12(2,1)+v1_mk*v2_kn/denom*berrycurv_unit
               Q12(2,2)=Q12(2,2)+v1_mk*v2_km/denom*berrycurv_unit

               Q21(1,1)=Q21(1,1)+v2_nk*v1_kn/denom*berrycurv_unit
               Q21(1,2)=Q21(1,2)+v2_nk*v1_km/denom*berrycurv_unit
               Q21(2,1)=Q21(2,1)+v2_mk*v1_kn/denom*berrycurv_unit
               Q21(2,2)=Q21(2,2)+v2_mk*v1_km/denom*berrycurv_unit
               endif
            endif
         enddo

         if(mingap.eq.huge(1.0d0)) mingap=0.0d0

         g11=0.5d0*real(Q11(1,1)+Q11(2,2))
         g22=0.5d0*real(Q22(1,1)+Q22(2,2))
         g12=0.25d0*real(Q12(1,1)+Q21(1,1)+ &
                         Q12(2,2)+Q21(2,2))

         Omega=zi*(Q12-Q21)
         trOmega=Omega(1,1)+Omega(2,2)
         Omega_su2=Omega
         Omega_su2(1,1)=Omega(1,1)-0.5d0*trOmega
         Omega_su2(2,2)=Omega(2,2)-0.5d0*trOmega

         omega_norm2=0.5d0*real( &
              Omega_su2(1,1)*Omega_su2(1,1) + &
              Omega_su2(1,2)*Omega_su2(2,1) + &
              Omega_su2(2,1)*Omega_su2(1,2) + &
              Omega_su2(2,2)*Omega_su2(2,2))

         detg=g11*g22-g12*g12
         if(detg.lt.0.0d0 .and. abs(detg).lt.1.0d-12) detg=0.0d0
         if(omega_norm2.lt.0.0d0 .and. abs(omega_norm2).lt.1.0d-12) &
              omega_norm2=0.0d0

         lhs=sqrt(max(detg,0.0d0))
         rhs=0.5d0*sqrt(max(omega_norm2,0.0d0))
         margin=lhs-rhs

         ineq2d(ik,1)=g11
         ineq2d(ik,2)=g22
         ineq2d(ik,3)=g12
         ineq2d(ik,4)=detg
         ineq2d(ik,5)=lhs
         ineq2d(ik,6)=omega_norm2
         ineq2d(ik,7)=rhs
         ineq2d(ik,8)=margin
         ineq2d(ik,9)=real(Omega_su2(1,1))
         ineq2d(ik,10)=aimag(Omega_su2(1,1))
         ineq2d(ik,11)=real(Omega_su2(1,2))
         ineq2d(ik,12)=aimag(Omega_su2(1,2))
         ineq2d(ik,13)=real(Omega_su2(2,1))
         ineq2d(ik,14)=aimag(Omega_su2(2,1))
         ineq2d(ik,15)=real(Omega_su2(2,2))
         ineq2d(ik,16)=aimag(Omega_su2(2,2))
         ineq2d(ik,17)=mingap
         ineq2d(ik,18)=dble(nskip_gap)

         write(6,'(A,I6,2I4,3F12.6,18ES16.6E3)') &
              '# QGT_INEQ_2D_DVEC ik frac, values: ', &
              ik,wk,ineq2d(ik,:)
      enddo

#ifdef MPI_USE
      call MPI_ALLREDUCE(ineq2d,ineq2d_,18*nk,MPI_REAL8, &
                         MPI_SUM,mpi_comm_earth,mpierr)
      ineq2d=ineq2d_
#endif

      return
      end subroutine kubo_qgt_ineq_2d_dvec_kramers_pair
