!#define MPI_USE
!#undef  MPI_USE        

    PROGRAM VASPBERRY

    implicit real*8 (a-h, o-z)
    complex*8, allocatable :: coeff(:)
    complex*16, allocatable :: Siju(:,:),Sijd(:,:),Sijt(:,:)
    complex*16, allocatable :: coeff1u(:),coeff1d(:)
    complex*16, allocatable :: coeff2u(:),coeff2d(:)
    complex*16,allocatable :: cener(:)
    real*8,    allocatable :: berrycurv(:),recivec(:,:),wklp(:,:,:)
    real*8,    allocatable :: berrycurv_tot(:)
    real*8,    allocatable :: rnfield(:),rnfield_tot(:)
    real*8,    allocatable :: recivec_tot(:,:)
    real*8,    allocatable :: recilat_tot(:,:)
    real*8,    allocatable :: xrecivec(:,:),xberrycurv(:),wklist(:,:)
    real*8,    allocatable :: wnklist(:,:),xrnfield(:)
    real*8,    allocatable :: recilat(:,:),xrecilat(:,:),occ(:)
    real*8,    allocatable :: selectivity(:),xselectivity(:)
    real*8,    allocatable :: selectivity_w(:)
    real*8,    allocatable :: spectrum(:,:),xspectrum(:,:)
    real*8,    allocatable :: berrycurv_kubo(:,:),xberrycurv_kubo(:,:)
    real*8,    allocatable :: berrycurv_kubo_tot(:,:)
    real*8,    allocatable :: e_range(:) 
    real*16,   allocatable :: ener(:)
    integer,   allocatable :: ig(:,:),nplist(:)
    integer*4, allocatable :: iatlist(:)! for proj band
    integer*4, allocatable :: iklist(:) ! for unfolding band
    integer*4                 nklist    ! for unfolding band
    integer*4                 jk        ! for unfolding band
    real*8,    allocatable :: pb(:,:,:) ! projected band 
    real*8,    allocatable :: sw(:,:,:) ! sw for unfolded band
    real*8,    allocatable :: sw_recivec(:,:) ! PBZ for unfold (cartesian)
    real*8,    allocatable :: sw_recilat(:,:) ! PBZ for unfold (fractional)
    logical                   flag_atom_project
    integer*4                 ie, natom, natlist
    dimension selectivitymax(4),selectivitymin(4)
    dimension a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),a2xa3(3),sumkg(3)
    dimension wk(3),wkk(3,5),ikk(5),isgg(2,5),npl(5),itr(5),itrim(5)
    dimension nbmax(3),xb(3),berrymax(4),berrymin(4),ng(3),rs(3)
    complex*16 csum1,csum2
    complex*16  detS(4),detA,detLOOP
    integer k, n, nkx, nky,nini,nmax,ns,ne,icd,ivel
    character*256 filename,foname,fonameo,fbz,ver_tag,vdirec
    character*256 klist_fname,sw_fname, atlist_fname, proj_fname
    character*20,external :: int2str
    character*20 dummy
    data c/0.262465831d0/ ! constant c = 2m/hbar**2 [1/eV Ang^2]
    real*8  rfield,rnnfield,rnnfield_bottom
    real*8  rnnfield_tot,rnnfield_bottom_tot
    real*8  init_e, fina_e
    real*8  theta, phi ! angle of incident light (theta: along z, phi: along x )
    real*8, allocatable:: w_half_klist(:,:)
    integer, allocatable:: i_half_klist(:),i_trim_klist(:)
    integer :: myrank, nprocs, ierr, mpierr
    integer ::  mpi_comm_earth
#ifdef MPI_USE
    include 'mpif.h'
    INTEGER         status(MPI_STATUS_SIZE)
    call MPI_INIT(ierr)
    mpi_comm_earth = MPI_COMM_WORLD
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
    if(myrank == 0)then
     write(6,*)"THIS IS ROOT:",myrank
    endif
    if(myrank == 0)then
     time_1=MPI_WTIME()
    endif
#else
    nprocs=1
    myrank=0
    mpi_comm_earth = 0
#endif

    ver_tag="# VASPBERRY (Ver 1.0), by Hyun-Jung Kim. 2018. Aug. 23."
    pi=4.*atan(1.)

    !default settings
    kperiod=2
    nkx=12 ;nky=12 
    ispinor=2  ! 2 for soc, 1 for non-soc
    itr=0;itrim=0

!!$*  reading general informations
    call parse(filename,foname,nkx,nky,ispinor,icd,ixt,fbz, &
       ivel,ikubo,iqm,iz,ihf,nini,nmax,nn,kperiod,it,iskp,ine,ver_tag,&
       iwf,ikwf,ng,rs,imag,init_e,fina_e,nediv,sigma,&
       klist_fname,sw_fname,atlist_fname,flag_atom_project,&
       theta,phi)
    !! if(myrank == 0)call creditinfo(ver_tag)

    if (ixt .ne. 0) call extendingBZ(fbz,ixt)
    if (it .eq. 1) call test
    if(myrank == 0)write(6,'(A,A)')"# File reading... : ",filename
    call inforead(irecl,ispin,nk,nband,ecut,a1,a2,a3,filename)
    if (ispin .eq. 2) ispinor=1
    allocate(cener(nband),ener(nband),occ(nband))
    if(myrank == 0)write(6,'(A,I9)')"# TOTAL RECORD LENGTH = ",irecl
    allocate(wklist(3,nk),nplist(nk))
    if(iz == 0)then
      iz2=1 ! enhancing factor for variable size define in subroutines
      allocate(berrycurv(nk))
      allocate(berrycurv_tot(nk))
      allocate(xberrycurv(kperiod*2*kperiod*2*nk))
    endif


    ! read wklist and nplist
    do ik=1,nk 
        ne_temp0=ne
        ne=0
        irec=3+(ik-1)*(nband+1)
        read(10,rec=irec) xnplane, (wk(i),i=1,3), &
                        (cener(inn),occ(inn),inn=1,nband)
        wklist(:,ik)=wk(:)
        nplist(ik)=nint(xnplane)
        
        do n=1,nband; ne=ne+nint(occ(n)); enddo
        
        if(ik .gt. 1 .and. ne_temp0 .ne. ne .and. ine .eq. 0) then
          write(6,*)"error. !!! ne(K) /= ne(K') !!!",ne,ne_temp0,ik ;stop
        endif
       
    enddo

    if(ine .ne. 0) ne=ine ! manually specified ne ; useful for the semimetal
    ! check whether multi or single band calculation is performed
    if((nini.eq.nmax))then
        nini=nmax
    else if(nmax .eq. 999999)then
        if(ispin .eq. 2 .and. icd .eq. 0) nmax=ne
        if(ispin .eq. 1 .and. ispinor .eq. 2 .and. icd .eq. 0) nmax=ne
        if(ispin .eq. 1 .and. ispinor .eq. 1 .and. icd .eq. 0) nmax=ne/2
        if(ispin .eq. 2 .and. icd .ge. 1)then;nini=ne;nmax=nini+1;endif
        if(ispin .eq. 1 .and. ispinor .eq. 1 .and. icd .eq. 1)then
         nini=ne/2
         nmax=nini+1
        elseif(ispin .eq. 1 .and. ispinor .eq. 1 .and. icd .gt. 1) then
         nini=ne/2
         nmax=nband
        endif
        
        if(ispinor .eq. 2 .and. icd .eq. 1)then
         nini=ne
         nmax=nini+1
        elseif(ispinor .eq. 2 .and. icd .gt. 1) then
         nini=ne
         nmax=nband
        endif
    endif ! check multi or single ?

    !>----write information
    ns=nmax-nini+1
    if(myrank == 0)then
        write(6,'(A,I6)')   "# NELECT     : ",ne*ispin
        if (ispinor .eq. 2)then
            write(6,'(A,I6,A)')"# ISPIN      : ",ispin," (LSORBIT =.TRUE.)"
        else
            write(6,'(A,I6,A)')"# ISPIN      : ",ispin," (LSORBIT =.FALSE.)"
        endif

        write(6,'(A,F11.4)')  "# ENCUT (eV) : ",ecut
        write(6,'(A,I6)')     "# NKPOINT    : ",nk
        write(6,'(A,I6,A,I4)')"#  K-GRID    : ",nkx,"   X",nky
        write(6,'(A,I6)')     "# NBANDS     : ",nband
        write(6,'(A,3F13.6)') "# LATTVEC A1 : ",(a1(i),i=1,3)
        write(6,'(A,3F13.6)') "# LATTVEC A2 : ",(a2(i),i=1,3)
        write(6,'(A,3F13.6)') "# LATTVEC A3 : ",(a3(i),i=1,3)
    endif
    call recilatt(b1,b2,b3,dSkxky, a1,a2,a3,nkx,nky)
    
    if(myrank == 0)then
        write(6,'(A,3F13.6)') "# RECIVEC B1 : ",(b1(i),i=1,3)
        write(6,'(A,3F13.6)') "# RECIVEC B2 : ",(b2(i),i=1,3)
        write(6,'(A,3F13.6)') "# RECIVEC B3 : ",(b3(i),i=1,3)
        if(icd.eq.0)write(6,'(A,F13.6)')  "#  dk^2 = |dk1xk2| = ",dSkxky
    endif
    call reciproperty(nbmax,npmax, b1,b2,b3,ecut,ispinor)
    if(myrank == 0)then
        write(6,'(A,I6)')     "# NPMAX      : ",npmax
        write(6,'(A,I6)')     "#  NB1MAX    : ",nbmax(1)
        write(6,'(A,I6)')     "#  NB2MAX    : ",nbmax(2)
        write(6,'(A,I6)')     "#  NB3MAX    : ",nbmax(3)
    endif
    !>----write information
 #ifdef MPI_USE
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 #endif


    allocate(ig(3,npmax))
    allocate(coeff(npmax))
    if(iz == 0)then
     allocate(wklp(3,5,nk))
     allocate(recivec(3,nk),recilat(3,nk))
     allocate(recivec_tot(3,nk),recilat_tot(3,nk))
     allocate(xrecivec(3,kperiod*2*kperiod*2*nk))
     allocate(xrecilat(3,kperiod*2*kperiod*2*nk))
    endif
    
    nbtot=(2*nbmax(1)+2)*(2*nbmax(2)+2)*(2*nbmax(3)+1)
    allocate(coeff1u(nbtot),coeff1d(nbtot))
    allocate(coeff2u(nbtot),coeff2d(nbtot))
    if(icd.ge.1) allocate(selectivity(nk))
    if(icd.ge.1) allocate(selectivity_w(nk))
    if(icd.ge.2) allocate(spectrum(nediv,nk))
    if(icd.ge.1) allocate(xselectivity(kperiod*2*kperiod*2*nk))
    if(icd.ge.2) allocate(xspectrum(nediv,kperiod*2*kperiod*2*nk))
    if(icd.ge.2) allocate(e_range(nediv))
    if(ikubo.ge.1) allocate(berrycurv_kubo(nk,3))
    if(ikubo.ge.1) allocate(xberrycurv_kubo(kperiod*2*kperiod*2*nk,3))
    if(ikubo.ge.1) allocate(berrycurv_kubo_tot(nk,3))
    !! if(ivel .eq. 1 .and. myrank==0) then
    !!     ni=nini 
    !!     nj=nini
    !!     call vel_expectation(a1,a2,a3,b1,b2,b3,&
    !!      kperiod,nbmax,npmax,ecut,ispinor,ispin,nband,ne,nk,wklist,&
    !!      nplist,ni,nj,irecl,filename,foname)
    !! endif


!!! ######### LOOP for Berry curvature using Kubo formula ###############################
    do isp=1,ispin   ! ispin start
!       elseif(ikubo.eq.1.and. myrank==0) then
        if(ikubo.ge.1               ) then
#ifdef MPI_USE
            if(myrank == 0)then
               time_2=MPI_WTIME()
            endif
#endif
            chernnumber_total=0d0
            chernnumber=0d0
            berrycurv_kubo_tot = 0
            do ie=nini, nmax
               if(myrank .eq. 0) then
                  write(6,'(A  )')" "
                  write(6,'(A,I)')"# BERRY CURVATURE FOR BAND INDEX n= ", ie
                  write(6,'(A,I)')"#                     SPIN INDEX s= ", isp
                  write(6,'(A  )')" "
               endif
#ifdef MPI_USE
                call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
! calculate BC for n-th band and write it into files 
                berrycurv_kubo = 0d0
                call kubo_berry_curvature(berrycurv_kubo,     &
                       b1,b2,b3,wklist,isp,                   &
                       nband,ecut,ispinor,nplist,nbmax,npmax, &
                       nk,ie,nmax,iqm,                           &
                       nprocs,myrank,mpi_comm_earth)


                if(myrank .eq. 0) then
                   do ik=1,nk
!                     write(6,'(A,3F12.4)') 'BC_kubo:',berrycurv_kubo(ik,:)
                      do j=1,3
                        recivec(j,ik)=wklist(1,ik)*b1(j)+   &
                                      wklist(2,ik)*b2(j)+   &
                                      wklist(3,ik)*b3(j)
                        recilat(j,ik)=wklist(j,ik)
                      enddo
                   enddo

                   call get_ext_variable(xrecivec,xrecilat,xberrycurv_kubo, &
                                   kext,recivec,recilat,berrycurv_kubo,     &
                                   nk,kperiod,nk,iz2,                       &
                                   b1,b2,b3) ! extending over exteded-BZ

                   !!xberrycurv_kubo = berrycurv_kubo
                   !if(ikubo .eq. 1) then
                   !   call get_sorted_xrvari(xrecivec,xrecilat,xberrycurv_kubo, &
                   !                          kext,kperiod,nk,iz2) ! sorting
                   !endif
                   call write_result(isp,ispin,ispinor,fonameo,foname,filename, &
                                  irecl,ecut,nk,nkx,nky,nband,b1,b2,b3,kperiod, &
                                  dSkxky,ie,ie    ,xrecivec,xrecilat,kext,      &
                                  xberrycurv_kubo,chernnumber,0.,0.,0.,         &
                                  iz2,icd,iz,ivel,ikubo,nprocs)
                endif ! end myrank
!            enddo ! loop ie

! for total berry curvature, this part is tested
!           do ie=nini, nmax
            berrycurv_kubo = 0d0
            call kubo_berry_curvature_tot(berrycurv_kubo,     &
                   b1,b2,b3,wklist,isp,                   &
                   nband,ecut,ispinor,nplist,nbmax,npmax, &
                   nk,ie,nmax,iqm,                        &
                   nprocs,myrank,mpi_comm_earth)


                if(myrank .eq. 0) then
!! chernnumber is not tested, it could be wrong, don't trust it
                   chernnumber=0d0
                   do ik=1,nk
                      do j=1,3
                        recivec(j,ik)=wklist(1,ik)*b1(j)+   &
                                      wklist(2,ik)*b2(j)+   &
                                      wklist(3,ik)*b3(j)
                        recilat(j,ik)=wklist(j,ik)
                      enddo
                      berrycurv_kubo_tot(ik,:) = berrycurv_kubo_tot(ik,:) +         &
                                             berrycurv_kubo(ik,:)
                   enddo
                endif !myrank
            enddo ! loop ie

            if(myrank .eq. 0 .and. nini .ne. nmax) then
                call get_ext_variable(xrecivec,xrecilat,xberrycurv_kubo,kext, &
                                recivec,recilat,berrycurv_kubo_tot,           &
                                nk,kperiod,nk,iz2,                            &
                                b1,b2,b3) ! extending over exteded-BZ
                !xberrycurv_kubo = berrycurv_kubo_tot
                !if(ikubo .eq. 1) then
                !   call get_sorted_xrvari(xrecivec,xrecilat,xberrycurv_kubo, &
                !                        kext,kperiod,nk,iz2) ! sorting
                !endif
                call write_result(isp,ispin,ispinor,fonameo,foname,filename, &
                               irecl,ecut,nk,nkx,nky,nband,b1,b2,b3,kperiod, &
                               dSkxky,nini,nmax,xrecivec,xrecilat,kext,      &
                               xberrycurv_kubo,chernnumber_total,0.,0.,0.,   &
                               iz2,icd,iz,ivel,ikubo,nprocs)
            endif ! myrank
#ifdef MPI_USE
            if(myrank == 0)then
               time_3=MPI_WTIME()
            endif
#endif
        endif !ikubo over
    enddo !ispin loop over
!!! ######### LOOP END for Berry curvature using Kubo formula ###########################


    !! if spin up and dn need to write seperately
    if(iskp .eq. 1 .and. myrank==0)call write_special_kpoint(b1,b2,b3)
#ifdef MPI_USE
    if(myrank==0)then
#endif
    write(6,'(A)')"# DONE! "
    do isp=1, ispin
       if(isp .eq. 1 .and. ispinor .eq. 1 .and. ispin .eq. 2) then
          write(6,'(A,A,A)')"#  Results are summarized in ",TRIM(foname),&
                           ".UP.dat for spin-1"                           
      else if (isp.eq.2 .and. ispinor.eq.1 .and. ispin.eq.2)then          
         write(6,'(A,A,A)')"#  Results are summarized in ",              &
                           TRIM(foname),".DN.dat for spin-2"              
      else if (isp .eq. 1 .and. ispinor .eq. 2) then                      
         write(6,'(A,A,A)')"#  Results are summarized in ",              &
                           TRIM(foname),".dat"                            
      else if (isp.eq.1 .and. ispinor.eq.1 .and. ispin.eq.1) then         
         write(6,'(A,A,A)')"#  Results are summarized in ",              &
                           TRIM(foname),".dat"
       endif
    enddo
#ifdef MPI_USE
    endif
#endif
 9999    if(myrank==0) write(6,*)"end of program"

#ifdef MPI_USE
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

    deallocate(ig)
    deallocate(coeff)
    deallocate(coeff1u)
    deallocate(coeff1d)
    deallocate(coeff2u)
    deallocate(coeff2d)
    deallocate(recilat)
    deallocate(recivec)
    deallocate(xrecilat)
    deallocate(xrecivec)
    deallocate(wklp)
    deallocate(wklist)
    deallocate(cener)
    deallocate(ener)
    if(iz == 0)then
       deallocate(berrycurv)
       deallocate(berrycurv_tot)
       deallocate(xberrycurv)
    elseif(iz == 1)then
       deallocate(rnfield)
       deallocate(rnfield_tot)
       deallocate(wnklist)
       deallocate(xrnfield)
       deallocate(i_trim_klist)
       deallocate(i_half_klist)
       deallocate(w_half_klist)
    endif
    deallocate(occ)
    deallocate(nplist)
    if(icd.ge.1)deallocate(selectivity)
    if(icd.ge.1)deallocate(xselectivity)
#ifdef MPI_USE
    if(myrank == 0 .and. iz+ivel+icd+iwf .eq. 0)then
       time_4=MPI_Wtime()
       write(6,*) "Data reading      : ",time_2-time_1
       write(6,*) "Parallel sequence : ",time_3-time_1
       write(6,*) "End sequence      : ",time_4-time_3
    endif

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_FINALIZE(ierr)
#endif
    end PROGRAM VASPBERRY

!!! >>>> end main program <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


!!$*  subroutine for reading basic information
    subroutine inforead(irecl,ispin,nk,nband,ecut,a1,a2,a3,filename)
        implicit real*8(a-h,o-z)
        character*256 filename
        dimension a1(3),a2(3),a3(3)

        irecl=24
        open(unit=10,file=filename,access='direct',recl=irecl,&
         iostat=iost,status='old')
        if (iost .ne. 0) write(6,*) '0.open error - iostat =',iost

        read(10,rec=1)xirecl,xispin,xiprec !RDUM,RISPIN,RTAG(in real type)
        close(10)
        irecl=nint(xirecl);ispin=nint(xispin);iprec=nint(xiprec) ! set to integer
        if(iprec.eq.45210) then
         write(0,*) '*** error - WAVECAR_double requires complex*16';stop
        endif
        open(unit=10,file=filename,access='direct',recl=irecl, &
            iostat=iost,status='old')
        if (iost.ne.0) write(6,*) '1.open error - iostat =',iost
        read(10,rec=2) xnk,xnband,ecut,  &               !RNKPTS,RNB_TOT,ENCUT
            (a1(j),j=1,3),(a2(j),j=1,3),(a3(j),j=1,3)       !A1(3),A2(3),A3(3)
        nk=nint(xnk)
        nband=nint(xnband)

        return
    end subroutine inforead


