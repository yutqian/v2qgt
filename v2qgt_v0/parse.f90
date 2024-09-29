!!$   parse command line arguments
      subroutine parse(filename,foname,nkx,nky,ispinor,icd,ixt,fbz,&
           ivel,ikubo,iqm,iz,ihf,nini,nmax,nn,kperiod,it,iskp,ine,&
           ver_tag,&
           iwf,ikwf,ng,rs,imag,init_e,fina_e,nediv,sigma,&
           klist_fname,sw_fname,atlist_fname,flag_atom_project,&
           theta,phi)
          implicit real*8(a-h,o-z)
          character*256 filename,foname,fbz,ver_tag,vdirec
          character*256 klist_fname, sw_fname, atlist_fname
          real*8 x,y
          character*256 option,value
          integer iarg,narg,ia,nkx,nky,ispinor,iskp,ine,ng(3)
          dimension rs(3)
          real*8    init_e, fina_e
          real*8    theta, phi
          logical   flag_atom_project
          nini=1;it=0;iskp=0;ine=0;icd=0;ixt=0;ivel=0;iz=0;ihf=0
          iwf=0;ikwf=1;ng=0;imag=0;rs=0.;ikubo=0;nn=0;nediv=1000
          init_e =  0.0d0;fina_e=10.0d0;sigma=0.01;iqm=0
          theta = 0.0d0 ; phi = 0.0d0
          nmax=999999
          iarg=iargc()
          nargs=iarg/2
          filename="WAVECAR"
          klist_fname="sw_klist.dat"
          atlist_fname="sw_atoms.dat"
          sw_fname="sw_spin0.dat" ! note that collinear case is not supported yet. 09.Oct. 2021, HJ Kim
          foname="BERRYCURV"
          fbz="BERRYCURV.tot.dat"
          flag_atom_project = .false.
          if(iarg.ne.2*nargs) then
          call help(ver_tag)
          endif
          do ia=1,nargs
          call getarg(2*ia-1,option)
          call getarg(2*ia,value)
          if(option == "-f") then
              !read(value,*) filename
              filename = trim(value)
              else if(option == "-o") then
              read(value,*) foname
              else if(option == "-kx") then
              read(value,*) nkx
              else if(option == "-ky") then
              read(value,*) nky
              else if(option == "-s") then
              read(value,*) ispinor
              else if(option == "-ii") then
              read(value,*) nini
              else if(option == "-if") then
              read(value,*) nmax
              else if(option == "-is") then
              read(value,*) nini
              nini=nini;nmax=nini
              else if(option == "-kp") then
              read(value,*) kperiod
              else if(option == "-t") then
              read(value,*) it
              else if(option == "-skp") then
              read(value,*) iskp
              else if(option == "-ne") then
              read(value,*) ine
              else if(option == "-ixt") then
              read(value,*) ixt
              else if(option == "-fbz") then
              read(value,*) fbz
              else if(option == "-cd") then
              read(value,*) icd
              else if(option =="-klist") then
              read(value,*) klist_fname
              else if(option =="-sw_file") then
              read(value,*) sw_fname
              else if(option =="-atlist") then
              flag_atom_project = .true.
              read(value,*) atlist_fname
              else if(option == "-ien") then
              read(value,*) init_e
              else if(option == "-fen") then
              read(value,*) fina_e
              else if(option == "-theta") then
              read(value,*) theta
              else if(option == "-phi") then
              read(value,*) phi
              else if(option == "-nediv") then
              read(value,*) nediv
              else if(option == "-sigma") then
              read(value,*) sigma
              else if(option == "-kubo") then
              read(value,*) ikubo
              else if(option == "-qm") then
              read(value,*) iqm
              else if(option == "-nn") then ! band index nn
              read(value,*) nn
              else if(option == "-vel") then
              read(value,*) ivel
              else if(option == "-z2") then
              read(value,*) iz
              else if(option == "-hf") then
              read(value,*) ihf
              else if(option == "-wf") then
              read(value,*) iwf
              else if(option == "-k") then
              read(value,*) ikwf
              else if(option == "-ng") then
              read(value(:),*)ng(1:3)
              else if(option == "-ishift") then
              read(value(:),*)rs(1:3)
              else if(option == "-im") then
              read(value,*) imag
              else if(option =="-h") then
                  call help(ver_tag)
              else
              call help(ver_tag)
              endif
          enddo
          if(icd.gt.1 .and. TRIM(foname) .ne. 'BERRYCURV' )then
          if(icd .eq.1) then
          write(foname,'(A,A)')"CIRC_DICHROISM.",TRIM(foname)
          elseif(icd.eq.2) then
          write(foname,'(A,A)')"CIRC_DICHROISM_W.",TRIM(foname)
          endif
          else if (icd .gt. 1 .and. TRIM(foname) .eq. 'BERRYCURV') then
          if(icd .eq.1) then
              foname="CIRC_DICHROISM"
          elseif(icd.eq.2) then
              foname="OPT_TRANS_RATE"
          elseif(icd.eq.3) then
              foname="OPT_TRANS_RATE_UNFOLD"
          endif
          else if (icd+ivel .eq. 0 .and. TRIM(foname) .ne. 'BERRYCURV')then
          write(foname,'(A,A)')"BERRYCURV.",TRIM(foname)
          else if (ivel .eq. 1 .and. TRIM(foname) .ne. 'BERRYCURV') then
          write(foname,'(A,A)')"VEL_EXPT.",TRIM(foname)
          else if (ivel .eq. 1 .and. TRIM(foname) .eq. 'BERRYCURV') then
          foname="VEL_EXPT"
          else if (iz .eq. 1 .and. TRIM(foname) .eq. 'BERRYCURV') then
          foname="NFIELD"
          else if (ikubo .ge. 1 .and. TRIM(foname) .eq. 'BERRYCURV') then
          foname="BERRYCURV_KUBO"
          endif
  
          if(iwf .ge. 1)then
          ine=iwf
          if(iwf .lt. 10) then
          if(ikwf .lt. 10) then
          write(foname,'(A,I1,A,I1)')"PARCHG-W-K00",ikwf,"-E00",iwf
          elseif(ikwf .ge. 10 .and. ikwf .lt. 100)then
          write(foname,'(A,I2,A,I1)')"PARCHG-W-K0",ikwf,"-E00",iwf
          elseif(ikwf .ge. 100 ) then
          write(foname,'(A,I3,A,I1)')"PARCHG-W-K",ikwf,"-E00",iwf
          endif
          elseif(iwf .ge. 10 .and. iwf .lt. 100)then
          if(ikwf .lt. 10) then
          write(foname,'(A,I1,A,I2)')"PARCHG-W-K00",ikwf,"-E0",iwf
          elseif(ikwf .ge. 10 .and. ikwf .lt. 100)then
          write(foname,'(A,I2,A,I2)')"PARCHG-W-K0",ikwf,"-E0",iwf
          elseif(ikwf .ge. 100 ) then
          write(foname,'(A,I3,A,I2)')"PARCHG-W-K",ikwf,"-E0",iwf
          endif
          elseif(iwf .ge. 100 ) then
          if(ikwf .lt. 10) then
          write(foname,'(A,I1,A,I3)')"PARCHG-W-K00",ikwf,"-E",iwf
          elseif(ikwf .ge. 10 .and. ikwf .lt. 100)then
          write(foname,'(A,I2,A,I3)')"PARCHG-W-K0",ikwf,"-E",iwf
          elseif(ikwf .ge. 100 ) then
          write(foname,'(A,I3,A,I3)')"PARCHG-W-K",ikwf,"-E",iwf
          endif
          endif
          endif
  
          return
      end subroutine parse
  
