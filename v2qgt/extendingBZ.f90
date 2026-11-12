      subroutine extendingBZ(fbz,ixt)
         implicit real*8 (a-h, o-z)
         character*256  fbz
         character*200  A,S,P
         dimension b1(3),b2(3),b3(3),xb(3)
         real*8, allocatable :: xrecivec(:,:),xrecilat(:,:),xdata_(:)
         real*8, allocatable ::recivec(:,:),recilat(:,:),data_(:)
         integer IOstatus    

         IOstatus=0;II=1;ik=0
         open(51,file=fbz,status='unknown')
         open(61,file='EXT.dat',status='unknown')
         read(51,'(A)',IOSTAT=IOstatus)A
         S=A(1:1)
         do while (TRIM(S) .eq. '#' .or. II .eq. 35)
          write(6,'(A)')TRIM(A)
          write(61,'(A)')TRIM(A)
          read(51,'(A)',IOSTAT=IOstatus)A
          if (TRIM(A(3:9)) .eq. 'NKPOINT') then
           P=A(23:27);read(P,*)nk
           else if(TRIM(A(3:12)) .eq. 'RECIVEC B1') then
            P=A(26:34);read(P,*)b1(1)
            P=A(39:47);read(P,*)b1(2)
            P=A(52:60);read(P,*)b1(3)
           else if (TRIM(A(3:12)) .eq. 'RECIVEC B2') then
            P=A(26:34);read(P,*)b2(1)
            P=A(39:47);read(P,*)b2(2)
            P=A(52:60);read(P,*)b2(3)
           else if (TRIM(A(3:12)) .eq. 'RECIVEC B3') then
            P=A(26:34);read(P,*)b3(1)
            P=A(39:47);read(P,*)b3(2)
            P=A(52:60);read(P,*)b3(3)
          endif
          S=TRIM(A(1:1))
          II=ICHAR(TRIM(A))
         enddo
         allocate(recivec(3,nk))
         allocate(recilat(3,nk))
         allocate(data_(nk))
         allocate(xrecivec(3,ixt*2*ixt*2*nk))
         allocate(xrecilat(3,ixt*2*ixt*2*nk))
         allocate(xdata_(ixt*2*ixt*2*nk))

         backspace(51)
         do while(IOstatus .eq. 0)
          ik=ik+1
          read(51,*,IOSTAT=IOstatus)(recivec(i,ik),i=1,3),data_(ik),&
                                     (recilat(i,ik),i=1,3)
         enddo

         do ik=1,nk
            write(6,'(3F11.6,A,F16.6,A,3F11.6)')(recivec(i,ik),i=1,3), &
                "     ",data_(ik),&
                "                  ",(recilat(i,ik),i=1,3)
         enddo

         kk=0   ! extend 
         do ib2=-1*(ixt-1)+1,ixt
          do ib1=-1*(ixt-1)+1,ixt 
           do ik=1,nk
            kk=kk+1
            xrecivec(1,kk)=recivec(1,ik) +(ib1-1)*b1(1)+(ib2-1)*b2(1)
            xrecivec(2,kk)=recivec(2,ik) +(ib1-1)*b1(2)+(ib2-1)*b2(2)
            xrecivec(3,kk)=recivec(3,ik) +(ib1-1)*b1(3)+(ib2-1)*b2(3)
            xrecilat(1,kk)=recilat(1,ik) +(ib1-1)
            xrecilat(2,kk)=recilat(2,ik) +(ib2-1)
            xrecilat(3,kk)=recilat(3,ik)
            xdata_(kk)=data_(ik)
            kext=kk
           enddo
          enddo
         enddo

         write(6,*)" "
         write(6,'(A,I1,A,I1)')"# EXTENDING DATA GRID by : ",ixt,' x ',ixt
         write(6,'(A)')"# SORTING K-grids..."
         do k=ixt-1,1,-1  ! sorting kx
          do j=1,k
           if(xrecivec(1,j+1) .gt. xrecivec(1,j))then
            xb(:)=xrecivec(:,j)
            xrecivec(:,j)=xrecivec(:,j+1)
            xrecivec(:,j+1)=xb(:)
            xb(:)=xrecilat(:,j)
            xrecilat(:,j)=xrecilat(:,j+1)
            xrecilat(:,j+1)=xb(:)
            xtemp=xdata_(j)
            xdata_(j)=xdata_(j+1)
            xdata_(j+1)=xtemp
           endif
          enddo
         enddo
         do k=ixt-1,1,-1  ! sorting ky
          do j=1,k
           if(xrecivec(1,j+1) .eq. xrecivec(1,j))then
            if(xrecivec(2,j+1) .gt. xrecivec(2,j))then
            xb(:)=xrecivec(:,j)
            xrecivec(:,j)=xrecivec(:,j+1)
            xrecivec(:,j+1)=xb(:)
            xb(:)=xrecilat(:,j)
            xrecilat(:,j)=xrecilat(:,j+1)
            xrecilat(:,j+1)=xb(:)
            xtemp=xdata_(j)
            xdata_(j)=xdata_(j+1)
            xdata_(j+1)=xtemp
            endif
           endif
          enddo
         enddo

         do ik=1,kext
             write(61,'(3F11.6,A,F16.6,A,3F11.6)')(xrecivec(i,ik),i=1,3),&
                "     ",xdata_(ik),"                  ", &
                (xrecilat(i,ik),i=1,3)
         enddo
         write(6,'(A)')"# DONE! result is in 'EXT.dat' "

         stop
      end subroutine extendingBZ
