!This function is to calculate determinant of the complex matrix
!The source is adoped from : https://dualm.wordpress.com/2012/01/06/computing-determinant-in-fortran/
      subroutine get_det(determinant, mat, N)
         implicit none
         integer*4, intent(in) :: N
         complex*16  mat(N,N)
         integer*4   i, info
         integer*4   ipiv(N)
         complex*16  determinant
         real*8      sgn

         ipiv = 0
         call zgetrf(N, N, mat, N, ipiv, info)

         determinant = (1d0,0d0)
         do i = 1, N
            determinant = determinant*mat(i, i)
         end do

         sgn = 1d0
         do i = 1, N
            if(ipiv(i) /= i) then
               sgn = -sgn
            end if
         end do
         determinant = sgn*determinant

      end subroutine get_det

!$$*  subroutine to get det(A) where A is n x n complex matrix
!   This subroutine is adopted from below,
!   "http://computer-programming-forum.com/49-fortran/9e079d718158c944.htm"
!   The author is : Che-Ping Su
!c****************************************************
!    Description:     Calculate the determinant of
!                     general complex matrix
!    Input:           A - matrix to be calculated
!                     N - the order of matrix A
!    Output:          DT - Determinant of matrix A
!    Last Updated:      03/17/1998
!****************************************************
      subroutine getdetA(DT, S,N)
         implicit real*8(a-h,o-z)
         complex*16  S(N,N),DT,TM,TC,W(N,1)
         complex*16  A(N,N)
         A=(0.,0.)
         A=S
         L=1
         K=2
         flag=(0.,0.)
         DT=cmplx(1.)
   10    TM=A(L,L)
         if(A(L,L) .eq. (0.0,0.0)) then
            do I= L+1,N
               if(A(L,I) .ne. (0.0,0.0)) then
                  do J=1,N
                     w(J,1)=A(J,L)
                     A(J,L)=A(J,I)
                     A(J,I)=W(J,1)
                  enddo
                  flag=flag+(1.,0.)
                  goto 10
               endif
            enddo
         endif

         do 20 J=L,N
            A(L,J)=A(L,J)/TM
   20    continue
         do 30 I=k,N
            TC=A(I,L)
            do 30 J=L,N
               A(I,J)=A(I,J)-A(L,J)*TC
   30    continue
         L=L+1
         k=k+1
         DT=DT*TM
         if(L-N) 10,40,40
   40    DT=(-1.,0.)**flag*DT*A(N,N)
         return
      end subroutine getdetA

      subroutine mpi_job_distribution_chain(njob, ourjob, ourjob_disp,&
                                           nprocs)
         implicit none
         integer*4    njob
         integer*4    mynjob, nprocs
         integer*4    cpuid, mpierr
         integer*4    nresidue
#ifdef MPI_USE
         integer*4    ourjob(nprocs)
         integer*4    ourjob_disp(0:nprocs-1)
#else
         integer*4    ourjob(1)
         integer*4    ourjob_disp(0)
#endif

         mynjob = floor ( real(njob)/real(nprocs) )
         nresidue = nint (real(njob) - real(mynjob) * real(nprocs))

         ourjob = 0

         do cpuid = 1, nprocs
            if( cpuid .le. nresidue ) then
               ourjob(cpuid) = mynjob + 1
            else
               ourjob(cpuid) = mynjob
            endif
         enddo

         ourjob_disp(0) = 0
#ifdef MPI_USE
         do cpuid = 1, nprocs-1
            ourjob_disp(cpuid)= ourjob_disp(cpuid - 1) + ourjob(cpuid)
         enddo
#endif

      end subroutine mpi_job_distribution_chain

      subroutine find_jk(nklist,iklist,ik, jk)
         implicit none
         integer*4    nklist, ik, jk, kk
         integer*4    iklist(nklist)
         do kk = 1, nklist
            if(iklist(kk) .eq. ik) then
               jk = kk
               return
            endif
         enddo
      end subroutine find_jk

      function int2str(w) result(string)
         implicit none
         character*20  string
         integer*4,    intent(in)  :: w
         write(string,*) w
         return
      end function int2str
