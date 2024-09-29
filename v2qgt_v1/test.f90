      subroutine test
          implicit real*8 (a-h, o-z)
          complex*8   a,b,c,d
          character*256  foname
          a=(2.2,-1.3)
          b=(3.2,-5.4)
      
          c=a+(0.,1.)*b       
          d=a-(0.,1.)*b       

          write(6,*)"AAA",a
          write(6,*)"BBB",b
          write(6,*)"CCC",c,abs(c)
          write(6,*)"DDD",d,abs(d)

          stop
      end subroutine test

