c * 
c * fsdeitm implements the algorithm of the paper 
c * "Fast Solution for the Diagonal Elements of the Inverse of a
c * Diagonal Matrix"
c * 
      subroutine fsdeitm(right, diag, left, lambda, m)
      integer m, i
      complex*16 left(1:m), diag(1:m), right(1:m), lambda(1:m)


c * inup and indn are source vectors that when multiplied with the tridiagonal
c * matrix give a vector that contains all zeros exept for the last/first position.
c * The last/first position will probably contain a non-zero value.

      complex*16 inup(1:m), indn(1:m)
c *   intensity comming from left or above / right or above
      complex*16 factor, tmp, nz, tmp1, tmp2, tmp3, tmp4

      inup(1) = 1
      do i = 1, m
        tmp = diag(i) * inup(i)
        if (i.GT.1) then
          tmp = tmp + right(i) * inup(i-1)
        end if
        if (i.LT.m) then
          inup(i+1) = -tmp / left(i)
          tmp3 =      -tmp / left(i)
        endif


        if (i.NE.m .AND. i.NE.1) then
          tmp1= diag(i) * inup(i) + right(i) * inup(i-1)
          tmp2= tmp1 + inup(i+1) * left(i)
        endif
      end do

      indn(m) = 1
      do i = m, 1, -1
        tmp = diag(i) * indn(i)
        if (i.LT.m) then
          tmp = tmp + left(i) * indn(i+1)
        end if
        if (i.GT.1) then
          indn(i-1) = -tmp / right(i)
        endif

        if (i.NE.m .AND. i.NE.1) then
          tmp1= diag(i) * indn(i) + right(i) * indn(i-1)
          tmp2= tmp1 + indn(i+1) * left(i)
        endif
      end do

c * inup and indn now contain solutions that (except for their respective end points,
c * consist entirely of zeros, except at their respective ends.
c * Now, splice the start of inup with the end of indn, each with some factor, so that
c * applying the spliced vector to the tridiagonal matrix results in a vector that has
c * only one non-zero element. Use that to calculate the value of its corresponding lambda
c * value. The non-zero value can be on the left or right of the splice. Calculate both,
c * and use the better one.
  
      do i = 1,m

        nz = diag(i) * inup(i)
        if (i.NE.1) then
          nz = nz + right(i) * inup(i-1)
        endif

c * find out by how much the from the right has to be multiplied to get a zerro at i+1
        factor = inup(i) / indn(i)

        if (i.NE.m) then
          nz = nz + factor * left(i) * indn(i+1)
        endif

        lambda(i) = inup(i)/nz

c       write(6,*) "lambda", i, lambda(i)
      enddo
      end
