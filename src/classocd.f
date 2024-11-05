      subroutine classocd(rj, xj, n, lambda, sxr)
      integer n, i
      double complex x, sxr, rj(1:n), xj(1:n)
      double precision lambda

      x = 0

      do 10 i=1,n
         x = x + conjg(xj(i)) * rj(i)
 10   continue
      
      x = x/n

      if ((abs(x) - lambda) .GT. 0) then
         sxr = (1 - (lambda / abs(x))) * x
      else
         sxr = 0
      endif

      return
      end



