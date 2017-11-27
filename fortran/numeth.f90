module numeth

   implicit none

   save

   private
   public :: factor, solve, bfactor, bsolve

contains

   subroutine factor(a)

   ! This subroutine factors a general matrix [a].
 
   integer :: i, j, neqn
   real(8), dimension(:, :), intent(inout) :: a

   neqn = size(a, 1)

   do j = 2, neqn
      do i = 1, j-1
         a(j, i) = a(j, i) - dot_product(a(j, 1:i-1), a(1:i-1, i))
         a(i, j) = a(i, j) - dot_product(a(i, 1:i-1), a(1:i-1, j))
         a(j, i) = a(j, i)/a(i, i)
      end do
      a(j, j) = a(j, j) - dot_product(a(j, 1:j-1), a(1:j-1, j))
   end do

end subroutine factor

subroutine solve(a, b)

   ! This subroutine solves [a]{x} = {b} using the previously factored
   ! coefficient matrix [a] from subroutine 'factor'.
   ! The subroutine returns the solution {x} by overwriting b.

   integer :: i, neqn
   real(8), dimension(:, :), intent(in) :: a
   real(8), dimension(:), intent(inout) :: b

   neqn = size(a, 1)

   ! Forward substitution
   do i = 2, neqn
      b(i) = b(i) - dot_product(a(i, 1:i-1), b(1:i-1))
   end do
 
   ! Backward substitution
   b(neqn) = b(neqn)/a(neqn, neqn)
   do i = neqn-1, 1, -1
      b(i) = (b(i) - dot_product(a(i, i+1:neqn), b(i+1:neqn)))/a(i, i)
   end do

end subroutine solve
 
subroutine bfactor(a)

   ! This subroutine factors a matrix [a] stored in banded form.

   integer :: i, j, n, l, k, bw, neqn
   real(8) :: c
   real(8), dimension(:, :), intent(inout) :: a

   bw = size(a, 1)
   neqn = size(a, 2)

   do n = 1, neqn
      do l = 2, bw
         if (a(l, n) == 0.) cycle
         i = n+l-1
         c = a(l, n)/a(1, n)
         j = 0
      do k = l, bw
         j = j+1
         a(j, i) = a(j, i) - c*a(k, n)
      end do
      a(l, n) = c
      end do
   end do 

end subroutine bfactor

subroutine bsolve(a, b)

   ! This subroutine solves [a]{x} = {b} using the previously factored
   ! coefficient matrix [a] from subroutine 'bfactor'.
   ! The subroutine returns the solution {x} by overwriting b.

   integer :: i, n, l, m, k, bw, neqn
   real(8), dimension(:, :), intent(inout) :: a
   real(8), dimension(:) :: b

   bw = size(a, 1)
   neqn = size(a, 2)

   do n = 1,neqn
      do l = 2, bw
         if (a(l, n) == 0.) cycle 
         i = n+l-1
         b(i) = b(i)-a(l, n)*b(n)
      end do
  !        print*,'b(n) = ', b(n)   
 !        print*,'a(1, n) = ', a(1, n)   
      b(n) = b(n)/a(1, n)
!                print*,'b(n) = ', b(n)
   end do

   do m = 2, neqn
      n = neqn+1-m
      do l = 2, bw
         if (a(l, n) == 0.) cycle 
         k = n+l-1
         b(n) = b(n)-a(l, n)*b(k)
      end do
   end do

end subroutine bsolve

end module numeth

