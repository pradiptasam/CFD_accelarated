module jacobi

  use cudafor
  implicit none

contains

attributes(global) subroutine jacobistep(psinew, psi, m, n)
!attributes(global) subroutine jacobistep(psinew, psi)
  integer, value :: m, n
  double precision, dimension(0:m+1, 0:n+1) :: psinew, psi
!  double precision, dimension(:,:) :: psinew, psi
  integer :: i, j 
  
  i = (blockIdx%x-1)*blockDim%x+threadIdx%x
  j = (blockIdx%y-1)*blockDim%y+threadIdx%y

  psinew(i, j) = 0.25d0*(psi(i+1, j) + psi(i-1, j) + &
                             psi(i, j+1) + psi(i, j-1)     )

end subroutine jacobistep

subroutine jacobistepvort(zetnew, psinew, zet, psi, m, n, re)

  integer :: m, n
  double precision :: re 
  double precision, dimension(0:m+1, 0:n+1) :: zetnew, zet, psinew, psi

  psinew(1:m, 1:n) = 0.25d0*(psi(2:m+1, 1:n) + psi(0:m-1, 1:n) + &
                             psi(1:m, 2:n+1) + psi(1:m, 0:n-1) - &
                             zet(1:m,   1:n))

  zetnew(1:m, 1:n) = 0.25d0*(zet(2:m+1, 1:n) + zet(0:m-1, 1:n) +     &
                             zet(1:m, 2:n+1) + zet(1:m, 0:n-1)   ) - &
                   re/16.0*((psi(1:m, 2:n+1) - psi(1:m, 0:n-1)) *    &
                            (zet(2:m+1, 1:n) - zet(0:m-1, 1:n)) -    &
                            (psi(2:m+1, 1:n) - psi(0:m-1, 1:n)) *    &
                            (zet(1:m, 2:n+1) - zet(1:m, 0:n-1))  )

end subroutine jacobistepvort

double precision function deltasq(new, old, m, n)

  integer :: m, n
  double precision, dimension(0:m+1, 0:n+1) :: new, old

  integer :: ierr

  deltasq =   sum((new(1:m,1:n)-old(1:m,1:n))**2)

end function deltasq

attributes(global) subroutine copyCuda2D(d_dst, d_src, m, n)

  integer, value :: m, n
  double precision, dimension(0:m+1, 0:n+1) :: d_dst, d_src
  integer :: i, j

  i = (blockIdx%x-1)*blockDim%x+threadIdx%x
  j = (blockIdx%y-1)*blockDim%y+threadIdx%y

  d_dst(i,j) = d_src(i,j)

end subroutine


end module jacobi
