program cfd

  use boundary
  use jacobi
  use cfdio
  use nvtx

  implicit none

! Output frequency
  
  integer, parameter :: printfreq = 1000

! Variables associated with convergence

  double precision :: error, bnorm

! Set tolerance for convergence; zero or negative means do not check

  double precision, parameter :: tolerance = 0.0d0

! Main arrays

  double precision, allocatable ::  psi(:,:), zet(:,:)
  double precision, allocatable ::  psitmp(:,:), zettmp(:,:)

! Command-line arguments

  integer :: scalefactor,  numiter

  double precision :: re  ! re = 3.7 seems to be stability limit with Jacobi

  integer, parameter :: maxline = 32
  character(len=maxline) :: tmparg

!  Basic sizes of simulation

  integer, parameter :: bbase = 10
  integer, parameter :: hbase = 15
  integer, parameter :: wbase =  5
  integer, parameter :: mbase = 32
  integer, parameter :: nbase = 32

  logical :: irrotational = .true., checkerr = .false.

!  Some auxiliary parameters and variables

  integer :: m, n, b, h, w, i, j
  integer :: iter

  double precision :: tstart, tstop, ttot, titer, modvsq, hue

!  Are we stopping based on tolerance?

  if (tolerance .gt. 0.0) checkerr = .true.

!  Read in parameters

  if (command_argument_count() /= 2 .and. command_argument_count() /= 3) then

     write(*,*) 'Usage: cfd <scale> <numiter> [reynolds]'
     stop

  end if

  call get_command_argument(1, tmparg)
  read(tmparg,*) scalefactor
  call get_command_argument(2, tmparg)
  read(tmparg,*) numiter

  if (command_argument_count() == 3) then

     irrotational = .false.
     call get_command_argument(3, tmparg)
     read(tmparg,*) re
        
  else

     re = -1.0
     
  end if

  if (.not. checkerr) then
     write(*,fmt='('' Scale factor = '',i3,'', iterations = '', i6)') &
           scalefactor, numiter
  else
     write(*,fmt='('' Scale factor = '',i3,'', iterations = '', i6, &
          &'', tolerance = '', g11.4)') scalefactor, numiter, tolerance
  end if

  if (irrotational) then
        
     write(*,*) 'Irrotational flow'
        
  else

     write(*,fmt='('' Reynolds number = '', f6.3)') re
        
  end if

!  Calculate b, h & w and m & n
        
  b = bbase*scalefactor 
  h = hbase*scalefactor
  w = wbase*scalefactor 
  m = mbase*scalefactor
  n = nbase*scalefactor

  re = re / dble(scalefactor)

  write(*,fmt='('' Running CFD on '', i4, '' x '', i4, &
       &'' grid in serial '')') m, n

!  Allocate arrays, including halos on psi and tmp

  allocate(psi(0:m+1, 0:n+1))
  allocate(zet(0:m+1, 0:n+1))

  allocate(psitmp(0:m+1, 0:n+1))

  if (.not. irrotational) then

     allocate(zettmp(0:m+1, 0:n+1))

  end if

!  Zero the psi array
  call nvtxStartRange("Initialization")
  psi(:,:) = 0.0
  zet(:,:) = 0.0
  call nvtxEndRange

!  Set the psi boundary condtions which are constant

   call nvtxStartRange("boundaryPSI")
   call boundarypsi(psi, m, n, b, h, w)
   call nvtxEndRange

!  Compute normalisation factor for error

   bnorm = sum(psi(:,:)**2)

   if (.not. irrotational) then

!    Update the zeta boundary condtions which depend on psi

     call boundaryzet(zet, psi, m, n)

!    Update the normalisation

     bnorm = bnorm + sum(zet(:,:)**2)

  end if

   bnorm = sqrt(bnorm)

!  Begin iterative Jacobi loop

   write(*,*)
   write(*,*) 'Starting main loop ...'
   write(*,*)

   tstart = gettime()

  call nvtxStartRange("Overall Iteration")
  do iter = 1, numiter

!  Compute the new psi based on the old one

     call nvtxStartRange("Jacobi Step")
     if (irrotational) then

!  Call function with no vorticity
        call jacobistep(psitmp, psi, m, n)

     else

!  Call function containing vorticity

        call jacobistepvort(zettmp, psitmp, zet, psi, m, n, re)

     end if
     call nvtxEndRange

!  Compute current error value if required
     
     call nvtxStartRange("Calculate Error")
     if (checkerr .or. iter == numiter) then

        error = deltasq(psitmp, psi, m, n)

        if (.not. irrotational) then

           error = error + deltasq(zettmp, zet, m, n)

        end if

        error = sqrt(error)
        
        error = error / bnorm

     end if
     call nvtxEndRange

!  Quit early if we have reached required tolerance

     if (checkerr) then
        if (error .lt. tolerance) then
           write(*,*) 'CONVERGED iteration ', iter, ': terminating'
           exit
        end if
     end if

!  Copy back

     call nvtxStartRange("Switch Array")
     !$acc parallel 
     !$acc loop
     do i = 1, m
       !$acc loop
       do j = 1, n
         psi(i, j) = psitmp(i, j)
       end do
     end do
     !$acc end parallel

     if (.not. irrotational) then

        zet(1:m, 1:n) = zettmp(1:m, 1:n)

     end if
     call nvtxEndRange

     if (.not. irrotational) then

!    Update the zeta boundary condtions which depend on psi

        call boundaryzet(zet, psi, m, n)
        
     end if

!  End iterative Jacobi loop

     if (mod(iter,printfreq) == 0) then

        if (.not. checkerr) then
           write(*,*) 'completed iteration ', iter
        else
           write(*,*) 'completed iteration ', iter, ', error = ', error
        end if

     end if

  end do
  call nvtxEndRange

  if (iter .gt. numiter) iter = numiter

  tstop = gettime()

  ttot  = tstop-tstart
  titer = ttot/dble(iter)

  write(*,*) 
  write(*,*) '... finished'
  write(*,*)
  write(*,fmt='('' After    '', i6, '' iterations, error is '', g11.4)') &
        iter, error
  write(*,fmt='('' Time for '', i6, '' iterations was '',&
        &g11.4, '' seconds'')') iter, ttot
  write(*,fmt='('' Each individual iteration took '', g11.4, '' seconds'')') &
        titer
  write(*,*)
  write(*,*) 'Writing output file ...'

!  Output results

  call writedatafiles(psi, m, n, scalefactor)

!  Output gnuplot file

  call writeplotfile(m, n, scalefactor)

! Finish

  write(*,*) ' ... finished'
  write(*,*)
  write(*,*) 'CFD completed'
  write(*,*)

end program cfd

