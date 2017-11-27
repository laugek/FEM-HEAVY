module fedata

! This module is used to define global values and matrices

   ! Variables
   integer :: ne, nn, nb, np, nm, neqn, bw, loadcases, loadc

   ! Coordinates:
   real(8), dimension(:,:), allocatable :: x

   ! Elements:
   type element_def
      integer, dimension(4) :: ix
      integer :: id, mat, numnode
   end type element_def
   type(element_def), dimension(:), allocatable :: element

   ! Material properties:
   type matprop
      real(8) :: young, nu, thk, youngy, shear, dens, area
   end type matprop
   type(matprop), dimension(:), allocatable :: mprop
   
   ! Boundary conditions
   real(8), dimension(:,:), allocatable :: bound, loads
   real(8) :: accel(2)
   
   ! Working arrays:
   real(8), dimension(:,:), allocatable :: k, strain, stress
   real(8), dimension(:),   allocatable ::d, p

   !working integers


   ! i/o
   character(len = 20) :: filename
   logical, parameter :: plot2screen = .true.

   character(len = 20) :: antype

   ! Constants
   real(8), parameter ::  pi = 3.1415927 

   ! Parameters
   real(8), parameter :: scale_def = 1.0
   real(8), parameter :: scale_vec = 1.0
   logical, parameter :: banded = .true.
   logical, parameter :: penalty = .false.
   logical, parameter :: globalComliance = .false.   
   logical, parameter :: nodalStressValue = .true.   

   !
   ! Timing
   real(8) start_time, finish_time

end module fedata
