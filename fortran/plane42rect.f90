MODULE plane42rect 

 ! This module contains subroutines specific to the PLANE42 element.

 IMPLICIT NONE

 PRIVATE
 PUBLIC :: plane42rect_ke, plane42rect_re, plane42rect_ss, shell41_ke, shell41_ss, shell41_re
 

CONTAINS
SUBROUTINE shell41_ke(xe, young, dens, nu, thk, ke, me)

  ! This subroutine constructs the stiffness matrix for
  ! a rectangular 4-noded quad element.

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: young, nu, thk, dens
  REAL(8), DIMENSION(:), INTENT(IN) :: xe
  REAL(8), DIMENSION(:,:), INTENT(OUT) :: ke, me
  REAL(8) :: aa, bb, t, EE
  INTEGER :: i

  ! Initialize 
  ke = 0.0d0
  me = 0.0d0
  t = thk
  EE = young

  ! aa and bb
  aa = (xe(3)-xe(1))/2.
  bb = (xe(6)-xe(4))/2.

! Ke
      ke(1,1) = -10*t**3*(aa**4-bb**2*(nu-7.D0/2.D0)*aa**2/5+bb**4)*EE/(120*&
     &nu**2-120)/bb**3/aa**3
      ke(1,2) = -t**3*EE/aa**2/bb*(4*nu*aa**2+aa**2+10*bb**2)/(nu**2-1)/120
      ke(1,3) = -t**3*EE/aa/bb**2*(4*nu*bb**2+10*aa**2+bb**2)/(nu**2-1)/120
      ke(1,4) = -5*t**3*(aa**4+2.D0/5.D0*bb**2*(nu-7.D0/2.D0)*aa**2-2*bb**4)&
     &*EE/(120*nu**2-120)/bb**3/aa**3
      ke(1,5) = t**3*EE*((-1+nu)*aa**2-10*bb**2)/(120*bb*nu**2-120*bb)/aa**2
      ke(1,6) = -t**3*EE/aa/bb**2*(-4*nu*bb**2+5*aa**2-bb**2)/(nu**2-1)/120
      ke(1,7) = 5*t**3*EE*(aa**4+2.D0/5.D0*bb**2*(nu-7.D0/2.D0)*aa**2+bb**4)&
     &/(120*nu**2-120)/bb**3/aa**3
      ke(1,8) = -((-1+nu)*aa**2+5*bb**2)*t**3*EE/(120*bb*nu**2-120*bb)/aa**2
      ke(1,9) = (-t**3*EE*(-1+nu)*bb**2-5*t**3*EE*aa**2)/aa/bb**2/(nu**2-1)/&
     &120
      ke(1,10) = t**3*EE*(-2*bb**2*aa**2*nu+10*aa**4+7*bb**2*aa**2-5*bb**4)/(1&
     &20*nu**2-120)/bb**3/aa**3
      ke(1,11) = t**3*EE*(4*nu*aa**2+aa**2-5*bb**2)/(120*bb*nu**2-120*bb)/aa**&
     &2
      ke(1,12) = -t**3*EE/aa/bb**2*(-nu*bb**2+10*aa**2+bb**2)/(nu**2-1)/120
      ke(2,1) = -t**3*EE/aa**2/bb*(4*nu*aa**2+aa**2+10*bb**2)/(nu**2-1)/120
      ke(2,2) = t**3*((-1+nu)*aa**2-5*bb**2)*EE/(45*bb*nu**2-45*bb)/aa
      ke(2,3) = -t**3*EE*nu/(12*nu**2-12)
      ke(2,4) = -t**3*EE*((-1+nu)*aa**2-10*bb**2)/(120*bb*nu**2-120*bb)/aa**2
      ke(2,5) = -t**3*((-1+nu)*aa**2+10*bb**2)*EE/(180*bb*nu**2-180*bb)/aa
      ke(2,6) = 0
      ke(2,7) = ((-1+nu)*aa**2+5*bb**2)*t**3*EE/(120*bb*nu**2-120*bb)/aa**2
      ke(2,8) = t**3*((-1+nu)*aa**2-5*bb**2)*EE/(180*bb*nu**2-180*bb)/aa
      ke(2,9) = 0
      ke(2,10) = t**3*EE*(4*nu*aa**2+aa**2-5*bb**2)/(120*bb*nu**2-120*bb)/aa**&
     &2
      ke(2,11) = -((-1+nu)*aa**2+5.D0/2.D0*bb**2)*t**3*EE/aa/bb/(nu**2-1)/45
      ke(2,12) = 0
      ke(3,1) = -t**3*EE/aa/bb**2*(4*nu*bb**2+10*aa**2+bb**2)/(nu**2-1)/120
      ke(3,2) = -t**3*EE*nu/(12*nu**2-12)
      ke(3,3) = -1/bb*t**3*EE*(-nu*bb**2+5*aa**2+bb**2)/(nu**2-1)/aa/45
      ke(3,4) = -t**3*EE/aa/bb**2*(-4*nu*bb**2+5*aa**2-bb**2)/(nu**2-1)/120
      ke(3,5) = 0
      ke(3,6) = (-2*t**3*EE*(-1+nu)*bb**2-5*t**3*EE*aa**2)/aa/bb/(nu**2-1)/9&
     &0
      ke(3,7) = t**3*EE*(nu*bb**2+5*aa**2-bb**2)/(120*aa*nu**2-120*aa)/bb**2
      ke(3,8) = 0
      ke(3,9) = -1/bb*t**3*EE*(-nu*bb**2+5*aa**2+bb**2)/(nu**2-1)/aa/180
      ke(3,10) = t**3*EE*(-nu*bb**2+10*aa**2+bb**2)/(120*aa*nu**2-120*aa)/bb**&
     &2
      ke(3,11) = 0
      ke(3,12) = (-t**3*EE*(-1+nu)*bb**2-10*t**3*EE*aa**2)/aa/bb/(nu**2-1)/1&
     &80
      ke(4,1) = -5*t**3*(aa**4+2.D0/5.D0*bb**2*(nu-7.D0/2.D0)*aa**2-2*bb**4)&
     &*EE/(120*nu**2-120)/bb**3/aa**3
      ke(4,2) = -t**3*EE*((-1+nu)*aa**2-10*bb**2)/(120*bb*nu**2-120*bb)/aa**2
      ke(4,3) = -t**3*EE/aa/bb**2*(-4*nu*bb**2+5*aa**2-bb**2)/(nu**2-1)/120
      ke(4,4) = -10*t**3*(aa**4-bb**2*(nu-7.D0/2.D0)*aa**2/5+bb**4)*EE/(120*&
     &nu**2-120)/bb**3/aa**3
      ke(4,5) = t**3*EE*(4*nu*aa**2+aa**2+10*bb**2)/(120*bb*nu**2-120*bb)/aa**&
     &2
      ke(4,6) = -t**3*EE/aa/bb**2*(4*nu*bb**2+10*aa**2+bb**2)/(nu**2-1)/120
      ke(4,7) = t**3*EE*(-2*bb**2*aa**2*nu+10*aa**4+7*bb**2*aa**2-5*bb**4)/(12&
     &0*nu**2-120)/bb**3/aa**3
      ke(4,8) = -t**3*EE/aa**2/bb*(4*nu*aa**2+aa**2-5*bb**2)/(nu**2-1)/120
      ke(4,9) = -t**3*EE/aa/bb**2*(-nu*bb**2+10*aa**2+bb**2)/(nu**2-1)/120
      ke(4,10) = 5*t**3*EE*(aa**4+2.D0/5.D0*bb**2*(nu-7.D0/2.D0)*aa**2+bb**4&
     &)/(120*nu**2-120)/bb**3/aa**3
      ke(4,11) = ((-1+nu)*aa**2+5*bb**2)*t**3*EE/(120*bb*nu**2-120*bb)/aa**2
      ke(4,12) = (-t**3*EE*(-1+nu)*bb**2-5*t**3*EE*aa**2)/aa/bb**2/(nu**2-1)&
     &/120
      ke(5,1) = t**3*EE*((-1+nu)*aa**2-10*bb**2)/(120*bb*nu**2-120*bb)/aa**2
      ke(5,2) = -t**3*((-1+nu)*aa**2+10*bb**2)*EE/(180*bb*nu**2-180*bb)/aa
      ke(5,3) = 0
      ke(5,4) = t**3*EE*(4*nu*aa**2+aa**2+10*bb**2)/(120*bb*nu**2-120*bb)/aa**&
     &2
      ke(5,5) = t**3*((-1+nu)*aa**2-5*bb**2)*EE/(45*bb*nu**2-45*bb)/aa
      ke(5,6) = t**3*EE*nu/(12*nu**2-12)
      ke(5,7) = -t**3*EE/aa**2/bb*(4*nu*aa**2+aa**2-5*bb**2)/(nu**2-1)/120
      ke(5,8) = -((-1+nu)*aa**2+5.D0/2.D0*bb**2)*t**3*EE/aa/bb/(nu**2-1)/45
      ke(5,9) = 0
      ke(5,10) = -((-1+nu)*aa**2+5*bb**2)*t**3*EE/(120*bb*nu**2-120*bb)/aa**2
      ke(5,11) = t**3*((-1+nu)*aa**2-5*bb**2)*EE/(180*bb*nu**2-180*bb)/aa
      ke(5,12) = 0
      ke(6,1) = -t**3*EE/aa/bb**2*(-4*nu*bb**2+5*aa**2-bb**2)/(nu**2-1)/120
      ke(6,2) = 0
      ke(6,3) = (-2*t**3*EE*(-1+nu)*bb**2-5*t**3*EE*aa**2)/aa/bb/(nu**2-1)/9&
     &0
      ke(6,4) = -t**3*EE/aa/bb**2*(4*nu*bb**2+10*aa**2+bb**2)/(nu**2-1)/120
      ke(6,5) = t**3*EE*nu/(12*nu**2-12)
      ke(6,6) = -1/bb*t**3*EE*(-nu*bb**2+5*aa**2+bb**2)/(nu**2-1)/aa/45
      ke(6,7) = t**3*EE*(-nu*bb**2+10*aa**2+bb**2)/(120*aa*nu**2-120*aa)/bb**2
      ke(6,8) = 0
      ke(6,9) = (-t**3*EE*(-1+nu)*bb**2-10*t**3*EE*aa**2)/aa/bb/(nu**2-1)/18&
     &0
      ke(6,10) = t**3*EE*(nu*bb**2+5*aa**2-bb**2)/(120*aa*nu**2-120*aa)/bb**2
      ke(6,11) = 0
      ke(6,12) = -1/bb*t**3*EE*(-nu*bb**2+5*aa**2+bb**2)/(nu**2-1)/aa/180
      ke(7,1) = 5*t**3*EE*(aa**4+2.D0/5.D0*bb**2*(nu-7.D0/2.D0)*aa**2+bb**4)&
     &/(120*nu**2-120)/bb**3/aa**3
      ke(7,2) = ((-1+nu)*aa**2+5*bb**2)*t**3*EE/(120*bb*nu**2-120*bb)/aa**2
      ke(7,3) = t**3*EE*(nu*bb**2+5*aa**2-bb**2)/(120*aa*nu**2-120*aa)/bb**2
      ke(7,4) = t**3*EE*(-2*bb**2*aa**2*nu+10*aa**4+7*bb**2*aa**2-5*bb**4)/(12&
     &0*nu**2-120)/bb**3/aa**3
      ke(7,5) = -t**3*EE/aa**2/bb*(4*nu*aa**2+aa**2-5*bb**2)/(nu**2-1)/120
      ke(7,6) = t**3*EE*(-nu*bb**2+10*aa**2+bb**2)/(120*aa*nu**2-120*aa)/bb**2
      ke(7,7) = -10*t**3*(aa**4-bb**2*(nu-7.D0/2.D0)*aa**2/5+bb**4)*EE/(120*&
     &nu**2-120)/bb**3/aa**3
      ke(7,8) = t**3*EE*(4*nu*aa**2+aa**2+10*bb**2)/(120*bb*nu**2-120*bb)/aa**&
     &2
      ke(7,9) = t**3*EE*(4*nu*bb**2+10*aa**2+bb**2)/(120*aa*nu**2-120*aa)/bb**&
     &2
      ke(7,10) = -5*t**3*(aa**4+2.D0/5.D0*bb**2*(nu-7.D0/2.D0)*aa**2-2*bb**4&
     &)*EE/(120*nu**2-120)/bb**3/aa**3
      ke(7,11) = -t**3*EE*((-1+nu)*aa**2-10*bb**2)/(120*bb*nu**2-120*bb)/aa**&
     &2
      ke(7,12) = t**3*EE*(-4*nu*bb**2+5*aa**2-bb**2)/(120*aa*nu**2-120*aa)/bb*&
     &*2
      ke(8,1) = -((-1+nu)*aa**2+5*bb**2)*t**3*EE/(120*bb*nu**2-120*bb)/aa**2
      ke(8,2) = t**3*((-1+nu)*aa**2-5*bb**2)*EE/(180*bb*nu**2-180*bb)/aa
      ke(8,3) = 0
      ke(8,4) = -t**3*EE/aa**2/bb*(4*nu*aa**2+aa**2-5*bb**2)/(nu**2-1)/120
      ke(8,5) = -((-1+nu)*aa**2+5.D0/2.D0*bb**2)*t**3*EE/aa/bb/(nu**2-1)/45
      ke(8,6) = 0
      ke(8,7) = t**3*EE*(4*nu*aa**2+aa**2+10*bb**2)/(120*bb*nu**2-120*bb)/aa**&
     &2
      ke(8,8) = t**3*((-1+nu)*aa**2-5*bb**2)*EE/(45*bb*nu**2-45*bb)/aa
      ke(8,9) = -t**3*EE*nu/(12*nu**2-12)
      ke(8,10) = t**3*EE*((-1+nu)*aa**2-10*bb**2)/(120*bb*nu**2-120*bb)/aa**2
      ke(8,11) = -t**3*((-1+nu)*aa**2+10*bb**2)*EE/(180*bb*nu**2-180*bb)/aa
      ke(8,12) = 0
      ke(9,1) = (-t**3*EE*(-1+nu)*bb**2-5*t**3*EE*aa**2)/aa/bb**2/(nu**2-1)/&
     &120
      ke(9,2) = 0
      ke(9,3) = -1/bb*t**3*EE*(-nu*bb**2+5*aa**2+bb**2)/(nu**2-1)/aa/180
      ke(9,4) = -t**3*EE/aa/bb**2*(-nu*bb**2+10*aa**2+bb**2)/(nu**2-1)/120
      ke(9,5) = 0
      ke(9,6) = (-t**3*EE*(-1+nu)*bb**2-10*t**3*EE*aa**2)/aa/bb/(nu**2-1)/18&
     &0
      ke(9,7) = t**3*EE*(4*nu*bb**2+10*aa**2+bb**2)/(120*aa*nu**2-120*aa)/bb**&
     &2
      ke(9,8) = -t**3*EE*nu/(12*nu**2-12)
      ke(9,9) = -1/bb*t**3*EE*(-nu*bb**2+5*aa**2+bb**2)/(nu**2-1)/aa/45
      ke(9,10) = t**3*EE*(-4*nu*bb**2+5*aa**2-bb**2)/(120*aa*nu**2-120*aa)/bb*&
     &*2
      ke(9,11) = 0
      ke(9,12) = (-2*t**3*EE*(-1+nu)*bb**2-5*t**3*EE*aa**2)/aa/bb/(nu**2-1)/&
     &90
      ke(10,1) = t**3*EE*(-2*bb**2*aa**2*nu+10*aa**4+7*bb**2*aa**2-5*bb**4)/(1&
     &20*nu**2-120)/bb**3/aa**3
      ke(10,2) = t**3*EE*(4*nu*aa**2+aa**2-5*bb**2)/(120*bb*nu**2-120*bb)/aa**&
     &2
      ke(10,3) = t**3*EE*(-nu*bb**2+10*aa**2+bb**2)/(120*aa*nu**2-120*aa)/bb**&
     &2
      ke(10,4) = 5*t**3*EE*(aa**4+2.D0/5.D0*bb**2*(nu-7.D0/2.D0)*aa**2+bb**4&
     &)/(120*nu**2-120)/bb**3/aa**3
      ke(10,5) = -((-1+nu)*aa**2+5*bb**2)*t**3*EE/(120*bb*nu**2-120*bb)/aa**2
      ke(10,6) = t**3*EE*(nu*bb**2+5*aa**2-bb**2)/(120*aa*nu**2-120*aa)/bb**2
      ke(10,7) = -5*t**3*(aa**4+2.D0/5.D0*bb**2*(nu-7.D0/2.D0)*aa**2-2*bb**4&
     &)*EE/(120*nu**2-120)/bb**3/aa**3
      ke(10,8) = t**3*EE*((-1+nu)*aa**2-10*bb**2)/(120*bb*nu**2-120*bb)/aa**2
      ke(10,9) = t**3*EE*(-4*nu*bb**2+5*aa**2-bb**2)/(120*aa*nu**2-120*aa)/bb*&
     &*2
      ke(10,10) = -10*t**3*(aa**4-bb**2*(nu-7.D0/2.D0)*aa**2/5+bb**4)*EE/(12&
     &0*nu**2-120)/bb**3/aa**3
      ke(10,11) = -t**3*EE/aa**2/bb*(4*nu*aa**2+aa**2+10*bb**2)/(nu**2-1)/120
      ke(10,12) = t**3*EE*(4*nu*bb**2+10*aa**2+bb**2)/(120*aa*nu**2-120*aa)/bb&
     &**2
      ke(11,1) = t**3*EE*(4*nu*aa**2+aa**2-5*bb**2)/(120*bb*nu**2-120*bb)/aa**&
     &2
      ke(11,2) = -((-1+nu)*aa**2+5.D0/2.D0*bb**2)*t**3*EE/aa/bb/(nu**2-1)/45
      ke(11,3) = 0
      ke(11,4) = ((-1+nu)*aa**2+5*bb**2)*t**3*EE/(120*bb*nu**2-120*bb)/aa**2
      ke(11,5) = t**3*((-1+nu)*aa**2-5*bb**2)*EE/(180*bb*nu**2-180*bb)/aa
      ke(11,6) = 0
      ke(11,7) = -t**3*EE*((-1+nu)*aa**2-10*bb**2)/(120*bb*nu**2-120*bb)/aa**&
     &2
      ke(11,8) = -t**3*((-1+nu)*aa**2+10*bb**2)*EE/(180*bb*nu**2-180*bb)/aa
      ke(11,9) = 0
      ke(11,10) = -t**3*EE/aa**2/bb*(4*nu*aa**2+aa**2+10*bb**2)/(nu**2-1)/120
      ke(11,11) = t**3*((-1+nu)*aa**2-5*bb**2)*EE/(45*bb*nu**2-45*bb)/aa
      ke(11,12) = t**3*EE*nu/(12*nu**2-12)
      ke(12,1) = -t**3*EE/aa/bb**2*(-nu*bb**2+10*aa**2+bb**2)/(nu**2-1)/120
      ke(12,2) = 0
      ke(12,3) = (-t**3*EE*(-1+nu)*bb**2-10*t**3*EE*aa**2)/aa/bb/(nu**2-1)/1&
     &80
      ke(12,4) = (-t**3*EE*(-1+nu)*bb**2-5*t**3*EE*aa**2)/aa/bb**2/(nu**2-1)&
     &/120
      ke(12,5) = 0
      ke(12,6) = -1/bb*t**3*EE*(-nu*bb**2+5*aa**2+bb**2)/(nu**2-1)/aa/180
      ke(12,7) = t**3*EE*(-4*nu*bb**2+5*aa**2-bb**2)/(120*aa*nu**2-120*aa)/bb*&
     &*2
      ke(12,8) = 0
      ke(12,9) = (-2*t**3*EE*(-1+nu)*bb**2-5*t**3*EE*aa**2)/aa/bb/(nu**2-1)/&
     &90
      ke(12,10) = t**3*EE*(4*nu*bb**2+10*aa**2+bb**2)/(120*aa*nu**2-120*aa)/bb&
     &**2
      ke(12,11) = t**3*EE*nu/(12*nu**2-12)
      ke(12,12) = -1/bb*t**3*EE*(-nu*bb**2+5*aa**2+bb**2)/(nu**2-1)/aa/45

! 	me matrix 
      me(1,1) = 1727.D0/3150.D0*dens*bb*aa*thk
      me(1,2) = 461.D0/3150.D0*dens*bb*aa**2*thk
      me(1,3) = 461.D0/3150.D0*dens*bb**2*aa*thk
      me(1,4) = 613.D0/3150.D0*dens*bb*aa*thk
      me(1,5) = -137.D0/1575.D0*dens*bb*aa**2*thk
      me(1,6) = 199.D0/3150.D0*dens*bb**2*aa*thk
      me(1,7) = 197.D0/3150.D0*dens*bb*aa*thk
      me(1,8) = -58.D0/1575.D0*dens*bb*aa**2*thk
      me(1,9) = -58.D0/1575.D0*dens*bb**2*aa*thk
      me(1,10) = 613.D0/3150.D0*dens*bb*aa*thk
      me(1,11) = 199.D0/3150.D0*dens*bb*aa**2*thk
      me(1,12) = -137.D0/1575.D0*dens*bb**2*aa*thk
      me(2,1) = 461.D0/3150.D0*dens*bb*aa**2*thk
      me(2,2) = 16.D0/315.D0*dens*bb*aa**3*thk
      me(2,3) = dens*bb**2*aa**2*thk/25
      me(2,4) = 137.D0/1575.D0*dens*bb*aa**2*thk
      me(2,5) = -4.D0/105.D0*dens*bb*aa**3*thk
      me(2,6) = 2.D0/75.D0*dens*bb**2*aa**2*thk
      me(2,7) = 58.D0/1575.D0*dens*bb*aa**2*thk
      me(2,8) = -2.D0/105.D0*dens*bb*aa**3*thk
      me(2,9) = -4.D0/225.D0*dens*bb**2*aa**2*thk
      me(2,10) = 199.D0/3150.D0*dens*bb*aa**2*thk
      me(2,11) = 8.D0/315.D0*dens*bb*aa**3*thk
      me(2,12) = -2.D0/75.D0*dens*bb**2*aa**2*thk
      me(3,1) = 461.D0/3150.D0*dens*bb**2*aa*thk
      me(3,2) = dens*bb**2*aa**2*thk/25
      me(3,3) = 16.D0/315.D0*dens*bb**3*aa*thk
      me(3,4) = 199.D0/3150.D0*dens*bb**2*aa*thk
      me(3,5) = -2.D0/75.D0*dens*bb**2*aa**2*thk
      me(3,6) = 8.D0/315.D0*dens*bb**3*aa*thk
      me(3,7) = 58.D0/1575.D0*dens*bb**2*aa*thk
      me(3,8) = -4.D0/225.D0*dens*bb**2*aa**2*thk
      me(3,9) = -2.D0/105.D0*dens*bb**3*aa*thk
      me(3,10) = 137.D0/1575.D0*dens*bb**2*aa*thk
      me(3,11) = 2.D0/75.D0*dens*bb**2*aa**2*thk
      me(3,12) = -4.D0/105.D0*dens*bb**3*aa*thk
      me(4,1) = 613.D0/3150.D0*dens*bb*aa*thk
      me(4,2) = 137.D0/1575.D0*dens*bb*aa**2*thk
      me(4,3) = 199.D0/3150.D0*dens*bb**2*aa*thk
      me(4,4) = 1727.D0/3150.D0*dens*bb*aa*thk
      me(4,5) = -461.D0/3150.D0*dens*bb*aa**2*thk
      me(4,6) = 461.D0/3150.D0*dens*bb**2*aa*thk
      me(4,7) = 613.D0/3150.D0*dens*bb*aa*thk
      me(4,8) = -199.D0/3150.D0*dens*bb*aa**2*thk
      me(4,9) = -137.D0/1575.D0*dens*bb**2*aa*thk
      me(4,10) = 197.D0/3150.D0*dens*bb*aa*thk
      me(4,11) = 58.D0/1575.D0*dens*bb*aa**2*thk
      me(4,12) = -58.D0/1575.D0*dens*bb**2*aa*thk
      me(5,1) = -137.D0/1575.D0*dens*bb*aa**2*thk
      me(5,2) = -4.D0/105.D0*dens*bb*aa**3*thk
      me(5,3) = -2.D0/75.D0*dens*bb**2*aa**2*thk
      me(5,4) = -461.D0/3150.D0*dens*bb*aa**2*thk
      me(5,5) = 16.D0/315.D0*dens*bb*aa**3*thk
      me(5,6) = -dens*bb**2*aa**2*thk/25
      me(5,7) = -199.D0/3150.D0*dens*bb*aa**2*thk
      me(5,8) = 8.D0/315.D0*dens*bb*aa**3*thk
      me(5,9) = 2.D0/75.D0*dens*bb**2*aa**2*thk
      me(5,10) = -58.D0/1575.D0*dens*bb*aa**2*thk
      me(5,11) = -2.D0/105.D0*dens*bb*aa**3*thk
      me(5,12) = 4.D0/225.D0*dens*bb**2*aa**2*thk
      me(6,1) = 199.D0/3150.D0*dens*bb**2*aa*thk
      me(6,2) = 2.D0/75.D0*dens*bb**2*aa**2*thk
      me(6,3) = 8.D0/315.D0*dens*bb**3*aa*thk
      me(6,4) = 461.D0/3150.D0*dens*bb**2*aa*thk
      me(6,5) = -dens*bb**2*aa**2*thk/25
      me(6,6) = 16.D0/315.D0*dens*bb**3*aa*thk
      me(6,7) = 137.D0/1575.D0*dens*bb**2*aa*thk
      me(6,8) = -2.D0/75.D0*dens*bb**2*aa**2*thk
      me(6,9) = -4.D0/105.D0*dens*bb**3*aa*thk
      me(6,10) = 58.D0/1575.D0*dens*bb**2*aa*thk
      me(6,11) = 4.D0/225.D0*dens*bb**2*aa**2*thk
      me(6,12) = -2.D0/105.D0*dens*bb**3*aa*thk
      me(7,1) = 197.D0/3150.D0*dens*bb*aa*thk
      me(7,2) = 58.D0/1575.D0*dens*bb*aa**2*thk
      me(7,3) = 58.D0/1575.D0*dens*bb**2*aa*thk
      me(7,4) = 613.D0/3150.D0*dens*bb*aa*thk
      me(7,5) = -199.D0/3150.D0*dens*bb*aa**2*thk
      me(7,6) = 137.D0/1575.D0*dens*bb**2*aa*thk
      me(7,7) = 1727.D0/3150.D0*dens*bb*aa*thk
      me(7,8) = -461.D0/3150.D0*dens*bb*aa**2*thk
      me(7,9) = -461.D0/3150.D0*dens*bb**2*aa*thk
      me(7,10) = 613.D0/3150.D0*dens*bb*aa*thk
      me(7,11) = 137.D0/1575.D0*dens*bb*aa**2*thk
      me(7,12) = -199.D0/3150.D0*dens*bb**2*aa*thk
      me(8,1) = -58.D0/1575.D0*dens*bb*aa**2*thk
      me(8,2) = -2.D0/105.D0*dens*bb*aa**3*thk
      me(8,3) = -4.D0/225.D0*dens*bb**2*aa**2*thk
      me(8,4) = -199.D0/3150.D0*dens*bb*aa**2*thk
      me(8,5) = 8.D0/315.D0*dens*bb*aa**3*thk
      me(8,6) = -2.D0/75.D0*dens*bb**2*aa**2*thk
      me(8,7) = -461.D0/3150.D0*dens*bb*aa**2*thk
      me(8,8) = 16.D0/315.D0*dens*bb*aa**3*thk
      me(8,9) = dens*bb**2*aa**2*thk/25
      me(8,10) = -137.D0/1575.D0*dens*bb*aa**2*thk
      me(8,11) = -4.D0/105.D0*dens*bb*aa**3*thk
      me(8,12) = 2.D0/75.D0*dens*bb**2*aa**2*thk
      me(9,1) = -58.D0/1575.D0*dens*bb**2*aa*thk
      me(9,2) = -4.D0/225.D0*dens*bb**2*aa**2*thk
      me(9,3) = -2.D0/105.D0*dens*bb**3*aa*thk
      me(9,4) = -137.D0/1575.D0*dens*bb**2*aa*thk
      me(9,5) = 2.D0/75.D0*dens*bb**2*aa**2*thk
      me(9,6) = -4.D0/105.D0*dens*bb**3*aa*thk
      me(9,7) = -461.D0/3150.D0*dens*bb**2*aa*thk
      me(9,8) = dens*bb**2*aa**2*thk/25
      me(9,9) = 16.D0/315.D0*dens*bb**3*aa*thk
      me(9,10) = -199.D0/3150.D0*dens*bb**2*aa*thk
      me(9,11) = -2.D0/75.D0*dens*bb**2*aa**2*thk
      me(9,12) = 8.D0/315.D0*dens*bb**3*aa*thk
      me(10,1) = 613.D0/3150.D0*dens*bb*aa*thk
      me(10,2) = 199.D0/3150.D0*dens*bb*aa**2*thk
      me(10,3) = 137.D0/1575.D0*dens*bb**2*aa*thk
      me(10,4) = 197.D0/3150.D0*dens*bb*aa*thk
      me(10,5) = -58.D0/1575.D0*dens*bb*aa**2*thk
      me(10,6) = 58.D0/1575.D0*dens*bb**2*aa*thk
      me(10,7) = 613.D0/3150.D0*dens*bb*aa*thk
      me(10,8) = -137.D0/1575.D0*dens*bb*aa**2*thk
      me(10,9) = -199.D0/3150.D0*dens*bb**2*aa*thk
      me(10,10) = 1727.D0/3150.D0*dens*bb*aa*thk
      me(10,11) = 461.D0/3150.D0*dens*bb*aa**2*thk
      me(10,12) = -461.D0/3150.D0*dens*bb**2*aa*thk
      me(11,1) = 199.D0/3150.D0*dens*bb*aa**2*thk
      me(11,2) = 8.D0/315.D0*dens*bb*aa**3*thk
      me(11,3) = 2.D0/75.D0*dens*bb**2*aa**2*thk
      me(11,4) = 58.D0/1575.D0*dens*bb*aa**2*thk
      me(11,5) = -2.D0/105.D0*dens*bb*aa**3*thk
      me(11,6) = 4.D0/225.D0*dens*bb**2*aa**2*thk
      me(11,7) = 137.D0/1575.D0*dens*bb*aa**2*thk
      me(11,8) = -4.D0/105.D0*dens*bb*aa**3*thk
      me(11,9) = -2.D0/75.D0*dens*bb**2*aa**2*thk
      me(11,10) = 461.D0/3150.D0*dens*bb*aa**2*thk
      me(11,11) = 16.D0/315.D0*dens*bb*aa**3*thk
      me(11,12) = -dens*bb**2*aa**2*thk/25
      me(12,1) = -137.D0/1575.D0*dens*bb**2*aa*thk
      me(12,2) = -2.D0/75.D0*dens*bb**2*aa**2*thk
      me(12,3) = -4.D0/105.D0*dens*bb**3*aa*thk
      me(12,4) = -58.D0/1575.D0*dens*bb**2*aa*thk
      me(12,5) = 4.D0/225.D0*dens*bb**2*aa**2*thk
      me(12,6) = -2.D0/105.D0*dens*bb**3*aa*thk
      me(12,7) = -199.D0/3150.D0*dens*bb**2*aa*thk
      me(12,8) = 2.D0/75.D0*dens*bb**2*aa**2*thk
      me(12,9) = 8.D0/315.D0*dens*bb**3*aa*thk
      me(12,10) = -461.D0/3150.D0*dens*bb**2*aa*thk
      me(12,11) = -dens*bb**2*aa**2*thk/25
      me(12,12) = 16.D0/315.D0*dens*bb**3*aa*thk
     

 END SUBROUTINE shell41_ke
 
 SUBROUTINE shell41_ss(xe, de, z, young, nu, estress, estrain)
   ! This subrotuine constructs the element stress and strain
 
  REAL(8), INTENT(IN) :: young, nu, z
  REAL(8), DIMENSION(:), INTENT(IN)  :: xe, de
  REAL(8), DIMENSION(:), INTENT(OUT) :: estress, estrain
  REAL(8) :: aa, bb, y, x, fact
  REAL(8) :: Bmat(3, 12), Cmat(3, 3)
  INTEGER :: i,j

  Bmat = 0.0d0

  ! Build constitutive matrix (plane stress)
  Cmat = 0.0d0
  fact = young/(1.-nu**2.)
  Cmat(1, 1) = fact
  Cmat(1, 2) = fact*nu
  Cmat(2, 1) = fact*nu
  Cmat(2, 2) = fact
  Cmat(3, 3) = fact*(1.-nu)/2.

  ! aa and bb
  aa = (xe(3)-xe(1))/2.
  bb = (xe(6)-xe(4))/2.


  ! x and y
  x = 0.0d0
  y = 0.0d0
!	x = -aa
 !   y = bb

   ! Bmat build
 	  Bmat(1,1) = 3.D0/4.D0*x/aa**3-3.D0/4.D0*x*y/aa**3/bb
      Bmat(1,2) = -1/aa/4+3.D0/4.D0*x/aa**2+y/bb/aa/4-3.D0/4.D0*x*y/a&
      &a**2/bb
      Bmat(1,3) = 0.0d0
      Bmat(1,4) = -3.D0/4.D0*x/aa**3+3.D0/4.D0*x*y/aa**3/bb
      Bmat(1,5) = 1/aa/4+3.D0/4.D0*x/aa**2-y/bb/aa/4-3.D0/4.D0*x*y/aa&
     &**2/bb
      Bmat(1,6) = 0.0d0
      Bmat(1,7) = -3.D0/4.D0*x/aa**3-3.D0/4.D0*x*y/aa**3/bb
      Bmat(1,8) = 1/aa/4+3.D0/4.D0*x/aa**2+y/bb/aa/4+3.D0/4.D0*x*y/aa&
    &**2/bb
      Bmat(1,9) = 0.0d0
      Bmat(1,10) = 3.D0/4.D0*x/aa**3+3.D0/4.D0*x*y/aa**3/bb
      Bmat(1,11) = -1/aa/4+3.D0/4.D0*x/aa**2-y/bb/aa/4+3.D0/4.D0*x*y/&
     &aa**2/bb
      Bmat(1,12) = 0.0d0
      Bmat(2,1) = 3.D0/4.D0*y/bb**3-3.D0/4.D0*x*y/bb**3/aa
      Bmat(2,2) = 0
      Bmat(2,3) = -1/bb/4+x/bb/aa/4+3.D0/4.D0*y/bb**2-3.D0/4.D0*x*y/a&
     &a/bb**2
      Bmat(2,4) = 3.D0/4.D0*y/bb**3+3.D0/4.D0*x*y/bb**3/aa
      Bmat(2,5) = 0.0d0
      Bmat(2,6) = -1/bb/4-x/bb/aa/4+3.D0/4.D0*y/bb**2+3.D0/4.D0*x*y/a&
     &a/bb**2
      Bmat(2,7) = -3.D0/4.D0*y/bb**3-3.D0/4.D0*x*y/bb**3/aa
      Bmat(2,8) = 0.0d0
      Bmat(2,9) = 1/bb/4+x/bb/aa/4+3.D0/4.D0*y/bb**2+3.D0/4.D0*x*y/aa&
     &/bb**2
      Bmat(2,10) = -3.D0/4.D0*y/bb**3+3.D0/4.D0*x*y/bb**3/aa
      Bmat(2,11) = 0.0d0
      Bmat(2,12) = 1/bb/4-x/bb/aa/4+3.D0/4.D0*y/bb**2-3.D0/4.D0*x*y/a&
     &a/bb**2
      Bmat(3,1) = 1/bb/aa/2-3.D0/8.D0*x**2/aa**3/bb-3.D0/8.D0*y**2/bb&
     &**3/aa
      Bmat(3,2) = 1/bb/8+x/bb/aa/4-3.D0/8.D0*x**2/aa**2/bb
      Bmat(3,3) = 1/aa/8+y/bb/aa/4-3.D0/8.D0*y**2/aa/bb**2
      Bmat(3,4) = -1/bb/aa/2+3.D0/8.D0*x**2/aa**3/bb+3.D0/8.D0*y**2/b&
     &b**3/aa
      Bmat(3,5) = 1/bb/8-x/bb/aa/4-3.D0/8.D0*x**2/aa**2/bb
      Bmat(3,6) = -1/aa/8-y/bb/aa/4+3.D0/8.D0*y**2/aa/bb**2
      Bmat(3,7) = 1/bb/aa/2-3.D0/8.D0*x**2/aa**3/bb-3.D0/8.D0*y**2/bb&
     &**3/aa
      Bmat(3,8) = -1/bb/8+x/bb/aa/4+3.D0/8.D0*x**2/aa**2/bb
      Bmat(3,9) = -1/aa/8+y/bb/aa/4+3.D0/8.D0*y**2/aa/bb**2
      Bmat(3,10) = -1/bb/aa/2+3.D0/8.D0*x**2/aa**3/bb+3.D0/8.D0*y**2/&
     &bb**3/aa
      Bmat(3,11) = -1/bb/8-x/bb/aa/4+3.D0/8.D0*x**2/aa**2/bb
      Bmat(3,12) = 1/aa/8-y/bb/aa/4-3.D0/8.D0*y**2/aa/bb**2

 !     print*,'de = ', de 
 !     do j = 1, 3
 !     	print*,'Bmat = ', Bmat(j,1:12)
 !     end do 
  !    print*,'z = ', z

  ! Compute element strain
  estrain = -z*MATMUL(Bmat, de)

 !       print*,'strain = ', estrain

  ! Compute element stress
  estress = MATMUL(Cmat, estrain) ! is this true?

 END SUBROUTINE shell41_ss
 
 SUBROUTINE plane42rect_ke(xe, young, nu, thk, dens, ke, me)

  ! This subroutine constructs the stiffness matrix for
  ! a rectangular 4-noded quad element.

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: young, nu, thk, dens
  REAL(8), DIMENSION(:), INTENT(IN) :: xe
  REAL(8), DIMENSION(:,:), INTENT(OUT) :: ke, me

  REAL(8) :: Cmat(3,3), fact, aa, bb, t
  REAL(8), DIMENSION(2) :: xiVec, etaVec
  INTEGER :: w,i,j,ng   			!shape
  Real(8) :: xi, eta 				!shape
  REAL(8), DIMENSION(2,2):: gamma	!shape
  REAL(8):: detjac					!shape
  REAL(8), DIMENSION(8,3) :: temp	!shape
  REAL(8):: N(2,8), bmat(3,8), jac(2,2)
  

  ! Initialize 
  ke = 0.
  me = 0.

  ! Build constitutive matrix (plane stress)
  Cmat = 0.0d0
  fact = young/(1.-nu**2)
  Cmat(1, 1) = fact
  Cmat(1, 2) = fact*nu
  Cmat(2, 1) = fact*nu
  Cmat(2, 2) = fact
  Cmat(3, 3) = fact*(1.-nu)/2.
  

  !Assuming that always: gausspoints = 2
  etaVec(1) = -1/sqrt(3.)
  etaVec(2) = +1/sqrt(3.)
  xiVec(1) = -1/sqrt(3.)
  xiVec(2) = +1/sqrt(3.)
  W = 1
  t = thk
  
  ! finding the element stifness matrix and mass matrix
  do i = 1, 2!ng
    do j = 1, 2!ng
      eta = etaVec(i)
      xi = xiVec(j)
      call shape(xe, xi, eta, N, bmat, jac, detjac)
      temp = matmul(transpose(bmat),cmat)
      ke = ke + w*w*thk*matmul(temp,bmat)*detjac
      me = me + w*w*thk*dens*matmul(transpose(N),N)*detjac
    end do
  end do


 END SUBROUTINE plane42rect_ke

 SUBROUTINE shell41_re(xe, fe, re)

  ! This subroutine assembles the element surface loads.

  REAL(8), INTENT(IN) :: fe
  REAL(8), DIMENSION(:), INTENT(IN) :: xe
  REAL(8), INTENT(OUT) :: re(12)
  REAL(8) :: N(12) 

  REAL(8) :: aa, bb
 
  !      3
  !  l ______k
  !   |     |
  ! 4 |     |  2
  !   |_____|
  !  i       j
  !      1
  !
  ! node numbers: i, j, k, l
  ! face numbers: 1(j->i), 2(k->j), 3(l->k), 4(i->l)

  aa = (xe(3)-xe(1))/2.0d0
  bb = (xe(8)-xe(2))/2.0d0

  re = 0.

      N(1) = aa*bb
      N(2) = aa**2.*bb/3.
      N(3) = aa*bb**2./3.
      N(4) = aa*bb
      N(5) = -aa**2.*bb/3.
      N(6) = aa*bb**2./3.
      N(7) = aa*bb
      N(8) = -aa**2.*bb/3.
      N(9) = -aa*bb**2./3.
      N(10) = aa*bb
      N(11) = aa**2.*bb/3.
      N(12) = -aa*bb**2./3.

   re = N*fe ! Shape * pressure    

 END SUBROUTINE shell41_re

 SUBROUTINE plane42rect_re(xe, eface, fe, thk, re)

  ! This subroutine assembles the element surface loads.



  INTEGER, INTENT(IN) :: eface
  REAL(8), INTENT(IN) :: fe, thk
  REAL(8), DIMENSION(:), INTENT(IN) :: xe
  REAL(8), INTENT(OUT) :: re(8)
  REAL(8) :: aa, bb, Nface(2, 8), f(8)
  REAL(8):: hJac(2,1)
  REAL(8):: n(2,8), bmat(3,8), jac(2,2)
  REAL(8)::detjac, eta, xi, w
  
  !      3
  !  l ______k
  !   |     |
  ! 4 |     |  2
  !   |_____|
  !  i       j
  !      1
  !
  ! node numbers: i, j, k, l
  ! face numbers: 1(j->i), 2(k->j), 3(l->k), 4(i->l)


  aa = (xe(3)-xe(1))/2.0d0
  bb = (xe(8)-xe(2))/2.0d0
  
  ! We have gauspoint = 1 so: 
  w = 2.	! always
  eta = 0.0d0 	!one will be changed below
  xi = 0.0d0 	!one will be changed below
  
  Nface = 0.0d0 
  f = 0.0d0
  
  IF (eface == 1) THEN
   Nface(1, 1) = aa
   Nface(1, 3) = aa
   Nface(2, 2) = aa
   Nface(2, 4) = aa
   
    eta = -1.
    
    call shape(xe,xi,eta,n,bmat,jac,detjac)
    
    hJac(1,1) = -jac(1,2)
    hJac(2,1) = jac(1,1)
    f = f + w*thk*fe*matmul((TRANSPOSE(n)),hJac)

	
  ELSEIF (eface == 2) THEN
   Nface(1, 3) = bb
   Nface(1, 5) = bb
   Nface(2, 4) = bb
   Nface(2, 6) = bb

    xi = 1.

    call shape(xe,xi,eta,n,bmat,jac,detjac)
    
    hJac(1,1) = -jac(2,2)
    hJac(2,1) = jac(2,1)
    f = f + w*thk*fe*matmul((TRANSPOSE(n)),hJac)  

  ELSEIF (eface == 3) THEN
   Nface(1, 5) = aa
   Nface(1, 7) = aa
   Nface(2, 6) = aa
   Nface(2, 8) = aa

    eta = 1.

    call shape(xe,xi,eta,n,bmat,jac,detjac)
    
    hJac(1,1) = jac(1,2)
    hJac(2,1) = -jac(1,1)
    f = f + w*thk*fe*matmul((TRANSPOSE(n)),hJac)

	
  ELSEIF (eface == 4) THEN
   Nface(1, 1) = bb
   Nface(1, 7) = bb
   Nface(2, 2) = bb
   Nface(2, 8) = bb

    xi = -1.

    call shape(xe,xi,eta,n,bmat,jac,detjac)
    
    hJac(1,1) = jac(2,2)
    hJac(2,1) = -jac(2,1)
    f = f + w*thk*fe*matmul((TRANSPOSE(n)),hJac)

  ENDIF

  re = f

 END SUBROUTINE plane42rect_re

 SUBROUTINE plane42rect_ss(xe, de, young, nu, estress, estrain)

  ! This subrotuine constructs the element stress and strain
 
  REAL(8), INTENT(IN) :: young, nu
  REAL(8), DIMENSION(:), INTENT(IN)  :: xe, de
  REAL(8), DIMENSION(:), INTENT(OUT) :: estress, estrain
  REAL(8) :: aa, bb, yy, xx, fact, eta, xi, detjac
  INTEGER :: kkk, kkkk
  REAL(8), DIMENSION(2) :: etaVec, xiVec
  REAL(8) :: Bmat(3, 8), Cmat(3, 3)
  REAL(8):: n(2,8), jac(2,2)

  Bmat = 0.0d0

  !Assuming that always: gausspoints = 1
  eta = 0.0d0
  xi= 0.0d0

  call shape(xe, xi, eta, n, bmat, jac, detjac)
  
  ! Compute element strain
  estrain = MATMUL(Bmat, de)

  ! Build constitutive matrix (plane stress)
  Cmat = 0.0d0
  fact = young/(1.-nu**2)
  Cmat(1, 1) = fact
  Cmat(1, 2) = fact*nu
  Cmat(2, 1) = fact*nu
  Cmat(2, 2) = fact
  Cmat(3, 3) = fact*(1.-nu)/2.


  ! Compute element stress
  estress = MATMUL(Cmat, estrain)

  ! Compute principal stress and direction
 ! print*,'Plane42rect_ss function not defined yet: xe = ', xe, ' young = ', young, ' nu = ', nu

 END SUBROUTINE plane42rect_ss

	SUBROUTINE shape(xe, xi, eta, n, bmat, jac, detjac)

        REAL(8), DIMENSION(:), INTENT(IN)  :: xe
        Real(8), INTENT(IN) :: xi, eta
        
        REAL(8), DIMENSION(:,:), INTENT(OUT) :: Bmat, n, jac
		REAL(8), INTENT(OUT) :: detjac        

        REAL(8), DIMENSION(2,2) :: gamma        
        REAL(8), DIMENSION(3,4) :: L, temp
        REAL(8) :: gammaTilde(4,4) 
        REAL(8) :: nTilde(4,8)        
		REAL(8), DIMENSION(8) :: dn 

        gamma = 0.0d0
        jac = 0.0d0
        gammaTilde = 0.0d0
        n = 0.0d0
        nTilde = 0.0d0
        
		! L is setup
        L = 0.0d0
		L(1,1) = 1.
        L(2,4) = 1.
        L(3,2) = 1.
        L(3,3) = 1.

        ! Shape functions 
        n = 0.0d0
        n(1,1) = 0.25*(1.-xi)*(1.-eta) 
   	    n(2,2) = 0.25*(1.-xi)*(1.-eta) 
        
   	    n(1,3) = 0.25*(1.+xi)*(1.-eta)
        n(2,4) = 0.25*(1.+xi)*(1.-eta)
        
        n(1,5) = 0.25*(1.+xi)*(1.+eta)
        n(2,6) = 0.25*(1.+xi)*(1.+eta) 
        
        n(1,7) = 0.25*(1.-xi)*(1.+eta)
        n(2,8) = 0.25*(1.-xi)*(1.+eta) 

        ! differentiated shape fundtions 1-4
        ! N1
	    dn(1) = -0.25+ 0.25*eta 	!diff wrt xi
    	dn(2) = -0.25 + 0.25*xi		!diff wrt eta
        ! N2  
        dn(3) = +0.25 - 0.25*eta
    	dn(4) = -0.25 - 0.25*xi        
        ! N3
        dn(5) = +0.25 + 0.25*eta
    	dn(6) = +0.25 + 0.25*xi 
        ! N4
        dn(7) = -0.25 - 0.25*eta
    	dn(8) = +0.25 - 0.25*xi 
        
		! Jacobian [J] 
        jac(1,1) = dot_product(dn(1::2),xe(1::2))
        jac(1,2) = dot_product(dn(1::2),xe(2::2))
        jac(2,1) = dot_product(dn(2::2),xe(1::2))   
        jac(2,2) = dot_product(dn(2::2),xe(2::2))
                        
		! Determinant of jacobian |[J]| 
        detjac = jac(1,1)*jac(2,2)-jac(2,1)*jac(1,2) 
		
		! Inverted jacobian [gamma] 
        gamma(1,1) = jac(2,2)/detjac
        gamma(1,2) = -jac(1,2)/detjac
        gamma(2,1) = -jac(2,1)/detjac
        gamma(2,2) = jac(1,1)/detjac

        ! Gamma tilde
		gammaTilde = 0.0d0 
        gammaTilde(1,1) = gamma(1,1) 
        gammaTilde(1,2) = gamma(1,2)
        gammaTilde(2,1) = gamma(2,1)
        gammaTilde(2,2) = gamma(2,2)
        gammaTilde(3,3) = gamma(1,1) 
        gammaTilde(3,4) = gamma(1,2)
        gammaTilde(4,3) = gamma(2,1)
        gammaTilde(4,4) = gamma(2,2)  

        ! N tilde    
    	nTilde = 0.0d0
        nTilde(1,1) = dn(1)  
        nTilde(2,1) = dn(2)         
        nTilde(3,2) = dn(1)
        nTilde(4,2) = dn(2)
        
        nTilde(1,3) = dn(3)  
        nTilde(2,3) = dn(4)         
        nTilde(3,4) = dn(3)
        nTilde(4,4) = dn(4)
                        
        nTilde(1,5) = dn(5)  
        nTilde(2,5) = dn(6)         
        nTilde(3,6) = dn(5)
        nTilde(4,6) = dn(6)

        nTilde(1,7) = dn(7)  
        nTilde(2,7) = dn(8)         
        nTilde(3,8) = dn(7)
        nTilde(4,8) = dn(8)

        ! strain-displacement matrix [B]
        temp = matmul(L,gammaTilde)
        bmat = matmul(temp,nTilde) 

    END SUBROUTINE


END MODULE plane42rect
