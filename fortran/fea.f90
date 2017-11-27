module fea

   implicit none

   private
   public :: displ, initial, buildload, buildstiff, enforce, recover, eigen, topopt, buildstiff_topopt

contains

subroutine initial

! This subroutine is mainly used to allocate vectors and matrices
   use fedata
   use link1
   use plane42rect
   
   integer :: e,h,nen,bw_temp
   integer :: mdim
   integer, dimension(:), allocatable :: edof
   real(8), dimension(:), allocatable :: xe
   
   if (element(1)%id == 3) then		! changed 
   		mdim = 12
   		allocate (edof(mdim))
   		allocate (xe(8))
        neqn = 3*nn
     else 
   		mdim = 8
   		allocate (edof(mdim))
   		allocate (xe(mdim))
        neqn = 2*nn
   end if
      
   print*,'Number of elements = ', ne
   print*,' '
   print*,'Number of equations = ', neqn
        
   if (.not. banded) then
      allocate (k(neqn, neqn))
   else	! now we are using banded..
      bw = 0 !sets this to zero, will be updated
	  do e = 1, ne ! over alle elemenets
 	     ! Find coordinates etc...
 	     nen = element(e)%numnode
   			   do h = 1, nen
    		    xe(2*h-1) = x(element(e)%ix(h),1)		! x coordinate for h'th node
     		    xe(2*h  ) = x(element(e)%ix(h),2)		! y coordinate for h'th node
        		edof(3*h-2) = 3 * element(e)%ix(h) - 2  ! changed
	    	    edof(3*h-1) = 3 * element(e)%ix(h) - 1  ! 
		        edof(3*h)   = 3 * element(e)%ix(h)		 !>
      		   end do
                 !bw_temp = 2*(maxval(element(e)%ix(1:nen))-minval(element(e)%ix(1:nen))+1) ! how frank did it, no differnece, but smarter
                 ! bw_temp = maxval(edof(1:2*nen)) - minval(edof(1:2*nen)) + 1 ! saves bw_temp
                 bw_temp = maxval(edof(1:3*nen)) - minval(edof(1:3*nen)) + 1 ! saves bw_temp
        		if (bw_temp > bw) then
                	bw = bw_temp			! changes if it was the biggest.. 
                end if 
      end do
	  print*,'The bandwith is, bw = ',bw
      print*,' '
	  allocate (k(bw, neqn))		   
   end if
 
   loadcases=maxval(loads(:,5))
   print*,'loadcases',loadcases
   allocate (p(neqn), d(neqn))
   allocate (strain(ne, 3), stress(ne, 3))
 ! MODIFIED FOR WINDOWS COMPILER
    do e=1,ne
     strain(e,1:3) = 0.0d0
     stress(e,1:3) = 0.0d0
   enddo
    

end subroutine initial

subroutine topopt

   use numeth
   use processor
   use fedata

   integer :: e,i, writeDisplacements, writeDesV, writeComp, ee, eee, q, h, topIt, nTopIt, nen, useFilter
   real(8) :: s, c, vm2, gComp, plotvalMax, aa, bb, thk, rmin, sume, fac
   real(8), dimension(:), allocatable :: plotval, p1, p2, angle
   real(8), dimension(neqn) :: pComp
   real(8), dimension(ne) :: desV, desVold, vol, dg, df, dfnew
   real(8), dimension(ne, 2) :: position
   real(8) :: dg, Vmax, desVmin, penal, lambda, f, loadsr
   real(8), dimension(12) :: xe
   real(8), dimension(:), allocatable :: compVec
   real(8), dimension(ne,loadcases) ::dfd

   ! setting penalty factor
   penal = 4.
   
   
   !loop around df calc for loadcases
 
  
   ! filter options
   useFilter = 0
   rmin = 3.

   ! setting minimum relative desinty 
   desVmin = 0.000001

   ! calculating the elements volume (could be moved to initial...)
   do q = 1, ne
     nen = element(q)%numnode
   		do h = 1, nen
    	xe(2*h-1) = x(element(q)%ix(h),1)		! x coordinate for h'th node
     	xe(2*h) = x(element(q)%ix(h),2)		! y coordinate for h'th node     	  
        end do
	    aa = (xe(3)-xe(1))/2.0d0
	    bb = (xe(8)-xe(2))/2.0d0
        thk  = mprop(element(q)%mat)%thk
        vol(q) = 2.*aa*2.*bb*thk
    end do

   ! setting the maximum volume 
   Vmax = sum(vol)*0.2
   
   ! Calculating dg
   dg = vol

   ! desV initial to satisfy Vmax
   desV = Vmax/sum(vol);

   ! Topology optimization starts here 
   nTopIt = 300
   allocate(compVec(nTopIt))
   compVec = 0.
do topIt = 1, nTopIt
  loadsr=0.
do loadc=1,loadcases
   print*,'loadc', loadc

   ! setting the old design variable
   desVold = desV    

   ! Build load-vector
   call buildload

   ! Build stiffness matrix
   call buildstiff_topopt(desV, penal)
   
!   print*,'k before enforce = '
!   do i = 1,neqn
!     print "(100(f7.2,tr1))",k(i,1:neqn)
!     print*,''
!   end do
!   pause

   ! Remove rigid body modes
   call enforce

!   print*,'k after enforce = '
!   do i = 1,neqn
!     print "(100(f7.2,tr1))",k(i,1:neqn)
!     print*,''
!   end do
!   pause


   ! saves the load vector
   ! pComp = p

   if (.not. banded) then
      if (sum(k) > nb) then
        call factor(k)
      ! solve for displacement vector
        call solve(k, p)
      else
        ! print*,'WARNING in DISPL: Stiffness equal to zero'
         p = 0.0d0
      end if
   else    
		call bfactor(k)
        call bsolve(k, p)
   endif

   ! Transfer results
   d(1:neqn) = p(1:neqn)

    !Print the displacements
!	print*,'Displacements, d ='
!    do i = 249*3, 249*3+2
!	print *,d(i)
!    end do

   ! dfdx calc called here
   call dfcalc(desV, penal, f, df)
   compVec(topIt) = f
   dfd(:,loadc)=df
   loadsr=loadsr+1.
   end do
   
   do e= 1,ne
     df(e)=sum(dfd(e,:))/loadsr
     !print*,'df',df(e)
     
     end do
   print*,'f = ', f   
!   print*, 'df before filter = ', d

   if (useFilter == 1) then
     ! initialize new df
     dfnew = 0.
     ! Finding element X/Y positions 
     do e = 1, ne
       do h = 1, nen
		xe(2*h-1) = x(element(e)%ix(h),1)
 		xe(2*h) = x(element(e)%ix(h),2)       
     end do
     ! Position(x,y)	note: yes it should be a '+'
     position(e,1) = (xe(3) + xe(1))/2.0d0 ! x-pos
     position(e,2) = (xe(8) + xe(2))/2.0d0 ! y-pos
     !print*,'position(x,y) =', position(e,1), position(e,2)
     end do

     ! The filter
     do e = 1, ne
       sume =0.                    
       do i = 1, ne 
         fac = rmin - sqrt((position(i,1)-position(e,1))**2+(position(i,2)-position(e,2))**2)       
         sume = sume + max(0,fac)
         dfnew(e) = dfnew(e) + max(0,fac)*desV(i)*df(i)     
       end do    
       if (sume > 0) then !if the radius is too small then no elements are within... 
         dfnew(e) = dfnew(e)/(desV(e)*sume) ! .. this can become div by zero..
       else  ! ... then we will set the new df = to the original df
         print*,'Error: Radius of filter "rmin" is too small'
         dfnew(e) = df(e) 
       end if
       !print*, ' dfnew = ', dfnew(e)  
     end do
     df = dfnew
   end if 
   
   !print*, 'df after filter = ', df 

   ! Bi-section goes here
   call bisect(Vmax, df, dg, desVmin, desVold, desV, lambda)  

   ! Break condition goes here
        !print*, 'norm (desVold - desV) = ', sqrt(dot_product((desVold-desV),(desVold-desV)))
        !print*, ' norm(desV)*0.0001 = ', 0.001 * sqrt(dot_product(desV,desV))
        
        ! norm (desVold - desV)							< 0.001 *   norm( desV)         
   if ( sqrt(dot_product((desVold-desV),(desVold-desV))) < 0.001 * sqrt(dot_product(desV,desV)) ) then
   print*,'Achieve convergence after itterations:', topIt
      print*,''   
  	exit                
   end if   


! ending topop loop	
 end do
    print*,'Topology ended after itterations:', topIt
      print*,''    


     ! write displacements to matlab file
	writeDisplacements = 1
	
	if (writeDisplacements == 1) then
	! write datafile
	  open (13, file = trim(filename)//'_disp.m')

	! write displacements
	  write(13,'("disp =  [")')
      do ee = 1, neqn
		  write (13,*) d(ee)
      end do
      write(13,'("];")')
	  close(13)
      print*, 'Displacement printet to "*filename*_disp.m"'
      print*,''      
	else
   end if

     ! write design variable to matlab file
	writeDesV = 1	
	if (writeDesV == 1) then
	! write datafile
	  open (13, file = trim(filename)//'_desV.m')

	! write design variable
	  write(13,'("desV =  [")')
      do ee = 1, ne
		  write (13,*) desV(ee)
      end do
      write(13,'("];")')
	  close(13)
      print*, 'Design variable printet to "*filename*_desV.m"'
      print*,''      
	else
   end if   

    ! write design variable to matlab file
	writeComp = 1
	
	if (writeComp == 1) then
	! write datafile
	  open (13, file = trim(filename)//'_compliance.m')

	! write design variable
	  write(13,'("compVec =  [")')
      do  eee = 1, topIt
		  write (13,*) compVec(eee)
      end do
      write(13,'("];")')
	  close(13)
      print*, 'Comlicance variable printet to "*filename*_desV.m"'
      print*,''
	else
   end if  
 

	! Recover stress
   call recover

   !print*,'stress = 8 og 9', stress(8,1:3), stress(9,1:3)

   ! Output results
   call output

   ! Plot deformed shape
   ! call plot('deformed', 'xwin', 'color', 'Deformed')
   ! See the subroutine 'plot' in the 'processor'-module for further options

   ! Finding the maxvalue of the stress
   plotvalMax = 0.0d0 

   ! Plot element values
   allocate(plotval(ne))
   do e = 1, ne
      if (element(e)%id == 1) then
         plotval(e) = stress(e,1)
      elseif (element(e)%id == 2) then
      	!plotval(e) = stress(e,3)
	   	plotval(e) = sqrt(stress(e,1)**2 + stress(e,2)**2 - stress(e,1)*stress(e,2) + 3*(stress(e,3)**2)) ! von mises..
        	if (plotval(e) > plotvalMax) then 
            	plotvalMax = plotval(e) 
            end if
      elseif (element(e)%id == 3) then 
      	  plotval(e) =  sqrt(stress(e,1)**2 + stress(e,2)**2 - stress(e,1)*stress(e,2) + 3*(stress(e,3)**2)) ! von mises..        
      end if
      
   end do
   
   !call plot('elements', 'xwin', 'color', 'Von Mises', plotval)
   call plot('elements', 'matlab', 'color', 'Von Mises', plotval)
   ! See the subroutine 'plot' in the 'processor'-module for further options

	! Plot principle stresses and directions       
   if (element(1)%id == 2) then 
	   allocate(p1(ne))
	   allocate(p2(ne))
	   allocate(angle(ne))   
	   do e = 1, ne 	
		   	p1(e) = 0.5 * (stress(e,1)+stress(e,2)) + sqrt(( (stress(e,1) - stress(e,2)) / 2 )**2 + stress(e,3)**2)
		   	p2(e) = 0.5 * (stress(e,1)+stress(e,2)) - sqrt(( (stress(e,1) - stress(e,2)) / 2 )**2 + stress(e,3)**2)   
    	    c = ( stress(e,1)-stress(e,2) )/(p1(e) - p2(e))
			s = - 2*stress(e,3) / (p1(e) - p2(e))                 
      	    angle(e) = atan2(s,c)/2                                              
	   end do

    !   call plot('vector', 'xwin', 'gray','Principle Stress',p1,p2,angle)
     !  call plot('vector', 'matlab', 'gray','Principle Stress',p1,p2,angle)
	end if 

end subroutine topopt


subroutine displ
   ! This subroutine calculates displacements

   use numeth
   use processor
   use fedata

   integer :: e,i, writeDisplacements, ee
   real(8) :: s, c, vm2, gComp, plotvalMax
   real(8), dimension(:), allocatable :: plotval, p1, p2, angle
   real(8), dimension(neqn) :: pComp

   ! Build load-vector
   call buildload

   ! Build stiffness matrix
   call buildstiff
   
!   print*,'k before enforce = '
!   do i = 1,neqn
!     print "(100(f7.2,tr1))",k(i,1:neqn)
!     print*,''
!   end do
!   pause

   ! Remove rigid body modes
   call enforce

!   print*,'k after enforce = '
!   do i = 1,neqn
!     print "(100(f7.2,tr1))",k(i,1:neqn)
!     print*,''
!   end do
!   pause

   if (.not. banded) then
      ! NOTE: the if statement below is inserted
      ! to make sure that the solver does not try
      ! to solve the zero system: 0 x = 0 - since this is singular !
      ! REMOVE THE
      ! IF STATEMENT BELOW WHEN 
      ! CONTINUUM ELEMENTS ARE IMPLEMENTED
      if (sum(k) > nb) then
        !print*,'sum(k)',sum(k)
        !print*,'nb',nb
        call factor(k)
      ! solve for displacement vector
        call solve(k, p)
      else
        !print*,'sum(k) = ',sum(k)
        !print*,'nb =', nb
        ! print*,'WARNING in DISPL: Stiffness equal to zero'
         p = 0.0d0
      end if
   else    
		call bfactor(k)
        call bsolve(k, p)
   endif

   ! Transfer results
   d(1:neqn) = p(1:neqn)

    !Print the displacements
!	print*,'Displacements, d ='
!    do i = 249*3, 249*3+2
!	print *,d(i)
!    end do

    
    ! write displacements to matlab file
	writeDisplacements = 1
	
	if (writeDisplacements == 1) then
	! write datafile
	  open (13, file = trim(filename)//'_disp.m')

	! write displacements
	  write(13,'("disp =  [")')
      do ee = 1, neqn
		  write (13,*) d(ee)
      end do
      write(13,'("];")')
	  close(13)
	else
    print*, 'Displacement printet to "*filename*_disp.m"'
   end if          

    
	! Recover stress
   call recover

   print*,'stress = 8 og 9', stress(8,1:3), stress(9,1:3)

   ! Output results
   call output

   ! Plot deformed shape
   ! call plot('deformed', 'xwin', 'color', 'Deformed')
   ! See the subroutine 'plot' in the 'processor'-module for further options

   ! Finding the maxvalue of the stress
   plotvalMax = 0.0d0 

   ! Plot element values
   allocate(plotval(ne))
   do e = 1, ne
      if (element(e)%id == 1) then
         plotval(e) = stress(e,1)
      elseif (element(e)%id == 2) then
      	!plotval(e) = stress(e,3)
	   	plotval(e) = sqrt(stress(e,1)**2 + stress(e,2)**2 - stress(e,1)*stress(e,2) + 3*(stress(e,3)**2)) ! von mises..
        	if (plotval(e) > plotvalMax) then 
            	plotvalMax = plotval(e) 
            end if
      elseif (element(e)%id == 3) then 
      	  plotval(e) =  sqrt(stress(e,1)**2 + stress(e,2)**2 - stress(e,1)*stress(e,2) + 3*(stress(e,3)**2)) ! von mises..        
      end if
      
   end do
   
   call plot('elements', 'xwin', 'color', 'Von Mises', plotval)
   call plot('elements', 'matlab', 'color', 'Von Mises', plotval)
   ! See the subroutine 'plot' in the 'processor'-module for further options

	! Plot principle stresses and directions       
   if (element(1)%id == 2) then 
	   allocate(p1(ne))
	   allocate(p2(ne))
	   allocate(angle(ne))   
	   do e = 1, ne 	
		   	p1(e) = 0.5 * (stress(e,1)+stress(e,2)) + sqrt(( (stress(e,1) - stress(e,2)) / 2 )**2 + stress(e,3)**2)
		   	p2(e) = 0.5 * (stress(e,1)+stress(e,2)) - sqrt(( (stress(e,1) - stress(e,2)) / 2 )**2 + stress(e,3)**2)   
    	    c = ( stress(e,1)-stress(e,2) )/(p1(e) - p2(e))
			s = - 2*stress(e,3) / (p1(e) - p2(e))                 
      	    angle(e) = atan2(s,c)/2                                              
	   end do

    !   call plot('vector', 'xwin', 'gray','Principle Stress',p1,p2,angle)
     !  call plot('vector', 'matlab', 'gray','Principle Stress',p1,p2,angle)
	end if 


end subroutine displ

subroutine dfcalc(desV, penal, f, df)

  use fedata
   use link1
   use plane42rect

   integer :: e, i, j
   integer :: nen, irow, icol
   real(8) :: young, area , nu, thk , dens, penal

   integer :: mdim
   integer, dimension(:), allocatable :: edof
   real(8), dimension(:), allocatable :: xe
   real(8), dimension(:,:), allocatable :: ke, me
   real(8), dimension(ne), intent(out) :: df
   real(8), dimension(ne), intent(in) :: desV 
   real(8), intent(in) :: penal
   real(8), intent(out) :: f
   real(8), dimension(12) :: temp
     
   if (element(1)%id == 3) then
   		mdim = 12
   		allocate (edof(mdim))
   		allocate (xe(8))
        allocate (ke(mdim,mdim))
        allocate (me(mdim,mdim))
     else 
   		mdim = 8
   		allocate (edof(mdim))
   		allocate (xe(mdim))        
        allocate (ke(mdim,mdim))
        allocate (me(mdim,mdim))
   end if  

   ! sets compliance to zero
   f = 0.

   ! Calc df
   do e = 1, ne
      nen = element(e)%numnode
      do i = 1, nen
         xe(2*i-1) = x(element(e)%ix(i),1)
         xe(2*i  ) = x(element(e)%ix(i),2)
         edof(3*i-2) = 3 * element(e)%ix(i) - 2  ! changed
         edof(3*i-1) = 3 * element(e)%ix(i) - 1  ! 
         edof(3*i)   = 3 * element(e)%ix(i)		 !>
      end do
         young = mprop(element(e)%mat)%young
         nu = mprop(element(e)%mat)%nu
         thk  = mprop(element(e)%mat)%thk
         dens = mprop(element(e)%mat)%dens            
         call shell41_ke(xe, young, dens, nu, thk, ke, me)  
         temp = MATMUL(d(edof),ke)         
     	 df(e) = -penal*desV(e)**(penal-1)*DOT_PRODUCT(temp, d(edof));
         f = f + desV(e)**penal*DOT_PRODUCT(temp, d(edof))
   end do    
   
end subroutine dfcalc

subroutine bisect(Vmax, df, dg, desVmin, desVold, desV,lambda)

   use fedata    

	real(8), intent(out) :: lambda
    real(8), dimension(:), intent(out) :: desV
    real(8), intent(in) :: Vmax, desVmin
    real(8), dimension(:), intent(in) :: df, dg, desVold
    real(8) :: lambda1, lambda2, lambda, eta, g, Be
    integer :: e

    lambda1 = 0.00000001 
    lambda2 = 100000000. 
    eta = 0.2
    ! resets design variable
    desV = desVold*0
    do while ( (lambda2-lambda1)/(lambda1+lambda2) > 0.000001 )    
        lambda = (lambda1+lambda2)/2;       
    	do e = 1, ne
	        Be = -df(e)/(lambda*dg(e));
            if (desVold(e) * Be**eta <= desVmin) then
              desV(e) = desVmin
            else                
	            if (desVold(e) * Be**eta >= 1) then           
    	          desv(e) = 1
        	    else 
            	  desV(e) = desVold(e)*Be**eta
	            end if
             end if               
        end do
        g = DOT_PRODUCT(desV, dg) - Vmax;
        if (g >0) then 
          lambda1 = lambda
        else
          lambda2 = lambda            
        end if
    end while 



end subroutine bisect

subroutine buildload

   ! This subroutine builds the global load vector.
 
   use fedata
   use plane42rect
 
   integer :: i,j,h
   integer :: nen 
   integer :: dof
   integer :: eface
   integer :: e
   real(8) :: thk 
   real(8) :: fe
   integer :: mdim
   integer, dimension(:), allocatable :: edof
   real(8), dimension(:), allocatable :: xe, re
   
   if (element(1)%id == 3) then		! changed 
   		mdim = 12
   		allocate (edof(mdim))
   		allocate (xe(8))
   		allocate (re(mdim))      
     else 
   		mdim = 8
   		allocate (edof(mdim))
   		allocate (xe(mdim))
   		allocate (re(mdim))      
   end if   
	
if (element(1)%id == 3) then
    ! Build load vector
   	p(1:neqn) = 0.0d0 
   	do i = 1, np	
      if(loads(i,5)==loadc) then														
      ! Build nodal load contribution
      ! INPUT FILE NOTATION:	F=1, node #, dof label, value
      if (loads(i, 1) == 1) then					
    	dof = 3*loads(i,2) - (3 - loads(i,3))   			! changed : FZ, mX, mY 
		p(dof) = loads(i,4)									
      ! Build uniformly distributed surface (pressure) load contribution
      ! INPUT FILE NOTATION: 	SFE, element #, face #, value
      else if (loads(i, 1) == 2) then
        e = loads(i,2)
        print*,'e = ',e 
               nen = element(e)%numnode					! nen number of nodes in element
   			   do h = 1, nen
    		    xe(2*h-1) = x(element(e)%ix(h),1)		! x coordinate for h'th node
     		    xe(2*h  ) = x(element(e)%ix(h),2)		! y coordinate for h'th node
     		    edof(3*h-2) = 3 * element(e)%ix(h) - 2                  
     		    edof(3*h-1) = 3 * element(e)%ix(h) - 1  
      		    edof(3*h)   = 3 * element(e)%ix(h)		
      		   end do
        print*,'edof = ',edof                 
        fe = loads(i,4)									! the value of the load        
        print*,'fe = ',fe        
		call shell41_re(xe, fe, re)
        p(edof) = p(edof) + re

      else
         print *, 'Error in fea/buildload'
         print *, 'Load type not known'
         stop
      end if
      end if
  end do
  

else ! it is not a shell element
  
 ! Build load vector
   p(1:neqn) = 0.0d0 		
   do i = 1, np	
     if(loads(i,5)==loadc) then														
      ! Build nodal load contribution
      ! INPUT FILE NOTATION:	F=1, node #, dof label, value
      if (loads(i, 1) == 1) then					
    	dof = 2*loads(i,2) - (2 - loads(i,3) )				! finds the dof from node
		p(dof) = loads(i,4)									! puts the force value in the load vector	
      ! Build uniformly distributed surface (pressure) load contribution
      ! INPUT FILE NOTATION: 	SFE, element #, face #, value
      else if (loads(i, 1) == 2) then
        e = loads(i,2)
               nen = element(e)%numnode					! nen number of nodes in element
   			   do h = 1, nen
    		    xe(2*h-1) = x(element(e)%ix(h),1)		! x coordinate for h'th node
     		    xe(2*h  ) = x(element(e)%ix(h),2)		! y coordinate for h'th node
     		    edof(2*h-1) = 2 * element(e)%ix(h) - 1  ! x dof for h'th node
      		    edof(2*h)   = 2 * element(e)%ix(h)		! y dof for h'th node
      		   end do
        eface = loads(i,3)								! face of which the surface load is on
        fe = loads(i,4)									! the value of the load
        thk  = mprop(element(e)%mat)%thk 				! thickness
		call plane42rect_re(xe, eface, fe, thk, re)
        p(edof) = p(edof) + re
      else
         print *, 'Error in fea/buildload'
         print *, 'Load type not known'
         stop
      end if
      end if
  end do  
  print*,'p',p(dof)
end if

	
end subroutine buildload

subroutine buildstiff_topopt(desV, penal)

   ! This subroutine builds the global stiffness matrix from
   ! the local element stiffness matrices.

   use fedata
   use link1
   use plane42rect

   integer :: e, i, j
   integer :: nen, irow, icol
   real(8) :: young, area , nu, thk , dens

   integer :: mdim
   integer, dimension(:), allocatable :: edof
   real(8), dimension(:), allocatable :: xe
   real(8), dimension(:,:), allocatable :: ke, me
   real(8), dimension(ne), intent(in) :: desV
   real(8), intent(in) :: penal
   
   if (element(1)%id == 3) then
   		mdim = 12
   		allocate (edof(mdim))
   		allocate (xe(8))
        allocate (ke(mdim,mdim))
        allocate (me(mdim,mdim))
     else 
   		mdim = 8
   		allocate (edof(mdim))
   		allocate (xe(mdim))        
        allocate (ke(mdim,mdim))
        allocate (me(mdim,mdim))
   end if   

   ! Reset stiffness matrix
   if (.not. banded) then
      k(1:neqn, 1:neqn) = 0.0d0
   else
		k(1:bw,1:neqn) = 0.0d0
   end if
   
   do e = 1, ne
      ! Find coordinates and degrees of freedom
      nen = element(e)%numnode
      do i = 1, nen
         xe(2*i-1) = x(element(e)%ix(i),1)
         xe(2*i  ) = x(element(e)%ix(i),2)
         edof(3*i-2) = 3 * element(e)%ix(i) - 2  ! changed
         edof(3*i-1) = 3 * element(e)%ix(i) - 1  ! 
         edof(3*i)   = 3 * element(e)%ix(i)		 !>
      end do 
               
      ! Gather material properties and find element stiffness matrix
      select case( element(e)%id )
      case( 1 )
         young = mprop(element(e)%mat)%young
         area  = mprop(element(e)%mat)%area
         call link1_ke(xe, young, area, ke)
      case( 2 )
         ! calculate plane42rect element stiffness matrix here
           	young = mprop(element(e)%mat)%young
            nu = mprop(element(e)%mat)%nu
            thk  = mprop(element(e)%mat)%thk
            dens = mprop(element(e)%mat)%dens
			call plane42rect_ke(xe, young, nu, thk, dens, ke, me) 
       case( 3 )
            young = mprop(element(e)%mat)%young
            nu = mprop(element(e)%mat)%nu
            thk  = mprop(element(e)%mat)%thk
            dens = mprop(element(e)%mat)%dens            
            call shell41_ke(xe, young, dens, nu, thk, ke, me)        
      end select

      ! Assemble into global matrix
      if (.not. banded) then
         do i = 1, (mdim/4)*nen ! (mdim/4) = 2 if plane, 3 if shell
            do j = 1, (mdim/4)*nen
               k(edof(i), edof(j)) = k(edof(i), edof(j)) + ke(i, j)*desV(e)**penal
            end do
         end do
      else
			do i = 1, 3*nen
    	    	do j = 1, 3*nen
					icol = edof(j)-edof(i)+1
					if (icol > 0) then
	        	    	k(edof(j)-edof(i)+1, edof(i)) = k(edof(j)-edof(i)+1, edof(i)) + ke(i,j)*desV(e)**penal
                    end if 
            	end do
			end do
        end if 

   end do


end subroutine buildstiff_topopt

subroutine buildstiff

   ! This subroutine builds the global stiffness matrix from
   ! the local element stiffness matrices.

   use fedata
   use link1
   use plane42rect

   integer :: e, i, j
   integer :: nen, irow, icol
   real(8) :: young, area , nu, thk , dens

   integer :: mdim
   integer, dimension(:), allocatable :: edof
   real(8), dimension(:), allocatable :: xe
   real(8), dimension(:,:), allocatable :: ke, me
   
   if (element(1)%id == 3) then
     print*,'shell41 is the type'
     print*,''
   		mdim = 12
   		allocate (edof(mdim))
   		allocate (xe(8))
        allocate (ke(mdim,mdim))
        allocate (me(mdim,mdim))
     else 
   		mdim = 8
   		allocate (edof(mdim))
   		allocate (xe(mdim))        
        allocate (ke(mdim,mdim))
        allocate (me(mdim,mdim))
   end if   

   ! Reset stiffness matrix
   if (.not. banded) then
      k(1:neqn, 1:neqn) = 0.0d0
   else
		k(1:bw,1:neqn) = 0.0d0
   end if
   
   do e = 1, ne
      ! Find coordinates and degrees of freedom
      nen = element(e)%numnode
      do i = 1, nen
         xe(2*i-1) = x(element(e)%ix(i),1)
         xe(2*i  ) = x(element(e)%ix(i),2)
         edof(3*i-2) = 3 * element(e)%ix(i) - 2  ! changed
         edof(3*i-1) = 3 * element(e)%ix(i) - 1  ! 
         edof(3*i)   = 3 * element(e)%ix(i)		 !>
      end do 
               
      ! Gather material properties and find element stiffness matrix
      select case( element(e)%id )
      case( 1 )
         young = mprop(element(e)%mat)%young
         area  = mprop(element(e)%mat)%area
         call link1_ke(xe, young, area, ke)
      case( 2 )
         ! calculate plane42rect element stiffness matrix here
           	young = mprop(element(e)%mat)%young
            nu = mprop(element(e)%mat)%nu
            thk  = mprop(element(e)%mat)%thk
            dens = mprop(element(e)%mat)%dens
			call plane42rect_ke(xe, young, nu, thk, dens, ke, me) 
       case( 3 )
            young = mprop(element(e)%mat)%young
            nu = mprop(element(e)%mat)%nu
            thk  = mprop(element(e)%mat)%thk
            dens = mprop(element(e)%mat)%dens            
            call shell41_ke(xe, young, dens, nu, thk, ke, me)
      end select

      ! Assemble into global matrix
      if (.not. banded) then
         do i = 1, (mdim/4)*nen ! (mdim/4) = 2 if plane, 3 if shell
            do j = 1, (mdim/4)*nen
               k(edof(i), edof(j)) = k(edof(i), edof(j)) + ke(i, j)
            end do
         end do
      else
			do i = 1, 3*nen
    	    	do j = 1, 3*nen
					icol = edof(j)-edof(i)+1
					if (icol > 0) then
	        	    	k(edof(j)-edof(i)+1, edof(i)) = k(edof(j)-edof(i)+1, edof(i)) + ke(i,j)
                    end if 
            	end do
			end do
        end if 

   end do


end subroutine buildstiff

subroutine enforce

   ! This subroutine enforces the support boundary conditions.

   use fedata

   integer :: i, idof, n, kkk 
   real(8) :: penal

   ! bound: node id, uw,ux,uy (1,2,3), displacement

   if (.not. banded) then
      if (.not. penalty) then
         do i = 1, nb
            idof = int(3*(bound(i,1)-1) + bound(i,2)) ! changed
            p(1:neqn) = p(1:neqn) - k(1:neqn, idof) * bound(i, 3)
            p(idof) = bound(i, 3)
            k(1:neqn, idof) = 0.0d0
            k(idof, 1:neqn) = 0.0d0
            k(idof, idof) = 1.
         end do
      else
         penal = 10.e9*maxval(K)			! not accounted for
         do i = 1, nb
            idof = int(3*(bound(i,1)-1) + bound(i,2)) ! changed
            k(idof, idof) = k(idof, idof) + penal
            p(idof) = penal * bound(i, 3)  
         end do  
      endif
   else     	 			     		    ! not accounted for
		 do i = 1, nb
            idof = int(3*(bound(i,1)-1) + bound(i,2)) ! changed
            p(1:neqn) = p(1:neqn) - k(1:neqn, idof) * bound(i, 3)
            p(idof) = bound(i, 3)
            k(1:bw, idof) = 0.0d0 	! sets the coloumn = 0
            	do n = 1, idof-1 	! sets the diagonal down left = 0
       	            if (idof - n + 1 < bw + 1 ) then ! to ensure we dont go outside matrix
	              		k(idof - n + 1, n) = 0.0d0 !virker ikke for matricer hvor bw < idof-n+1
                    end if
            	end do             
            k(1, idof) = 1. 	! sets the top row = 1            
         end do
       
   endif
     
       
        

!	print*,'Global stiffnes matrix with BC, k ='
!	   do kkk = 2, 2
!	   !"(4(f4.2,tr1))" = "(4 coloumns(f=real number		4 spaces to write in including the decimal point 	2 spaces after the decimal point )
!	   print "(100(f7.2,tr1))", k(kkk,1:neqn)
!	   print *,' '
!       end do
!	   print *,' '


end subroutine enforce

subroutine recover

   ! This subroutine recovers the element stress, element strain, 
   ! and nodal reaction forces

   use fedata
   use link1
   use plane42rect

   integer :: e, i, nen,kkk
!   integer, parameter :: mdim = 8
!   integer :: edof(mdim)
!   real(8), dimension(mdim) :: xe, de
!   real(8), dimension(mdim, mdim) :: ke
   real(8) :: young, area , nu, dens, thk, z
   real(8), dimension(3) :: estrain, estress

   integer :: mdim
   integer, dimension(:), allocatable :: edof
   real(8), dimension(:), allocatable :: xe, de
   real(8), dimension(:,:), allocatable :: ke
   
   if (element(1)%id == 3) then		! changed 
   		mdim = 12
   		allocate (edof(mdim))
   		allocate (xe(8))
        allocate (de(mdim))
        allocate (ke(mdim,mdim))
     else 
   		mdim = 8
   		allocate (edof(mdim))
   		allocate (xe(mdim))
        allocate (de(mdim))        
        allocate (ke(mdim,mdim))        
   end if

   p = 0.0d0

   do e = 1, ne
      
      ! Find coordinates etc...
      nen = element(e)%numnode
      do i = 1,nen
         xe(2*i-1) = x(element(e)%ix(i), 1)
         xe(2*i)   = x(element(e)%ix(i), 2)
         edof(3*i-2) = 3 * element(e)%ix(i) - 2  ! changed       
         edof(3*i-1) = 3 * element(e)%ix(i) - 1  ! changed  
         edof(3*i)   = 3 * element(e)%ix(i)  ! changed  
         de(3*i-2) = d(edof(3*i-2))           ! changed  
         de(3*i-1) = d(edof(3*i-1))  ! changed  
         de(3*i)   = d(edof(3*i))  ! changed  
      end do
   
      ! Find stress and strain
      select case( element(e)%id )
      case( 1 )
         young = mprop(element(e)%mat)%young
         area  = mprop(element(e)%mat)%area
         call link1_ke(xe, young, area, ke)
         p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))
         call link1_ss(xe, de, young, estress, estrain)
         stress(e, 1:3) = estress
         strain(e, 1:3) = estrain
      case( 2 )
        young = mprop(element(e)%mat)%young
        nu = mprop(element(e)%mat)%nu
        call plane42rect_ss(xe, de, young, nu, estress, estrain)
        stress(e, 1:3) = estress
        strain(e, 1:3) = estrain

      case( 3 )
      	thk = mprop(element(e)%mat)%thk
        z = 1./2.*thk
        young = mprop(element(e)%mat)%young
        nu = mprop(element(e)%mat)%nu      	
      	call shell41_ss(xe, de, z, young, nu, estress, estrain)  
        stress(e, 1:3) = estress
        strain(e, 1:3) = estrain            	
      end select

 
   end do
   
!        print*,'stress = '
 !       do kkk = 1,4 
  !       print*, stress(kkk,1:3)
   !     end do
   
!        print*,'strain = '
!        do kkk = 1,4 
 !         print "(100(f7.2,tr1))", strain(kkk,1:3)
  !      end do
end subroutine recover

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eigen(neig)


   use numeth
   use processor
   use fedata

   integer :: e,i, pp, ppmax, neig, j, kk, b, idof, writeEigenvec, ii,jj, ee, eee, jjj, h, nen, writeXYval
   integer, intent(in) :: neig
   real(8) :: s, c, vm2, gComp, rp,error
   real(8), dimension(:), allocatable :: plotval, p1, p2, angle
   real(8), dimension(neqn) :: pComp
   real(8), dimension(neqn) :: x0,y0, xp, yp, temp,  lambda, omega
   real(8), dimension(neqn,neig) :: Zj, eigenvec
   real(8), dimension(8) :: xe
   
!	call stopwatch('star')

   ! Build load-vector
   call buildload
     
   ! Build stiffness matrix
   call buildstiff

   ! Remove rigid body modes
   call enforce

   ! First factor k, because it cannot be done in loop (multiple factor)
   	   if (.not. banded) then
	      if (sum(k) > nb) then
	        call factor(k)
	      else
	         print*,'WARNING in DISPL: Stiffness equal to zero'
	         p = 0.0d0
	      end if
	   else    
			call bfactor(k)
	   endif
	         
	do i  = 1, neig
    
	   ! Initial guess for eigenvector
	   x0 = 1.     
   
	   ! Compute Y0
	   call mmul(x0,y0)   

       ! Campute Zj = M D
       do j = 1, i-1
       		call mmul(eigenvec(:,j), Zj(:,j))
       end do	

   	   ! set max itterations
	   ppmax = 1000

	   ! initialize 
	   xp = x0 
	   yp = 0.0d0 

	   do pp = 1, ppmax
    	 ! solve K Xp = Yp-1 (we want xp from y0)
         ! BOUNDARY CONDITIONS.
      
         do b = 1, nb
    	     idof = int(3*(bound(b,1)-1) + bound(b,2))             
	         y0(idof) = 0.0d0
         end do

	     ! sets temp = y0, because wee need y0 later
    	 temp = y0
         
		   if (.not. banded) then
		      if (sum(k) > nb) then
           		!print*,'sum(k), nb',sum(k),nb                 
	    	    call solve(k, temp)
		      else
		         print*,'WARNING in DISPL: Stiffness equal to zero'
	    	     p = 0.0d0
		      end if
		   else
           		!print*,'sum(k), nb',sum(k),nb   
		        call bsolve(k, temp)
		   endif

        	 ! sets the xp to temp
	         xp = temp
             
             ! Boundary conditions
             do b = 1, nb
    	     	idof = int(3*(bound(b,1)-1) + bound(b,2))
	         	xp(idof) = 0.0d0
	         end do             

           ! Computes c
           ! and orthogonalize x
	       do j = 1, i-1
           		c = dot_product(xp,Zj(:,j))
                xp = xp - c*eigenvec(:,j)                         
    	   end do


    	! Compute Yp = M Xp
		call mmul(xp,yp)
		! !boundary conditions
             do b = 1, nb
    	     	idof = int(3*(bound(b,1)-1) + bound(b,2))
	         	yp(idof) = 0.0d0
	         end do                

    	!Compute rp
	    rp = sqrt(dot_product(xp,yp))

    	! Compute Yp
	   	yp = yp/rp
    
	      error = sqrt(dot_product((xp-x0),(xp-x0)))  /sqrt(dot_product(xp,xp))

    	! if statement goes here  
	      if (error < 1./100000000000000.) then
   	            print*,'Achieved convergence after itterations:',pp
    	      	exit                
	      end if      
	
		! sets the y0 and x0
		x0 = xp
	    y0 = yp

	   end do


   ! Compute the eugenvector 
   eigenvec(:,i) = x0/rp
   
   ! compute the eigenvalue
   lambda(i) = dot_product(x0,y0)/rp**2 
   omega(i) = sqrt(lambda(i))

   print*,'The eigenvalue(i):', i

   !print*,'Eigenvec =', eigenvec
   !print*,'Lambda(i) =', lambda(i)
   print*,'Omega(i) =', omega(i)

   print*,' '

     !call plot('eigenmode', 'xwin', 'color', ' Eigenmode ', (/lambda(i)/),eigenvec(1:neqn, i),(/1.d2, 1./30.d0/))
     !it works, but never exits last plot so have to shut down plato

end do 

	call plot('eigenmode', 'matlab', 'color', ' Eigenmode ', (/lambda(1)/),eigenvec(1:neqn, 1),(/10.d2, 2.d0/))     

!print*,'eigenvec',eigenvec(1:10,1:neig)

writeEigenvec = 1
if (writeEigenvec == 1) then
! write datafile
  open (13, file = trim(filename)//'_eigenval.m')

! write neqn
  write(13,'("neqn = ... ")')
  write (13,*) neqn
  write(13,'(";")')  

! write omega
  write(13,'("omega = [")')
  do ee = 1, neig
    write (13,'(f32.15)') omega(ee)
  end do
  write(13,'("];")')

  close(13)


! write datafile
  open (13, file = trim(filename)//'_eigendisp.m')
  
! write neqn
  write(13,'("neig = ... ")')
  write (13,*) neig
  write(13,'(";")') 

! write eigendisp
  write(13,'("eigendisp = [ ")')

! write eigen dispo
  do eee = 1,neig
    do jjj = 1,neqn
      write (13,*) eigenvec(jjj,eee)
    end do
  end do
  write(13,'("];")')
  close(13)

end if

end subroutine eigen


subroutine mmul(invector, outvector)

   ! This subroutine builds the global stiffness matrix from
   ! the local element stiffness matrices.

   use fedata
   use link1
   use plane42rect

   integer :: e, i, j
   integer :: nen, irow, icol
   integer, parameter :: mdim = 12 !8
   integer, dimension(mdim) :: edof
!   real(8), dimension(mdim) :: xe
   real(8), dimension(8) :: xe   
   real(8), dimension(mdim, mdim) :: ke , me
   real(8) :: young, area , nu, thk , dens
   real(8), dimension(neqn), intent(in) :: invector 
   real(8), dimension(neqn), intent(out) :: outvector
   
   ! Initialize the output.. 
   outvector = 0.0d0 
   
   do e = 1, ne
     !print*,'element =',e
      ! Find coordinates and degrees of freedom
      nen = element(e)%numnode
      do i = 1, nen
         xe(2*i-1) = x(element(e)%ix(i),1)
         xe(2*i  ) = x(element(e)%ix(i),2)
         edof(3*i-2) = 3 * element(e)%ix(i) - 2
         edof(3*i-1) = 3 * element(e)%ix(i) - 1           
         edof(3*i)   = 3 * element(e)%ix(i)
      end do
      
      ! Gather material properties and find element stiffness matrix
      select case( element(e)%id )
      case( 1 )
         young = mprop(element(e)%mat)%young
         area  = mprop(element(e)%mat)%area
         call link1_ke(xe, young, area, ke)
      case( 2 )
         ! calculate plane42rect element stiffness matrix here
           	young = mprop(element(e)%mat)%young
            nu = mprop(element(e)%mat)%nu
            thk  = mprop(element(e)%mat)%thk
            !if (thk > 0.9) then
             ! 	if (thk < 1.1) then
	          !  print*,'thk',thk
               ! end if
            !end if
            dens = mprop(element(e)%mat)%dens
			call plane42rect_ke(xe, young, nu, thk, dens, ke, me) 
       case( 3 )
            young = mprop(element(e)%mat)%young
            nu = mprop(element(e)%mat)%nu
            thk  = mprop(element(e)%mat)%thk
            dens = mprop(element(e)%mat)%dens            
            call shell41_ke(xe, young, dens, nu, thk, ke, me)                	
      end select

      ! Calculates the outvector       
      outvector(edof) = outvector(edof) + matmul(me,invector(edof))
      
   end do

end subroutine mmul


end module fea
