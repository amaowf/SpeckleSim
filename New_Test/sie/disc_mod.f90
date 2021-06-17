!  Rune_Code.f90 
!
!  FUNCTIONS:
!  Rune_Code - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Rune_Code
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

	module disc_mod
!
! Module defining the discretization data structure and
! subroutines for filling the structure with data using
! different discretization algorithms.

    use input_tri
    use Gauss_Quadrature_Fast   
    use libmath
    use tri_type_function
	    
    implicit none
	public  Data_Import	  
    type(structure), public :: struc
	 
	contains
	
	subroutine Data_Import()    !
      implicit none

      integer :: status_open            ! Rueckgabewert aus iostat bei open
      integer :: status_read
      integer :: j
      real(kind = dp) :: wert
      integer, dimension(:, :), allocatable :: Elements_data
      real(kind = dp), dimension(:, :), allocatable :: Node_Vector  
      
      !type(Structure), intent(out) :: struc
      print*, 'Object type: S (sphere) or R (rough_surface)'
      read(*, *) object_type

      print*, 'Illumination type: Plane (Plane, plane) or Gauss (Gauss) or focused beam (FB, fb)'
      read(*, *) illumination

      select case (object_type)
      case('S')
         D_factor = 0.5e-6 !sphere diameter
         n_disc = 258!66! (rec:6, 18, 66, 258, 1026, 4098, 16386, 65538, 262146, ...)16386!
         call recursive_triangulation(struc, n_disc, D_factor*0.5)
      case ('R')
         D_factor = 1.0
         write(*,*) ' Input the p-file name for the rough surface'
         write(*,'(A$)') ' (Not with quotation marks): '
         read(*,'(A)') Name_NodeVector
         write(*,*)

         Nc = 0
         open(unit=25, file = Name_NodeVector, status='old', action='read', iostat=status_open)
         oeffnen_1: if ( status_open == 0 ) then

         einlese_schleife_1: do
                  read (25, *, iostat=status_read) wert
                  if ( status_read /= 0 ) exit
                     Nc = Nc + 1
                  end do einlese_schleife_1!
         readif_1: if ( status_read > 0 ) then
                  write(*,*) 'Beim Lesen von Zeile ', Nc + 1, &
                  ' ist ein Fehler aufgetreten'
                  else ! status_read < 0
                     write(*,*)
                     write(*,*) ' Hinweis: das Dateiende wurde erreicht'
                     write(*,*) ' => Nc = ', Nc, 'nodes'
                  end if readif_1
               else oeffnen_1
                  write(*,*) 'Beim OEffnen der Datei trat &
                              &  Systemfehler Nr. ', status_open,' auf'
               end if oeffnen_1
         close( unit=25 )  !Datei schliessen

         Ne = 0
         write(*,*) ' Input the t-file name for the rough surface'
         write(*,'(A$)') ' (Not with quotation marks): '
         read(*,'(A)') Name_Element
         write(*,*)
         open(unit=25, file = Name_Element, status='old', action='read', iostat=status_open)
         oeffnen_2: if ( status_open == 0 ) then

               einlese_schleife_2: do
                  read (25, *, iostat=status_read) wert

                  if ( status_read /= 0 ) exit
                     Ne = Ne + 1
                  end do einlese_schleife_2
                  readif_2: if ( status_read > 0 ) then
                              write(*,*) 'Beim Lesen von Zeile ', Ne +1, &
                              ' ist ein Fehler aufgetreten'
                           else ! status_read < 0
                              write(*,*)
                              write(*,*) ' Hinweis: das Dateiende wurde erreicht'
                              write(*,*) ' => Ne = ', Ne, 'elements'
                           end if readif_2
               else oeffnen_2
                  write(*,*) 'Beim OEffnen der Datei trat &
                              &  Systemfehler Nr. ', status_open,' auf'
               end if oeffnen_2
         close( unit=25 )  !Datei schliessen

    !------------------------------------
        allocate (Node_Vector(1:3, 1:Nc))
        allocate (Elements_data(1:3, 1:Ne))  
        allocate(struc%elements(Ne))   
        allocate(struc%points(Nc))


        !pp: array of cartesian node coordinates x, y and z
        open(unit = 80, file = Name_NodeVector, status = 'old', action='read') !Small_
        read(80, *) Node_Vector
        close(80)
        
        !Array of node numberings for each element
        open(unit = 81, file = Name_Element, status = 'old', action='read') !Small_
        read(81, *) Elements_data
        close(81)        

    !-----------------------------------------------------------------------
        do j = 1, Nc !pp
            struc%points(j)%point = Node_Vector(:, j)*D_factor !for sphere, here is 2:4
        end do

        do j = 1, Ne !tt
            allocate(struc%elements(j)%corners(3)) ! 8-node element with node numbering
            struc%elements(j)%corners = Elements_data(1:3, j) !for sphere, the column range is 3:10
        end do

        call find_neighbours(struc)
        call find_midpoints(struc)
        call find_quadpoints(struc)
        
        case default
            call exit            
        end select

        Result_SEH = 'Result_SEH_tri.txt'
        return
    end subroutine Data_Import
	 
	subroutine recursive_triangulation(struct, n_disc, R_factor)
     ! Creating symmetric triangular discretization recursively
     ! n=6,18,58,258,1026, ...
         type(Structure) :: struct
         integer :: n_disc, i, ne, np_curr, ne_curr, lim!, split_num
         real(kind = dp) :: R_factor

         print *, 'Creating discretization recursively with #points = ', n_disc
         allocate(struct%points(n_disc))
         struct%points(1)%point=(/ 1.0, 0.0, 0.0 /)*R_factor
         struct%points(2)%point=(/-1.0, 0.0, 0.0 /)*R_factor
         struct%points(3)%point=(/ 0.0, 1.0, 0.0 /)*R_factor
         struct%points(4)%point=(/ 0.0,-1.0, 0.0 /)*R_factor
         struct%points(5)%point=(/ 0.0, 0.0, 1.0 /)*R_factor
         struct%points(6)%point=(/ 0.0, 0.0,-1.0 /)*R_factor

         ne=2*(n_disc-2) !int(8*4**(split_num-1))
 
         allocate(struct%elements(ne))
         do i=1,ne
         struct%elements(i)%nr_of_corners=3
         allocate(struct%elements(i)%corners(3))
         end do
         struct%elements(1)%corners=(/ 1, 5, 3/)
         struct%elements(2)%corners=(/ 3, 5, 2/)
         struct%elements(3)%corners=(/ 2, 5, 4/)
         struct%elements(4)%corners=(/ 4, 5, 1/)
         struct%elements(5)%corners=(/ 1, 3, 6/)
         struct%elements(6)%corners=(/ 3, 2, 6/)
         struct%elements(7)%corners=(/ 2, 4, 6/)
         struct%elements(8)%corners=(/ 4, 1, 6/)

         np_curr=6
         ne_curr=8

         lim=int(log10(ne/6.0)/log10(4.0))
         do i=1, lim
         call recursive_split(struct, np_curr, ne_curr, R_factor)
         ne_curr=ne_curr*4
         np_curr=(ne_curr/2)+2
         end do
         call find_neighbours(struct)
         call find_midpoints(struct)
         call find_quadpoints(struct)
    end subroutine
    
    !subroutine recursive_triangulation_surface(struct, n_disc, R_factor)
    ! ! Creating symmetric triangular discretization recursively
    ! ! n=6,18,58,258,1026, ...
    !     type(Structure) :: struct
    !     integer :: n_disc, i, ne, np_curr, ne_curr, lim, split_num
    !     real(kind = dp) :: R_factor
 !
    !     print *, 'Creating discretization recursively with #points = ', n_disc
    !     allocate(struct%points(n_disc))
    !     struct%points(1)%point=(/ 1.0, 0.0, 0.0 /)*R_factor
    !     struct%points(2)%point=(/-1.0, 0.0, 0.0 /)*R_factor
    !     struct%points(3)%point=(/ 0.0, 1.0, 0.0 /)*R_factor
    !     struct%points(4)%point=(/ 0.0,-1.0, 0.0 /)*R_factor
    !     struct%points(5)%point=(/ 0.0, 0.0, 1.0 /)*R_factor
    !     struct%points(6)%point=(/ 0.0, 0.0,-1.0 /)*R_factor
    !
    !     ne=2*(n_disc-2) !int(8*4**(split_num-1))
 !
    !     allocate(struct%elements(ne))
    !     do i=1,ne
    !     struct%elements(i)%nr_of_corners=3
    !     allocate(struct%elements(i)%corners(3))
    !     end do
    !     struct%elements(1)%corners=(/ 1, 5, 3/)
    !     struct%elements(2)%corners=(/ 3, 5, 2/)
    !     struct%elements(3)%corners=(/ 2, 5, 4/)
    !     struct%elements(4)%corners=(/ 4, 5, 1/)
    !     struct%elements(5)%corners=(/ 1, 3, 6/)
    !     struct%elements(6)%corners=(/ 3, 2, 6/)
    !     struct%elements(7)%corners=(/ 2, 4, 6/)
    !     struct%elements(8)%corners=(/ 4, 1, 6/)
 !
    !     np_curr=6
    !     ne_curr=8
 !
    !     lim=int(log10(ne/6.0)/log10(4.0))
    !     do i=1, lim
    !     call recursive_split(struct, np_curr, ne_curr, R_factor)
    !     ne_curr=ne_curr*4
    !     np_curr=(ne_curr/2)+2
    !     end do
    !     call find_neighbours(struct)
    !     call find_midpoints(struct)
    !     call find_quadpoints(struct)
 !
    !end subroutine recursive_triangulation_surface

    subroutine recursive_split(struct, n, ne, Rf)
        integer :: n, ne, i, counter, &
            ind_a, ind_b, ind_c
        type(Structure) :: struct
        type(Element), dimension(:), allocatable :: el_list
        real(kind = dp), intent(in) :: Rf !Scaling factor
        real(kind = dp) :: a(3),b(3),c(3)!, v1(3), v2(3)
        logical :: lo(3)
        allocate(el_list(4*ne))
        counter=n+1

        do i=1,ne
            lo=(/ .false., .false., .false. /)
            a=normalize(0.5*(struct%points(struct%elements(i)%corners(1))&
            %point + struct%points(struct%elements(i)%corners(3))%point ))*Rf
            b=normalize(0.5*(struct%points(struct%elements(i)%corners(1))&
            %point+struct%points(struct%elements(i)%corners(2))%point ))*Rf
            c=normalize(0.5*(struct%points(struct%elements(i)%corners(2))&
            %point+struct%points(struct%elements(i)%corners(3))%point ))*Rf
            ind_a=point_in_list(a,struct%points)
            lo(1)=ind_a.eq. -1
            if (lo(1))then ! a is not present in list
                struct%points(counter)%point=a
                ind_a=counter
                counter=counter+1
            end if

            ind_b=point_in_list(b,struct%points)
            lo(1)=ind_b.eq. -1
            if (lo(1))then ! a is not present in list
                struct%points(counter)%point=b
                ind_b=counter
                counter=counter+1
            end if

            ind_c=point_in_list(c,struct%points)
            lo(1)=ind_c.eq. -1
            if (lo(1))then ! a is not present in list
                struct%points(counter)%point=c
                ind_c=counter
                counter=counter+1
            end if
            struct%elements(ne+(i-1)*3+1)%corners=&
            (/ struct%elements(i)%corners(1), ind_b, ind_a/)
            struct%elements(ne+(i-1)*3+2)%corners=&
            (/ ind_b, struct%elements(i)%corners(2), ind_c/)
            struct%elements(ne+(i-1)*3+3)%corners=&
            (/ ind_a, ind_c, struct%elements(i)%corners(3)/)
            struct%elements(i)%corners=(/ ind_a, ind_b, ind_c /)
        end do
    end subroutine

    function point_in_list(p, list)
     ! Function returning the index of the point p in list
     ! and returns -1 if it is not in the list
        type(Point), dimension(:), allocatable, intent(in) :: list
        real(kind = dp) :: temp(3)
        real(kind = dp), intent(in) :: p(3)
        integer :: point_in_list
        integer :: i
        intrinsic sqrt
  
        point_in_list=-1
        do i=1,size(list)
            temp=list(i)%point-p
            temp(1)=sqrt(temp(1)**2+temp(2)**2+temp(3)**2)
            if (temp(1)<1e-30)then
                point_in_list=i
            return
            end if
        end do
        return
     end function

     subroutine find_quadpoints(struct)
         ! Subroutine finding 3 points in each element
         ! for gauss quad integration and saving the
         ! resulting points in struct%quadpoints
         ! Compatible with only triangular elements
         type(Structure) :: struct
         integer :: m, ne

         ne=size(struct%elements)
         allocate(struct%quadpoints(ne, 3))
    
        do m=1,ne
         struct%quadpoints(m,1)%point=(2.0/3.0) &
         *struct%points(struct%elements(m)%corners(1))%point &
         +(1.0/6.0)*struct%points(struct%elements(m)%corners(2))%point&
         +(1.0/6.0)*struct%points(struct%elements(m)%corners(3))%point
         struct%quadpoints(m,2)%point=(1.0/6.0) &
         *struct%points(struct%elements(m)%corners(1))%point &
         +(2.0/3.0)*struct%points(struct%elements(m)%corners(2))%point &
         +(1.0/6.0)*struct%points(struct%elements(m)%corners(3))%point
         struct%quadpoints(m,3)%point=(1.0/6.0) &
         *struct%points(struct%elements(m)%corners(1))%point &
         +(1.0/6.0)*struct%points(struct%elements(m)%corners(2))%point &
         +(2.0/3.0)*struct%points(struct%elements(m)%corners(3))%point
         end do
         !print*, 'struct%quadpoints(m,1)%point =', struct%quadpoints(5,1)%point
    end subroutine

    subroutine find_midpoints(struct)
    ! Subroutine finding all centroids of the elements in struct
    ! and saving the resulting points in struct%midpoints.
    ! Compatible with only triangular elements.
    integer :: ne, i
    type(Structure), intent(inout) :: struct
    
    real(kind = dp) :: v1(3),v2(3),v3(3), mid(3)

    ne=size(struct%elements)
    if (allocated(struct%midpoints))then
    deallocate(struct%midpoints)
    endif
    allocate(struct%midpoints(ne))
    do i=1,ne
    v1=struct%points(struct%elements(i)%corners(1))%point
    v2=struct%points(struct%elements(i)%corners(2))%point
    v3=struct%points(struct%elements(i)%corners(3))%point
    mid=(v1+v2+v3)/3.0
    struct%midpoints(i)%point=mid
    end do
    end subroutine

    subroutine find_neighbours(struct)
        ! Subroutine finding all neighbours in struct and
        ! saving the result in struct%neighbours.
        ! Compatible with arbitrary number of corners.
        integer :: ne, i, j, d1,d2, nn, &
        out_arr(2)
        type(Pair), dimension(:), allocatable :: temp_neighbours
        type(Structure), intent(inout) :: struct
        type(Structure) :: temp_struct

        nn=0
        ne=size(struct%elements)
        if (allocated(struct%neighbours))then
            deallocate(struct%neighbours)
        endif

        if (allocated(temp_struct%neighbours))then
            deallocate(temp_struct%neighbours)
        endif

        allocate(temp_neighbours(3*ne/2))
        allocate(temp_struct%neighbours(3*ne/2))

        do i = 1, ne
            do j = i+1,ne 
                d1 = size(struct%elements(i)%corners)
                d2 = size(struct%elements(j)%corners)
                out_arr = common_2(struct%elements(i)%corners, &
                        struct%elements(j)%corners, d1, d2 )
                if ( out_arr(2) /= -1 )then
                    nn = nn + 1
                    temp_struct%neighbours(nn)%elements=(/ i, j /)
                    temp_struct%neighbours(nn)%corners=(/out_arr(1:2)/)
                end if
            end do
        end do
        allocate(struct%neighbours(nn))
        struct%neighbours(1:nn) = temp_struct%neighbours(1:nn)
        deallocate(temp_struct%neighbours)
        return
    end subroutine

    function common_2(arr1, arr2, dim1, dim2)
    ! function checking whether arr1 and arr2 have 2 common elements.
    ! If true, the function returns the two common elements
    ! If false, the function returns (-1, -1)
    ! Compatible with arbitrary number of dimensions
    integer, intent(in) :: dim1, dim2
    integer :: i, j, n_equal, common_2(2), n_different
    integer, intent(in) :: arr1(dim1), arr2(dim2)
    common_2=(/-1,-1/)
    !print*, 'arr1=', arr1
    !print*, 'arr2=', arr2

    n_equal=0
    do i=1,dim1
        n_different=0
        do j=1,dim2
            if (arr1(i) == arr2 (j))then
                n_equal=n_equal+1
                common_2(n_equal)=arr1(i)
            endif
        end do
    end do
    return
    end function common_2

    function normalize(v)
        real(kind = dp), intent(in) :: v(3)
        real(kind = dp) :: normalize(3), l
        l = sqrt(v(1)**2+v(2)**2+v(3)**2)
        normalize(1:3)=v(1:3)/l
        return
    end function normalize


 
	
	subroutine Observation_points(p_t, n1, n2, ax1_m, ax2_m, ref_pl, r_local)
		integer, intent(in) :: p_t, n1, n2
        real(kind = dp ), intent(in) :: ax1_m, ax2_m, ref_pl		
		real(kind = dp ) :: de_1, de_2, ax1_min, ax2_min, rho, phi
		type(Point), dimension(:), allocatable :: r_local
		integer :: M, i, j
		
		!ref_pl = 0.0
		
		M = n1*n2
		allocate(r_local(M))		
		
		if (p_t .eq. 1)then ! xy-plot
            ax1_min = -ax1_m
		    ax2_min = -ax2_m
		    de_1 = (ax1_m - ax1_min)/(n1-1)
		    de_2 = (ax2_m - ax2_min)/(n2-1)		
			do i = 1, n2
				do j = 1, n1
					r_local((i-1)*n1+j)%point= &
					(/ ax1_m - de_1*(j-1), ax2_m - de_2*(i-1), ref_pl/)
				end do
			end do
        else if (p_t .eq. 2)then ! yz-plot
            ax1_min = -ax1_m
		    ax2_min = -ax2_m
		    de_1 = (ax1_m - ax1_min)/(n1-1)
		    de_2 = (ax2_m - ax2_min)/(n2-1)		            
			do i=1, n2
				do j=1, n1
					r_local((i-1)*n1 + j)%point=(/ref_pl, ax1_m - de_1*(j-1), &
					ax2_m - de_2*(i-1) /)
				end do
			end do
        
        else if (p_t .eq. 3)then ! xz-plot
            ax1_min = -ax1_m
		    ax2_min = -ax2_m
		    de_1 = (ax1_m - ax1_min)/(n1-1)
		    de_2 = (ax2_m - ax2_min)/(n2-1)
			do i=1, n2
				do j=1, n1
					r_local((i-1)*n1+j)%point=(/ax1_min + de_1*(j-1), ref_pl, &
					ax2_min + de_2*(i-1)/)
				end do
            end do
        else if (p_t .eq. 4)then ! circular - plot
            ! ax1_m: rho_max, 
            ! ax2_m: phi_max in xy plane 
		    de_1 = ax1_m/(n1-1)
		    de_2 = ax2_m/(n2-1)            		
			do i = 1, n2
                phi = 0 + de_2*(i-1)
				do j = 1, n1
                    rho = 0 + de_1*(j-1)
					r_local((i-1)*n1+j)%point=(/rho,  &
					phi, ref_pl/)
				end do
            end do
        else if (p_t .eq. 5)then ! circular - plot
            ax1_min = -ax1_m
            ax2_min = -ax2_m
            de_1 = (ax1_m-ax1_min)/(n1-1)
		    de_2 = (ax2_m-ax2_min)/(n2-1)              
			do i=1, n1                                
				do j=1, n2                    
					r_local((i-1)*n2 + j)%point =(/ax1_min + de_1*(i-1),  &
					ref_pl, ax2_min + de_2*(j-1)/)
				end do
            end do
		end if				
		return		
	end subroutine Observation_points
	
	subroutine Focused_Gaussian_beam(w0, r_local, E_00)	
		implicit none		
		real (kind = dp) :: f0, f, theta_max, a, x_j
		integer(kind = 4) :: m, nm, i, ngp!, pt, N1, N2, Nop    
		real(kind = dp), dimension(:), allocatable :: jn, djn, xx, ww
		real(kind = dp) :: c1, c2, rho, phi_r, fw, j0, j1, j2
		real(kind = dp) :: E0, r_max, phi_m, z_max
		real(kind = dp) :: indx_2, NA
		real(kind = dp) :: w0
		real(kind = dp) :: r0(3) !on of the focused beam
		complex(kind = dp) :: p1, p2, p3, I_00, I_01, I_02, trf    
		real(kind = dp) :: r_local(3)
		complex(kind = dp), intent(out) :: E_00(3)		
		intrinsic :: sin, sqrt    
		
		ngp = 7
		nm = 5
		m = 3
		r_max = 2e-6 !radius of the field to be calculated
		phi_m = 2*PI
		z_max = 5e-6 !half of the calculation range along the z-axis

		r0 = (/0.0*D_factor, 0.0*D_factor, 0.0*D_factor/) !center position of the focused point		
		E0 = 1.0
		indx_2 = 1.0
		NA = 0.8
		a = 0.0
		theta_max = asin(NA/indx_2)    
		f = 300e-6
		f0 = w0/(f*sin(theta_max))    
		
		allocate(jn(1:m))
		allocate(djn(1:m))    
		allocate(ww(1:ngp))
		allocate(xx(1:ngp))
		call legendre_handle(ngp, a, theta_max, xx, ww)
		!Ex ploarized light focused onto z = 0
		!Eq. (3.49 - 3.66) Novotny, Nano-Optics
		rho = sqrt((r_local(1) - r0(1))**2 + (r_local(2) - r0(2))**2)
		phi_r = acos((r_local(1) - r0(1))/rho)		
		I_01 = (0.0, 0.0)
		I_02 = (0.0, 0.0)
		I_00 = (0.0, 0.0)
		do i = 1, ngp        
			c1 = sin(xx(i))
			c2 = cos(xx(i))
			fw = exp(-1/f0**2*(c1/sin(theta_max))**2)             
			x_j = c1*k0*rho
			jn = lib_math_bessel_spherical_first_kind(x_j, 0, nm)
			j0 = jn(1)
			j1 = jn(3)
			j2 = jn(3)              
			trf = fw*sqrt(c2)*exp(im*k0*(r_local(3)-r0(3))*c2)*ww(i)
			I_00 = I_00 + c1*(1 + c2)*trf*j0!
			I_01 = I_01 + c1**2*j1*trf        
			I_02 = I_02 + c1*(1-c2)*j2*trf
		end do        
		p1 = I_00 + I_02*cos(2*phi_r) 
		p2 = I_02*sin(2*phi_r)
		p3 = -2*im*I_01*cos(phi_r)
		E_00 = im*k0*f/2*sqrt(sqrt(eps_r1)/indx_2)*E0*exp(-im*k0*f)*(/p1, p2, p3/) 
		return
	end subroutine
        
    subroutine Gaussian_beam_theta(w0, r_local, Es)!, theta, pol, lambda, r0
        real(kind = dp), intent(in) :: w0, r_local(3)
        real(kind = dp) :: zR, ra(3), Rz, wz, phi, A(3, 3), op, expf
        complex(kind = dp), intent(out) :: Es(3)
        intrinsic :: exp, sin, cos
        integer :: i

        zR = PI*w0**2/lambda
        op = 0.0
		! A: rotation matrix
        A(1, :) = (/cos(theta_in), op, -sin(theta_in)/)
        A(2, :) = (/op, cos(theta_in)/cos(theta_in), op/)
        A(3, :) = (/sin(theta_in), op, cos(theta_in)/)

        do i = 1, 3
            ra(i) = dot_product(A(i, :), r_local(:))
        end do

        if (ra(3) == 0.0) then
            Rz = 0.0
            phi = 0.0
            wz = w0
        else
            Rz = ra(3)*(1 + (zR/ra(3))**2)
            phi = atan(ra(3)/zR)
            wz = w0*sqrt(1+(ra(3)/zR)**2)
        end if

        expf = -(ra(1)**2+ra(2)**2)/wz**2
        if (expf < -100) then
            Es = 0.0
        else if (ra(3) == 0.0) then
            Es = E_in_hat*exp(-(ra(1)**2+ra(2)**2)/wz**2)*exp(im*phi)!
		else
			Es = E_in_hat*w0/wz*exp(-(ra(1)**2+ra(2)**2)/wz**2)*exp(-im*(k0*ra(3) + k0*(ra(1)**2+ra(2)**2)/(2*Rz)-phi))!
        end if	
	return
    end subroutine
end module disc_mod
