! File: lib_tree_helper_functions.F90
!
! Created on Tue Feb  5 16:33:05 2019
!
! @author: itodaiber
!
!
! Notes
! ----
!
! Preprocessor
! -----
!   standard value: _FMM_DIMENSION_ = 3
!
! Eclipse settings
! ------
!   Project properties -> Fortran General -> Paths and Symbols -> Symbols
!
! Limits
! ----
! Universal index
! -----
!   2D
!   single precision:
!    coordinate(kind=4) ->
!
!   3D
!
!
!
!


! spatial dimension, value = [2,3]
#define _FMM_DIMENSION_ 2

! 1: true, 0: false (-> spatial point is real)
#define _SPATIAL_POINT_IS_DOUBLE_ 1

! number of bytes of the universal index, value = [4,8,16]
! standard value: 8
!
! Constraint
! ----
!          |  _FMM_DIMENSION_  |
!   value  |    2     |   3    |
!   -----------------------------
!   single | [4,8]    | [8]    |
!   double | [8,16]   | [8,16] |
!
#define _UINDEX_BYTES_ 8



module lib_tree_helper_functions
    !$  use omp_lib
    use lib_tree_type
    use file_io
    implicit none

    private

    interface reallocate
        module procedure lib_tree_hf_reallocate_1d_uindex
        module procedure lib_tree_hf_reallocate_1d_data_element_list
    end interface

    interface concatenate
        module procedure lib_tree_hf_concatenate_1d_data_element_list
        module procedure lib_tree_hf_concatenate_1d_data_element_list_single
    end interface

    ! parameter
    integer(kind=1), public, parameter :: NUMBER_OF_BITS_PER_BYTE = 8
    integer(kind=UINDEX_BYTES), public, parameter :: TREE_BOX_IGNORE_ENTRY = -1
    integer(kind=1), public, parameter :: UNIVERSAL_INDEX_OVERFLOW = 0

    ! integer kind of the bit interleaving process, default = 1
    !
    ! TODO: bug fix value > 1
    integer(kind=1), private, parameter :: INTERLEAVE_BITS_INTEGER_KIND = 1

    ! module global variable
#if (_FMM_DIMENSION_ == 2)
    logical :: lib_tree_interleave_bits_lut_initialised = .false.
    integer(kind=INTERLEAVE_BITS_INTEGER_KIND), dimension (:,:,:) &
                                                , allocatable :: lib_tree_interleave_bits_lut

    logical :: lib_tree_deinterleave_bits_lut_initialised = .false.
    integer(kind=INTERLEAVE_BITS_INTEGER_KIND), dimension (:,:,:) &
                                                , allocatable :: lib_tree_deinterleave_bits_lut
#elif (_FMM_DIMENSION_ == 3)
    integer(kind=INTERLEAVE_BITS_INTEGER_KIND), dimension (:,:,:,:) &
                                                , allocatable :: lib_tree_interleave_bits_lut
    logical :: lib_tree_interleave_bits_lut_initialised = .false.

    integer(kind=INTERLEAVE_BITS_INTEGER_KIND), dimension (:,:,:,:) &
                                                , allocatable :: lib_tree_deinterleave_bits_lut
    logical :: lib_tree_deinterleave_bits_lut_initialised = .false.
#endif
    ! ~ module global variable ~

    ! public member functions
    public :: reallocate
    public :: concatenate
    public :: lib_tree_hf_destructor

    public :: lib_tree_spatial_point
    public :: lib_tree_universal_index

    public :: lib_tree_hf_get_universal_index
    public :: lib_tree_hf_get_parent

    public :: lib_tree_hf_get_children_all
    public :: lib_tree_hf_get_centre_of_box
    public :: lib_tree_hf_get_neighbour_all_xD

!    public :: lib_tree_hf_get_neighbourhood_size

    ! test functions
    public :: lib_tree_hf_test_functions
    public :: lib_tree_hf_benchmark

contains

    ! Cleans up the memory
    subroutine lib_tree_hf_destructor()
        implicit none

        !$OMP SINGLE
#if (_FMM_DIMENSION_ == 2)
        if (lib_tree_interleave_bits_lut_initialised) then
           lib_tree_interleave_bits_lut_initialised = .false.
           deallocate (lib_tree_interleave_bits_lut)
        end if
        if (lib_tree_deinterleave_bits_lut_initialised) then
            lib_tree_deinterleave_bits_lut_initialised = .false.
           deallocate (lib_tree_deinterleave_bits_lut)
        end if
#elif (_FMM_DIMENSION_ == 3)
        if (lib_tree_interleave_bits_lut_initialised) then
            lib_tree_interleave_bits_lut_initialised = .false.
            deallocate (lib_tree_interleave_bits_lut)
        end if

        if (lib_tree_deinterleave_bits_lut_initialised) then
            lib_tree_deinterleave_bits_lut_initialised = .false.
            deallocate (lib_tree_deinterleave_bits_lut)
        end if
#endif
        !$OMP END SINGLE

    end subroutine lib_tree_hf_destructor

    ! Calculates the universal index *n* of a given normalised floating point *point_x*.
    ! Related to the level *l*.
    !
    !   n = (2**d)**(l-1)*N_1 + (2**d)**(l-2)*N_2 + ... + (2**d)*N_l-1 + N_l      (49)
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
    !
    !
    ! Example:
    !   Point |  x_10    |  x_2    |  n(l=2)
    !   --------------------------------------------
    !      1  |  0.125   | 0.001   |  2*0 + 1*0 = 0
    !      2  |  0.3125  | 0.0101  |  2*0 + 1*1 = 1
    !      3  |  0.375   | 0.011   |  2*0 + 1*1 = 1
    !      4  |  0.625   | 0.101   |  2*1 + 1*0 = 2
    !      5  |  0.825   | 0.111   |  2*1 + 1*1 = 3
    !
    !
    ! Arguments
    ! ----
    !   point_x: type(lib_tree_spatial_point)
    !       normalised floating point number (0.0 .. 1.0)
    !       HINT: datatype real is also possible, use *_1D_float() instead
    !
    !   l: integer(kind=1)
    !       number of layers
    !
    ! Returns
    ! ----
    !   the universal index *uindex*.
    !
    !   uindex: type(lib_tree_universal_index)
    !
    function lib_tree_hf_get_universal_index(point_x, l) result(uindex)
        implicit none

        ! dummy arguments
        type(lib_tree_spatial_point), intent(in) :: point_x
        integer(kind=1), intent (in) :: l
        type(lib_tree_universal_index) :: uindex

        ! auxiliary
        integer(kind=COORDINATE_BINARY_BYTES), dimension(TREE_DIMENSIONS) :: coordinate_binary

        ! Example
        ! ----
        !
        ! coordinate_binary: kind=4, dimension=3                                            2**-1
        ! element   byte representation                          base(2)                     v        base(10)
        ! 1:       |---04---|---03---|---02---|---01---|   e.g.: |00000000|00000000|00000000|11000000|  0.75
        ! 2:       |---04---|---03---|---02---|---01---|   e.g.: |00000000|00000000|00000000|10000000|  0.5
        ! 3:       |---04---|---03---|---02---|---01---|   e.g.: |00000000|00000000|10000000|00000001|  0.005859375
        !
        ! cb_buffer: kind=1, dimension=3*4=12
        ! element   byte representation                          base(2)                                base(10)
        ! 1-4:     |---01---|---02---|---03---|---04---|   e.g.: |00000000|00000000|00000000|11000000|  0.75
        ! 5-8:     |---05---|---06---|---07---|---08---|   e.g.: |00000000|00000000|00000000|10000000|  0.5
        ! 9-12:    |---09---|---10---|---11---|---12---|   e.g.: |00000000|00000000|10000000|00000001|  0.005859375
        !                                     <-- *
        !                    interleave column wise
        !
        ! equivalence (coordinate_binary, cb_buffer)
        !
        ! Interleave bits
        ! -----
        !
        ! interleaved_bits: kind=1, dimension=3*4=12
        ! element   byte representation
        ! 1-4:     |04-08-12|04-08-12|04-08-12|03-07-11|   e.g.: |11010000|00000000|00000001|00100000|
        ! 5-8:     |03-07-11|03-07-11|02-06-10|02-06-10|   e.g.: |00000000|00000000|00000000|00000000|
        ! 9-12:    |02-06-10|01-05-09|01-05-09|01-05-09|   e.g.: |00000000|00000000|00000000|00000000|
        !
        !   xx-yy-zz describes the bytes which were interleaved
        !
        ! ib_buffer: kind=1, dimension=3
        ! element   byte representation
        ! 1-3:     |---03---|---02---|---01---|
        !
        !
        ! ib_buffer = interleave(cb_buffer(9), cb_buffer(5), cb_buffer(1))
        ! iterleaved_bits(1:3) = ib_buffer
        !
        ! ...
        !
        ! ib_buffer = interleave(cb_buffer(12), cb_buffer(8), cb_buffer(4))
        ! interleaved_bits(10:12) = ib_buffer
        !
        integer(kind=INTERLEAVE_BITS_INTEGER_KIND) &
            ,dimension(TREE_DIMENSIONS * COORDINATE_BINARY_BYTES/INTERLEAVE_BITS_INTEGER_KIND) &
            :: interleaved_bits
        !
        ! doubel precision
        !   Number of decimal digits: ca. 16
        !   bit precision: 53
        !
        ! 3d example
        ! ----
        ! smalest cube
        !   edge length: 10**-16
        !   volume: 10**-48
        !
        ! largest cube
        !   edge length: 1
        !   volume: 1
        !
        ! number of smalest cubes in the biggest cube
        !   number = 1 / 10**-48 = 10**48
        !
        ! determination of the integer kind
        ! -----
        ! 32 bit word range: −(2**31) to 2**31 − 1
        ! 64 bit word range: −(2**63) to 2**63 − 1
        !
        !
        integer(kind=UINDEX_BYTES) :: interleaved_bits_dimension_0

        !        interleaved_bits number of bytes/INTERLEAVE_BITS_INTEGER_KIND     -   output number of bytes/INTERLEAVE_BITS_INTEGER_KIND + 1
        !   TREE_DIMENSIONS * COORDINATE_BINARY_BYTES/INTERLEAVE_BITS_INTEGER_KIND - UINDEX_BYTES /INTERLEAVE_BITS_INTEGER_KIND + 1
        ! = ( TREE_DIMENSIONS*COORDINATE_BINARY_BYTES - UINDEX_BYTES)/INTERLEAVE_BITS_INTEGER_KIND + 1
        equivalence (interleaved_bits(( TREE_DIMENSIONS*COORDINATE_BINARY_BYTES - UINDEX_BYTES)/INTERLEAVE_BITS_INTEGER_KIND + 1), &
                     interleaved_bits_dimension_0)

        integer(kind=INTERLEAVE_BITS_INTEGER_KIND) &
            ,dimension(TREE_DIMENSIONS * COORDINATE_BINARY_BYTES/INTERLEAVE_BITS_INTEGER_KIND) &
            :: cb_buffer

        equivalence (coordinate_binary, cb_buffer)

        integer(kind=INTERLEAVE_BITS_INTEGER_KIND), dimension(TREE_DIMENSIONS) :: ob_buffer_DIM
        integer(kind=INTERLEAVE_BITS_INTEGER_KIND), dimension(TREE_DIMENSIONS) :: ib_buffer

        integer(kind=1) :: i
        integer(kind=1) :: ii
#if (_UINDEX_BYTES_ == 16)
        integer(kind=16) :: buffer_diff
        integer(kind=8), dimension(2) :: buffer_diff_a
        equivalence(buffer_diff, buffer_diff_a)
#else
        integer(kind=8) :: buffer_diff
#endif

        ! calculate the x-dimensional binary coordinate
        coordinate_binary = lib_tree_hf_get_coordinate_binary_number_xD(point_x%x)

        ! interleave bits
        do i=COORDINATE_BINARY_BYTES/INTERLEAVE_BITS_INTEGER_KIND, 1, -1 ! interleave column wise
            do ii=1, TREE_DIMENSIONS  ! get column entries
                ob_buffer_DIM(ii) = cb_buffer(i + (ii-1)*COORDINATE_BINARY_BYTES/INTERLEAVE_BITS_INTEGER_KIND)
            end do
            ib_buffer = lib_tree_hf_interleave_bits_use_lut(ob_buffer_DIM)

            ! e.g.    12                4       (4..1)        3
            ! ii = total length - (total columns - i + 1) * length(ib_buffer) + 1
            ! ii = TREE_DIMENSIONS * COORDINATE_BINARY_BYTES/INTERLEAVE_BITS_INTEGER_KIND - (COORDINATE_BINARY_BYTES/INTERLEAVE_BITS_INTEGER_KIND - i + 1) * TREE_DIMENSIONS + 1
            ! ii = TREE_DIMENSIONS * (COORDINATE_BINARY_BYTES/INTERLEAVE_BITS_INTEGER_KIND - COORDINATE_BINARY_BYTES/INTERLEAVE_BITS_INTEGER_KIND + i - 1) + 1
            ! ii = TREE_DIMENSIONS * (i - 1) + 1
            ii = int(TREE_DIMENSIONS * (i-1) + 1, 1)

            interleaved_bits(ii:ii+TREE_DIMENSIONS-1) = ib_buffer(:)
        end do

        ! calculate the universal index n
#if (_UINDEX_BYTES_ == 16)
        buffer_diff = int(UINDEX_BYTES,16)*int(NUMBER_OF_BITS_PER_BYTE,16) - int(l,16)*int(TREE_DIMENSIONS,16)
        i = int(ishft(buffer_diff,-15),1)
        if ( i .eq. 0 ) then ! not a negativ number
            uindex%n = ibits(interleaved_bits_dimension_0, &
                            UINDEX_BYTES*NUMBER_OF_BITS_PER_BYTE-l*TREE_DIMENSIONS, &
                            l*TREE_DIMENSIONS)
        else if ( (buffer_diff_a(1) .eq. 0) .and. (buffer_diff_a(2) .eq. 0))  then
            uindex%n = interleaved_bits_dimension_0
        else
            uindex%n = UNIVERSAL_INDEX_OVERFLOW
            print *, "lib_tree_hf_get_universal_index ERROR: universal index overflow"
            print *, "  l*TREE_DIMENSIONS = ", l*TREE_DIMENSIONS
        end if
#else
        buffer_diff = UINDEX_BYTES*NUMBER_OF_BITS_PER_BYTE - l*TREE_DIMENSIONS
        if (buffer_diff .gt. 0) then
            uindex%n = ibits(interleaved_bits_dimension_0, &
                            UINDEX_BYTES*NUMBER_OF_BITS_PER_BYTE-l*TREE_DIMENSIONS, &
                            l*TREE_DIMENSIONS)
        else if (buffer_diff .eq. 0) then
            uindex%n = interleaved_bits_dimension_0
        else
            uindex%n = UNIVERSAL_INDEX_OVERFLOW
            print *, "lib_tree_hf_get_universal_index ERROR: universal index overflow"
            print *, "  l*TREE_DIMENSIONS = ", l*TREE_DIMENSIONS
        end if
#endif
!        uindex%interleaved_bits = interleaved_bits

        uindex%l = l

    end function lib_tree_hf_get_universal_index

    ! Calculates the binary coordinate of a normalised floating point (0..1).
    !
    !   String(n, l)=(N_1, N_2, ..., N_l), N_j=0, ..., 2**d−1, j=1, ..., l (48)
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
    !
    ! Example:
    !   Point |  x_10    |  x_2
    !   --------------------------
    !      1  |  0.125   | 0.001    -> String(n,l)=(0,0,1)
    !      2  |  0.3125  | 0.0101   -> String(n,l)=(0,1,0,1)
    !      3  |  0.375   | 0.011    -> String(n,l)=(0,1,1)
    !      4  |  0.625   | 0.101    -> String(n,l)=(1,0,1)
    !      5  |  0.825   | 0.111    -> String(n,l)=(1,1,1)
    !
    !
    !
    ! Floting point bitwise structure
    !   Bit no.: 7 6 5 4   3 2 1 0
    !   Byte 4: |S|E|E|E| |E|E|E|E|
    !   Byte 3: |E|M|M|M| |M|M|M|M|
    !   Byte 2: |M|M|M|M| |M|M|M|M|
    !   Byte 1: |M|M|M|M| |M|M|M|M|
    !
    ! legend:
    !   S: algebraic sign
    !   E: Exponent
    !   M: Mantissa
    !   Point x_a: (base a)
    !
    !
    ! Arguments
    ! ----
    !   f: float
    !       normalised single precision floating point number (0.0 .. 1.0)
    !
    ! Returns
    ! ----
    !   the binary representation of the floating point number (only the decimal place).
    !
    !   coordinate_binary: 4 bytes
    !
    function lib_tree_hf_get_coordinate_binary_number_1D_float(f) result (coordinate_binary)
        implicit none

        ! dummy arguments
        real, intent(in) :: f
        integer(kind=4) :: coordinate_binary

        ! variable for the binary access
        real :: f_buffer
        byte, dimension(4) :: f_byte
        equivalence (f_buffer, f_byte)

        byte :: f_exponent
        integer(kind=4) :: f_mantissa
        integer(kind=4) :: f_integer_buffer
        equivalence (f_integer_buffer, f_byte)

        ! auxiiary variables
        integer(kind=4) :: shift

        ! parametres
        integer(kind=1), parameter :: INTEGER_SIGNED_MIN_ABS = 127
        integer(kind=1), parameter :: BITS_SIGN = 1
        integer(kind=1), parameter :: BITS_EXPONENT = 8
        integer(kind=1), parameter :: BITS_MANTISSA = 23


        f_buffer = f

        ! --- extract the exponent from byte 4 and 3 ---
        ! Bit no.: 7 6 5 4   3 2 1 0
        ! Byte 4: |S|E|E|E| |E|E|E|E|
        ! Byte 3: |E|M|M|M| |M|M|M|M|
        !
        ! legend:
        !   S: algebraic sign
        !   E: Exponent
        !   M: Mantissa
        f_exponent = ishft(f_byte(4), BITS_SIGN)

        if (btest(f_byte(3), 7)) then
            f_exponent = ibset(f_exponent, 0)
        else
            f_exponent = ibclr(f_exponent, 0)
        end if

        ! --- extract the mantissa from byte 3, 2 and 1 ---
        ! f_integer_buffer:
        !   Bit no.: 7 6 5 4   3 2 1 0
        !   Byte 4: |S|E|E|E| |E|E|E|E|
        !   Byte 3: |E|M|M|M| |M|M|M|M|
        !   Byte 2: |M|M|M|M| |M|M|M|M|
        !   Byte 1: |M|M|M|M| |M|M|M|M|
        !
        ! f_mantissa
        !   Bit no.: 7 6 5 4   3 2 1 0
        !   Byte 4: |1|M|M|M| |M|M|M|M|   bit 7: "virtual 1 of mantissa: 1.MMMMM"
        !   Byte 3: |M|M|M|M| |M|M|M|M|
        !   Byte 2: |M|M|M|M| |M|M|M|M|
        !   Byte 2: |0|0|0|0| |0|0|0|0|

        f_mantissa = ishft(f_integer_buffer, BITS_SIGN + BITS_EXPONENT - 1)
        f_mantissa = ibset(f_mantissa, 31)  ! set virtual 1 of mantissa

        ! --- convert exponent ---
        ! binary representation: sigend integer, but interpreted as unsigned integer
        ! e.g.: -2 (base 10) = 0111 1101 (base 2) displayed as 125 (base 10)
        !
        ! -> f_exponent = INT_MIN_ABS - f_exponent
        !               = 127 - f_exponent
        !

        ! --- generate binary coordinate (only the decimal place) ---
        shift = INTEGER_SIGNED_MIN_ABS-f_exponent-1    ! -1: to respect the virtual 1 of mantissa
        coordinate_binary = ishft(f_mantissa, -shift)

    end function lib_tree_hf_get_coordinate_binary_number_1D_float

    ! Calculates the binary coordinate of a normalised floating point (0..1).
    !
    !   String(n, l)=(N_1, N_2, ..., N_l), N_j=0, ..., 2**d−1, j=1, ..., l (48)
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
    !
    ! Example:
    !   Point |  x_10    |  x_2
    !   --------------------------
    !      1  |  0.125   | 0.001    -> String(n,l)=(0,0,1)
    !      2  |  0.3125  | 0.0101   -> String(n,l)=(0,1,0,1)
    !      3  |  0.375   | 0.011    -> String(n,l)=(0,1,1)
    !      4  |  0.625   | 0.101    -> String(n,l)=(1,0,1)
    !      5  |  0.825   | 0.111    -> String(n,l)=(1,1,1)
    !
    !
    !
    ! Floting point bitwise structure
    !   Bit no.: 7 6 5 4   3 2 1 0
    !   Byte 7: |S|E|E|E| |E|E|E|E|
    !   Byte 6: |E|E|E|E| |M|M|M|M|
    !   Byte 5: |M|M|M|M| |M|M|M|M|
    !   Byte 4: |M|M|M|M| |M|M|M|M|
    !   Byte 3: |M|M|M|M| |M|M|M|M|
    !   Byte 2: |M|M|M|M| |M|M|M|M|
    !   Byte 1: |M|M|M|M| |M|M|M|M|
    !   Byte 0: |M|M|M|M| |M|M|M|M|
    !
    ! legend:
    !   S: algebraic sign
    !   E: Exponent
    !   M: Mantissa
    !   Point x_a: (base a)
    !
    !
    ! Arguments
    ! ----
    !   f: float
    !       normalised single precision floating point number (0.0 .. 1.0)
    !
    ! Returns
    ! ----
    !   the binary representation of the floating point number (only the decimal place).
    !
    !   coordinate_binary: 8 bytes
    !
    function lib_tree_hf_get_coordinate_binary_number_1D_double(f) result (coordinate_binary)
        implicit none

        ! dummy arguments
        double precision, intent(in) :: f
        integer(kind=8) :: coordinate_binary

        ! variable for the binary access
        double precision :: f_buffer
        integer(kind=8) :: f_integer_buffer
        equivalence (f_integer_buffer, f_buffer)

        integer(kind=8) :: f_exponent
        integer(kind=8) :: f_mantissa


        ! auxiiary variables
        integer(kind=8) :: shift

        ! parametres
        integer(kind=2), parameter :: EXPONENT_CONVENTION = 1023
        integer(kind=1), parameter :: BITS_SIGN = 1
        integer(kind=1), parameter :: BITS_EXPONENT = 11
        integer(kind=1), parameter :: BITS_MANTISSA = 52


        f_buffer = f

        ! --- extract the exponent from byte 7 and 6 ---
        !  Bit no.: 7 6 5 4   3 2 1 0
        !   Byte 7: |S|E|E|E| |E|E|E|E|
        !   Byte 6: |E|E|E|E| |M|M|M|M|
        !
        ! Result: f_exponent
        !  Bit no.: 7 6 5 4   3 2 1 0
        !   Byte 1: |0|0|0|0| |0|E|E|E|
        !   Byte 0: |E|E|E|E| |E|E|E|E|
        !
        ! legend:
        !   S: algebraic sign
        !   E: Exponent
        !   M: Mantissa
        f_exponent = ishft(f_integer_buffer, -BITS_MANTISSA)    ! shift right

        ! --- extract the mantissa from byte 0-6 ---
        ! f_integer_buffer:
        !   Bit no.: 7 6 5 4   3 2 1 0
        !   Byte 7: |S|E|E|E| |E|E|E|E|
        !   Byte 6: |E|E|E|E| |M|M|M|M|
        !   Byte 5: |M|M|M|M| |M|M|M|M|
        !   Byte 4: |M|M|M|M| |M|M|M|M|
        !   Byte 3: |M|M|M|M| |M|M|M|M|
        !   Byte 2: |M|M|M|M| |M|M|M|M|
        !   Byte 1: |M|M|M|M| |M|M|M|M|
        !   Byte 0: |M|M|M|M| |M|M|M|M|
        !
        ! Result: f_mantissa
        !   Bit no.: 7 6 5 4   3 2 1 0
        !   Byte 7: |1|M|M|M| |M|M|M|M|   bit 7: "virtual 1 of mantissa: 1.MMMMM"
        !   Byte 6: |M|M|M|M| |M|M|M|M|
        !   Byte 5: |M|M|M|M| |M|M|M|M|
        !   Byte 4: |M|M|M|M| |M|M|M|M|
        !   Byte 3: |M|M|M|M| |M|M|M|M|
        !   Byte 2: |M|M|M|M| |M|M|M|M|
        !   Byte 1: |M|M|M|M| |M|0|0|0|
        !   Byte 0: |0|0|0|0| |0|0|0|0|
        f_mantissa = ishft(f_integer_buffer, BITS_SIGN + BITS_EXPONENT - 1)
        f_mantissa = ibset(f_mantissa, 63)


        ! --- convert exponent ---
        ! binary representation: sigend integer, but interpreted as unsigned integer
        ! e.g.: -2 (base 10) = 0000 0011 1111 1101 (base 2) displayed as 1021 (base 10)
        !
        ! -> f_exponent = EXPONENT_CONVENTION - f_exponent
        !               = 1023 - f_exponent
        !

        ! --- generate binary coordinate (only the decimal place) ---
        shift = EXPONENT_CONVENTION - int(f_exponent) - 1
        coordinate_binary = ishft(f_mantissa, -shift)   ! right shift

    end function lib_tree_hf_get_coordinate_binary_number_1D_double
    
    ! Calculates the binary coordinate of a normalised floating point vector (0..1, 0..1, 0..1).
    ! For each element separately.
    !
    ! Arguments
    ! ----
    !   f: double precision, dimension(TREE_DIMENSIONS)
    !       normalised single precision floating point number (0.0..1.0, 0.0..1.0, 0.0..1.0)
    !
    ! Returns
    ! ----
    !   the binary representation of the floating point vector (only the decimal place).
    !
    !   coordinate_binary_D: vector<int(kind=8)>
    !
    function lib_tree_hf_get_coordinate_binary_number_xD(f) result (coordinate_binary_xD)
        implicit none

        ! dummy arguments
#if(_SPATIAL_POINT_IS_DOUBLE_ == 1)
        double precision, dimension(:), intent(in) :: f
#elif(_SPATIAL_POINT_IS_DOUBLE_ == 0)
        real, dimension(:), intent(in) :: f
#endif
        integer(kind=COORDINATE_BINARY_BYTES), dimension(size(f)) :: coordinate_binary_xD

        ! auxiliary variables
#if(_SPATIAL_POINT_IS_DOUBLE_ == 1)
        double precision :: f_buffer
#elif(_SPATIAL_POINT_IS_DOUBLE_ == 0)
        real :: f_buffer
#endif
        integer(kind=1) :: i

        do i = 1, int(size(f), 1)
            f_buffer = f(i)
#if(_SPATIAL_POINT_IS_DOUBLE_ == 1)
            coordinate_binary_xD(i) = lib_tree_hf_get_coordinate_binary_number_1D_double(f_buffer)
#elif(_SPATIAL_POINT_IS_DOUBLE_ == 0)
            coordinate_binary_xD(i) = lib_tree_hf_get_coordinate_binary_number_1D_float(f_buffer)
#endif
        end do

    end function lib_tree_hf_get_coordinate_binary_number_xD

    ! Calculates the universal index of the parent box.
    !
    !   Parent(n) = [n/(2^d)]     (57)
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
    !
    ! case: one-dimensional tree
    ! l=1
    !  |     0     |     1     |
    !  -------------------------
    ! l=2             ^
    !                 |
    !  |  0  |  1  |  2  |  3  |
    !  -------------------------
    !
    ! example:
    !   Box(n,l) = (2,2)
    !   Parent(2,2) = (1,1)
    !
    !
    ! Arguments
    ! ----
    !   uindex: type(lib_tree_universal_index)
    !       universal index of a box
    !   step: integer, optional (standard value: 1)
    !       value = 1 : function returns the universal index of the parent box
    !       value >= 1: function returns the universal index of the (grand*step)parent box
    !
    ! Returns
    ! ----
    !   the universal index of the parent box.
    !
    !   parent_uindex: type(lib_tree_universal_index)
    !
    function lib_tree_hf_get_parent(n, step) result (parent_n)
        implicit none

        ! dummy arguments
        integer(kind=UINDEX_BYTES), intent (in) :: n
        integer(kind=1), intent (in), optional :: step
        integer(kind=UINDEX_BYTES) :: parent_n

        ! auxiliary
        integer(kind=1) :: m_step

        m_step = 1
        if(present(step))m_step=step

        parent_n = ishft(n,-(TREE_DIMENSIONS*m_step))

    end function

    ! Calculates the universal index of all children's boxes
    !
    !   ChildrenAll(n) = {2^d * n + j}, j=0, ..., 2^d − 1   (58)
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
    !
    ! example: one-dimensional tree
    ! ----
    !    l=1
    !     |     0     |     1     |
    !     -------------------------
    !                      / \
    !    l=2              v   v
    !     |  0  |  1  |  2  |  3  |
    !     -------------------------
    !
    !      Box(n,l) = (1,1)
    !      ChildrenAll(1,1) = {(2,2), (3,2)}
    !
    ! Arguments
    ! ----
    !   n: integer(kind=4)
    !       universal index of a box
    !
    ! Returns
    ! ----
    !   the universal indexes of all children boxes of the given box.
    !
    !   children_n: Integer(kind=4), dimension(2^d)
    function lib_tree_hf_get_children_all(n) result (children_n)
        implicit none

        ! dummy arguments
        integer(kind=UINDEX_BYTES), intent (in) :: n
        integer(kind=UINDEX_BYTES), dimension(2**TREE_DIMENSIONS) :: children_n

        ! auxiliary variables
        integer(kind=1) :: j

        do j = 0, 2**TREE_DIMENSIONS - 1
            children_n(j+1) = ishft(1,TREE_DIMENSIONS) * n + j
        end  do

    end function

    ! Calculates the centre of a box(n,l).
    !
    !   x_m(n, l) = 2^(−l) ( n + 2^(−1))       (87)
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
    !
    ! Arguments
    ! ----
    !   n: integer
    !       number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
    !   l: integer
    !       number of the level
    ! Returns
    ! ----
    !   x_m: [float, double]
    !       counting system undependet value
    !
    function lib_tree_hf_get_centre_of_box(n,l) result (point)
        implicit none

        ! dummy arguments
        integer(kind=UINDEX_BYTES), intent (in) :: n
        integer(kind=1), intent (in) :: l
        type(lib_tree_spatial_point) :: point

        ! auxiliary variables
        integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS) :: n_treeD
        integer(kind=1) :: i
        integer(kind=UINDEX_BYTES) :: buffer

        buffer = n

        n_treeD = lib_tree_hf_deinterleave_bits_1D_to_treeD(buffer)

        do i=1, TREE_DIMENSIONS
#if (_SPATIAL_POINT_IS_DOUBLE_ == 1)
            point%x(i) = 2.0D0**(-l) * (n_treeD(i) + 0.5D0)
#elif (_SPATIAL_POINT_IS_DOUBLE_ == 0)
            point%x(i) = 2.0**(-l) * (n_treeD(i) + 0.5)
#endif
        end do

    end function lib_tree_hf_get_centre_of_box

    ! Calculates the universal index of all k-th neigbour's boxes
    !
    !   N_(min)^(Neighbours)(d) = 2^d − 1                (19)
    !   N_(max)^(Neighbours)(d) = 3^d − 1                (19)
    !
    !   NeighborAll^(k)(n, l) = {(n−k, l), (n+k, l)}     (74)
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
    !
    ! example: one-dimensional tree
    ! ----
    !    l=1
    !     |     0     |     1     |
    !     -------------------------
    !    l=2
    !     |  0  |  1  |  2  |  3  |
    !     -------------------------
    !             -k <-- * --> +k
    !
    !      Box(n,l) = (2,2)
    !      NeighbourAll(1,2,2) = {(1,2), (3,2)}
    !
    ! Arguments
    ! ----
    !   k: integer(kind=UINDEX_BYTES)
    !       k-th neighbour
    !   n: integer(kind=UINDEX_BYTES)
    !       universal index of a box
    !   l: integer(kind=UINDEX_BYTES)
    !       number of the level
    !
    ! Returns
    ! ----
    !   the universal indexes of all neigbour boxes of the given box.
    !
    !   neighbour_all: Integer(kind=UINDEX_BYTES), dimension(2^d)
    function lib_tree_hf_get_neighbour_all_1D(k,n,l) result (neighbour_1d)
        implicit none

        ! dummy arguments
        integer(kind=UINDEX_BYTES), intent (in) :: k
        integer(kind=UINDEX_BYTES), intent (in) :: n
        integer(kind=1), intent (in) :: l
        integer(kind=UINDEX_BYTES), dimension(3) :: neighbour_1d

        ! auxiliary variables
        integer(kind=UINDEX_BYTES), parameter :: ignore_entry = -1
        integer(kind=UINDEX_BYTES), parameter :: lower_boundary = 0
        integer(kind=UINDEX_BYTES) :: upper_boundary
        integer(kind=UINDEX_BYTES) :: buffer_n

        upper_boundary = 2**l - 1

        ! calculate left neighbour
        buffer_n = n-k
        if (buffer_n < lower_boundary) then
            buffer_n = ignore_entry
        end if

        neighbour_1d(1) = buffer_n

        ! own box
        neighbour_1d(2) = n

        ! calculate right neigbour
        buffer_n = n+k
        if (buffer_n > upper_boundary) then
            buffer_n = ignore_entry
        end if

        neighbour_1d(3) = buffer_n

    end function

    ! Calculates the universal index of all k-th neigbour's boxes
    !
    !   N_(min)^(Neighbours)(d) = 2^d − 1                (19)
    !   N_(max)^(Neighbours)(d) = 3^d − 1                (19)
    !
    !   NeighborAll^(k)(n, l) = {(n−k, l), (n+k, l)}     (74)
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
    !
    ! example: one-dimensional tree
    ! ----
    !    l=1
    !     |     0     |     1     |
    !     -------------------------
    !    l=2
    !     |  0  |  1  |  2  |  3  |
    !     -------------------------
    !             -k <-- * --> +k
    !
    !      Box(n,l) = (2,2)
    !      NeighbourAll(1,2,2) = {(1,2), (3,2)}
    !
    ! Arguments
    ! ----
    !   k: integer(kind=UINDEX_BYTES)
    !       k-th neighbour
    !   n: integer(kind=COORDINATE_BINARY_BYTES)
    !       universal index of a box
    !   l: integer(kind=1)
    !       number of the level
    !
    ! Returns
    ! ----
    !   the universal indexes of all neigbour boxes of the given box.
    !
    !   neighbour_all: Integer(kind=COORDINATE_BINARY_BYTES), dimension(3**TREE_DIMENSIONS-1)
    !
    !
    !
    function lib_tree_hf_get_neighbour_all_xD(k,n,l) result (neighbour_all)
        implicit none

        ! parameter
        integer(kind=UINDEX_BYTES), parameter :: lower_boundary = 0
        integer(kind=1), parameter :: NUMBER_OF_1D_NEIGHBOURS = 3   ! includes its own box number
        integer(kind=1), parameter :: NUMBER_OF_NEIGHBOURS = 3**TREE_DIMENSIONS-1

        ! dummy arguments
        integer(kind=UINDEX_BYTES), intent (in) :: k
        integer(kind=UINDEX_BYTES), intent (in) :: n
        integer(kind=1), intent (in) :: l
        integer(kind=UINDEX_BYTES), dimension(NUMBER_OF_NEIGHBOURS) :: neighbour_all

        ! auxiliary variables
        integer(kind=UINDEX_BYTES) :: upper_boundary

        integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS) :: n_deinterleaved

        integer(kind=UINDEX_BYTES), dimension(NUMBER_OF_1D_NEIGHBOURS, TREE_DIMENSIONS) :: buffer_n_1d
        integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS) :: buffer_n
        integer(kind=UINDEX_BYTES) :: buffer_neighbour

        integer(kind=1) :: i
        integer(kind=1) :: ii
#if (_FMM_DIMENSION_ == 3)
        integer(kind=1) :: iii
#endif
        integer(kind=1) :: neighbour_counter


        upper_boundary = 2**l - 1

        n_deinterleaved = lib_tree_hf_deinterleave_bits_1D_to_treeD(n)

        do i=1, TREE_DIMENSIONS
            buffer_n_1d(:,i) = lib_tree_hf_get_neighbour_all_1D(k, n_deinterleaved(i), l)
        end do

        do i=1, NUMBER_OF_NEIGHBOURS
            neighbour_all(i) = TREE_BOX_IGNORE_ENTRY
        end do

#if (_FMM_DIMENSION_ == 2)
        neighbour_counter = 0
        do i=1, NUMBER_OF_1D_NEIGHBOURS
            buffer_n(1) = buffer_n_1d(i,1)
            if ( buffer_n(1) .ne. TREE_BOX_IGNORE_ENTRY) then
                do ii=1, NUMBER_OF_1D_NEIGHBOURS
                    buffer_n(2) = buffer_n_1d(ii,2)
                    if (buffer_n(2) .ne. TREE_BOX_IGNORE_ENTRY) then
                        buffer_neighbour = lib_tree_hf_interleave_bits_treeD_to_1D(buffer_n)
                        if (buffer_neighbour .ne. n) then
                            neighbour_counter = neighbour_counter + int(1,1)
                            neighbour_all(neighbour_counter) = buffer_neighbour
                        end if
                    end if
                end do
            end if
        end do
#elif (_FMM_DIMENSION_ == 3)
        neighbour_counter = 0
        do i=1, NUMBER_OF_1D_NEIGHBOURS
            buffer_n(1) = buffer_n_1d(i,1)
            if ( buffer_n(1) .ne. TREE_BOX_IGNORE_ENTRY) then
                do ii=1, NUMBER_OF_1D_NEIGHBOURS
                    buffer_n(2) = buffer_n_1d(ii,2)
                    if (buffer_n(2) .ne. TREE_BOX_IGNORE_ENTRY) then
                        do iii=1, NUMBER_OF_1D_NEIGHBOURS
                            buffer_n(3) = buffer_n_1d(iii,3)
                            if (buffer_n(3) .ne. TREE_BOX_IGNORE_ENTRY) then
                                buffer_neighbour = lib_tree_hf_interleave_bits_treeD_to_1D(buffer_n)
                                if (buffer_neighbour .ne. n) then
                                    neighbour_counter = neighbour_counter + int(1,1)
                                    neighbour_all(neighbour_counter) = buffer_neighbour
                                end if
                            end if
                        end do
                    end if
                end do
            end if
        end do
#endif
    end function

    ! Returns the look-up table (LUT) for interleaved bits. If necessary, the LUT is recalculated.
    !
    ! *Hint*
    !   This function has only global dependencies.
    !
    ! Dependences
    ! ----
    !   FMM_DIMENSION
    !       number of dimensions
    !   INTERLEAVE_BITS_INTEGER_KIND
    !       number of bytes of the integer
    !
    ! Returns
    ! ----
    !   rv: integer, dimension(:,:,:[,,:])
    !       the look-up tabel
    function lib_tree_hf_get_interleave_bits_lut() result(rv)
        implicit none
        ! parameter
        integer(kind=1), parameter :: integer_kind = INTERLEAVE_BITS_INTEGER_KIND
        integer(kind=integer_kind), parameter :: integer_range_high = huge(integer_range_high)
        integer(kind=integer_kind), parameter :: integer_range_low = -integer_range_high-1

#if (_FMM_DIMENSION_ == 2)
        ! parameter
        character (len = *), parameter :: file_lut="pre_calc/bit_interleaving_LUT_2d.dat"

        ! dummy
        integer(kind=integer_kind), dimension (integer_range_low:integer_range_high, &
                                               integer_range_low:integer_range_high, &
                                               1:2) :: rv

        ! check if LUT is already calculated
        ! if yes: load
        ! if not: calculate
        if (file_exists(file_lut)) then
            OPEN(UNIT=14, FILE=file_lut, ACTION="read", STATUS="old", &
                 FORM='unformatted')
            READ(14) rv
            CLOSE(UNIT=14)
        else
            rv = lib_tree_hf_creat_interleave_bits_lut()
            OPEN(UNIT=13, FILE=file_lut, ACTION="write", STATUS="replace", &
                 FORM="unformatted")
            WRITE(13) rv
            CLOSE(UNIT=13)
        end if
#elif (_FMM_DIMENSION_ == 3)
        ! parameter
        character (len = *), parameter :: file_lut="pre_calc/bit_interleaving_LUT_3d.dat"

        ! dummy
        integer(kind=integer_kind), dimension (integer_range_low:integer_range_high, &
                                           integer_range_low:integer_range_high, &
                                           integer_range_low:integer_range_high, &
                                           1:3) :: rv

        ! check if LUT is already calculated
        ! if yes: load
        ! if not: calculate
        if (file_exists(file_lut)) then
            OPEN(UNIT=14, FILE=file_lut, ACTION="read", STATUS="old", &
                 FORM='unformatted')
            READ(14) rv
            CLOSE(UNIT=14)
        else
            rv = lib_tree_hf_creat_interleave_bits_lut()
            OPEN(UNIT=13, FILE=file_lut, ACTION="write", STATUS="replace", &
                 FORM="unformatted")
            WRITE(13) rv
            CLOSE(UNIT=13)
        end if
#endif

    end function lib_tree_hf_get_interleave_bits_lut

    ! Creates the look-up table (LUT) for the interleaved bits.
    !
    ! *Hint*
    !   This function has only global dependencies.
    !
    ! Dependences
    ! ----
    !   FMM_DIMENSION
    !       number of dimensions
    !   INTERLEAVE_BITS_INTEGER_KIND
    !       number of bytes of the integer
    !
    ! Returns
    ! ----
    !   rv: integer, dimension(:,:,:[,,:])
    !       the look-up tabel
    !
    function lib_tree_hf_creat_interleave_bits_lut() result(rv)
        implicit none
        ! parameter
        integer(kind=1), parameter :: integer_kind = INTERLEAVE_BITS_INTEGER_KIND
        integer(kind=integer_kind), parameter :: integer_range_high = huge(integer_range_high)
        integer(kind=integer_kind), parameter :: integer_range_low = -integer_range_high-1

#if (_FMM_DIMENSION_ == 2)
        ! dummy
        integer(kind=integer_kind), dimension (integer_range_low:integer_range_high, &
                                               integer_range_low:integer_range_high, &
                                               1:2) :: rv

        ! auxiliary
        integer(kind=integer_kind*2) :: i
        integer(kind=integer_kind*2) :: ii

        integer(kind=integer_kind), dimension(2) :: buffer
        integer(kind=1) :: p
        integer(kind=1) :: p_old

        ! ---- OMP semaphore ----
        !$  logical :: mutex_rv
        !$  mutex_rv = .true.
        ! ~~~~ OMP semaphore ~~~~

        p_old = 0
!        !$OMP PARALLEL DO PRIVATE(i, ii, buffer)
        do i = integer_range_low, integer_range_high
            do ii = integer_range_low, integer_range_high
                buffer(1) = int(i,1)
                buffer(2) = int(ii,1)
                buffer = lib_tree_hf_interleave_bits(buffer)
                !$  if (mutex_rv .eqv. .false.) then
                !$      do
                !$          if (mutex_rv .eqv. .true.) then
                !$              exit
                !$          end if
                !$      end do
                !$  end if
                !$  mutex_rv = .false.
                rv(i,ii, 1) = buffer(1)
                rv(i,ii, 2) = buffer(2)
                !$  mutex_rv = .true.
            end do
            p = int(100.0*(i-integer_range_low)/(integer_range_high-integer_range_low), 1)
            if (int(p/10, 1) .ne. int(p_old/10, 1)) then
                print *, "Interleave bits: create LUT: ", p, "%"
            end if
            p_old = p
        end do
!        !$OMP END PARALLEL DO

#elif (_FMM_DIMENSION_ == 3)
        ! dummy
        integer(kind=integer_kind), dimension (integer_range_low:integer_range_high, &
                                               integer_range_low:integer_range_high, &
                                               integer_range_low:integer_range_high, &
                                               1:3) :: rv

        ! auxiliary
        integer(kind=integer_kind*2) :: i
        integer(kind=integer_kind*2) :: ii
        integer(kind=integer_kind*2) :: iii

        integer(kind=integer_kind), dimension(3) :: buffer
        integer(kind=1) :: p
        integer(kind=1) :: p_old
        ! ---- OMP semaphore ----
        !$  logical :: mutex_rv
        !$  mutex_rv = .true.
        ! ~~~~ OMP semaphore ~~~~

        p_old = 0
!        !$OMP PARALLEL DO PRIVATE(i, ii, iii, buffer)
        do i = integer_range_low, integer_range_high
            do ii = integer_range_low, integer_range_high
                do iii = integer_range_low, integer_range_high
                    buffer(1) = int(i,1)
                    buffer(2) = int(ii,1)
                    buffer(3) = int(iii,1)
                    buffer = lib_tree_hf_interleave_bits(buffer)
                    !$  if (mutex_rv .eqv. .false.) then
                    !$      do
                    !$          if (mutex_rv .eqv. .true.) then
                    !$              exit
                    !$          end if
                    !$      end do
                    !$  end if
                    !$  mutex_rv = .false.
                    rv(i,ii, iii,1) = buffer(1)
                    rv(i,ii, iii,2) = buffer(2)
                    rv(i,ii, iii,3) = buffer(3)
                    !$  mutex_rv = .true.
                end do
            end do
            p = int(100.0*(i-integer_range_low)/(integer_range_high-integer_range_low), 1)
            if (int(p/10, 1) .ne. int(p_old/10, 1)) then
                print *, "Interleave bits: create LUT: ", p, "%"
            end if
            p_old = p
        end do
!        !$OMP END PARALLEL DO
#endif
    end function lib_tree_hf_creat_interleave_bits_lut

    ! Calculates the bit interleaving of x-dimensional integers.
    ! The kind of the integer is defined with INTERLEAVE_BITS_INTEGER_KIND.
    ! The dafault value is 1 (1 byte);
    !
    ! Argument
    ! ----
    !   x: integer, x-dimensional
    !
    !
    ! Example
    ! ----
    !   x1: 0.100   => 0.5   (base 10)
    !   x2: 0.010   => 0.25  (base 10)
    !
    !   x1: 0.1 |0 |0
    !   x2: 0. 0| 1| 0
    !  ------------------
    !  x2D: 0.10|01|00
    function lib_tree_hf_interleave_bits(x) result(rv)
        implicit none
        ! parameter
        integer(kind=1), parameter :: x_kind = INTERLEAVE_BITS_INTEGER_KIND

        ! dummy
        integer(kind=x_kind), dimension(:), intent(in) :: x
        integer(kind=x_kind), dimension(size(x)) :: rv

        ! auxiliary
        integer(kind=1) :: i
        integer(kind=1) :: ii
        integer(kind=1) :: x_dimension
        integer(kind=1) :: bit_number
        integer(kind=1) :: x_element

        x_dimension = int(size(x), 1)

        do ii = 1, x_dimension
            rv(ii) = 0
        end do

        do ii = 1, x_dimension
            do i = 0, x_kind * NUMBER_OF_BITS_PER_BYTE - 1
                ! calculates the "global" bit number
                ! two element example (1 byte / element):
                ! bit_number=9; element 2: |15 ... 8| element 1: |7 ... 0|
                bit_number = i*x_dimension+(x_dimension - ii)
                ! calculates the element; bit_number=9 => element=2
                x_element = int(bit_number / (x_kind * NUMBER_OF_BITS_PER_BYTE), 1) + int(1, 1)
                ! calculates the "local" bit_number; bit_number=9, element=2 => bit_number=1
                bit_number = bit_number - (x_element-int(1,1))*x_kind*NUMBER_OF_BITS_PER_BYTE
                if (btest(x(ii), i)) then       ! first bit number = 0
                    rv(x_element) = ibset(rv(x_element), bit_number)
                end if
            end do
        end do
    end function lib_tree_hf_interleave_bits

    ! Calculates the interleaved bits form a x-dimensional integer.
    !
    ! Arguments
    ! ----
    !   x: integer, dimension(x)
    !       integers to interleave
    !
    ! Returns
    ! ----
    !   rv: integer, dimension(:,:,:[,,:]) // 2- or 3-dimensional case
    !       interleaved integers
    !
    ! Example
    ! ----
    !   Interleaving of two one-byte integers
    !
    !   x1: 0 0 0 0| 0 0 1 0   => 2 (base 10)
    !   x2:  0 0 0 0 |0 0 0 0    => 0 (base 10)
    !      -------------------
    !       00000000|00001000   => 0|4 (base 10)
    !
    !   rv1 = 4
    !   rv2 = 0
    !
    function lib_tree_hf_interleave_bits_use_lut(x) result(rv)
        implicit none
        ! parameter
        integer(kind=1), parameter :: x_kind = INTERLEAVE_BITS_INTEGER_KIND

        integer(kind=x_kind), parameter :: integer_range_high = huge(integer_range_high)
        integer(kind=x_kind), parameter :: integer_range_low = -integer_range_high-1

        ! dummy
        integer(kind=x_kind), dimension(TREE_DIMENSIONS), intent(in) :: x
        integer(kind=x_kind), dimension(TREE_DIMENSIONS) :: rv

#if (_FMM_DIMENSION_ == 2)
        ! allocate memory
!        !$OMP SINGLE
        if (.NOT. lib_tree_interleave_bits_lut_initialised) then
            allocate( lib_tree_interleave_bits_lut(integer_range_low:integer_range_high, &
                                                     integer_range_low:integer_range_high, &
                                                     1:2) )
            lib_tree_interleave_bits_lut = lib_tree_hf_creat_interleave_bits_lut()
            lib_tree_interleave_bits_lut_initialised = .true.
        end if
!        !$OMP END SINGLE
        rv(1) = lib_tree_interleave_bits_lut(x(1), x(2), 1)
        rv(2) = lib_tree_interleave_bits_lut(x(1), x(2), 2)
#elif (_FMM_DIMENSION_ == 3)
        ! allocate memory
!        !$OMP SINGLE
        if (.NOT. lib_tree_interleave_bits_lut_initialised) then
            allocate( lib_tree_interleave_bits_lut(integer_range_low:integer_range_high, &
                                                 integer_range_low:integer_range_high, &
                                                 integer_range_low:integer_range_high, &
                                                 1:3) )

            lib_tree_interleave_bits_lut = lib_tree_hf_creat_interleave_bits_lut()
            lib_tree_interleave_bits_lut_initialised = .true.
        end if
!        !$OMP END SINGLE
        rv(1) = lib_tree_interleave_bits_lut(x(1), x(2), x(3), 1)
        rv(2) = lib_tree_interleave_bits_lut(x(1), x(2), x(3), 2)
        rv(3) = lib_tree_interleave_bits_lut(x(1), x(2), x(3), 3)
#else
        rv = lib_tree_hf_interleave_bits(x)
#endif
    end function lib_tree_hf_interleave_bits_use_lut

    ! Returns the look-up table (LUT) for deinterleaved bits. If necessary, the LUT is recalculated.
    !
    ! *Hint*
    !   This function has only global dependencies.
    !
    ! Dependences
    ! ----
    !   FMM_DIMENSION
    !       number of dimensions
    !   INTERLEAVE_BITS_INTEGER_KIND
    !       number of bytes of the integer
    !
    ! Returns
    ! ----
    !   rv: integer, dimension(:,:,:[,,:])
    !       the look-up tabel
    function lib_tree_hf_get_deinterleave_bits_lut() result(rv)
        implicit none
        ! parameter
        integer(kind=1), parameter :: integer_kind = INTERLEAVE_BITS_INTEGER_KIND
        integer(kind=integer_kind), parameter :: integer_range_high = huge(integer_range_high)
        integer(kind=integer_kind), parameter :: integer_range_low = -integer_range_high-1

#if (_FMM_DIMENSION_ == 2)
        ! parameter
        character (len = *), parameter :: file_lut="pre_calc/bit_deinterleaving_LUT_2d.dat"

        ! dummy
        integer(kind=integer_kind), dimension (integer_range_low:integer_range_high, &
                                               integer_range_low:integer_range_high, &
                                               1:2) :: rv

        ! check if LUT is already calculated
        ! if yes: load
        ! if not: calculate
        if (file_exists(file_lut)) then
            OPEN(UNIT=14, FILE=file_lut, ACTION="read", STATUS="old", &
                 FORM='unformatted')
            READ(14) rv
            CLOSE(UNIT=14)
        else
            rv = lib_tree_hf_creat_deinterleave_bits_lut()
            OPEN(UNIT=13, FILE=file_lut, ACTION="write", STATUS="replace", &
                 FORM="unformatted")
            WRITE(13) rv
            CLOSE(UNIT=13)
        end if
#elif (_FMM_DIMENSION_ == 3)
        ! parameter
        character (len = *), parameter :: file_lut="pre_calc/bit_deinterleaving_LUT_3d.dat"

        ! dummy
        integer(kind=integer_kind), dimension (integer_range_low:integer_range_high, &
                                           integer_range_low:integer_range_high, &
                                           integer_range_low:integer_range_high, &
                                           1:3) :: rv

        ! check if LUT is already calculated
        ! if yes: load
        ! if not: calculate
        if (file_exists(file_lut)) then
            OPEN(UNIT=14, FILE=file_lut, ACTION="read", STATUS="old", &
                 FORM='unformatted')
            READ(14) rv
            CLOSE(UNIT=14)
        else
            rv = lib_tree_hf_creat_deinterleave_bits_lut()
            OPEN(UNIT=13, FILE=file_lut, ACTION="write", STATUS="replace", &
                 FORM="unformatted")
            WRITE(13) rv
            CLOSE(UNIT=13)
        end if
#endif

    end function lib_tree_hf_get_deinterleave_bits_lut

    ! Creates the look-up table (LUT) for the deinterleaved bits.
    !
    ! *Hint*
    !   This function has only global dependencies.
    !
    ! Dependences
    ! ----
    !   FMM_DIMENSION
    !       number of dimensions
    !   INTERLEAVE_BITS_INTEGER_KIND
    !       number of bytes of the integer
    !
    ! Returns
    ! ----
    !   rv: integer, dimension(:,:,:[,,:])
    !       the look-up tabel
    !
    function lib_tree_hf_creat_deinterleave_bits_lut() result(rv)
        implicit none
        ! parameter
        integer(kind=1), parameter :: integer_kind = INTERLEAVE_BITS_INTEGER_KIND
        integer(kind=integer_kind), parameter :: integer_range_high = huge(integer_range_high)
        integer(kind=integer_kind), parameter :: integer_range_low = -integer_range_high-1

#if (_FMM_DIMENSION_ == 2)
        ! dummy
        integer(kind=integer_kind), dimension (integer_range_low:integer_range_high, &
                                               integer_range_low:integer_range_high, &
                                               1:2) :: rv

        ! auxiliary
        integer(kind=integer_kind*2) :: i
        integer(kind=integer_kind*2) :: ii

        integer(kind=integer_kind), dimension(2) :: buffer
        integer(kind=1) :: p
        integer(kind=1) :: p_old

        p_old = 0
        do i = integer_range_low, integer_range_high
            do ii = integer_range_low, integer_range_high
                buffer(1) = int(i,1)
                buffer(2) = int(ii,1)
                buffer = lib_tree_hf_deinterleave_bits(buffer)
                rv(i,ii, 1) = buffer(1)
                rv(i,ii, 2) = buffer(2)
            end do
            p = int(100.0*(i-integer_range_low)/(integer_range_high-integer_range_low), 1)
            if (int(p/10, 1) .ne. int(p_old/10, 1)) then
                print *, "Deinterleave bits: create LUT: ", p, "%"
            end if
            p_old = p
        end do

#elif (_FMM_DIMENSION_ == 3)
        ! dummy
        integer(kind=integer_kind), dimension (integer_range_low:integer_range_high, &
                                               integer_range_low:integer_range_high, &
                                               integer_range_low:integer_range_high, &
                                               1:3) :: rv

        ! auxiliary
        integer(kind=integer_kind*2) :: i
        integer(kind=integer_kind*2) :: ii
        integer(kind=integer_kind*2) :: iii

        integer(kind=integer_kind), dimension(3) :: buffer
        integer(kind=1) :: p
        integer(kind=1) :: p_old

        p_old = 0
        do i = integer_range_low, integer_range_high
            do ii = integer_range_low, integer_range_high
                do iii = integer_range_low, integer_range_high
                    buffer(1) = int(i,1)
                    buffer(2) = int(ii,1)
                    buffer(3) = int(iii,1)
                    buffer = lib_tree_hf_deinterleave_bits(buffer)
                    rv(i,ii, iii,1) = buffer(1)
                    rv(i,ii, iii,2) = buffer(2)
                    rv(i,ii, iii,3) = buffer(3)
                end do
            end do
            p = int(100.0*(i-integer_range_low)/(integer_range_high-integer_range_low), 1)
            if (int(p/10, 1) .ne. int(p_old/10, 1)) then
                print *, "Deinterleave bits: create LUT: ", p, "%"
            end if
            p_old = p
        end do
#endif
    end function lib_tree_hf_creat_deinterleave_bits_lut

    ! deinterleavs the
    ! The kind of the integer is defined with INTERLEAVE_BITS_INTEGER_KIND.
    ! The dafault value is 1 (1 byte);
    !
    ! Argument
    ! ----
    !   x:
    !       counting system undependet value
    !
    !
    ! Example
    ! ----
    !  x2D: 0.10|01|00
    !  ------------------
    !   x1: 0.1 |0 |0
    !   x2: 0. 0| 1| 0
    !
    !   x1: 0.100   => 0.5   (base 10)
    !   x2: 0.010   => 0.25  (base 10)
    !
    function lib_tree_hf_deinterleave_bits(x) result(rv)
        implicit none
        ! parameter
        integer(kind=1), parameter :: x_kind = INTERLEAVE_BITS_INTEGER_KIND

        ! dummy
        integer(kind=x_kind), dimension(TREE_DIMENSIONS), intent(in) :: x
        integer(kind=x_kind), dimension(TREE_DIMENSIONS) :: rv

        ! auxiliary
        integer(kind=1) :: i
        integer(kind=1) :: ii
        integer(kind=1) :: x_dimension
        integer(kind=1) :: bit_number
        integer(kind=1) :: target_element

        x_dimension = int(size(x), 1)

        do ii = 1, x_dimension
            rv(ii) = 0
        end do

        do ii = 1, x_dimension
            do i = 0, x_kind * NUMBER_OF_BITS_PER_BYTE - 1                                      ! e.g.: 16 bit =  2 byte * 8 bit / byte
                ! global bit number = (ii-1)*x_kind * NUMBER_OF_BITS_PER_BYTE + i
                ! target_element = TREE_DIMENSIONS - mod(gloabel bit number, TREE_DIMENSIONS)
                ! target local bit number = int(global bit number / TREE_DIMENSIONS)
                target_element = int((ii-1)*x_kind * NUMBER_OF_BITS_PER_BYTE + i, 1)
                bit_number = int(target_element / TREE_DIMENSIONS, 1)
                target_element = int(TREE_DIMENSIONS - mod(target_element, TREE_DIMENSIONS), 1)
                if (btest(x(ii), i)) then       ! first bit number = 0
                    rv(target_element) = ibset(rv(target_element), bit_number)
                end if
            end do
        end do
    end function lib_tree_hf_deinterleave_bits

    ! deinterleavs the
    ! The kind of the integer is defined with INTERLEAVE_BITS_INTEGER_KIND.
    ! The dafault value is 1 (1 byte);
    !
    ! Argument
    ! ----
    !   x:
    !       counting system undependet value
    !
    !
    ! Example
    ! ----
    !  x2D: 0.01|01|00
    !  ------------------
    !   x1: 0. 1| 0| 0
    !   x2: 0.0 |1 |0
    !
    !   x1: 0.100   => 0.5   (base 10)
    !   x2: 0.010   => 0.25  (base 10)
    !
    function lib_tree_hf_deinterleave_bits_use_lut(x) result(rv)
        implicit none
        ! parameter
        integer(kind=1), parameter :: x_kind = INTERLEAVE_BITS_INTEGER_KIND

        integer(kind=x_kind), parameter :: integer_range_high = huge(integer_range_high)
        integer(kind=x_kind), parameter :: integer_range_low = -integer_range_high-1

        ! dummy
        integer(kind=x_kind), dimension(TREE_DIMENSIONS), intent(in) :: x
        integer(kind=x_kind), dimension(TREE_DIMENSIONS) :: rv

#if (_FMM_DIMENSION_ == 2)
        ! allocate memory
!        !$OMP SINGLE
        if (.NOT. lib_tree_deinterleave_bits_lut_initialised) then
            allocate( lib_tree_deinterleave_bits_lut(integer_range_low:integer_range_high, &
                                                       integer_range_low:integer_range_high, &
                                                       1:2) )
            lib_tree_deinterleave_bits_lut = lib_tree_hf_creat_deinterleave_bits_lut()

            lib_tree_deinterleave_bits_lut_initialised = .true.
        end if
!        !$OMP END SINGLE
        rv(1) = lib_tree_deinterleave_bits_lut(x(1), x(2), 1)
        rv(2) = lib_tree_deinterleave_bits_lut(x(1), x(2), 2)
#elif (_FMM_DIMENSION_ == 3)
        ! allocate memory
!        !$OMP SINGLE
        if (.NOT. lib_tree_deinterleave_bits_lut_initialised) then
            allocate( lib_tree_deinterleave_bits_lut(integer_range_low:integer_range_high, &
                                                   integer_range_low:integer_range_high, &
                                                   integer_range_low:integer_range_high, &
                                                   1:3) )
            lib_tree_deinterleave_bits_lut = lib_tree_hf_creat_deinterleave_bits_lut()

            lib_tree_deinterleave_bits_lut_initialised = .true.
        end if
!        !$OMP END SINGLE
        rv(1) = lib_tree_deinterleave_bits_lut(x(1), x(2), x(3), 1)
        rv(2) = lib_tree_deinterleave_bits_lut(x(1), x(2), x(3), 2)
        rv(3) = lib_tree_deinterleave_bits_lut(x(1), x(2), x(3), 3)
#else
        rv = lib_tree_hf_deinterleave_bits(x)
#endif
    end function lib_tree_hf_deinterleave_bits_use_lut

    ! interleavs the
    ! The kind of the integer is defined with INTERLEAVE_BITS_INTEGER_KIND.
    ! The dafault value is 1 (1 byte);
    !
    ! Argument
    ! ----
    !   x:
    !       counting system undependet value
    !
    !
    ! Example
    ! ----
    !   x1: 0.100   => 0.5   (base 10)
    !   x2: 0.010   => 0.25  (base 10)
    !
    !   x1: 0.1 |0 |0
    !   x2: 0. 0| 1| 0
    !  ------------------
    !  x1D: 0.10|01|00
    !
    function lib_tree_hf_interleave_bits_treeD_to_1D(n) result(rv)
        implicit none

        ! dummy
        integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS), intent(in) :: n
        integer(kind=UINDEX_BYTES) :: rv

        ! auxiliary
        integer(kind=1) :: i
#if (_UINDEX_BYTES_ == 16)
        integer(kind=2) :: ii
        integer(kind=2) :: bit_number
#else
        integer(kind=1) :: ii
        integer(kind=1) :: bit_number
#endif

        rv = 0
        do i = 1, TREE_DIMENSIONS
            do ii=0, UINDEX_BYTES*NUMBER_OF_BITS_PER_BYTE
                ! gloabal bit number
                bit_number = ii*TREE_DIMENSIONS+(TREE_DIMENSIONS - i)

                if ( bit_number >= UINDEX_BYTES*NUMBER_OF_BITS_PER_BYTE ) then
                    exit
                end if

                if (btest(n(i), ii)) then       ! first bit number = 0
                    rv = ibset(rv, bit_number)
                end if
             end do
        end do
    end function lib_tree_hf_interleave_bits_treeD_to_1D

    function lib_tree_hf_interleave_bits_treeD_to_1D_use_lut(x) result(rv)
        implicit none
        ! dummy
        integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS), intent(in) :: x
        integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS) :: rv

        ! auxiliar
        integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS) :: coordinate_binary
        integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS) :: interleaved_coordinate_binary

        ! Example
        ! ----
        !
        ! coordinate_binary: kind=4, dimension=3                                            2**-1
        ! element   byte representation                          base(2)                     v        base(10)
        ! 1:       |---04---|---03---|---02---|---01---|   e.g.: |00000000|00000000|00000000|11000000|  0.75
        ! 2:       |---04---|---03---|---02---|---01---|   e.g.: |00000000|00000000|00000000|10000000|  0.5
        ! 3:       |---04---|---03---|---02---|---01---|   e.g.: |00000000|00000000|10000000|00000001|  0.005859375
        !
        ! cb_buffer: kind=1, dimension=3*4=12
        ! element   byte representation                          base(2)                                base(10)
        ! 1-4:     |---01---|---02---|---03---|---04---|   e.g.: |00000000|00000000|00000000|11000000|  0.75
        ! 5-8:     |---05---|---06---|---07---|---08---|   e.g.: |00000000|00000000|00000000|10000000|  0.5
        ! 9-12:    |---09---|---10---|---11---|---12---|   e.g.: |00000000|00000000|10000000|00000001|  0.005859375
        !                                     <-- *
        !                    interleave column wise
        !
        ! equivalence (coordinate_binary, cb_buffer)
        !
        ! Interleave bits
        ! -----
        !
        ! interleaved_bits: kind=1, dimension=3*4=12
        ! element   byte representation
        ! 1-4:     |04-08-12|04-08-12|04-08-12|03-07-11|   e.g.: |11010000|00000000|00000001|00100000|
        ! 5-8:     |03-07-11|03-07-11|02-06-10|02-06-10|   e.g.: |00000000|00000000|00000000|00000000|
        ! 9-12:    |02-06-10|01-05-09|01-05-09|01-05-09|   e.g.: |00000000|00000000|00000000|00000000|
        !
        !   xx-yy-zz describes the bytes which were interleaved
        !
        ! ib_buffer: kind=1, dimension=3
        ! element   byte representation
        ! 1-3:     |---03---|---02---|---01---|
        !
        !
        ! ib_buffer = interleave(cb_buffer(9), cb_buffer(5), cb_buffer(1))
        ! iterleaved_bits(1:3) = ib_buffer
        !
        ! ...
        !
        ! ib_buffer = interleave(cb_buffer(12), cb_buffer(8), cb_buffer(4))
        ! interleaved_bits(10:12) = ib_buffer
        !
        integer(kind=INTERLEAVE_BITS_INTEGER_KIND) &
            ,dimension(TREE_DIMENSIONS * UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND) &
            :: interleaved_bits

        equivalence (interleaved_bits, interleaved_coordinate_binary)

        !
        ! doubel precision
        !   Number of decimal digits: ca. 16
        !   bit precision: 53
        !
        ! 3d example
        ! ----
        ! smalest cube
        !   edge length: 10**-16
        !   volume: 10**-48
        !
        ! largest cube
        !   edge length: 1
        !   volume: 1
        !
        ! number of smalest cubes in the biggest cube
        !   number = 1 / 10**-48 = 10**48
        !
        ! determination of the integer kind
        ! -----
        ! 32 bit word range: −(2**31) to 2**31 − 1
        ! 64 bit word range: −(2**63) to 2**63 − 1
        !
        !
        integer(kind=UINDEX_BYTES) :: interleaved_bits_dimension_0

        !   TREE_DIMENSIONS * UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND - UINDEX_BYTES /INTERLEAVE_BITS_INTEGER_KIND + 1
        ! = UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND * ( TREE_DIMENSIONS - 1) + 1
        equivalence (interleaved_bits(UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND * ( TREE_DIMENSIONS - 1) + 1), &
                     interleaved_bits_dimension_0)

        integer(kind=INTERLEAVE_BITS_INTEGER_KIND) &
            ,dimension(TREE_DIMENSIONS * UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND) &
            :: cb_buffer

        equivalence (coordinate_binary, cb_buffer)

        integer(kind=INTERLEAVE_BITS_INTEGER_KIND), dimension(TREE_DIMENSIONS) :: ob_buffer_DIM
        integer(kind=INTERLEAVE_BITS_INTEGER_KIND), dimension(TREE_DIMENSIONS) :: ib_buffer

        integer(kind=1) :: i
        integer(kind=1) :: ii

        coordinate_binary = x

        ! interleave bits
        do i=UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND, 1, -1 ! interleave column wise
            do ii=1, TREE_DIMENSIONS  ! get column entries
                ob_buffer_DIM(ii) = cb_buffer(i + (ii-1)*UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND)
            end do
            ib_buffer = lib_tree_hf_interleave_bits_use_lut(ob_buffer_DIM)

            ! e.g.    12                4       (4..1)        3
            ! ii = total length - (total columns - i + 1) * length(ib_buffer) + 1
            ! ii = TREE_DIMENSIONS * UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND - (UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND - i + 1) * TREE_DIMENSIONS + 1
            ! ii = TREE_DIMENSIONS * (UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND - UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND + i - 1) + 1
            ! ii = TREE_DIMENSIONS * (i - 1) + 1
            ii = int(TREE_DIMENSIONS * (i-1) + 1, 1)

            interleaved_bits(ii:ii+TREE_DIMENSIONS-1) = ib_buffer(:)
        end do

        rv = interleaved_coordinate_binary

    end function lib_tree_hf_interleave_bits_treeD_to_1D_use_lut

    ! deinterleavs the
    ! The kind of the integer is defined with INTERLEAVE_BITS_INTEGER_KIND.
    ! The dafault value is 1 (1 byte);
    !
    ! Argument
    ! ----
    !   x:
    !       counting system undependet value
    !
    !
    ! Example
    ! ----
    !  x1D: 0.10|01|00
    !  ------------------
    !   x1: 0.1 |0 |0
    !   x2: 0. 0| 1| 0
    !
    !   x1: 0.100   => 0.5   (base 10)
    !   x2: 0.010   => 0.25  (base 10)
    !
    function lib_tree_hf_deinterleave_bits_1D_to_treeD(x) result(rv)
        implicit none

        ! dummy
        integer(kind=UINDEX_BYTES), intent(in) :: x
        integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS) :: rv

        ! auxiliary
#if (_UINDEX_BYTES_ == 16)
        integer(kind=2) :: i
#else
        integer(kind=1) :: i
#endif
        integer(kind=1) :: bit_number
        integer(kind=1) :: target_element


        do i = 1, TREE_DIMENSIONS
            rv(i) = 0
        end do

        do i = 0, UINDEX_BYTES * NUMBER_OF_BITS_PER_BYTE - 1
            ! global bit number = i
            ! target_element = TREE_DIMENSIONS - mod( gloabel bit number, TREE_DIMENSIONS)
            ! target local bit number = int(global bit number / TREE_DIMENSIONS)
            target_element = int(i, 1)
            bit_number = int(target_element / TREE_DIMENSIONS, 1)
            target_element = int(TREE_DIMENSIONS - (mod(target_element, TREE_DIMENSIONS)), 1)
            if (btest(x, i)) then       ! first bit number = 0
                rv(target_element) = ibset(rv(target_element), bit_number)
            end if
        end do
    end function lib_tree_hf_deinterleave_bits_1D_to_treeD

    ! deinterleavs the
    ! The kind of the integer is defined with INTERLEAVE_BITS_INTEGER_KIND.
    ! The dafault value is 1 (1 byte);
    !
    ! Argument
    ! ----
    !   x:
    !       counting system undependet value
    !
    !
    ! Example
    ! ----
    !  x1D: 0.01|01|00
    !  ------------------
    !   x1: 0. 1| 0| 0
    !   x2: 0.0 |1 |0
    !
    !   x1: 0.100   => 0.5   (base 10)
    !   x2: 0.010   => 0.25  (base 10)
    !
    function lib_tree_hf_deinterleave_bits_1D_to_treeD_use_lut(x) result(rv)
        implicit none
        ! dummy
        integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS), intent(in) :: x
        integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS) :: rv

        ! auxiliar
        integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS) :: coordinate_binary
        integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS) :: deinterleaved_coordinate_binary

        ! Example
        ! ----
        !
        ! coordinate_binary: kind=4, dimension=3                                   2**-9    2**-1
        ! element   byte representation                          base(2)            v        v        base(10)
        ! 1:       |---04---|---03---|---02---|---01---|   e.g.: |00000000|00000000|00000000|11000000|  0.75
        ! 2:       |---04---|---03---|---02---|---01---|   e.g.: |00000000|00000000|00000000|10000000|  0.5
        ! 3:       |---04---|---03---|---02---|---01---|   e.g.: |00000000|00000000|10000000|00000001|  0.005859375
        !
        ! cb_buffer: kind=1, dimension=3*4=12
        ! element   byte representation                          base(2)                                base(10)
        ! 1-4:     |---01---|---02---|---03---|---04---|   e.g.: |00000000|00000000|00000000|11000000|  0.75
        ! 5-8:     |---05---|---06---|---07---|---08---|   e.g.: |00000000|00000000|00000000|10000000|  0.5
        ! 9-12:    |---09---|---10---|---11---|---12---|   e.g.: |00000000|00000000|10000000|00000001|  0.005859375
        !                                     <-- *
        !                    deinterleave element wise
        !
        ! equivalence (coordinate_binary, cb_buffer)
        !
        ! Deinterleave bits
        ! -----
        !
        ! deinterleaved_bits: kind=1, dimension=3*4=12
        ! element   byte representation
        ! 1-4:     |01-02-03|04-05-06|07-08-09|10-11-12|   e.g.: |00000000|10000000|00000000|00000000|
        ! 5-8:     |01-02-03|04-05-06|07-08-09|10-11-12|   e.g.: |00000000|10000000|00000000|00000000|
        ! 9-12:    |01-02-03|04-05-06|07-08-09|10-11-12|   e.g.: |00000000|00000000|00100000|00100001|
        !
        !   xx-yy-zz describes the bytes which were interleaved
        !
        ! ib_buffer: kind=1, dimension=3
        ! element   byte representation
        ! 1-3:     |---03---|---02---|---01---|
        !
        !
        ! dib_buffer = interleave(cb_buffer(9), cb_buffer(5), cb_buffer(1))
        ! deiterleaved_bits(12) = dib_buffer(1)
        ! deiterleaved_bits(8) = dib_buffer(2)
        ! deiterleaved_bits(4) = dib_buffer(3)
        !
        ! ...
        !
        ! dib_buffer = interleave(cb_buffer(12), cb_buffer(8), cb_buffer(4))
        ! deiterleaved_bits(1) = dib_buffer(1)
        ! deiterleaved_bits(5) = dib_buffer(2)
        ! deiterleaved_bits(9) = dib_buffer(3)
        !
        integer(kind=INTERLEAVE_BITS_INTEGER_KIND) &
            ,dimension(TREE_DIMENSIONS * UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND) &
            :: deinterleaved_bits

        equivalence (deinterleaved_bits, deinterleaved_coordinate_binary)

        !
        ! doubel precision
        !   Number of decimal digits: ca. 16
        !   bit precision: 53
        !
        ! 3d example
        ! ----
        ! smalest cube
        !   edge length: 10**-16
        !   volume: 10**-48
        !
        ! largest cube
        !   edge length: 1
        !   volume: 1
        !
        ! number of smalest cubes in the biggest cube
        !   number = 1 / 10**-48 = 10**48
        !
        ! determination of the integer kind
        ! -----
        ! 32 bit word range: −(2**31) to 2**31 − 1
        ! 64 bit word range: −(2**63) to 2**63 − 1
        !
        !
        integer(kind=UINDEX_BYTES) :: deinterleaved_bits_dimension_0

        !   TREE_DIMENSIONS * UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND - UINDEX_BYTES /INTERLEAVE_BITS_INTEGER_KIND + 1
        ! = UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND * ( TREE_DIMENSIONS - 1) + 1
        equivalence (deinterleaved_bits(UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND * ( TREE_DIMENSIONS - 1) + 1), &
                     deinterleaved_bits_dimension_0)

        integer(kind=INTERLEAVE_BITS_INTEGER_KIND) &
            ,dimension(TREE_DIMENSIONS * UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND) &
            :: cb_buffer

        equivalence (coordinate_binary, cb_buffer)

        integer(kind=INTERLEAVE_BITS_INTEGER_KIND), dimension(TREE_DIMENSIONS) :: ob_buffer_DIM
        integer(kind=INTERLEAVE_BITS_INTEGER_KIND), dimension(TREE_DIMENSIONS) :: dib_buffer        ! deinterleaved binary buffer

        integer(kind=1) :: i
        integer(kind=1) :: ii

        coordinate_binary = x

        ! deinterleave bits
        do i=UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND, 1, -1 ! interleave element wise
            ! get entries
            ! e.g.    12                4       (4..1)        3
            ! ii = total length - (total columns - i + 1) * length(ib_buffer) + 1
            ! ii = TREE_DIMENSIONS * UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND - (UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND - i + 1) * TREE_DIMENSIONS + 1
            ! ii = TREE_DIMENSIONS * (UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND - UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND + i - 1) + 1
            ! ii = TREE_DIMENSIONS * (i - 1) + 1
            ii = int(TREE_DIMENSIONS * (i-1) + 1, 1)
            ob_buffer_DIM(:) = cb_buffer(ii:ii+TREE_DIMENSIONS-1)

            dib_buffer = lib_tree_hf_deinterleave_bits_use_lut(ob_buffer_DIM)

            do ii=1, TREE_DIMENSIONS  ! set column entries
                deinterleaved_bits(i + (ii-1)*UINDEX_BYTES/INTERLEAVE_BITS_INTEGER_KIND) = dib_buffer(ii)
            end do
        end do

        rv = deinterleaved_coordinate_binary
    end function lib_tree_hf_deinterleave_bits_1D_to_treeD_use_lut

    ! reallocate a 1-dimensional array
    !
    ! Arguments
    ! ----
    !   a: 1-dimensional uindex array
    !       original array, will be replaced with the resized array
    !   n: positiv integer
    !       number of additional array elements
    !
    subroutine lib_tree_hf_reallocate_1d_uindex(a, n)
        use lib_tree_type
        implicit none
        type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: a
        type(lib_tree_universal_index), dimension(:), allocatable :: temp
        integer,intent(in) :: n
        integer :: ni_old

        if ( allocated(a) ) then
            ni_old = size(a)

            allocate(temp(ni_old+n))

            temp(1:ni_old) = a

            call move_alloc(temp,a)

        else
            allocate(a(n))
        end if
    end subroutine lib_tree_hf_reallocate_1d_uindex

    ! reallocates the data element list with additional n elements
    !
    ! Arguments
    ! ----
    !   a: 1-dimensional lib_tree_data_element array
    !       original array, will be replaced with the resized array
    !   n: pos. integer
    !       number of additional array elements
    !
    ! copy of toolbox.reallocate_1d
    subroutine lib_tree_hf_reallocate_1d_data_element_list(a,n)
        implicit none
        type(lib_tree_data_element),dimension(:),allocatable,intent(inout) :: a
        type(lib_tree_data_element),dimension(:),allocatable :: temp
        integer,intent(in) :: n
        integer :: ni_old

        if ( allocated(a) ) then
            ni_old = size(a)

            allocate(temp(ni_old+n))

            temp(1:ni_old) = a

            call move_alloc(temp,a)
        else
            allocate(a(n))
        end if

    end subroutine lib_tree_hf_reallocate_1d_data_element_list

    subroutine lib_tree_hf_concatenate_1d_data_element_list(a, b)
        implicit none
        ! dummy
        type(lib_tree_data_element),dimension(:),allocatable,intent(inout) :: a
        type(lib_tree_data_element),dimension(:),allocatable,intent(in) :: b

        ! auxiliaray
        integer :: size_a_org

        if ( allocated(a) ) then
            size_a_org = size(a)
        else
            size_a_org = 0
        end if
        call reallocate(a, size(b))

        a(size_a_org+1:) = b

    end subroutine lib_tree_hf_concatenate_1d_data_element_list

    subroutine lib_tree_hf_concatenate_1d_data_element_list_single(a, b)
        implicit none
        ! dummy
        type(lib_tree_data_element), dimension(:), allocatable, intent(inout) :: a
        type(lib_tree_data_element), intent(in) :: b

        ! auxiliaray
        type(lib_tree_data_element),dimension(:),allocatable :: buffer_b

        allocate(buffer_b(1))
        buffer_b(1) = b

        call concatenate(a, buffer_b)

    end subroutine lib_tree_hf_concatenate_1d_data_element_list_single

    ! ----------------- test functions -----------------
    function lib_tree_hf_test_functions() result(error_counter)
        implicit none

        integer :: error_counter

        error_counter = 0

        if (.not. test_lib_tree_hf_get_universal_index()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_get_universal_index_l_max()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_get_parent()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_get_children_all()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_get_centre_of_box()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_get_coordinate_binary_number_xD()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_get_neighbour_all_xD()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_get_neighbour_all_xD_2()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_interleave_bits()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_interleave_bits_2()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_interleave_bits_use_lut()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_interleave_bits_use_lut_2()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_deinterleave_bits()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_deinterleave_bits_2()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_deinterleave_bits_use_lut()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_deinterleave_bits_use_lut_2()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_interleave_bits_treeD_to_1D()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_interleave_bits_treeD_to_1D_use_lut()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_deinterleave_bits_1D_to_treeD()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_hf_deinterleave_bits_1D_to_treeD_use_lut()) then
            error_counter = error_counter + 1
        end if


        print *, "-------------lib_tree_hf_test_functions-------------"
        if (error_counter == 0) then
            print *, "lib_tree_hf_test_functions tests: OK"
        else
            print *, error_counter,"lib_tree_hf_test_functions test(s) FAILED"
        end if
        print *, "----------------------------------------------------"

        contains

        function test_lib_tree_hf_get_universal_index() result(rv)
            implicit none

            ! dummy
            logical :: rv

            type(lib_tree_spatial_point) :: point
            type(lib_tree_universal_index) :: universal_index

            integer(kind=1) :: l
            type(lib_tree_universal_index) :: universal_index_ground_trouth

            l = 1
            universal_index_ground_trouth%l = l

            point%x(1) = 0.75
            point%x(2) = 0.5
#if (_FMM_DIMENSION_ == 2)
            universal_index_ground_trouth%n = 3
#elif (_FMM_DIMENSION_ == 3)
            point%x(3) = 2.0**(-9.0) + 2.0**(-8)

            universal_index_ground_trouth%n = 6
#else
            print *, "test_lib_tree_hf_get_universal_index: Dimension not defines: ", FMM_DIMENSION
#endif

            universal_index = lib_tree_hf_get_universal_index(point, l)

            if (universal_index%n == universal_index_ground_trouth%n) then
                rv = .true.
                print *, "test_lib_tree_hf_get_universal_index: ", "OK"
            else
                rv = .false.
                print *, "test_lib_tree_hf_get_universal_index: ", "FAILED"
            end if
        end function test_lib_tree_hf_get_universal_index

        function test_lib_tree_hf_get_universal_index_l_max() result(rv)
            implicit none

            ! dummy
            logical :: rv

            type(lib_tree_spatial_point) :: point
            type(lib_tree_universal_index) :: universal_index

            integer(kind=1) :: l
            type(lib_tree_universal_index) :: universal_index_ground_trouth

#if (_FMM_DIMENSION_ == 2)
            l = int(COORDINATE_BINARY_BYTES * 1.0*NUMBER_OF_BITS_PER_BYTE / TREE_DIMENSIONS)
            universal_index_ground_trouth%l = l
#if (_SPATIAL_POINT_IS_DOUBLE_ == 0)
            point%x(1) = 1.0-2.0**(-l)
            point%x(2) = 1.0-2.0**(-l)
#if (_UINDEX_BYTES_ == 4)
            universal_index_ground_trouth%n = -1            ! 2**(dl)-1:
#elif (_UINDEX_BYTES_ == 8)
            universal_index_ground_trouth%n = X'FFFFFFFF'   ! 2**(dl)-1: 4294967295
#elif (_UINDEX_BYTES_ == 16)
            universal_index_ground_trouth%n = X'FFFFFFFFFFFFFFFF'   ! 2**(dl)-1:
#endif

#elif (_SPATIAL_POINT_IS_DOUBLE_ == 1)
            point%x(1) = 1.0d+0-2.0d+0**(-l)
            point%x(2) = 1.0d+0-2.0d+0**(-l)
            universal_index_ground_trouth%n = -1
#if (_UINDEX_BYTES_ == 4)
            universal_index_ground_trouth%n = -1            ! 2**(dl)-1:
#elif (_UINDEX_BYTES_ == 8)
            universal_index_ground_trouth%n = -1   ! 2**(dl)-1:
#elif (_UINDEX_BYTES_ == 16)
            universal_index_ground_trouth%n = X'FFFFFFFFFFFFFFFF'   ! 2**(dl)-1:
#endif

#endif

#elif (_FMM_DIMENSION_ == 3)
            l = int(COORDINATE_BINARY_BYTES * 1.0*NUMBER_OF_BITS_PER_BYTE / TREE_DIMENSIONS)
            universal_index_ground_trouth%l = l
#if (_SPATIAL_POINT_IS_DOUBLE_ == 0)
            point%x(1) = 1.0-2.0**(-l-1)
            point%x(2) = 1.0-2.0**(-l-1)
            point%x(3) = 1.0-2.0**(-l-1)
            universal_index_ground_trouth%n = X'3FFFFFFF' ! 2**(dl)-1
#elif (_SPATIAL_POINT_IS_DOUBLE_ == 1)
            point%x(1) = 1.0d+0-2.0d+0**(-l)
            point%x(2) = 1.0d+0-2.0d+0**(-l)
            point%x(3) = 1.0d+0-2.0d+0**(-l)
            universal_index_ground_trouth%n = X'7FFFFFFFFFFFFFFF' ! 2**(dl)-1   ! todo: 9223372036854775807
#endif


#else
            print *, "test_lib_tree_hf_get_universal_index_l_max: Dimension not defines: ", FMM_DIMENSION
#endif

            universal_index = lib_tree_hf_get_universal_index(point, l)

            if (universal_index%n == universal_index_ground_trouth%n) then
                rv = .true.
                print *, "test_lib_tree_hf_get_universal_index_l_max: ", "OK"
            else
                rv = .false.
                print *, "test_lib_tree_hf_get_universal_index_l_max: ", "FAILED"
                print *, "  universal_index%l = ", universal_index%l
                print *, "  universal_index%n = ", universal_index%n
            end if
        end function test_lib_tree_hf_get_universal_index_l_max

        function test_lib_tree_hf_get_parent() result (rv)
            implicit none

            ! dummy
            logical :: rv

            integer(kind=UINDEX_BYTES) :: n
            integer(kind=UINDEX_BYTES) :: n_parent

            integer(kind=1) :: l
            integer(kind=UINDEX_BYTES) :: n_parent_ground_trouth

            l = 1
#if (_FMM_DIMENSION_ == 2)
            n = 2
            n_parent_ground_trouth = 0
#elif (_FMM_DIMENSION_ == 3)
            n = 6
            n_parent_ground_trouth = 0
#else
            print *, "test_lib_tree_hf_get_parent: Dimension not defines: ", FMM_DIMENSION
#endif

            n_parent = lib_tree_hf_get_parent(n)

            if (n_parent == n_parent_ground_trouth) then
                print *, "test_lib_tree_hf_get_parent: ", "ok"
                rv = .true.
            else
                print *, "test_lib_tree_hf_get_parent: ", "FAILED"
                rv = .false.
            end if

        end function test_lib_tree_hf_get_parent

        function test_lib_tree_hf_get_children_all() result (rv)
            implicit none

            ! dummy
            logical :: rv

            integer(kind=UINDEX_BYTES) :: n
            integer(kind=UINDEX_BYTES), dimension(2**TREE_DIMENSIONS) :: children_n
            integer(kind=UINDEX_BYTES), dimension(2**TREE_DIMENSIONS) :: children_n_ground_truth

            integer(kind=1) :: i

            n = 1
#if (_FMM_DIMENSION_ == 2)
            children_n_ground_truth(1) = 4
            children_n_ground_truth(2) = 5
            children_n_ground_truth(3) = 6
            children_n_ground_truth(4) = 7
#elif (_FMM_DIMENSION_ == 3)
            children_n_ground_truth(1) = 8
            children_n_ground_truth(2) = 9
            children_n_ground_truth(3) = 10
            children_n_ground_truth(4) = 11
            children_n_ground_truth(5) = 12
            children_n_ground_truth(6) = 13
            children_n_ground_truth(7) = 14
            children_n_ground_truth(8) = 15
#else
            print *, "test_lib_tree_hf_get_children_all: Dimension not defines: ", FMM_DIMENSION
#endif
            children_n = lib_tree_hf_get_children_all(n)

            print *, "test_lib_tree_hf_get_children_all: "
            rv = .true.
            do i=1, 2**TREE_DIMENSIONS
                if (children_n(i) == children_n_ground_truth(i)) then
                    print *, " ", i, ": ok"

                else
                    print *, " ", i, ": FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_tree_hf_get_children_all

        function test_lib_tree_hf_get_centre_of_box() result (rv)
            implicit none

            ! dummy
            logical :: rv

            integer(kind=UINDEX_BYTES) :: n
            integer(kind=1) :: l
            type(lib_tree_spatial_point) :: point

            type(lib_tree_spatial_point) :: point_ground_trouth

            integer :: i

            l=1
            n=1
#if (_FMM_DIMENSION_ == 2)
            point_ground_trouth%x(1) = 0.25
            point_ground_trouth%x(2) = 0.75
#elif (_FMM_DIMENSION_ == 3)
            point_ground_trouth%x(1) = 0.25
            point_ground_trouth%x(2) = 0.25
            point_ground_trouth%x(3) = 0.75
#else
            print *, "test_lib_tree_hf_get_centre_of_box: Dimension not defines: ", FMM_DIMENSION
#endif
            point = lib_tree_hf_get_centre_of_box(n,l)

            print *, "test_lib_tree_hf_get_centre_of_box:"
            rv = .true.
            do i=1,TREE_DIMENSIONS
                if (point%x(i) == point_ground_trouth%x(i)) then
                    print *, "  dim: ", i,": ", "ok"
                else
                    print *, "  dim: ", i,": ", "FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_tree_hf_get_centre_of_box

        function test_lib_tree_hf_get_neighbour_all_xD() result (rv)
            implicit none
            !dummy
            logical :: rv

            ! parameter
            integer(kind=1), parameter :: number_of_neighbours = 3**TREE_DIMENSIONS-1

            ! auxiliary
            integer(kind=UINDEX_BYTES) :: k
            integer(kind=1) :: l
            integer(kind=UINDEX_BYTES) :: n
            integer(kind=UINDEX_BYTES), dimension(number_of_neighbours) :: neighbour_all
            integer(kind=UINDEX_BYTES), dimension(number_of_neighbours) :: neighbour_all_ground_truth

            integer(kind=1) :: i

#if (_FMM_DIMENSION_ == 2)
            k = 1
            l = 2
            n = 3
            ! ------------------
            ! |   |   ||   |   |
            ! | 5 | 7 || 13| 15|
            ! |---1--------3----
            ! |   |   ||   |   |
            ! | 4 | 6 || 12| 14|
            ! ==================
            ! |   | X ||   |   |
            ! | 1 | 3 || 9 | 11|
            ! |---0--------2----
            ! |   |   ||   |   |
            ! | 0 | 2 || 8 | 10|
            ! ------------------

            neighbour_all_ground_truth = (/ 0,1,4, &
                                            2, 6, &
                                            8,9,12 /)

#elif (_FMM_DIMENSION_ == 3)
            k = 1
            l = 2
            n = 35
            !      __________________
            !     /   3    //   7   /|
            !    /        //       / |
            !   /--------//-------/  |
            !  /    1   //   5   /|  |
            ! /        //       / | 7|
            ! ------------------  |  |
            ! |   |   ||   |   | 5| /|
            ! | 9 | 13|| 41| 45|  |//|
            ! |---1--------5----  // |
            ! |   |   ||   |   | //  |
            ! | 8 | 12|| 40| 44|//| 6|
            ! ==================/ |  |
            ! |   |   ||   |   |  |  /
            ! | 1 | 5 || 33| 37| 4| /
            ! |---0--------4----  |/
            ! |   |   ||   |   |  /
            ! | 0 | 4 || 32| 36| /
            ! ------------------/

            neighbour_all_ground_truth = (/ 4,5,12, &
                                            6,7,14, &
                                            20,21,28, &
                                            32,33,40, &
                                            34, 42, &
                                            48,49,56, &
                                            36,37,44, &
                                            38,39,46, &
                                            52,53,60 /)
#else
            print *, "test_lib_tree_hf_get_neighbour_all_xD: Dimension not defines: ", FMM_DIMENSION
#endif

            neighbour_all = lib_tree_hf_get_neighbour_all_xD(k, n, l)

            print *, "test_lib_tree_hf_get_neighbour_all_xD:"
            rv = .true.
            do i=1,number_of_neighbours
                if (neighbour_all_ground_truth(i) == neighbour_all(i)) then
                    print *, "  dim ",i,": ", "ok"
                else
                    print *, "  dim ",i,": ", "FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_tree_hf_get_neighbour_all_xD

        function test_lib_tree_hf_get_neighbour_all_xD_2() result (rv)
            implicit none
            !dummy
            logical :: rv

            ! parameter
            integer(kind=1), parameter :: number_of_neighbours = 3**TREE_DIMENSIONS-1

            ! auxiliary
            integer(kind=UINDEX_BYTES) :: k
            integer(kind=1) :: l
            integer(kind=UINDEX_BYTES) :: n
            integer(kind=UINDEX_BYTES), dimension(number_of_neighbours) :: neighbour_all
            integer(kind=UINDEX_BYTES), dimension(number_of_neighbours) :: neighbour_all_ground_truth

            integer(kind=1) :: i

#if (_FMM_DIMENSION_ == 2)
            k = 1
            l = 2
            n = 11
            ! ------------------
            ! |   |   ||   |   |
            ! | 5 | 7 || 13| 15|
            ! |---1--------3----
            ! |   |   ||   |   |
            ! | 4 | 6 || 12| 14|
            ! ==================
            ! |   |   ||   | X |
            ! | 1 | 3 || 9 | 11|
            ! |---0--------2----
            ! |   |   ||   |   |
            ! | 0 | 2 || 8 | 10|
            ! ------------------

            neighbour_all_ground_truth = (/ 8,9,12, &
                                            10, 14, &
                                            -1,-1,-1 /)

#elif (_FMM_DIMENSION_ == 3)
            k = 1
            l = 1
            n = 5
            !      __________________
            !     /   3    //   7   /|
            !    /        //       / |
            !   /--------//-------/  |
            !  /    1   //   5   /|  |
            ! /        //       / | 7|
            ! ------------------  |  |
            ! |   |   ||   |   | 5| /|
            ! | 9 | 13|| 41| 45|  |//|
            ! |---1--------5----  // |
            ! |   |   ||   |   | //  |
            ! | 8 | 12|| 40| 44|//| 6|
            ! ==================/ |  |
            ! |   |   ||   |   |  |  /
            ! | 1 | 5 || 33| 37| 4| /
            ! |---0--------4----  |/
            ! |   |   ||   |   |  /
            ! | 0 | 4 || 32| 36| /
            ! ------------------/

            neighbour_all_ground_truth = (/ 0,1,2, &
                                            3,4,6, &
                                            7,-1,-1, &
                                            -1,-1,-1, &
                                            -1,-1,-1, &
                                            -1,-1,-1, &
                                            -1,-1,-1, &
                                            -1,-1,-1, &
                                            -1,-1 /)
#else
            print *, "test_lib_tree_hf_get_neighbour_all_xD_2: Dimension not defines: ", FMM_DIMENSION
#endif

            neighbour_all = lib_tree_hf_get_neighbour_all_xD(k, n, l)

            print *, "test_lib_tree_hf_get_neighbour_all_xD_2:"
            rv = .true.
            do i=1,number_of_neighbours
                if (neighbour_all_ground_truth(i) == neighbour_all(i)) then
                    print *, "  dim ",i,": ", "ok"
                else
                    print *, "  dim ",i,": ", "FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_tree_hf_get_neighbour_all_xD_2

        function test_lib_tree_hf_get_coordinate_binary_number_xD() result (rv)
            implicit none

            ! dummy
            logical :: rv

#if(_SPATIAL_POINT_IS_DOUBLE_ == 1)
            double precision, dimension(TREE_DIMENSIONS) :: f
#elif(_SPATIAL_POINT_IS_DOUBLE_ == 0)
            real, dimension(TREE_DIMENSIONS) :: f
#endif
            integer(kind=COORDINATE_BINARY_BYTES), dimension(TREE_DIMENSIONS) :: coordinate_binary_xD
            integer(kind=COORDINATE_BINARY_BYTES), dimension(TREE_DIMENSIONS) :: coordinate_binary_xD_ground_trouth

            integer :: i

#if (_FMM_DIMENSION_ == 2)
            f(1) = 0.25
            f(2) = 0.375
            if (COORDINATE_BINARY_BYTES == 4) then
                coordinate_binary_xD_ground_trouth(1) = 1073741824      ! |0100 0000|0000 0000| ... |0000 0000| byte 3-0
                coordinate_binary_xD_ground_trouth(2) = 1610612736      ! |0110 0000|0000 0000| ... |0000 0000| byte 3-0
            else if (COORDINATE_BINARY_BYTES == 8) then
                coordinate_binary_xD_ground_trouth(1) = 2**62           ! |0100 0000|0000 0000| ... |0000 0000| byte 7-0
                coordinate_binary_xD_ground_trouth(2) = 2**62 + 2**61   ! |0110 0000|0000 0000| ... |0000 0000| byte 7-0
            end if
#elif (_FMM_DIMENSION_ == 3)
            f(1) = 0.25
            f(2) = 0.375
            f(3) = 0.4375
            if (COORDINATE_BINARY_BYTES == 4) then
                coordinate_binary_xD_ground_trouth(1) = 1073741824          ! |0100 0000|0000 0000| ... |0000 0000| byte 3-0
                coordinate_binary_xD_ground_trouth(2) = 1610612736          ! |0110 0000|0000 0000| ... |0000 0000| byte 3-0
                coordinate_binary_xD_ground_trouth(3) = 1879048192          ! |0111 0000|0000 0000| ... |0000 0000| byte 3-0
            else if (COORDINATE_BINARY_BYTES == 8) then
                coordinate_binary_xD_ground_trouth(1) = 2**62               ! |0100 0000|0000 0000| ... |0000 0000| byte 7-0
                coordinate_binary_xD_ground_trouth(2) = 2**62+2**61         ! |0110 0000|0000 0000| ... |0000 0000| byte 7-0
                coordinate_binary_xD_ground_trouth(3) = 2**62+2**61+2**60   ! |0111 0000|0000 0000| ... |0000 0000| byte 7-0
            end if
#else
            print *, "test_lib_tree_hf_get_coordinate_binary_number_xD: Dimension not defines: ", FMM_DIMENSION
#endif

            coordinate_binary_xD = lib_tree_hf_get_coordinate_binary_number_xD(f)

            print *, "test_lib_tree_hf_get_coordinate_binary_number_xD:"
            rv = .true.
            do i=1,TREE_DIMENSIONS
                if (coordinate_binary_xD(i) == coordinate_binary_xD_ground_trouth(i)) then
                    print *, "  dim ",i,": ", "ok"
                else
                    print *, "  dim ",i,": ", "FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_tree_hf_get_coordinate_binary_number_xD

        function test_lib_tree_hf_interleave_bits() result (rv)
            implicit none

            ! dummy
            logical :: rv

            ! parameter
            integer(kind=1), parameter :: x_kind = INTERLEAVE_BITS_INTEGER_KIND

            integer(kind=x_kind), dimension(TREE_DIMENSIONS) :: x
            integer(kind=x_kind), dimension(size(x)) :: interleaved_bits

            integer(kind=x_kind), dimension(size(x)) :: interleaved_bits_ground_trouth

            integer :: i

            x(1) = 4  ! 0100
            x(2) = 1  ! 0001
#if (_FMM_DIMENSION_ == 2)
            interleaved_bits_ground_trouth(1) = 2**0 + 2**5 !           |0010 0001|
            interleaved_bits_ground_trouth(2) = 0           ! |0000 0000|
#elif (_FMM_DIMENSION_ == 3)
            x(3) = 2  ! 0010
            interleaved_bits_ground_trouth(1) = 2**1 + 2**3        !                     |0000 1010|
            interleaved_bits_ground_trouth(2) = 2**0               !           |0000 0001|
            interleaved_bits_ground_trouth(3) = 0                  ! |0000 0000|
#else
            print *, "test_lib_tree_hf_interleave_bits: Dimension not defines: ", FMM_DIMENSION
#endif

            interleaved_bits = lib_tree_hf_interleave_bits(x)

            print *, "test_lib_tree_hf_interleave_bits:"
            rv = .true.
            do i=1, TREE_DIMENSIONS
                if (interleaved_bits(i) == interleaved_bits_ground_trouth(i)) then
                    print *, "  dim ",i,": ", "ok"
                else
                    print *, "  dim ",i,": ", "FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_tree_hf_interleave_bits

        function test_lib_tree_hf_interleave_bits_2() result (rv)
            implicit none

            ! dummy
            logical :: rv

            ! parameter
            integer(kind=1), parameter :: x_kind = INTERLEAVE_BITS_INTEGER_KIND

            integer(kind=x_kind), dimension(TREE_DIMENSIONS) :: x
            integer(kind=x_kind), dimension(size(x)) :: interleaved_bits

            integer(kind=x_kind), dimension(size(x)) :: interleaved_bits_ground_trouth

            integer :: i

            x(1) = 20  ! 0001 0100
            x(2) = 65  ! 0100 0001
#if (_FMM_DIMENSION_ == 2)
            interleaved_bits_ground_trouth(1) = 2**0 + 2**5  !           |0010 0001|
            interleaved_bits_ground_trouth(2) = 2**1 + 2**4  ! |0001 0010|
#elif (_FMM_DIMENSION_ == 3)
            x(3) = 40  ! 0010 1000
            interleaved_bits_ground_trouth(1) = 2**1                       !                     |0000 0010|
            interleaved_bits_ground_trouth(2) = 2**0 + 2**1 + 2**6 - 2**7  !           |1100 0011|
            interleaved_bits_ground_trouth(3) = 2**3                       ! |0000 1000|
#else
            print *, "test_lib_tree_hf_interleave_bits_2: Dimension not defines: ", FMM_DIMENSION
#endif

            interleaved_bits = lib_tree_hf_interleave_bits(x)

            print *, "test_lib_tree_hf_interleave_bits_2:"
            rv = .true.
            do i=1, TREE_DIMENSIONS
                if (interleaved_bits(i) == interleaved_bits_ground_trouth(i)) then
                     print *, "  dim ",i,": ", "ok"
                else
                    print *, "  dim ",i,": ", "FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_tree_hf_interleave_bits_2

        function test_lib_tree_hf_interleave_bits_use_lut() result (rv)
            implicit none

            ! dummy
            logical :: rv

            ! parameter
            integer(kind=1), parameter :: x_kind = INTERLEAVE_BITS_INTEGER_KIND

            integer(kind=x_kind), dimension(TREE_DIMENSIONS) :: x
            integer(kind=x_kind), dimension(size(x)) :: interleaved_bits

            integer(kind=x_kind), dimension(size(x)) :: interleaved_bits_ground_trouth

            integer :: i

            x(1) = 4  ! 0100
            x(2) = 1  ! 0001
#if (_FMM_DIMENSION_ == 2)
            interleaved_bits_ground_trouth(1) = 2**0 + 2**5 !           |0010 0001|
            interleaved_bits_ground_trouth(2) = 0           ! |0000 0000|
#elif (_FMM_DIMENSION_ == 3)
            x(3) = 2  ! 0010
            interleaved_bits_ground_trouth(1) = 2**1 + 2**3        !                     |0000 1010|
            interleaved_bits_ground_trouth(2) = 2**0               !           |0000 0001|
            interleaved_bits_ground_trouth(3) = 0                  ! |0000 0000|
#else
            print *, "test_lib_tree_hf_interleave_bits_use_lut: Dimension not defines: ", FMM_DIMENSION
#endif

            interleaved_bits = lib_tree_hf_interleave_bits_use_lut(x)

            print *, "test_lib_tree_hf_interleave_bits_use_lut:"
            rv = .true.
            do i=1, TREE_DIMENSIONS
                if (interleaved_bits(i) == interleaved_bits_ground_trouth(i)) then
                    print *, "  dim ",i,": ", "ok"
                else
                    print *, "  dim ",i,": ", "FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_tree_hf_interleave_bits_use_lut

        function test_lib_tree_hf_interleave_bits_use_lut_2() result (rv)
            implicit none

            ! dummy
            logical :: rv

            ! parameter
            integer(kind=1), parameter :: x_kind = INTERLEAVE_BITS_INTEGER_KIND

            integer(kind=x_kind), dimension(TREE_DIMENSIONS) :: x
            integer(kind=x_kind), dimension(size(x)) :: interleaved_bits

            integer(kind=x_kind), dimension(size(x)) :: interleaved_bits_ground_trouth

            integer :: i

            x(1) = 20  ! 0001 0100
            x(2) = 65  ! 0100 0001
#if (_FMM_DIMENSION_ == 2)
            interleaved_bits_ground_trouth(1) = 2**0 + 2**5  !           |0010 0001|
            interleaved_bits_ground_trouth(2) = 2**1 + 2**4  ! |0001 0010|
#elif (_FMM_DIMENSION_ == 3)
            x(3) = 40  ! 0010 1000
            interleaved_bits_ground_trouth(1) = 2**1                       !                     |0000 0010|
            interleaved_bits_ground_trouth(2) = 2**0 + 2**1 + 2**6 - 2**7  !           |1100 0011|
            interleaved_bits_ground_trouth(3) = 2**3                       ! |0000 1000|
#else
            print *, "test_lib_tree_hf_interleave_bits_use_lut_2: Dimension not defines: ", FMM_DIMENSION
#endif

            interleaved_bits = lib_tree_hf_interleave_bits_use_lut(x)

            print *, "test_lib_tree_hf_interleave_bits_use_lut_2:"
            rv = .true.
            do i=1, TREE_DIMENSIONS
                if (interleaved_bits(i) == interleaved_bits_ground_trouth(i)) then
                     print *, "  dim ",i,": ", "ok"
                else
                    print *, "  dim ",i,": ", "FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_tree_hf_interleave_bits_use_lut_2

        function test_lib_tree_hf_deinterleave_bits() result (rv)
            implicit none

            ! dummy
            logical :: rv

            ! parameter
            integer(kind=1), parameter :: x_kind = INTERLEAVE_BITS_INTEGER_KIND

            ! dummy
            integer(kind=x_kind), dimension(TREE_DIMENSIONS) :: x
            integer(kind=x_kind), dimension(size(x)) :: deinterleaved_bits

            integer(kind=x_kind), dimension(size(x)) :: deinterleaved_bits_ground_trouth

            integer :: i

            deinterleaved_bits_ground_trouth(1) = 4  ! 0100
            deinterleaved_bits_ground_trouth(2) = 1  ! 0001

#if (_FMM_DIMENSION_ == 2)
            x(1) = 2**0 + 2**5 !           |0010 0001|
            x(2) = 0           ! |0000 0000|
#elif (_FMM_DIMENSION_ == 3)
            deinterleaved_bits_ground_trouth(3) = 2  ! 0010
            x(1) = 2**1 + 2**3        !                     |0000 1010|
            x(2) = 2**0               !           |0000 0001|
            x(3) = 0                  ! |0000 0000|
#else
            print *, "test_lib_tree_hf_deinterleave_bits: Dimension not defines: ", FMM_DIMENSION
#endif

            deinterleaved_bits = lib_tree_hf_deinterleave_bits(x)

            print *, "test_lib_tree_hf_deinterleave_bits:"
            rv = .true.
            do i=1, TREE_DIMENSIONS
                if (deinterleaved_bits(i) == deinterleaved_bits_ground_trouth(i)) then
                    print *, "  dim ",i,": ", "ok"
                else
                    print *, "  dim ",i,": ", "FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_tree_hf_deinterleave_bits

        function test_lib_tree_hf_deinterleave_bits_2() result (rv)
            implicit none

            ! dummy
            logical :: rv

            ! parameter
            integer(kind=1), parameter :: x_kind = INTERLEAVE_BITS_INTEGER_KIND

            ! dummy
            integer(kind=x_kind), dimension(TREE_DIMENSIONS) :: x
            integer(kind=x_kind), dimension(size(x)) :: deinterleaved_bits

            integer(kind=x_kind), dimension(size(x)) :: deinterleaved_bits_ground_trouth

            integer :: i

            deinterleaved_bits_ground_trouth(1) = 2**2 + 2**6  ! 0100 0100
            deinterleaved_bits_ground_trouth(2) = 2**0 + 2**5  ! 0010 0001

#if (_FMM_DIMENSION_ == 2)
            x(1) = 2**0 + 2**5 !           |0010 0001|
            x(2) = 2**2 + 2**5 ! |0010 0100|
#elif (_FMM_DIMENSION_ == 3)
            deinterleaved_bits_ground_trouth(3) = 2**1 + 2**4  ! 0001 0010
            x(1) = 2**1 + 2**3        !                     |0000 1010|
            x(2) = 2**0 + 2**4        !           |0001 0001|
            x(3) = 2**0 + 2**4        ! |0001 0001|
#else
            print *, "test_lib_tree_hf_deinterleave_bits_2: Dimension not defines: ", FMM_DIMENSION
#endif

            deinterleaved_bits = lib_tree_hf_deinterleave_bits(x)

            print *, "test_lib_tree_hf_deinterleave_bits_2:"
            rv = .true.
            do i=1, TREE_DIMENSIONS
                if (deinterleaved_bits(i) == deinterleaved_bits_ground_trouth(i)) then
                    print *, "  dim ",i,": ", "ok"
                else
                    print *, "  dim ",i,": ", "FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_tree_hf_deinterleave_bits_2

        function test_lib_tree_hf_deinterleave_bits_use_lut() result (rv)
            implicit none

            ! dummy
            logical :: rv

            ! parameter
            integer(kind=1), parameter :: x_kind = INTERLEAVE_BITS_INTEGER_KIND

            ! dummy
            integer(kind=x_kind), dimension(TREE_DIMENSIONS) :: x
            integer(kind=x_kind), dimension(size(x)) :: deinterleaved_bits

            integer(kind=x_kind), dimension(size(x)) :: deinterleaved_bits_ground_trouth

            integer :: i

            deinterleaved_bits_ground_trouth(1) = 4  ! 0100
            deinterleaved_bits_ground_trouth(2) = 1  ! 0001

#if (_FMM_DIMENSION_ == 2)
            x(1) = 2**0 + 2**5 !           |0010 0001|
            x(2) = 0           ! |0000 0000|
#elif (_FMM_DIMENSION_ == 3)
            deinterleaved_bits_ground_trouth(3) = 2  ! 0010
            x(1) = 2**1 + 2**3        !                     |0000 1010|
            x(2) = 2**0               !           |0000 0001|
            x(3) = 0                  ! |0000 0000|
#else
            print *, "test_lib_tree_hf_deinterleave_bits_use_lut: Dimension not defines: ", FMM_DIMENSION
#endif

            deinterleaved_bits = lib_tree_hf_deinterleave_bits_use_lut(x)

            print *, "test_lib_tree_hf_deinterleave_bits_use_lut:"
            rv = .true.
            do i=1, TREE_DIMENSIONS
                if (deinterleaved_bits(i) == deinterleaved_bits_ground_trouth(i)) then
                    print *, "  dim ",i,": ", "ok"
                else
                    print *, "  dim ",i,": ", "FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_tree_hf_deinterleave_bits_use_lut

        function test_lib_tree_hf_deinterleave_bits_use_lut_2() result (rv)
            implicit none

            ! dummy
            logical :: rv

            ! parameter
            integer(kind=1), parameter :: x_kind = INTERLEAVE_BITS_INTEGER_KIND

            ! dummy
            integer(kind=x_kind), dimension(TREE_DIMENSIONS) :: x
            integer(kind=x_kind), dimension(size(x)) :: deinterleaved_bits

            integer(kind=x_kind), dimension(size(x)) :: deinterleaved_bits_ground_trouth

            integer :: i

            deinterleaved_bits_ground_trouth(1) = 2**2 + 2**6  ! 0100 0100
            deinterleaved_bits_ground_trouth(2) = 2**0 + 2**5  ! 0010 0001

#if (_FMM_DIMENSION_ == 2)
            x(1) = 2**0 + 2**5 !           |0010 0001|
            x(2) = 2**2 + 2**5 ! |0010 0100|
#elif (_FMM_DIMENSION_ == 3)
            deinterleaved_bits_ground_trouth(3) = 2**1 + 2**4  ! 0001 0010
            x(1) = 2**1 + 2**3        !                     |0000 1010|
            x(2) = 2**0 + 2**4        !           |0001 0001|
            x(3) = 2**0 + 2**4        ! |0001 0001|
#else
            print *, "test_lib_tree_hf_deinterleave_bits_use_lut_2: Dimension not defines: ", FMM_DIMENSION
#endif

            deinterleaved_bits = lib_tree_hf_deinterleave_bits_use_lut(x)

            print *, "test_lib_tree_hf_deinterleave_bits_use_lut_2:"
            rv = .true.
            do i=1, TREE_DIMENSIONS
                if (deinterleaved_bits(i) == deinterleaved_bits_ground_trouth(i)) then
                    print *, "  dim ",i,": ", "ok"
                else
                    print *, "  dim ",i,": ", "FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_tree_hf_deinterleave_bits_use_lut_2

        function test_lib_tree_hf_interleave_bits_treeD_to_1D() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer(kind=1), parameter :: x_kind = UINDEX_BYTES

            ! dummy
            integer(kind=x_kind), dimension(TREE_DIMENSIONS) :: x
            integer(kind=x_kind) :: interleaved_bits

            integer(kind=x_kind) :: interleaved_bits_ground_trouth

#if (_FMM_DIMENSION_ == 2)
            x(1) = 2**0 + 2**5                          ! |0000 0000|0010 0001|
            x(2) = 2**2 + 2**4                          ! |0000 0000|0001 0100|
            interleaved_bits_ground_trouth = 2**1 + 2**4 + 2**8 + 2**11 ! |0000 1001|0001 0010|
#elif (_FMM_DIMENSION_ == 3)
            x(1) = 2**0 + 2**5                          ! |0000 0000|0010 0001|
            x(2) = 2**2 + 2**4                          ! |0000 0000|0001 0100|
            x(3) = 2**1 + 2**6                          ! |0000 0000|0100 0010|
            interleaved_bits_ground_trouth = 2**2 + 2**3 + 2**7 + 2**13 + 2**17 + 2**18   ! |0000 0110|0010 0000|1000 1100|
#else
            print *, "test_lib_tree_hf_interleave_bits_1D_to_treeD: Dimension not defines: ", FMM_DIMENSION
#endif

            interleaved_bits = lib_tree_hf_interleave_bits_treeD_to_1D(x)

            rv = .true.
            if (interleaved_bits == interleaved_bits_ground_trouth) then
                print *, "test_lib_tree_hf_interleave_bits_1D_to_treeD: ", "ok"
            else
                print *, "test_lib_tree_hf_interleave_bits_1D_to_treeD: ", "FAILED"
                rv = .false.
            end if

        end function test_lib_tree_hf_interleave_bits_treeD_to_1D

        function test_lib_tree_hf_interleave_bits_treeD_to_1D_use_lut() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliar
            integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS) :: x
            integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS) :: interleaved_x
            integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS) :: interleaved_x_ground_truth

            integer(kind=1) :: i

#if (_FMM_DIMENSION_ == 2)
            x(1) = 2**2 + 2**5  ! |0010 0100|
            x(2) = 2**0 + 2**7  ! |1000 0001|

            interleaved_x_ground_truth(1) = 2**0 + 2**5 + 2**11 + 2**14 !             .. |0100 1000|0010 0001|
            interleaved_x_ground_truth(2) = 0                           ! |0000 0000|
#elif (_FMM_DIMENSION_ == 3)
            x(1) = 2**2 + 2**5  ! |0010 0100|
            x(2) = 2**0 + 2**7  ! |1000 0001|
            x(3) = 2**1 + 2**6  ! |0100 0010|

            interleaved_x_ground_truth(1) = 2**1 + 2**3 + 2**8 &
                                            + 2**17 + 2**18 + 2**22     !             .. |0100 0110|0000 0001|0000 1010|
            interleaved_x_ground_truth(2) = 0                           ! |0000 0000|
            interleaved_x_ground_truth(3) = 0                           ! |0000 0000|
#else
            print *, "lib_tree_hf_interleave_bits_treeD_to_1D_use_lut: Dimension not defines: ", FMM_DIMENSION
#endif

            interleaved_x = lib_tree_hf_interleave_bits_treeD_to_1D_use_lut(x)

            print *, "lib_tree_hf_interleave_bits_treeD_to_1D_use_lut:"
            rv = .true.
            do i=1, TREE_DIMENSIONS
                if (interleaved_x(i) == interleaved_x_ground_truth(i)) then
                    print *, "  dim ",i,": ", "ok"
                else
                    print *, "  dim ",i,": ", "FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_tree_hf_interleave_bits_treeD_to_1D_use_lut

        function test_lib_tree_hf_deinterleave_bits_1D_to_treeD() result (rv)
            implicit none

            ! dummy
            logical :: rv

            ! parameter
            integer(kind=1), parameter :: x_kind = UINDEX_BYTES

            ! dummy
            integer(kind=x_kind) :: x
            integer(kind=x_kind), dimension(TREE_DIMENSIONS) :: deinterleaved_bits

            integer(kind=x_kind), dimension(TREE_DIMENSIONS) :: deinterleaved_bits_ground_trouth

            integer :: i

#if (_FMM_DIMENSION_ == 2)
            x = 2**0 + 2**5 !           |0010 0001|
            deinterleaved_bits_ground_trouth(1) = 2**2  ! 0000 0100
            deinterleaved_bits_ground_trouth(2) = 2**0  ! 0000 0001
#elif (_FMM_DIMENSION_ == 3)
            x = 2**0 + 2**5 !           |0010 0001|
            deinterleaved_bits_ground_trouth(1) = 2**1  ! 0000 0010
            deinterleaved_bits_ground_trouth(2) = 0     ! 0000 0000
            deinterleaved_bits_ground_trouth(3) = 2**0  ! 0000 0001
#else
            print *, "test_lib_tree_hf_deinterleave_bits_1D_to_treeD: Dimension not defines: ", FMM_DIMENSION
#endif

            deinterleaved_bits = lib_tree_hf_deinterleave_bits_1D_to_treeD(x)

            print *, "test_lib_tree_hf_deinterleave_bits_1D_to_treeD:"
            rv = .true.
            do i=1, TREE_DIMENSIONS
                if (deinterleaved_bits(i) == deinterleaved_bits_ground_trouth(i)) then
                    print *, "  dim ",i,": ", "ok"
                else
                    print *, "  dim ",i,": ", "FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_tree_hf_deinterleave_bits_1D_to_treeD

        function test_lib_tree_hf_deinterleave_bits_1D_to_treeD_use_lut() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer(kind=1), parameter :: x_kind = UINDEX_BYTES

            ! dummy
            integer(kind=x_kind), dimension(TREE_DIMENSIONS) :: x
            integer(kind=x_kind), dimension(TREE_DIMENSIONS) :: deinterleaved_bits

            integer(kind=x_kind), dimension(TREE_DIMENSIONS) :: deinterleaved_bits_ground_trouth

            integer :: i

#if (_FMM_DIMENSION_ == 2)
            x(1) = 2**0 + 2**5 !           |0010 0001|
            x(2) = 0
            deinterleaved_bits_ground_trouth(1) = 2**2  ! 0000 0100
            deinterleaved_bits_ground_trouth(2) = 2**0  ! 0000 0001
#elif (_FMM_DIMENSION_ == 3)
            x(1) = 2**0 + 2**5 !           |0010 0001|
            x(2) = 0
            x(3) = 0
            deinterleaved_bits_ground_trouth(1) = 2**1  ! 0000 0010
            deinterleaved_bits_ground_trouth(2) = 0     ! 0000 0000
            deinterleaved_bits_ground_trouth(3) = 2**0  ! 0000 0001
#else
            print *, "test_lib_tree_hf_deinterleave_bits_1D_to_treeD_use_lut: Dimension not defines: ", FMM_DIMENSION
#endif

            deinterleaved_bits = lib_tree_hf_deinterleave_bits_1D_to_treeD_use_lut(x)

            print *, "test_lib_tree_hf_deinterleave_bits_1D_to_treeD_use_lut:"
            rv = .true.
            do i=1, TREE_DIMENSIONS
                if (deinterleaved_bits(i) == deinterleaved_bits_ground_trouth(i)) then
                    print *, "  dim ",i,": ", "ok"
                else
                    print *, "  dim ",i,": ", "FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_tree_hf_deinterleave_bits_1D_to_treeD_use_lut

    end function lib_tree_hf_test_functions

    ! ----------------- benchmark functions -----------------
    subroutine lib_tree_hf_benchmark
        implicit none

        call benchmark_lib_tree_hf_interleave_bits_use_lut()
        call benchmark_lib_tree_hf_deinterleave_bits_use_lut()

        call benchmark_lib_tree_hf_interleave_bits_treeD_to_1D_use_lut()
        call benchmark_lib_tree_hf_deinterleave_bits_1D_to_treeD_use_lut()

        contains

        subroutine benchmark_lib_tree_hf_interleave_bits_use_lut()
            implicit none

            integer(kind=INTERLEAVE_BITS_INTEGER_KIND), dimension(TREE_DIMENSIONS) :: x
            integer(kind=INTERLEAVE_BITS_INTEGER_KIND), dimension(TREE_DIMENSIONS) :: buffer

            integer(kind=8) :: number_of_runs = 10**10
            integer(kind=8) :: i
            real :: start, finish
            double precision :: delta

            x(1) = 2
            x(2) = 0
#if (_FMM_DIMENSION_ == 3)
            x(3) = 0
#endif
            print *, "benchmark_lib_tree_hf_interleave_bits_use_lut"
            call cpu_time(start)
            buffer = lib_tree_hf_interleave_bits_use_lut(x)
            call cpu_time(finish)
            print *, "  Interleave + LUT Time = ", finish-start, " seconds."

            call cpu_time(start)
            do i=1, number_of_runs
                buffer = lib_tree_hf_interleave_bits_use_lut(x)
            end do
            call cpu_time(finish)
            delta = finish-start
            print *, "  Interleave + LUT Time (second run) = ", delta/number_of_runs, " seconds."

            number_of_runs = 10**17
            call cpu_time(start)
            do i=1, number_of_runs
                buffer = lib_tree_hf_interleave_bits(x)
            end do
            call cpu_time(finish)
            delta = finish-start
            print *, "  Interleave Time = ", delta/number_of_runs, " seconds."
            print *, ""
        end subroutine benchmark_lib_tree_hf_interleave_bits_use_lut

        subroutine benchmark_lib_tree_hf_deinterleave_bits_use_lut()
            implicit none

            integer(kind=INTERLEAVE_BITS_INTEGER_KIND), dimension(TREE_DIMENSIONS) :: x
            integer(kind=INTERLEAVE_BITS_INTEGER_KIND), dimension(TREE_DIMENSIONS) :: buffer

            integer(kind=8) :: number_of_runs = 10**10
            integer(kind=8) :: i
            real :: start, finish
            double precision :: delta

            x(1) = 2
            x(2) = 0
#if (_FMM_DIMENSION_ == 3)
            x(3) = 0
#endif
            print *, "benchmark_lib_tree_hf_interleave_bits_use_lut"
            call cpu_time(start)
            buffer = lib_tree_hf_deinterleave_bits_use_lut(x)
            call cpu_time(finish)
            print *, "  Deinterleave + LUT Time = ", finish-start, " seconds."

            call cpu_time(start)
            do i=1, number_of_runs
                buffer = lib_tree_hf_deinterleave_bits_use_lut(x)
            end do
            call cpu_time(finish)
            delta = finish-start
            print *, "  Deinterleave + LUT Time (second run) = ", delta/number_of_runs, " seconds."

            number_of_runs = 10**17
            call cpu_time(start)
            do i=1, number_of_runs
                buffer = lib_tree_hf_deinterleave_bits(x)
            end do
            call cpu_time(finish)
            delta = finish-start
            print *, "  Deinterleave Time = ", delta/number_of_runs, " seconds."
            print *, ""
        end subroutine benchmark_lib_tree_hf_deinterleave_bits_use_lut

        subroutine benchmark_lib_tree_hf_interleave_bits_treeD_to_1D_use_lut()
            implicit none

            integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS) :: x
            integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS) :: buffer

            integer(kind=8) :: number_of_runs = 10**10
            integer(kind=8) :: i
            real :: start, finish
            DOUBLE PRECISION :: delta

            x(1) = 2
            x(2) = 0
#if (_FMM_DIMENSION_ == 3)
            x(3) = 0
#endif
            print *, "benchmark_lib_tree_hf_interleave_bits_treeD_to_1D_use_lut"
            call cpu_time(start)
            buffer = lib_tree_hf_interleave_bits_treeD_to_1D_use_lut(x)
            call cpu_time(finish)
            print *, "  Interleave + LUT Time = ", finish-start, " seconds."

            call cpu_time(start)
            do i=1, number_of_runs
                buffer = lib_tree_hf_interleave_bits_treeD_to_1D_use_lut(x)
            end do
            call cpu_time(finish)
            delta = finish-start
            print *, "  Interleave + LUT Time (second run) = ", delta/number_of_runs, " seconds."

            number_of_runs = 10**18
            call cpu_time(start)
            do i=1, number_of_runs
                buffer = lib_tree_hf_interleave_bits_treeD_to_1D(x)
            end do
            call cpu_time(finish)
            delta = finish-start
            print *, "  Interleave Time = ", delta/number_of_runs, " seconds."
            print *, ""
        end subroutine benchmark_lib_tree_hf_interleave_bits_treeD_to_1D_use_lut

        subroutine benchmark_lib_tree_hf_deinterleave_bits_1D_to_treeD_use_lut()
            implicit none

            integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS) :: x
            integer(kind=UINDEX_BYTES), dimension(TREE_DIMENSIONS) :: buffer

            integer(kind=8) :: number_of_runs = 10**10
            integer(kind=8) :: i
            real :: start, finish
            DOUBLE PRECISION :: delta

            x(1) = 2
            x(2) = 0
#if (_FMM_DIMENSION_ == 3)
            x(3) = 0
#endif
            print *, "benchmark_lib_tree_hf_deinterleave_bits_1D_to_treeD_use_lut"
            call cpu_time(start)
            buffer = lib_tree_hf_deinterleave_bits_1D_to_treeD_use_lut(x)
            call cpu_time(finish)
            print *, "  Deinterleave + LUT Time = ", finish-start, " seconds."

            call cpu_time(start)
            do i=1, number_of_runs
                buffer = lib_tree_hf_deinterleave_bits_1D_to_treeD_use_lut(x)
            end do
            call cpu_time(finish)
            delta = finish-start
            print *, "  Deinterleave + LUT Time (second run) = ", delta/number_of_runs, " seconds."

            number_of_runs = 10**18
            call cpu_time(start)
            do i=1, number_of_runs
                buffer = lib_tree_hf_deinterleave_bits_1D_to_treeD(x(1))
            end do
            call cpu_time(finish)
            delta = finish-start
            print *, "  Deinterleave Time = ", delta/number_of_runs, " seconds."
            print *, ""
        end subroutine benchmark_lib_tree_hf_deinterleave_bits_1D_to_treeD_use_lut

    end subroutine lib_tree_hf_benchmark

end module lib_tree_helper_functions
