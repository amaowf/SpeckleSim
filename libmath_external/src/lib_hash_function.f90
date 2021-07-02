module lib_hash_function
#ifndef __GFORTRAN__
    USE IFPORT          ! ifort function: rand
#endif
    implicit none
    private

    public :: hash_fnv1a

    ! test functions
    public :: lib_test_hash_function

    ! interface
    interface hash_fnv1a
        module procedure hash_fnv1a32_2_byte
        module procedure hash_fnv1a32_4_byte
        module procedure hash_fnv1a32_8_byte
#ifdef __GFORTRAN__
        module procedure hash_fnv1a32_16_byte
#endif
        module procedure hash_fnv1a64_4_byte
        module procedure hash_fnv1a64_8_byte
#ifdef __GFORTRAN__
        module procedure hash_fnv1a64_16_byte
#endif
    end interface

    contains

    ! FNV 1a 32 bit
    !
    ! Reference: http://www.isthe.com/chongo/tech/comp/fnv/index.html#FNV-1a
    ! hash = offset_basis
    ! for each octet_of_data to be hashed
    !         hash = hash xor octet_of_data
    !         hash = hash * FNV_prime
    ! return hash
    !
    ! NOTE
    ! ----
    !   Replaced the bitwise copy process with a equivalence (performance reason).
    !
    ! Source code source:
    !   http://www.isthe.com/chongo/tech/comp/fnv/fnv32.f
    !
    !    *#######################################################################
    !    * 32 bit Fowler/Noll/Vo FNV-1a hash code
    !    * Public Domain
    !    *#######################################################################
    !    * fixed 32 -> 8 bit conversion bug            2013/10/12 AJA
    !    * translated to FORTRAN                       2013/03/11 Andy Allinger
    !    *                                             andy_a@users.sourceforge.net
    !    *
    !    * Revision: 5.1                               2009/06/30 09:13:32
    !    *#######################################################################
    !    * Fowler/Noll/Vo hash
    !    *
    !    * The basis of this hash algorithm was taken from an idea sent
    !    * as reviewer comments to the IEEE POSIX P1003.2 committee by:
    !    *
    !    *      Phong Vo (http://www.research.att.com/info/kpv/)
    !    *      Glenn Fowler (http://www.research.att.com/~gsf/)
    !    *
    !    * In a subsequent ballot round:
    !    *
    !    *      Landon Curt Noll (http://www.isthe.com/chongo/)
    !    *
    !    * improved on their algorithm.  Some people tried this hash
    !    * and found that it worked rather well.  In an EMail message
    !    * to Landon, they named it the ``Fowler/Noll/Vo'' or FNV hash.
    !    *
    !    * FNV hashes are designed to be fast while maintaining a low
    !    * collision rate. The FNV speed allows one to quickly hash lots
    !    * of data while maintaining a reasonable collision rate.  See:
    !    *
    !    *      http://www.isthe.com/chongo/tech/comp/fnv/index.html
    !    *
    !    * for more details as well as other forms of the FNV hash.
    !    *
    !    *#######################################################################
    !    * Please do not copyright this code.  This code is in the public domain.
    !    *
    !    * LANDON CURT NOLL DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
    !    * INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO
    !    * EVENT SHALL LANDON CURT NOLL BE LIABLE FOR ANY SPECIAL, INDIRECT OR
    !    * CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
    !    * USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
    !    * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
    !    * PERFORMANCE OF THIS SOFTWARE.
    !    *
    !    * By:
    !    *   chongo <Landon Curt Noll> /\oo/\
    !    *      http://www.isthe.com/chongo/
    !    *
    !    * Share and Enjoy!  :-)
    !    *
    !    *#######################################################################
    !    * 32 bit FNV-1 and FNV-1a non-zero initial basis
    !    *
    !    * The FNV-1 initial basis is the FNV-0 hash of the following 32 octets:
    !    *
    !    *              chongo <Landon Curt Noll> /\../\
    !    *
    !    * NOTE: The \'s above are not back-slashing escape characters.
    !    * They are literal ASCII  backslash 0x5c characters.
    !    *
    !    * NOTE: The FNV-1a initial basis is the same value as FNV-1 by definition.
    !    *
    !    *#######################################################################
    !    * perform a 32 bit Fowler/Noll/Vo FNV-1a hash on a buffer
    !    *
    !    * input:
    !    *   BUFFER  - start of buffer to hash
    !    *   LENGTH  - length of buffer in octets [bytes]
    !    *   HASH    - previous hash value or [standard initial value]
    !    *
    !    * assumption:  INTEGER's are at least 32 bits
    !    *
    !    * NOTE:  To use the recommended 32 bit FNV-1a hash, use the initial value:
    !    *               HASH = X'811C9DC5'
    !    *        as the argument on the first call to FNV32
    !    *
    !    * !NOTE:  Only pass the previous hash value if you always hash the same
    !    *         data in the same order (or else you will get different answers!)
    !    *
    !    * returns:
    !    *   32 bit HASH as default INTEGER
    !    *
    !    *#######################################################################
          SUBROUTINE FNV32 (BUFFER, LENGTH, HASH)
           IMPLICIT NONE
           INTEGER LENGTH, HASH
           INTEGER*1 BUFFER
           DIMENSION BUFFER(LENGTH)

           INTEGER PRIME ; PARAMETER (PRIME = 16777619)
           INTEGER J, K
           integer(kind=1), dimension(4) :: K_buffer
           INTEGER*1 B

           equivalence (K, K_buffer)                                    ! 190417 itodaiber: bitwise copy process replaced by the equivalence
           equivalence (B, K_buffer(1))                                 ! 190417 itodaiber: bitwise copy process replaced by the equivalence

    !    *#######################################################################
    !    *                begin
    !    *#######################################################################
    !    *          FNV-1a hash each octet in the buffer
           DO 90 J = 1, LENGTH
             B = BUFFER(J)
             K_buffer(2:4) = 0                                          ! 190417 itodaiber: bitwise copy process replaced by the equivalence
!             K = 0                                                     ! 190417 itodaiber: bitwise copy process replaced by the equivalence
!             DO 80 I = 0, 7           ! copy each bit from B to K      ! 190417 itodaiber: bitwise copy process replaced by the equivalence
!               IF (BTEST(B, I)) K = IBSET(K, I)                        ! 190417 itodaiber: bitwise copy process replaced by the equivalence
!      80     CONTINUE ! next i                                         ! 190417 itodaiber: bitwise copy process replaced by the equivalence

    !    *          xor the bottom with the current octet
             HASH = IEOR(HASH, K)

    !    *          multiply by the 32 bit FNV magic prime mod 2^32
             HASH = HASH * PRIME
             HASH = IAND(HASH, X'FFFFFFFF')      ! discard > 32 bits
      90   CONTINUE ! next j

      END SUBROUTINE!############## of file fnv32.f ##############################

      SUBROUTINE FNV64 (BUFFER, LENGTH, HASH)
           IMPLICIT NONE
           INTEGER(kind=4) LENGTH
           INTEGER(kind=8) HASH
           INTEGER*1 BUFFER
           DIMENSION BUFFER(LENGTH)

           INTEGER(kind=8), parameter :: PRIME = x'100000001b3'
           INTEGER(kind=8) J, K
           integer(kind=1), dimension(8) :: K_buffer
           INTEGER*1 B

           equivalence (K, K_buffer)                                    ! 190417 itodaiber: bitwise copy process replaced by the equivalence
           equivalence (B, K_buffer(1))                                 ! 190417 itodaiber: bitwise copy process replaced by the equivalence

    !    *#######################################################################
    !    *                begin
    !    *#######################################################################
    !    *          FNV-1a hash each octet in the buffer
           DO 90 J = 1, LENGTH
             B = BUFFER(J)
             K_buffer(2:8) = 0                                          ! 190417 itodaiber: bitwise copy process replaced by the equivalence
!             K = 0                                                     ! 190417 itodaiber: bitwise copy process replaced by the equivalence
!             DO 80 I = 0, 7           ! copy each bit from B to K      ! 190417 itodaiber: bitwise copy process replaced by the equivalence
!               IF (BTEST(B, I)) K = IBSET(K, I)                        ! 190417 itodaiber: bitwise copy process replaced by the equivalence
!      80     CONTINUE ! next i                                         ! 190417 itodaiber: bitwise copy process replaced by the equivalence

    !    *          xor the bottom with the current octet
             HASH = IEOR(HASH, K)

    !    *          multiply by the 64 bit FNV magic prime mod 2^64
             HASH = HASH * PRIME
             HASH = IAND(HASH, X'FFFFFFFFFFFFFFFF')      ! discard > 64 bits
      90   CONTINUE ! next j

      END SUBROUTINE!############## of file fnv32.f ##############################

!    function hash_fnv1a(buffer) result(hash)
!        ! dummy
!        integer(kind=4), intent(in) :: buffer
!        integer(kind=4) :: hash
!
!        ! auxiliary
!        integer(kind=4) :: buffer_buffer
!        integer(kind=1), dimension(4) :: buffer_list
!
!        equivalence (buffer_buffer, buffer_list)
!
!        buffer_buffer = buffer
!
!        hash = -2128831035
!        call FNV32(buffer_list, size(buffer_list), hash)
!
!    end function hash_fnv1a

    ! fnv1a hash function
    !
    ! Arguments
    ! ----
    !   buffer: integer(kind=4)
    !       data to be hashed
    !   max_value: integer(kind=4)
    !       maximum value of the hash
    !
    ! Returns
    ! ----
    !   a 32-bit hash value
    !
    function hash_fnv1a32_2_byte(buffer, max_value) result(hash)
        implicit none
        ! dummy
        integer(kind=2), intent(in) :: buffer
        integer(kind=4) :: max_value
        integer(kind=4) :: hash

        ! auxiliary
        double precision :: buffer_hash
        integer(kind=2) :: buffer_buffer
        integer(kind=1), dimension(2) :: buffer_list

        equivalence (buffer_buffer, buffer_list)

        buffer_buffer = buffer

        hash = -2128831035
        call FNV32(buffer_list, size(buffer_list), hash)

!        hash = int(real(hash,8) /  4294967296.0D0 * int(max_value-1,8) + max_value/2,4)
        !        <--     hash     -->
        !       |--------------------|
        !    -2**31               2**31-1        <-- integer(kind=8)
        !     -0.5                  0.5          <-- if hash < 0: / 2.0D0**32; else / (2*(2**31-1))
        ! -max_value/2+1        max_value/2-1    <-- * (max_value-2)
        !      1                  max_value      <-- + max_value / 2
        if (hash .lt. 0) then
            buffer_hash = real(hash,8) / 2.0D0**32
        else if (hash .gt. 0) then
            buffer_hash = real(hash,8) / (2.0D0**31 - 1.0D0) / 2.0D0
        else
            buffer_hash = 0.0D0
        end if
        hash = int(buffer_hash * int(max_value-2,8) + max_value/2.0D0,4)

    end function hash_fnv1a32_2_byte

    ! fnv1a hash function
    !
    ! Arguments
    ! ----
    !   buffer: integer(kind=4)
    !       data to be hashed
    !   max_value: integer(kind=4)
    !       maximum value of the hash
    !
    ! Returns
    ! ----
    !   a 32-bit hash value
    !
    function hash_fnv1a32_4_byte(buffer, max_value) result(hash)
        implicit none
        ! dummy
        integer(kind=4), intent(in) :: buffer
        integer(kind=4) :: max_value
        integer(kind=4) :: hash

        ! auxiliary
        double precision :: buffer_hash
        integer(kind=4) :: buffer_buffer
        integer(kind=1), dimension(4) :: buffer_list

        equivalence (buffer_buffer, buffer_list)

        buffer_buffer = buffer

        hash = -2128831035
        call FNV32(buffer_list, size(buffer_list), hash)

!        hash = int(real(hash,8) /  4294967296.0D0 * int(max_value-1,8) + max_value/2,4)
        !        <--     hash     -->
        !       |--------------------|
        !    -2**31               2**31-1        <-- integer(kind=8)
        !     -0.5                  0.5          <-- if hash < 0: / 2.0D0**32; else / (2*(2**31-1))
        ! -max_value/2+1        max_value/2-1    <-- * (max_value-2)
        !      1                  max_value      <-- + max_value / 2
        if (hash .lt. 0) then
            buffer_hash = real(hash,8) / 2.0D0**32
        else if (hash .gt. 0) then
            buffer_hash = real(hash,8) / (2.0D0**31 - 1.0D0) / 2.0D0
        else
            buffer_hash = 0.0D0
        end if
        hash = int(buffer_hash * int(max_value-2,8) + max_value/2.0D0,4)

    end function hash_fnv1a32_4_byte

    ! fnv1a hash function
    !
    ! Arguments
    ! ----
    !   buffer: integer(kind=8)
    !       data to be hashed
    !   max_value: integer(kind=4)
    !       maximum value of the hash
    !
    ! Returns
    ! ----
    !   a 32-bit hash value
    !
    function hash_fnv1a32_8_byte(buffer, max_value) result(hash)
        implicit none
        ! dummy
        integer(kind=8), intent(in) :: buffer
        integer(kind=4) :: max_value
        integer(kind=4) :: hash

        ! auxiliary
!        integer(kind=4) :: m_hash
        double precision :: buffer_hash
        integer(kind=8) :: buffer_buffer
        integer(kind=1), dimension(8) :: buffer_list

        equivalence (buffer_buffer, buffer_list)

        buffer_buffer = buffer

        hash = -2128831035
        call FNV32(buffer_list, size(buffer_list), hash)

!        hash = int(real(hash,8) /  4294967296.0D0 * int(max_value-1,8) + max_value/2,4)
        !        <--     hash     -->
        !       |--------------------|
        !    -2**31               2**31-1        <-- integer(kind=8)
        !     -0.5                  0.5          <-- if hash < 0: / 2.0D0**32; else / (2*(2**31-1))
        ! -max_value/2+1        max_value/2-1    <-- * (max_value-2)
        !      1                  max_value      <-- + max_value / 2
        if (hash .lt. 0) then
            buffer_hash = real(hash,8) / 2.0D0**32
        else if (hash .gt. 0) then
            buffer_hash = real(hash,8) / (2.0D0**31 - 1.0D0) / 2.0D0
        else
            buffer_hash = 0.0D0
        end if
        hash = int(buffer_hash * int(max_value-2,8) + max_value/2.0D0,4)

    end function hash_fnv1a32_8_byte

#ifdef __GFORTRAN__
    ! fnv1a hash function
    !
    ! Arguments
    ! ----
    !   buffer: integer(kind=16)
    !       data to be hashed
    !   max_value: integer(kind=4)
    !       maximum value of the hash
    !
    ! Returns
    ! ----
    !   a 32-bit hash value
    !
    function hash_fnv1a32_16_byte(buffer, max_value) result(hash)
        implicit none
        ! dummy
        integer(kind=16), intent(in) :: buffer
        integer(kind=4) :: max_value
        integer(kind=4) :: hash

        ! auxiliary
        double precision :: buffer_hash
        integer(kind=16) :: buffer_buffer
        integer(kind=1), dimension(16) :: buffer_list

        equivalence (buffer_buffer, buffer_list)

        buffer_buffer = buffer

        hash = -2128831035
        call FNV32(buffer_list, size(buffer_list), hash)

!        hash = int(real(hash,8) /  4294967296.0D0 * int(max_value-1,8) + max_value/2,4)

        !        <--     hash     -->
        !       |--------------------|
        !    -2**31               2**31-1        <-- integer(kind=8)
        !     -0.5                  0.5          <-- if hash < 0: / 2.0D0**32; else / (2*(2**31-1))
        ! -max_value/2+1        max_value/2-1    <-- * (max_value-2)
        !      1                  max_value      <-- + max_value / 2
        if (hash .lt. 0) then
            buffer_hash = real(hash,8) / 2.0D0**32
        else if (hash .gt. 0) then
            buffer_hash = real(hash,8) / (2.0D0**31 - 1.0D0) / 2.0D0
        else
            buffer_hash = 0.0D0
        end if
        hash = int(buffer_hash * int(max_value-2,8) + max_value/2.0D0,4)

    end function hash_fnv1a32_16_byte
#endif

    ! fnv1a hash function
    !
    ! Arguments
    ! ----
    !   buffer: integer(kind=4)
    !       data to be hashed
    !   max_value: integer(kind=8)
    !       maximum value of the hash
    !
    ! Returns
    ! ----
    !   a 64-bit hash value
    !
    function hash_fnv1a64_4_byte(buffer, max_value) result(hash)
        implicit none
        ! dummy
        integer(kind=4), intent(in) :: buffer
        integer(kind=8) :: max_value
        integer(kind=8) :: hash

        ! auxiliary
!        integer(kind=4) :: m_hash
        integer(kind=4) :: buffer_buffer
        integer(kind=1), dimension(4) :: buffer_list
        double precision :: buffer_hash

        equivalence (buffer_buffer, buffer_list)

        buffer_buffer = buffer

        hash = -5472609002491880230_8
        call FNV64(buffer_list, size(buffer_list), hash)

        !        <--     hash     -->
        !       |--------------------|
        !    -2**63               2**63-1        <-- integer(kind=8)
        !     -0.5                  0.5          <-- if hash < 0: / 2.0D0**64; else / (2*(2**63-1))
        ! -max_value/2+1        max_value/2-1    <-- * (max_value-2)
        !      1                  max_value      <-- + max_value / 2
        if (hash .lt. 0) then
            buffer_hash = real(hash,8) / 2.0D0**64
        else if (hash .gt. 0) then
            buffer_hash = real(hash,8) / (2.0D0**63 - 1.0D0) / 2.0D0
        else
            buffer_hash = 0.0D0
        end if
        hash = int(buffer_hash * int(max_value-2,8) + max_value/2.0D0,8)

    end function hash_fnv1a64_4_byte

    ! fnv1a hash function
    !
    ! Arguments
    ! ----
    !   buffer: integer(kind=8)
    !       data to be hashed
    !   max_value: integer(kind=8)
    !       maximum value of the hash
    !
    ! Returns
    ! ----
    !   a 64-bit hash value
    !
    function hash_fnv1a64_8_byte(buffer, max_value) result(hash)
        implicit none
        ! dummy
        integer(kind=8), intent(in) :: buffer
        integer(kind=8) :: max_value
        integer(kind=8) :: hash

        ! auxiliary
!        integer(kind=4) :: m_hash
        integer(kind=8) :: buffer_buffer
        integer(kind=1), dimension(8) :: buffer_list
        double precision :: buffer_hash

        equivalence (buffer_buffer, buffer_list)

        buffer_buffer = buffer

        hash = -5472609002491880230_8
        call FNV64(buffer_list, size(buffer_list), hash)

        !        <--     hash     -->
        !       |--------------------|
        !    -2**63               2**63-1        <-- integer(kind=8)
        !     -0.5                  0.5          <-- if hash < 0: / 2.0D0**64; else / (2*(2**63-1))
        ! -max_value/2+1        max_value/2-1    <-- * (max_value-2)
        !      1                  max_value      <-- + max_value / 2
        if (hash .lt. 0) then
            buffer_hash = real(hash,8) / 2.0D0**64
        else if (hash .gt. 0) then
            buffer_hash = real(hash,8) / (2.0D0**63 - 1.0D0) / 2.0D0
        else
            buffer_hash = 0.0D0
        end if
        hash = int(buffer_hash * int(max_value-2,8) + max_value/2.0D0,8)

    end function hash_fnv1a64_8_byte

#ifdef __GFORTRAN__
    ! fnv1a hash function
    !
    ! Arguments
    ! ----
    !   buffer: integer(kind=16)
    !       data to be hashed
    !   max_value: integer(kind=8)
    !       maximum value of the hash
    !
    ! Returns
    ! ----
    !   a 64-bit hash value
    !
    function hash_fnv1a64_16_byte(buffer, max_value) result(hash)
        implicit none
        ! dummy
        integer(kind=16), intent(in) :: buffer
        integer(kind=8) :: max_value
        integer(kind=8) :: hash

        ! auxiliary
!        integer(kind=4) :: m_hash
        integer(kind=16) :: buffer_buffer
        integer(kind=1), dimension(16) :: buffer_list
        double precision :: buffer_hash

        equivalence (buffer_buffer, buffer_list)

        buffer_buffer = buffer

        hash = -5472609002491880230_8
        call FNV64(buffer_list, size(buffer_list), hash)

        !        <--     hash     -->
        !       |--------------------|
        !    -2**63               2**63-1        <-- integer(kind=8)
        !     -0.5                  0.5          <-- if hash < 0: / 2.0D0**64; else / (2*(2**63-1))
        ! -max_value/2+1        max_value/2-1    <-- * (max_value-2)
        !      1                  max_value      <-- + max_value / 2
        if (hash .lt. 0) then
            buffer_hash = real(hash,8) / 2.0D0**64
        else if (hash .gt. 0) then
            buffer_hash = real(hash,8) / (2.0D0**63 - 1.0D0) / 2.0D0
        else
            buffer_hash = 0.0D0
        end if
        hash = int(buffer_hash * int(max_value-2,8) + max_value/2.0D0,8)

    end function hash_fnv1a64_16_byte
#endif

!
!    ! fnv1a hash function
!    !
!    ! Arguments
!    ! ----
!    !   buffer: integer(kind=16)
!    !       data to be hashed
!    !   max_value: integer(kind=4)
!    !       maximum value of the hash
!    !
!    ! Returns
!    ! ----
!    !   a 32-bit hash value
!    !
!    function hash_fnv1a64_16_byte(buffer, max_value) result(hash)
!        implicit none
!        ! dummy
!        integer(kind=16), intent(in) :: buffer
!        integer(kind=8) :: max_value
!        integer(kind=8) :: hash
!
!        ! auxiliary
!        integer(kind=16) :: buffer_buffer
!        integer(kind=1), dimension(16) :: buffer_list
!
!        equivalence (buffer_buffer, buffer_list)
!
!        buffer_buffer = buffer
!
!        hash = (/0b1100101111110010100111001110010010000100001000100010001100100101_8/)
!        call FNV32(buffer_list, size(buffer_list), hash)
!
!        hash = int(real(hash,8) /  4294967296.0D0 * int(max_value-1,8) + max_value/2,4)
!
!    end function hash_fnv1a64_16_byte

!    ! ********************************************************* hash
!    !berechnet ndigit-Hash wert der Koordinaten Matrix a,b,max,n ,versuch i
!    ! Source: K. Frenner: fem_sparse.F90
!    integer(kind=8) function hash(a,b,max,i,idum)
!    implicit none
!    integer,intent(in)::a,b,max,i
!    integer(kind=8),intent(inout):: idum
!      integer i2
!      idum=int(a+int(b,kind=8)*int(max,kind=8),kind=8) !als start
!      do i2=1,i+2,1  !2=Offset zum einschwingen !
!        idum=int(modulo(real(16807.0D0*idum,kind=16),2147483647.0D0),kind=8)
!      end do
!      hash=int(real(idum,kind=8)/2147483647.0D0*real(max-1,kind=8),kind=8)
!    end function

!    ! ********************************************************* hash
!    !berechnet ndigit-Hash wert der Koordinaten Matrix a,b,max,n ,versuch i
!    ! Source: K. Frenner: fem_sparse.F90
!    integer(kind=8) function hash_kf(a,b,max,i)
!    implicit none
!    integer(kind=8),intent(in)::a,b,max
!    integer, intent(in)::i
!    integer(kind=8):: idum
!      integer i2
!      idum=int(a+int(b,kind=8)*int(max,kind=8),kind=8) !als start
!      do i2=1,i+2,1  !2=Offset zum einschwingen !
!        idum=int(modulo(real(16807.0D0*idum,kind=16),2147483647.0D0),kind=8)
!      end do
!      hash_kf=int(real(idum,kind=8)/2147483647.0D0*real(max-1,kind=8),kind=8)
!    end function
!
!
!    ! ********************************************************* hashpp
!    !beschleunigte Fkt
!    ! Source: K. Frenner: fem_sparse.F90
!    integer(kind=8) function hashpp_kf(max,idum)
!    implicit none
!    integer(kind=8),intent(inout):: idum
!    integer(kind=8),intent(in)::max
!      idum=int(modulo(real(16807.0D0*idum,kind=16),2147483647.0D0),kind=8)
!      hashpp_kf=int(real(idum,kind=8)/2147483647.0D0*real(max-1,kind=8),kind=8)
!    end function

    ! ********************************************************* hash
    !berechnet ndigit-Hash wert der Koordinaten Matrix a,b,max,n ,versuch i
    integer(kind=8) function hash_kf(a,b,max,i,idum) result(hash)
    implicit none
    integer(kind=8),intent(in)::a,b,max,i
    integer(kind=8),intent(inout):: idum
      integer i2
      idum=int(a+int(b,kind=8)*int(max,kind=8),kind=8) !als start
      do i2=1,i+2,1  !2=Offset zum einschwingen !
        idum=int(modulo(real(16807.0D0*idum,kind=16),2147483647.0D0),kind=8)
      end do
      hash=int(real(idum,kind=8)/2147483647.0D0*real(max-1,kind=8),kind=8)
    end function

    ! ********************************************************* hashpp
    !beschleunigte Fkt
    integer(kind=8) function hashpp_kf(max,idum) result (hashpp)
    implicit none
    integer(kind=8),intent(inout):: idum
    integer(kind=8),intent(in)::max
      idum=int(modulo(real(16807.0D0*idum,kind=16),2147483647.0D0),kind=8)
      hashpp=int(real(idum,kind=8)/2147483647.0D0*real(max-1,kind=8),kind=8)
    end function

    ! ********************************************************* hash
    !berechnet ndigit-Hash wert der Koordinaten Matrix a,b,max,n ,versuch i
    integer(kind=8) function hash_kf_8_byte(a,max,i,idum) result(hash)
    implicit none
        ! parameter
        integer(kind=1), parameter :: HASH_DIMENSION = 4

        ! dummy
        integer(kind=8), intent(in) :: a
        integer(kind=4), intent(in)::i, max
        integer(kind=8), dimension(HASH_DIMENSION), intent(inout):: idum

        ! auxiliary
        integer(kind=4) i2
        integer(kind=8) :: aa
        integer(kind=2), dimension(HASH_DIMENSION) :: buffer_a
        real(kind=8) :: buffer
        integer(kind=1) :: ii

        equivalence(aa, buffer_a)

        aa = a

        if (a .lt. 0) then
            print *, "hash_kf_8_byte: a < 0"
            print *, "  ", a
        end if

        do ii=1, HASH_DIMENSION
            idum(ii)=int(ishft(int(buffer_a(ii), kind=8),-8)+int(4,kind=8)*int(max,kind=8),kind=8) !als start
        end do
        do i2=1,i+2,1  !2=Offset zum einschwingen !
            do ii=1, HASH_DIMENSION
                idum(ii)=int(modulo(real(16807.0D0*idum(ii),kind=16),2147483647.0D0),kind=8)
            end do
        end do
        buffer=1.0
        do ii=1, HASH_DIMENSION
            buffer = buffer * real(idum(i),kind=16)/2147483647.0D0
        end do
        hash=int(buffer*real(max-1,kind=16),kind=8)

        if (hash .lt. 0) then
            print *, "hash_kf_8_byte: hash < 0"
            print *, "  ", hash
        end if
    end function

    ! ********************************************************* hashpp
    !beschleunigte Fkt
    integer(kind=8) function hashpp_kf_8_byte(max,idum) result (hashpp)
    implicit none
        ! parameter
        integer(kind=1), parameter :: HASH_DIMENSION = 4

        integer(kind=8), dimension(2), intent(inout):: idum
        integer(kind=4), intent(in)::max

        !auxiliary
        real(kind=8) :: buffer
        integer(kind=1) :: i

        if (idum(1) .lt. 0) then
            print *, "hashpp_kf_8_byte: idum(1) <= 0"
            print *, "  ", idum
        end if

        do i=1, HASH_DIMENSION
            idum(i)=int(modulo(real(16807.0D0*idum(i),kind=16),2147483647.0D0),kind=8)
        end do

        buffer=1.0
        do i=1, HASH_DIMENSION
            buffer = buffer * real(idum(i),kind=16)/2147483647.0D0
        end do
        hashpp=int(buffer*real(max-1,kind=16),kind=8)

        if (hashpp .lt. 0) then
            print *, "hashpp_kf_8_byte: hashpp < 0"
            print *, "  ", hashpp
        end if
    end function

#ifdef __GFORTRAN__
    ! ********************************************************* hash
    !berechnet ndigit-Hash wert der Koordinaten Matrix a,b,max,n ,versuch i
    integer(kind=8) function hash_kf_16_byte(a,b,max,i,idum) result(hash)
    implicit none
        integer(kind=16), intent(in) :: a
        integer(kind=8),intent(in)::b,max,i
        integer(kind=8), dimension(2), intent(inout):: idum
        integer(kind=8) i2

        ! auxiliary
        integer(kind=16) :: aa
        integer(kind=8), dimension(2) :: buffer_a
        real(kind=8) :: buffer

        equivalence(aa, buffer_a)

        aa = a

        idum(1)=int(buffer_a(1)+int(b,kind=8)*int(max,kind=8),kind=8) !als start
        idum(2)=int(buffer_a(2)+int(b,kind=8)*int(max,kind=8),kind=8) !als start
        do i2=1,i+2,1  !2=Offset zum einschwingen !
            idum(1)=int(modulo(real(16807.0D0*idum(1),kind=16),2147483647.0D0),kind=8)
            idum(2)=int(modulo(real(16807.0D0*idum(2),kind=16),2147483647.0D0),kind=8)
        end do
        buffer = real(idum(1),kind=8)/2147483647.0D0
        buffer = buffer * real(idum(2),kind=8)/2147483647.0D0
        hash=int(buffer*real(max-1,kind=8),kind=8)
    end function

    ! ********************************************************* hashpp
    !beschleunigte Fkt
    integer(kind=8) function hashpp_kf_16_byte(max,idum) result (hashpp)
    implicit none
        integer(kind=8), dimension(2),intent(inout):: idum
        integer(kind=8), intent(in)::max

        !auxiliary
        real(kind=8) :: buffer
        idum(1)=int(modulo(real(16807.0D0*idum(1),kind=16),2147483647.0D0),kind=8)
        idum(2)=int(modulo(real(16807.0D0*idum(2),kind=16),2147483647.0D0),kind=8)

        buffer = real(idum(1),kind=8)/2147483647.0D0
        buffer = buffer * real(idum(2),kind=8)/2147483647.0D0
        hashpp=int(buffer*real(max-1,kind=8),kind=8)
    end function
#endif

!    function hash_8_byte(a, max) result(hash)
!        implicit none
!        ! dummy
!        integer(kind=8) :: a
!        integer(kind=4) :: max
!        integer(kind=4) :: hash
!
!        ! auxiliary
!        integer(kind=8) :: idum
!        idum = a
!
!        hash = int(numerical_recipes_ran(idum) * real((max-1)), 4)
!
!    end function hash_8_byte

!    ! Numerical Recipes in F90
!    ! p. 1142
!    FUNCTION numerical_recipes_ran(idum)
!        IMPLICIT NONE
!        INTEGER, PARAMETER :: K4B=8!selected_int_kind(9)
!        INTEGER(K4B), INTENT(INOUT) :: idum
!        REAL :: numerical_recipes_ran
!        INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
!        REAL :: am
!        INTEGER(K4B):: ix=-1,iy=-1,k
!        if (idum <= 0 .or. iy < 0) then
!            am=nearest(1.0,-1.0)/IM
!            iy=ior(ieor(888889999,abs(idum)),1)
!            ix=ieor(777755555,abs(idum))
!            idum=abs(idum)+1
!        end if
!        ix=ieor(ix,ishft(ix,13))
!        ix=ieor(ix,ishft(ix,-17))
!        ix=ieor(ix,ishft(ix,5))
!        k=iy/IQ
!        iy=IA*(iy-k*IQ)-IR*k
!        if (iy < 0) iy=iy+IM
!        numerical_recipes_ran=am*ior(iand(IM,ieor(ix,iy)),1)
!    END FUNCTION numerical_recipes_ran


    ! ----------- test functions -----------
    function lib_test_hash_function() result(error_counter)
        implicit none
        ! dummy
        integer :: error_counter

        error_counter = 0

        if (.not. test_hash_fnv1a_8_byte()) then
            error_counter = error_counter + 1
        end if
!        if (.not. test_hash_8_byte()) then
!            error_counter = error_counter + 1
!        end if
!        if (.not. test_hash_kf_8_byte()) then
!            error_counter = error_counter + 1
!        end if
!        if (.not. test_hashpp_kf_8_byte()) then
!            error_counter = error_counter + 1
!        end if

!        if (.not. test_hashpp()) then
!            print *, "test_hashpp: error"
!        end if

        contains

!        function test_hash() result (rv)
!            implicit none
!            ! dummy
!            logical :: rv
!
!            integer ::i
!            integer(kind=8) :: a, start, max, rv_hash
!
!            a = 40
!            max = 10**7
!            i=2
!            start = 1
!
!            rv_hash = hash_kf(start, a, max, i)
!            rv_hash = hash_kf(start, a, max, i)
!            rv_hash = hash_kf(start, a, max, i)
!            rv_hash = hash_kf(start, rv_hash, max, i)
!            rv_hash = hash_kf(start, rv_hash, max, i)
!            rv_hash = hash_kf(start, rv_hash, max, i)
!
!            rv = .true.
!
!        end function

!        function test_hashpp() result (rv)
!            implicit none
!
!            ! dummy
!            logical :: rv
!
!            integer(kind=8) :: idum
!            integer(kind=8) :: max
!            integer(kind=8) :: hash
!
!            max = 2
!            idum = 1
!
!            hash = hashpp_kf(max, idum)
!            hash = hashpp_kf(max, idum)
!            hash = hashpp_kf(max, idum)
!            hash = hashpp_kf(max, idum)
!            hash = hashpp_kf(max, idum)
!            hash = hashpp_kf(max, idum)
!            hash = hashpp_kf(max, idum)
!            hash = hashpp_kf(max, idum)
!
!            rv = .true.
!
!        end function test_hashpp

        function test_hash_kf_8_byte() result (rv)
            implicit none
            ! dummy
            logical :: rv

            integer(kind=8), dimension(3) :: a
            integer(kind=4) :: max
            integer(kind=8), dimension(4) :: idum
            integer(kind=4) :: hash_i

            integer(kind=8), dimension(3) :: hash

            integer :: i

            integer,parameter :: seed = 86456
            call srand(seed)

            max = 10
            hash_i = 10
            a(2) = int(rand()*1D7, kind=8)
            a(1) = a(2) - 1
            a(3) = a(2) + 1

            do i=1, 3
                hash(i) = hash_kf_8_byte(a(i), max, i, idum)
            end do

            if ((hash(1) .ne. hash(2)) .and. &
                (hash(1) .ne. hash(3)) .and. &
                (hash(2) .ne. hash(3))) then
                print *, "test_hash_kf_8_byte: OK"
                rv = .true.
            else
                print *, "test_hash_kf_8_byte: FAILED"
                print *, "     a:", a(1), a(2), a(3)
                print *, "  hash:", hash(1), hash(2), hash(3)
                rv = .false.
            end if

        end function test_hash_kf_8_byte

        function test_hashpp_kf_8_byte() result (rv)
            implicit none
            ! dummy
            logical :: rv

            integer(kind=8), dimension(3) :: a
            integer(kind=4) :: max
            integer(kind=8), dimension(4) :: idum
            integer(kind=4) :: hash_i

            integer(kind=8), dimension(3) :: hash

            integer :: i
            integer :: ii

            integer,parameter :: seed = 86456
            call srand(seed)

            max = 10
            hash_i = 10
            a(2) = int(rand()*1D7, kind=8)
            a(1) = a(2) - 1
            a(3) = a(2) + 1

            do i=1, 3
                hash(i) = hash_kf_8_byte(a(i), max, i, idum)
                do ii=1, 5
                    hash(i) = hashpp_kf_8_byte(max, idum)
                end do
            end do

            if ((hash(1) .ne. hash(2)) .and. &
                (hash(1) .ne. hash(3)) .and. &
                (hash(2) .ne. hash(3))) then
                print *, "test_hashpp_kf_8_byte: OK"
                rv = .true.
            else
                print *, "test_hashpp_kf_8_byte: FAILED"
                print *, "     a:", a(1), a(2), a(3)
                print *, "  hash:", hash(1), hash(2), hash(3)
                rv = .false.
            end if

        end function test_hashpp_kf_8_byte

!        function test_hash_8_byte() result(rv)
!            implicit none
!            ! dummy
!            logical :: rv
!
!            ! auxiliary
!            integer(kind=4) :: max
!            integer(kind=8) :: n
!            integer(kind=4) :: hash
!
!            max = 10
!            n = 74
!            hash = hash_8_byte(n,max)
!
!            n = -74
!            hash = hash_8_byte(n,max)
!
!!            n = 0
!!            hash = hash_8_byte(n,max)
!
!            n = 7
!            hash = hash_8_byte(int(-7,8),max)
!            hash = hash_8_byte(n,max)
!            n = 7
!            hash = hash_8_byte(int(-7,8),max)
!            hash = hash_8_byte(n,max)
!
!            n = 6
!            hash = hash_8_byte(int(-6,8),max)
!            hash = hash_8_byte(n,max)
!
!            n = 70
!            hash = hash_8_byte(n,max)
!
!            n = 7
!            hash = hash_8_byte(n,max)
!
!
!            rv = .true.
!        end function test_hash_8_byte

        function test_hash_fnv1a_8_byte() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer(kind=4) :: max
            integer(kind=8) :: n
            integer(kind=4) :: hash

            max = 10
            n = 74
            hash = hash_fnv1a32_8_byte(n,max)

            n = -74
            hash = hash_fnv1a32_8_byte(n,max)

!            n = 0
!            hash = hash_8_byte(n,max)

            n = 7
            hash = hash_fnv1a32_8_byte(n,max)
            hash = hash_fnv1a32_8_byte(n,max)

            n = 6
            hash = hash_fnv1a32_8_byte(n,max)

            n = 70
            hash = hash_fnv1a32_8_byte(n,max)
            hash = hash_fnv1a32_8_byte(int(hash,8),max)

            n = 7
            hash = hash_fnv1a32_8_byte(n,max)


            rv = .true.
        end function test_hash_fnv1a_8_byte

    end function lib_test_hash_function
end module lib_hash_function
