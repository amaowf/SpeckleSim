!    Copyright (C) 2020  Max Daiber-Huppert <max_daiber-huppert@gmx.de>
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
! Created on Thu Jan 30 13:07:51 2020
! 
! @author: Max Daiber-Huppert
!

module lib_field_polarisation_type
    implicit none

    type jones_vector_type
        double complex :: x
        double complex :: y
    end type


    !
    !     m_11  m_12
    ! M =
    !     m_21  m_22
    !
    type jones_matrix_type
        double complex :: m_11
        double complex :: m_12
        double complex :: m_21
        double complex :: m_22
    end type

    contains

        ! Retruns
        ! ----
        !   J: type(jones_vector_type)
        !       horizontal / x polarised light
        !
        ! Reference: Experimentalphysik 2, W. Demtröder
        !            p. 289
        function lib_field_polarisation_jones_vector_get_linear_h() result(J)
            implicit none
            ! dummy
            type(jones_vector_type) :: J

            J%x = dcmplx(1, 0)
            J%y = dcmplx(0, 0)

        end function

        ! Retruns
        ! ----
        !   J: type(jones_vector_type)
        !       vertical / y polarised light
        !
        ! Reference: Experimentalphysik 2, W. Demtröder
        !            p. 289
        function lib_field_polarisation_jones_vector_get_linear_v() result(J)
            implicit none
            ! dummy
            type(jones_vector_type) :: J

            J%x = dcmplx(0, 0)
            J%y = dcmplx(1, 0)

        end function

        ! Argument
        ! ----
        !   phi: double precision
        !       angel to x-axis [0, PI/2] [rad]
        !
        ! Retruns
        ! ----
        !   J: type(jones_vector_type)
        !       linear polarisated light with the angle phi to the x-axis
        !
        ! Reference: Experimentalphysik 2, W. Demtröder
        !            eq. 9.52
        function lib_field_polarisation_jones_vector_get_linear_rot(phi) result(J)
            implicit none
            ! dummy
            double precision, intent(in) :: phi

            type(jones_vector_type) :: J

            J%x = dcmplx(cos(phi), 0)
            J%y = dcmplx(sin(phi), 0)

        end function

        ! Retruns
        ! ----
        !   J: type(jones_vector_type)
        !       circular polaricated light (sigma_plus)
        !
        ! Reference: Experimentalphysik 2, W. Demtröder
        !            eq. 9.53a
        function lib_field_polarisation_jones_vector_get_circular_plus() result(J)
            implicit none
            ! dummy
            type(jones_vector_type) :: J

            J%x = dcmplx(1d0 / sqrt(2d0), 0)
            J%y = dcmplx(0, 1d0 / sqrt(2d0))

        end function

        ! Retruns
        ! ----
        !   J: type(jones_vector_type)
        !       circular polaricated light (sigma_plus)
        !
        ! Reference: Experimentalphysik 2, W. Demtröder
        !            eq. 9.53b
        function lib_field_polarisation_jones_vector_get_circular_minus() result(J)
            implicit none
            ! dummy
            type(jones_vector_type) :: J

            J%x = dcmplx(1d0 / sqrt(2d0), 0)
            J%y = dcmplx(0, -1d0 / sqrt(2d0))

        end function

        ! Argument
        ! ----
        !   phi: double precision
        !       angel [rad]
        !
        ! Retruns
        ! ----
        !   J: type(jones_vector_type)
        !       linear polarisated light with the angle phi to the x-axis
        !
        ! Reference: Experimentalphysik 2, W. Demtröder
        !            p. 289
        function lib_field_polarisation_jones_vector_get_elliptical(phi) result(J)
            implicit none
            ! dummy
            double precision, intent(in) :: phi

            type(jones_vector_type) :: J

            J%x = dcmplx(1, 0)
            J%y = dcmplx(-cos(phi), sin(phi))

        end function

        ! Returns
        ! ----
        !   M: type(jones_matrix_type)
        !       x polariser Jones matrix
        !
        ! Reference: Experimentalphysik 2, W. Demtröder
        !            eq. 9.54a
        function lib_field_polarisation_jones_matrix_get_x_polariser() result(M)
            implicit none
            ! dummy
            type(jones_matrix_type) :: M

            M%m_11 = dcmplx(1d0, 0)
            M%m_12 = dcmplx(0, 0)
            M%m_21 = dcmplx(0, 0)
            M%m_22 = dcmplx(0, 0)
        end function

        ! Argument
        ! ----
        !   phi: double precision
        !       anglte to the x-axis [rad]
        !
        ! Returns
        ! ----
        !   M: type(jones_matrix_type)
        !       rotation matrix of component (angle phi to the x-axis)
        !
        ! Reference: Experimentalphysik 2, W. Demtröder
        !            eq. 9.54c
        function lib_field_polarisation_jones_matrix_rotate_component(phi) result(M)
            implicit none
            ! dummy
            double precision, intent(in) :: phi

            type(jones_matrix_type) :: M

            ! auxiliary
            double precision :: cos_phi
            double precision :: sin_phi

            cos_phi = cos(phi)
            sin_phi = sin(phi)

            M%m_11 = dcmplx(cos_phi**2, 0)
            M%m_12 = dcmplx(sin_phi * cos_phi, 0)
            M%m_21 = dcmplx(sin_phi * cos_phi, 0)
            M%m_22 = dcmplx(sin_phi**2, 0)
        end function
end module
