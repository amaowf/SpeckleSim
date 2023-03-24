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

module lib_math_public
    use lib_math_auxiliaries
    use lib_math_bessel
    use lib_math_convergence
    use lib_math_factorial
    use lib_math_legendre
    use lib_math_hermite
    use lib_math_solver
    use lib_math_solver_GMRES
    use lib_math_type
    use lib_math_type_operator
#ifdef __GFORTRAN__
!    use lib_math_wigner
!    use lib_math_vector_spherical_harmonics
#endif
    use lib_math_constants
    use lib_hash_function
    use lib_math_probability_distribution
    use lib_math_set_theory
    use lib_math_triangle
    implicit none

end module lib_math_public
