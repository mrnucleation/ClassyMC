!================================================================================
! Spherical Harmonics Module for Fast Multipole Method
!
! This module provides utilities for computing spherical harmonics and
! related quantities needed for FMM multipole expansions and translations.
!
! Key functions:
!   - Associated Legendre polynomials P_l^m(x)
!   - Spherical harmonics Y_l^m(theta, phi)
!   - Solid harmonics for multipole/local expansions
!   - Translation coefficients for M2M, M2L, L2L operators
!
! Conventions:
!   - Uses real-valued spherical harmonics where possible
!   - Indexing: for order p, we have (p+1)^2 coefficients
!   - Linear index: idx(l,m) = l^2 + l + m + 1 for l=0..p, m=-l..l
!
! References:
!   - Greengard & Rokhlin (1987) - A Fast Algorithm for Particle Simulations
!   - Cheng, Greengard, Rokhlin (1999) - A Fast Adaptive Multipole Algorithm
!================================================================================
module SphericalHarmonics
  use VarPrecision
  implicit none
  
  private
  
  ! Public procedures
  public :: SH_Init
  public :: SH_Cleanup
  public :: SH_Index
  public :: SH_AssocLegendre
  public :: SH_ComputeYlm
  public :: SH_ComputeSolidHarmonic
  public :: SH_ComputeIrregularSolid
  public :: SH_M2M_Coeff
  public :: SH_M2L_Coeff
  public :: SH_L2L_Coeff
  public :: SH_CartesianToSpherical
  public :: SH_GetExpansionSize
  
  ! Mathematical constants
  real(dp), parameter :: PI_SH = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: SQRT_PI = 1.7724538509055159_dp
  complex(dp), parameter :: IMAG_UNIT = (0.0_dp, 1.0_dp)
  
  ! Precomputed factorials and coefficients
  integer :: maxOrder = -1
  real(dp), allocatable :: factorial(:)      ! n! for n = 0..2*maxOrder
  real(dp), allocatable :: invFactorial(:)   ! 1/n!
  real(dp), allocatable :: sqrtFact(:)       ! sqrt(n!)
  real(dp), allocatable :: A_lm(:,:)         ! Normalization coefficients
  
  ! Wigner rotation matrix elements (for rotations, optional)
  logical :: coefficientsInitialized = .false.
  
contains

!================================================================================
! Initialize spherical harmonics module with given maximum order
!================================================================================
subroutine SH_Init(p)
  implicit none
  integer, intent(in) :: p
  integer :: n, l, m
  real(dp) :: norm
  
  if (coefficientsInitialized .and. p <= maxOrder) return
  
  ! Cleanup any existing allocations
  call SH_Cleanup()
  
  maxOrder = p
  
  ! Allocate factorial arrays (need up to (2p)!)
  allocate(factorial(0:2*p+1))
  allocate(invFactorial(0:2*p+1))
  allocate(sqrtFact(0:2*p+1))
  allocate(A_lm(0:p, -p:p))
  
  ! Compute factorials
  factorial(0) = 1.0_dp
  do n = 1, 2*p+1
    factorial(n) = factorial(n-1) * real(n, dp)
  end do
  
  ! Compute inverse factorials
  do n = 0, 2*p+1
    invFactorial(n) = 1.0_dp / factorial(n)
  end do
  
  ! Compute square root of factorials
  do n = 0, 2*p+1
    sqrtFact(n) = sqrt(factorial(n))
  end do
  
  ! Compute normalization coefficients A_l^m
  ! A_l^m = sqrt((2l+1)/(4*pi) * (l-|m|)!/(l+|m|)!)
  A_lm = 0.0_dp
  do l = 0, p
    do m = -l, l
      norm = real(2*l + 1, dp) / (4.0_dp * PI_SH)
      norm = norm * factorial(l - abs(m)) / factorial(l + abs(m))
      A_lm(l, m) = sqrt(norm)
    end do
  end do
  
  coefficientsInitialized = .true.
  
end subroutine SH_Init

!================================================================================
! Cleanup allocated arrays
!================================================================================
subroutine SH_Cleanup()
  implicit none
  
  if (allocated(factorial)) deallocate(factorial)
  if (allocated(invFactorial)) deallocate(invFactorial)
  if (allocated(sqrtFact)) deallocate(sqrtFact)
  if (allocated(A_lm)) deallocate(A_lm)
  
  maxOrder = -1
  coefficientsInitialized = .false.
  
end subroutine SH_Cleanup

!================================================================================
! Get linear index for (l, m) pair
! Returns index in range 1 to (p+1)^2
!================================================================================
pure function SH_Index(l, m) result(idx)
  implicit none
  integer, intent(in) :: l, m
  integer :: idx
  
  idx = l*l + l + m + 1
  
end function SH_Index

!================================================================================
! Get total number of expansion coefficients for order p
!================================================================================
pure function SH_GetExpansionSize(p) result(nCoeff)
  implicit none
  integer, intent(in) :: p
  integer :: nCoeff
  
  nCoeff = (p + 1) * (p + 1)
  
end function SH_GetExpansionSize

!================================================================================
! Compute Associated Legendre polynomials P_l^m(x) for all l, m up to order p
! Uses stable recurrence relation
! Output: Plm(0:p, 0:p) where Plm(l,m) = P_l^m(x) for m >= 0
!         For m < 0: P_l^(-m) = (-1)^m * (l-m)!/(l+m)! * P_l^m
!================================================================================
subroutine SH_AssocLegendre(x, p, Plm)
  implicit none
  real(dp), intent(in) :: x        ! cos(theta), must be in [-1, 1]
  integer, intent(in) :: p         ! Maximum order
  real(dp), intent(out) :: Plm(0:p, 0:p)  ! P_l^m(x)
  
  integer :: l, m
  real(dp) :: somx2, fact
  
  Plm = 0.0_dp
  
  ! P_0^0 = 1
  Plm(0, 0) = 1.0_dp
  
  if (p == 0) return
  
  ! Compute sqrt(1 - x^2) = sin(theta)
  somx2 = sqrt(max(0.0_dp, 1.0_dp - x*x))
  
  ! Compute diagonal elements P_m^m using:
  ! P_m^m = -(2m-1) * sqrt(1-x^2) * P_{m-1}^{m-1}
  fact = 1.0_dp
  do m = 1, p
    Plm(m, m) = -fact * somx2 * Plm(m-1, m-1)
    fact = fact + 2.0_dp
  end do
  
  ! Compute P_{m+1}^m using:
  ! P_{m+1}^m = x * (2m+1) * P_m^m
  do m = 0, p-1
    Plm(m+1, m) = x * real(2*m + 1, dp) * Plm(m, m)
  end do
  
  ! Compute remaining elements using recurrence:
  ! (l-m) * P_l^m = x * (2l-1) * P_{l-1}^m - (l+m-1) * P_{l-2}^m
  do m = 0, p-2
    do l = m+2, p
      Plm(l, m) = (x * real(2*l - 1, dp) * Plm(l-1, m) - &
                  real(l + m - 1, dp) * Plm(l-2, m)) / real(l - m, dp)
    end do
  end do
  
end subroutine SH_AssocLegendre

!================================================================================
! Compute complex spherical harmonic Y_l^m(theta, phi)
! Y_l^m = A_l^m * P_l^|m|(cos(theta)) * exp(i*m*phi)
!================================================================================
subroutine SH_ComputeYlm(theta, phi, p, Ylm)
  implicit none
  real(dp), intent(in) :: theta, phi  ! Spherical angles
  integer, intent(in) :: p            ! Maximum order
  complex(dp), intent(out) :: Ylm((p+1)*(p+1))  ! Y_l^m values
  
  real(dp) :: Plm(0:p, 0:p)
  real(dp) :: cosTheta
  complex(dp) :: expPhi
  integer :: l, m, idx
  
  cosTheta = cos(theta)
  
  ! Compute Legendre polynomials
  call SH_AssocLegendre(cosTheta, p, Plm)
  
  ! Compute spherical harmonics
  do l = 0, p
    do m = -l, l
      idx = SH_Index(l, m)
      expPhi = exp(IMAG_UNIT * real(m, dp) * phi)
      
      if (m >= 0) then
        Ylm(idx) = A_lm(l, m) * Plm(l, m) * expPhi
      else
        ! Y_l^{-m} = (-1)^m * conj(Y_l^m)
        Ylm(idx) = A_lm(l, m) * Plm(l, abs(m)) * expPhi
        if (mod(abs(m), 2) == 1) Ylm(idx) = -Ylm(idx)
      end if
    end do
  end do
  
end subroutine SH_ComputeYlm

!================================================================================
! Convert Cartesian to spherical coordinates
!================================================================================
pure subroutine SH_CartesianToSpherical(x, y, z, r, theta, phi)
  implicit none
  real(dp), intent(in) :: x, y, z
  real(dp), intent(out) :: r, theta, phi
  
  real(dp) :: rxy
  
  r = sqrt(x*x + y*y + z*z)
  
  if (r < 1.0e-15_dp) then
    theta = 0.0_dp
    phi = 0.0_dp
    return
  end if
  
  theta = acos(z / r)
  
  rxy = sqrt(x*x + y*y)
  if (rxy < 1.0e-15_dp) then
    phi = 0.0_dp
  else
    phi = atan2(y, x)
  end if
  
end subroutine SH_CartesianToSpherical

!================================================================================
! Compute regular solid harmonic R_l^m(r, theta, phi) = r^l * Y_l^m(theta, phi)
! These are used for multipole moments
!================================================================================
subroutine SH_ComputeSolidHarmonic(x, y, z, p, Rlm)
  implicit none
  real(dp), intent(in) :: x, y, z    ! Cartesian coordinates
  integer, intent(in) :: p           ! Maximum order
  complex(dp), intent(out) :: Rlm((p+1)*(p+1))  ! Regular solid harmonics
  
  real(dp) :: r, theta, phi
  complex(dp) :: Ylm((p+1)*(p+1))
  real(dp) :: r_power
  integer :: l, m, idx
  
  call SH_CartesianToSpherical(x, y, z, r, theta, phi)
  call SH_ComputeYlm(theta, phi, p, Ylm)
  
  r_power = 1.0_dp
  do l = 0, p
    do m = -l, l
      idx = SH_Index(l, m)
      Rlm(idx) = r_power * Ylm(idx)
    end do
    r_power = r_power * r
  end do
  
end subroutine SH_ComputeSolidHarmonic

!================================================================================
! Compute irregular solid harmonic S_l^m = r^{-(l+1)} * Y_l^m(theta, phi)
! These are used for local expansions
!================================================================================
subroutine SH_ComputeIrregularSolid(x, y, z, p, Slm)
  implicit none
  real(dp), intent(in) :: x, y, z    ! Cartesian coordinates
  integer, intent(in) :: p           ! Maximum order
  complex(dp), intent(out) :: Slm((p+1)*(p+1))  ! Irregular solid harmonics
  
  real(dp) :: r, theta, phi
  complex(dp) :: Ylm((p+1)*(p+1))
  real(dp) :: r_inv_power
  integer :: l, m, idx
  
  call SH_CartesianToSpherical(x, y, z, r, theta, phi)
  
  if (r < 1.0e-15_dp) then
    Slm = (0.0_dp, 0.0_dp)
    return
  end if
  
  call SH_ComputeYlm(theta, phi, p, Ylm)
  
  r_inv_power = 1.0_dp / r
  do l = 0, p
    do m = -l, l
      idx = SH_Index(l, m)
      Slm(idx) = r_inv_power * Ylm(idx)
    end do
    r_inv_power = r_inv_power / r
  end do
  
end subroutine SH_ComputeIrregularSolid

!================================================================================
! M2M Translation: Translate multipole expansion from child to parent
!
! Given multipole moments M_child centered at origin, compute contribution
! to parent multipole M_parent centered at (dx, dy, dz)
!
! M_parent(j,k) += sum_{l,m} M2M_coeff(j,k,l,m) * M_child(l,m)
!================================================================================
subroutine SH_M2M_Coeff(dx, dy, dz, p, M_child, M_parent)
  implicit none
  real(dp), intent(in) :: dx, dy, dz        ! Translation vector (child to parent)
  integer, intent(in) :: p                   ! Maximum order
  complex(dp), intent(in) :: M_child((p+1)*(p+1))   ! Child multipole
  complex(dp), intent(inout) :: M_parent((p+1)*(p+1)) ! Parent multipole (accumulated)
  
  complex(dp) :: Rlm((p+1)*(p+1))
  integer :: j, k, l, m, idx_jk, idx_lm
  complex(dp) :: coeff
  real(dp) :: sign_factor
  
  ! Compute regular solid harmonics at translation point
  call SH_ComputeSolidHarmonic(dx, dy, dz, p, Rlm)
  
  ! M2M translation using addition theorem
  ! M_j^k(parent) = sum_{l=0}^{j} sum_{m=-l}^{l} 
  !                 C(j,k,l,m) * R_{j-l}^{k-m}(d) * M_l^m(child)
  do j = 0, p
    do k = -j, j
      idx_jk = SH_Index(j, k)
      
      do l = 0, j
        do m = max(-l, k-(j-l)), min(l, k+(j-l))
          if (abs(k-m) > j-l) cycle
          
          idx_lm = SH_Index(l, m)
          
          ! Translation coefficient with proper normalization
          coeff = ComputeM2MCoeff(j, k, l, m) * Rlm(SH_Index(j-l, k-m))
          M_parent(idx_jk) = M_parent(idx_jk) + coeff * M_child(idx_lm)
        end do
      end do
    end do
  end do
  
contains
  
  pure function ComputeM2MCoeff(j, k, l, m) result(c)
    integer, intent(in) :: j, k, l, m
    complex(dp) :: c
    real(dp) :: num, den
    integer :: jml, kmm
    
    jml = j - l
    kmm = k - m
    
    ! Coefficient from addition theorem
    num = factorial(j) * factorial(l)
    den = factorial(jml) * factorial(l - abs(m)) * factorial(jml - abs(kmm))
    
    if (den > 0.0_dp) then
      c = cmplx(sqrt(num/den), 0.0_dp, dp)
    else
      c = (0.0_dp, 0.0_dp)
    end if
    
    ! Sign correction
    if (mod(abs(m), 2) == 1 .and. m < 0) c = -c
    
  end function ComputeM2MCoeff
  
end subroutine SH_M2M_Coeff

!================================================================================
! M2L Translation: Convert multipole expansion to local expansion
!
! Given multipole moments M centered at origin, compute local expansion L
! centered at (dx, dy, dz) in the far field
!
! Uses the simplified addition theorem:
! L_j^k(target) = sum_{l=0}^{p} sum_{m=-l}^{l}
!                 (-1)^l * S_{j+l}^{m-k}(d) * M_l^m(source)
!
! where S is the irregular solid harmonic (1/r^{n+1} * Y_n^m)
!
! Note: The irregular solid harmonics already contain the spherical harmonic
! normalization, so no additional normalization is needed here.
!================================================================================
subroutine SH_M2L_Coeff(dx, dy, dz, p, M_source, L_target)
  implicit none
  real(dp), intent(in) :: dx, dy, dz        ! Translation vector (source to target)
  integer, intent(in) :: p                   ! Maximum order
  complex(dp), intent(in) :: M_source((p+1)*(p+1))   ! Source multipole
  complex(dp), intent(inout) :: L_target((p+1)*(p+1)) ! Target local (accumulated)
  
  complex(dp) :: Slm((2*p+1)*(2*p+1))
  integer :: j, k, l, m, idx_jk, idx_lm, idx_jpl
  complex(dp) :: coeff
  real(dp) :: sign_factor
  integer :: jpl, mmk
  
  ! Compute irregular solid harmonics at translation point
  ! Need order 2p for M2L
  call SH_ComputeIrregularSolid(dx, dy, dz, 2*p, Slm)
  
  ! M2L translation using addition theorem
  ! L_j^k(target) = sum_{l=0}^{p} sum_{m=-l}^{l}
  !                 (-1)^l * S_{j+l}^{m-k}(d) * M_l^m(source)
  do j = 0, p
    do k = -j, j
      idx_jk = SH_Index(j, k)
      
      do l = 0, p
        do m = -l, l
          jpl = j + l
          mmk = m - k  ! Note: m-k, not k-m
          if (abs(mmk) > jpl) cycle
          
          idx_lm = SH_Index(l, m)
          idx_jpl = SH_Index(jpl, mmk)
          
          ! Sign factor: (-1)^l
          sign_factor = 1.0_dp
          if (mod(l, 2) == 1) sign_factor = -1.0_dp
          
          coeff = cmplx(sign_factor, 0.0_dp, dp) * Slm(idx_jpl)
          L_target(idx_jk) = L_target(idx_jk) + coeff * M_source(idx_lm)
        end do
      end do
    end do
  end do
  
end subroutine SH_M2L_Coeff

!================================================================================
! L2L Translation: Translate local expansion from parent to child
!
! Given local expansion L_parent centered at origin, compute contribution
! to child local L_child centered at (dx, dy, dz)
!
! L_child(j,k) = sum_{l,m} L2L_coeff(j,k,l,m) * L_parent(l,m)
!================================================================================
subroutine SH_L2L_Coeff(dx, dy, dz, p, L_parent, L_child)
  implicit none
  real(dp), intent(in) :: dx, dy, dz        ! Translation vector (parent to child)
  integer, intent(in) :: p                   ! Maximum order
  complex(dp), intent(in) :: L_parent((p+1)*(p+1))   ! Parent local
  complex(dp), intent(inout) :: L_child((p+1)*(p+1)) ! Child local (accumulated)
  
  complex(dp) :: Rlm((p+1)*(p+1))
  integer :: j, k, l, m, idx_jk, idx_lm
  complex(dp) :: coeff
  
  ! Compute regular solid harmonics at translation point
  call SH_ComputeSolidHarmonic(dx, dy, dz, p, Rlm)
  
  ! L2L translation using addition theorem
  ! L_j^k(child) = sum_{l=j}^{p} sum_{m} 
  !                C(j,k,l,m) * R_{l-j}^{m-k}(d) * L_l^m(parent)
  do j = 0, p
    do k = -j, j
      idx_jk = SH_Index(j, k)
      
      do l = j, p
        do m = max(-l, k-(l-j)), min(l, k+(l-j))
          if (abs(m-k) > l-j) cycle
          
          idx_lm = SH_Index(l, m)
          
          ! Translation coefficient
          coeff = ComputeL2LCoeff(j, k, l, m) * Rlm(SH_Index(l-j, m-k))
          L_child(idx_jk) = L_child(idx_jk) + coeff * L_parent(idx_lm)
        end do
      end do
    end do
  end do
  
contains
  
  pure function ComputeL2LCoeff(j, k, l, m) result(c)
    integer, intent(in) :: j, k, l, m
    complex(dp) :: c
    real(dp) :: num, den
    integer :: lmj, mmk
    
    lmj = l - j
    mmk = m - k
    
    ! Coefficient from addition theorem
    num = factorial(l) * factorial(j)
    den = factorial(lmj) * factorial(j - abs(k)) * factorial(lmj - abs(mmk))
    
    if (den > 0.0_dp) then
      c = cmplx(sqrt(num/den), 0.0_dp, dp)
    else
      c = (0.0_dp, 0.0_dp)
    end if
    
  end function ComputeL2LCoeff
  
end subroutine SH_L2L_Coeff

!================================================================================
end module SphericalHarmonics
!================================================================================
