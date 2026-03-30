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
! Compute regular solid harmonic (unnormalized convention):
!   R_l^m(r,theta,phi) = r^l * P_l^|m|(cos theta) * exp(i*m*phi) / (l+|m|)!
!
! This convention gives the simplest translation operators:
!   R_n^m(a+b) = sum_j sum_k R_j^k(a) * R_{n-j}^{m-k}(b)
! with NO additional coefficients needed.
!
! Properties: conj(R_l^m) = R_l^{-m}
!================================================================================
subroutine SH_ComputeSolidHarmonic(x, y, z, p, Rlm)
  implicit none
  real(dp), intent(in) :: x, y, z    ! Cartesian coordinates
  integer, intent(in) :: p           ! Maximum order
  complex(dp), intent(out) :: Rlm((p+1)*(p+1))  ! Regular solid harmonics
  
  real(dp) :: r, theta, phi
  real(dp) :: Plm(0:p, 0:p)
  real(dp) :: cosTheta
  real(dp) :: r_power
  complex(dp) :: expPhi
  integer :: l, m, idx
  
  call SH_CartesianToSpherical(x, y, z, r, theta, phi)
  
  Rlm = (0.0_dp, 0.0_dp)
  
  if (r < 1.0e-15_dp) then
    ! At origin, only R_0^0 = 1/0! = 1 is nonzero
    Rlm(SH_Index(0,0)) = (1.0_dp, 0.0_dp)
    return
  endif
  
  cosTheta = cos(theta)
  call SH_AssocLegendre(cosTheta, p, Plm)
  
  r_power = 1.0_dp
  do l = 0, p
    do m = -l, l
      idx = SH_Index(l, m)
      expPhi = exp(IMAG_UNIT * real(m, dp) * phi)
      ! R_l^m = r^l * P_l^|m|(cos theta) * exp(i*m*phi) / (l + |m|)!
      Rlm(idx) = r_power * Plm(l, abs(m)) * expPhi * invFactorial(l + abs(m))
    end do
    r_power = r_power * r
  end do
  
end subroutine SH_ComputeSolidHarmonic

!================================================================================
! Compute irregular solid harmonic (unnormalized convention):
!   S_l^m(r,theta,phi) = (l-|m|)! * P_l^|m|(cos theta) * exp(i*m*phi) / r^{l+1}
!
! The addition theorem for the Laplace kernel states:
!   1/|a-b| = sum_{l,m} (-1)^m * R_l^{-m}(b) * S_l^m(a)   for |b| < |a|
!
! Properties: conj(S_l^m) = S_l^{-m}
!================================================================================
subroutine SH_ComputeIrregularSolid(x, y, z, p, Slm)
  implicit none
  real(dp), intent(in) :: x, y, z    ! Cartesian coordinates
  integer, intent(in) :: p           ! Maximum order
  complex(dp), intent(out) :: Slm((p+1)*(p+1))  ! Irregular solid harmonics
  
  real(dp) :: r, theta, phi
  real(dp) :: Plm(0:p, 0:p)
  real(dp) :: cosTheta
  real(dp) :: r_inv_power
  complex(dp) :: expPhi
  integer :: l, m, idx
  
  call SH_CartesianToSpherical(x, y, z, r, theta, phi)
  
  Slm = (0.0_dp, 0.0_dp)
  if (r < 1.0e-15_dp) return
  
  cosTheta = cos(theta)
  call SH_AssocLegendre(cosTheta, p, Plm)
  
  r_inv_power = 1.0_dp / r
  do l = 0, p
    do m = -l, l
      idx = SH_Index(l, m)
      expPhi = exp(IMAG_UNIT * real(m, dp) * phi)
      ! S_l^m = (l-|m|)! * P_l^|m|(cos theta) * exp(i*m*phi) / r^{l+1}
      Slm(idx) = factorial(l - abs(m)) * Plm(l, abs(m)) * expPhi * r_inv_power
    end do
    r_inv_power = r_inv_power / r
  end do
  
end subroutine SH_ComputeIrregularSolid

!================================================================================
! M2M Translation: Translate multipole expansion from child to parent
!
! Given multipole moments M_child centered at child, compute contribution
! to parent multipole M_parent centered at parent.
! Translation vector t = (dx,dy,dz) = child_center - parent_center.
!
! Formula (derived from regular solid harmonic addition theorem):
!   M'_n^m += sum_{j=0}^{n} sum_k M_j^k * conj(R_{n-j}^{m-k}(t))
!
! Since conj(R_l^m) = R_l^{-m} for our convention:
!   M'_n^m += sum_{j=0}^{n} sum_k M_j^k * R_{n-j}^{k-m}(t)
!
! Constraint: |k-m| <= n-j
!================================================================================
subroutine SH_M2M_Coeff(dx, dy, dz, p, M_child, M_parent)
  implicit none
  real(dp), intent(in) :: dx, dy, dz        ! Translation vector (child to parent)
  integer, intent(in) :: p                   ! Maximum order
  complex(dp), intent(in) :: M_child((p+1)*(p+1))   ! Child multipole
  complex(dp), intent(inout) :: M_parent((p+1)*(p+1)) ! Parent multipole (accumulated)
  
  complex(dp) :: Rlm((p+1)*(p+1))
  integer :: n, m_n, j, k, nj, km
  integer :: idx_nm, idx_jk, idx_tr
  
  ! Compute regular solid harmonics at translation point
  call SH_ComputeSolidHarmonic(dx, dy, dz, p, Rlm)
  
  ! M2M translation
  do n = 0, p
    do m_n = -n, n
      idx_nm = SH_Index(n, m_n)
      do j = 0, n
        do k = -j, j
          nj = n - j
          km = k - m_n
          if (abs(km) > nj) cycle
          idx_jk = SH_Index(j, k)
          idx_tr = SH_Index(nj, km)
          M_parent(idx_nm) = M_parent(idx_nm) + M_child(idx_jk) * Rlm(idx_tr)
        end do
      end do
    end do
  end do
  
end subroutine SH_M2M_Coeff

!================================================================================
! M2L Translation: Convert multipole expansion to local expansion
!
! Given multipole moments M at source center, compute local expansion L
! at target center. Translation vector d = (dx,dy,dz) = target - source.
!
! Derived from the re-expansion of irregular solid harmonics:
!   S_n^m(delta + d) = sum_{j,k} (-1)^{j-k} * S_{n+j}^{m-k}(d) * R_j^k(delta)
!
! Therefore:
!   L_j^k += sum_{n=0}^{p} sum_{m=-n}^{n} M_n^m * (-1)^{j-k} * S_{n+j}^{m-k}(d)
!
! Requires S up to order 2p.
! Verified: monopole gives L_j^k = (-1)^{j-k} S_j^{-k}(d), consistent with
! Taylor expansion of 1/|delta+d|.
!================================================================================
subroutine SH_M2L_Coeff(dx, dy, dz, p, M_source, L_target)
  implicit none
  real(dp), intent(in) :: dx, dy, dz        ! Translation vector (source to target)
  integer, intent(in) :: p                   ! Maximum order
  complex(dp), intent(in) :: M_source((p+1)*(p+1))   ! Source multipole
  complex(dp), intent(inout) :: L_target((p+1)*(p+1)) ! Target local (accumulated)
  
  complex(dp) :: Slm((2*p+1)*(2*p+1))
  integer :: j, k, n, m_n, njp, mmk
  integer :: idx_jk, idx_nm, idx_s
  real(dp) :: sign_factor
  
  ! Compute irregular solid harmonics at translation point
  ! Need order 2p for M2L since n+j can be up to 2p
  call SH_ComputeIrregularSolid(dx, dy, dz, 2*p, Slm)
  
  ! M2L translation
  do j = 0, p
    do k = -j, j
      idx_jk = SH_Index(j, k)
      
      ! Sign factor: (-1)^{j-k}
      if (mod(abs(j - k), 2) == 0) then
        sign_factor = 1.0_dp
      else
        sign_factor = -1.0_dp
      endif
      
      do n = 0, p
        do m_n = -n, n
          njp = n + j
          mmk = m_n - k
          if (abs(mmk) > njp) cycle
          if (njp > 2*p) cycle
          
          idx_nm = SH_Index(n, m_n)
          idx_s = SH_Index(njp, mmk)
          
          L_target(idx_jk) = L_target(idx_jk) + &
            cmplx(sign_factor, 0.0_dp, dp) * M_source(idx_nm) * Slm(idx_s)
        end do
      end do
    end do
  end do
  
end subroutine SH_M2L_Coeff

!================================================================================
! L2L Translation: Translate local expansion from parent to child
!
! Given local expansion L_parent centered at parent, compute contribution
! to child local L_child centered at child.
! Translation vector tau = (dx,dy,dz) = child_center - parent_center.
!
! Derived from regular solid harmonic addition theorem:
!   R_n^m(delta' + tau) = sum_{j=0}^{n} sum_k R_j^k(delta') * R_{n-j}^{m-k}(tau)
!
! Therefore:
!   L'_j^k = sum_{n=j}^{p} sum_m L_n^m * R_{n-j}^{m-k}(tau)
!
! Constraint: |m-k| <= n-j
! Verified: correctly reproduces Taylor expansion shift for p=2 on-axis.
!================================================================================
subroutine SH_L2L_Coeff(dx, dy, dz, p, L_parent, L_child)
  implicit none
  real(dp), intent(in) :: dx, dy, dz        ! Translation vector (parent to child)
  integer, intent(in) :: p                   ! Maximum order
  complex(dp), intent(in) :: L_parent((p+1)*(p+1))   ! Parent local
  complex(dp), intent(inout) :: L_child((p+1)*(p+1)) ! Child local (accumulated)
  
  complex(dp) :: Rlm((p+1)*(p+1))
  integer :: j, k, n, m_n, nj, mk
  integer :: idx_jk, idx_nm, idx_r
  
  ! Compute regular solid harmonics at translation point
  call SH_ComputeSolidHarmonic(dx, dy, dz, p, Rlm)
  
  ! L2L translation
  do j = 0, p
    do k = -j, j
      idx_jk = SH_Index(j, k)
      do n = j, p
        do m_n = -n, n
          nj = n - j
          mk = m_n - k
          if (abs(mk) > nj) cycle
          idx_nm = SH_Index(n, m_n)
          idx_r = SH_Index(nj, mk)
          L_child(idx_jk) = L_child(idx_jk) + L_parent(idx_nm) * Rlm(idx_r)
        end do
      end do
    end do
  end do
  
end subroutine SH_L2L_Coeff

!================================================================================
end module SphericalHarmonics
!================================================================================
