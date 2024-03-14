module mc_data_mod
  use mc_precision_mod, only : dp
  use curand_device
  implicit none

  type pac
    integer(8) :: id                          ! Photon number unique id
    integer(8) :: seed                        ! Unique random seed of packet
    type(curandStateMRG32k3a) :: iseed
    real(dp) :: wl                         ! Wavelength of packet

    real(dp) :: cost, sint, cosp, sinp, phi
    real(dp) :: nzp, nxp, nyp
    real(dp) :: zp, xp, yp
    real(dp) :: w
    real(dp) :: tau, tau_p
    real(dp) :: e0dt
    integer :: zc, iscat, b, inc
    integer :: flag

    integer :: nscat

  end type pac


  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi, fourpi = 4.0_dp * pi
  real(dp), parameter :: third = 1.0_dp/3.0_dp
  real(dp), parameter :: sb_c = 5.670374419e-5_dp

  !! Host moments and absorption
  real(dp), allocatable, dimension(:) :: Jdot, Hdot, Kdot, Jdot_s, Hdot_s, Kdot_s
  real(dp), allocatable, dimension(:) ::  Adot, Adot_s

  !! Device moments and absorption
  real(dp), allocatable, dimension(:), device :: Jdot_d, Hdot_d, Kdot_d, Jdot_s_d, Hdot_s_d, Kdot_s_d
  real(dp), allocatable, dimension(:), device :: Adot_d, Adot_s_d

  !! Device grid variables
  integer, device :: nlay_d, nlev_d, nb_d, iscat_d
  real(dp), device :: mu_z_d
  real(dp), allocatable, dimension(:), device :: z_d, rho_d
  real(dp), allocatable, dimension(:,:), device :: rhokap_d, alb_d, g_d, k_abs_d 

  !! OLR measurements
  real(dp) :: OLR
  real(dp), device :: OLR_d


end module mc_data_mod
