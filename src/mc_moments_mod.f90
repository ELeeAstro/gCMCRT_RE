module mc_moments_mod
  use mc_precision_mod
  use mc_data_mod
  implicit none

contains

  attributes(device) subroutine estimators(d,ph)
    implicit none

    real(dp), intent(in) :: d
    type(pac), intent(in) :: ph

    integer :: istat


    istat = atomicadd(Jdot_d(ph%zc), d * ph%w * ph%e0dt)
    istat = atomicadd(Hdot_d(ph%zc), d*ph%nzp * ph%w * ph%e0dt)
    istat = atomicadd(Kdot_d(ph%zc), d*ph%nzp**2 * ph%w * ph%e0dt)

    istat = atomicadd(Adot_d(ph%zc), rho_d(ph%zc) * k_abs_d(ph%b,ph%zc)*d * ph%w * ph%e0dt)

  end subroutine estimators

end module mc_moments_mod
