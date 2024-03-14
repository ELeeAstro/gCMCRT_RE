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

    if (ph%inc == 1) then
      istat = atomicadd(Jdot_s_d(ph%zc), d * ph%w * ph%e0dt)
      istat = atomicadd(Hdot_s_d(ph%zc), d*ph%nzp * ph%w * ph%e0dt)
      istat = atomicadd(Kdot_s_d(ph%zc), d*ph%nzp**2 * ph%w * ph%e0dt)
      istat = atomicadd(Adot_s_d(ph%zc), rho_d(ph%zc) * k_abs_d(ph%b,ph%zc)*d * ph%w * ph%e0dt)
    else
      istat = atomicadd(Jdot_d(ph%zc), d * ph%w * ph%e0dt)
      istat = atomicadd(Hdot_d(ph%zc), d*ph%nzp * ph%w * ph%e0dt)
      istat = atomicadd(Kdot_d(ph%zc), d*ph%nzp**2 * ph%w * ph%e0dt)
      istat = atomicadd(Adot_d(ph%zc), rho_d(ph%zc) * k_abs_d(ph%b,ph%zc)*d * ph%w * ph%e0dt)
    end if



  end subroutine estimators

end module mc_moments_mod
