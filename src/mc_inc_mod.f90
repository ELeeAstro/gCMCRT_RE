module mc_inc_mod
  use mc_precision_mod, only : dp
  use mc_data_mod
  use mc_aux_mod
  implicit none

contains

  attributes(device) subroutine inc_stellar(ph)
    implicit none

    type(pac), intent(inout) :: ph

    ph%cost = -mu_z_d
    ph%nzp = ph%cost
    ph%zc = nlay_d
    ph%zp = z_d(nlev_d) - 1.0e-6_dp

  end subroutine inc_stellar

end module mc_inc_mod
