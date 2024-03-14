module mc_emit_mod
  use mc_precision_mod, only : dp
  use mc_data_mod
  use curand_device
  use cudafor
  implicit none


contains

  attributes(device) subroutine emit_iso(ph)
    implicit none

    type(pac), intent(inout) :: ph

    ph%cost = 2.0_dp*curand_uniform(ph%iseed) - 1.0_dp
    ph%nzp = ph%cost

  end subroutine emit_iso

  attributes(device) subroutine emit_iso_surf(ph)
    implicit none

    type(pac), intent(inout) :: ph

    ph%cost = sqrt(curand_uniform(ph%iseed))
    ph%nzp = ph%cost

    ph%zp = 1.0e-6_dp
    ph%zc = 1

  end subroutine emit_iso_surf

end module mc_emit_mod
