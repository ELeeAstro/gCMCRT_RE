module mc_scat_mod
  use mc_precision_mod, only : dp
  use mc_data_mod
  use curand_device
  use cudafor
  implicit none

contains

  attributes(device) subroutine scatt(ph)
    implicit none

    type(pac), intent(inout) :: ph


    select case(ph%iscat)
    case(1)

      ! Isotropic scattering
      ph%cost = 2.0_dp * curand_uniform(ph%iseed) - 1.0_dp

    case default
      print*, 'Invalid scattering phase function number: ', ph%iscat
      stop
    end select

    ! Give nzp the cost value
    ph%nzp = ph%cost

  end subroutine scatt

end module mc_scat_mod
