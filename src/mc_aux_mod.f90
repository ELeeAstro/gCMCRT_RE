module mc_aux_mod
  use mc_precision_mod, only : dp
  implicit none

contains

  subroutine hypsometric(nlay,Rd_air,grav,Tl,pe,z)
    implicit none

    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) :: Rd_air, grav, Tl
    real(dp), dimension(nlay+1), intent(in) :: pe

    real(dp), dimension(nlay+1), intent(out) :: z

    integer :: i

    z(1) = 0.0_dp
    do i = 1, nlay
      z(i+1) = z(i) + Rd_air(i)/grav(i) * Tl(i) * log(pe(i)/pe(i+1))
    end do

  end subroutine hypsometric

end module mc_aux_mod
