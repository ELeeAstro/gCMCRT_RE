module mc_tauint_mod
  use mc_precision_mod, only : dp
  use mc_data_mod
  use mc_moments_mod
  use cudafor
  implicit none

contains

  attributes(device) subroutine tauint_1D_pp(ph)
    implicit none

    type(pac), intent(inout) :: ph

    real(dp) :: dsz, taucell, d1
    integer :: zoffset, i

    !! Begin the tau integration
    ph%tau = 0.0_dp

    do while (ph%tau < ph%tau_p)

      !! Calculate dsz, the distance to the next vertical level
      if (ph%nzp > 0.0_dp) then
        ! Packet travelling upward, find distance to upper level
        dsz = (z_d(ph%zc+1)-ph%zp)/ph%nzp
        zoffset = 1
      else if (ph%nzp < 0.0_dp) then
        ! Packet travelling downward, find distance to lower level
        dsz = (z_d(ph%zc)-ph%zp)/ph%nzp
        zoffset = -1
      else
        ! Packet travelling directly in z plane
        ! Return, packet does not move in z direction
        return
      end if

      !! Calculate optical depth to level
      taucell = dsz * rhokap_d(ph%b,ph%zc)

      !! Check if packet ends path in this layer
      if ((ph%tau + taucell) >= ph%tau_p) then
        ! Packet stops in this cell - move distance then exit loop
        d1 = (ph%tau_p-ph%tau)/rhokap_d(ph%b,ph%zc)
        ph%zp = ph%zp + d1 * ph%nzp

        call estimators(d1,ph)

        ph%tau = ph%tau_p
      else
        ! Packet continues to level edge - update position, cell index and tau
        ph%zp = ph%zp + (dsz + 1.0e-6_dp) * ph%nzp

        call estimators(dsz,ph)

        ph%zc = ph%zc + zoffset

        ph%tau = ph%tau + taucell

      end if


      ! Check is packet has exited the domain
      if ((ph%zc > nlay_d) .or. (ph%zp >= z_d(nlev_d))) then
        ph%flag = 1
        exit
      else if ((ph%zc < 1) .or. (ph%zp <= z_d(1))) then
        ph%flag = -2
        exit
      end if

    end do
 
  end subroutine tauint_1D_pp

end module mc_tauint_mod
