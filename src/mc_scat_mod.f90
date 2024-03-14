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

    real(dp) :: bmu, costp, sintp, sinbt, ri1, cosi1, cosi3
    real(dp) :: q, u, hgg, g1, g2
    real(dp) :: gf1, gf2, gb1, gb2, alph
    real(dp) :: zeta1, zeta2

    select case(ph%iscat)
    case(1)

      ! Isotropic scattering
      ph%cost = 2.0_dp * curand_uniform(ph%iseed) - 1.0_dp
      ph%nzp = ph%cost
      return

    case(2)
      ! Rayleigh scattering via direct spherical coordinate sampling
      ! Assumes non-polarised incident packet
      q = 4.0_dp*curand_uniform(ph%iseed) - 2.0_dp
      u = (-q + sqrt(1.0_dp + q**2))**(third)
      bmu = u - 1.0_dp/u

    case(3)
      ! Sample from single HG function
      if (g_d(ph%b,ph%zc) /= 0.0_dp) then
        hgg = g_d(ph%b,ph%zc)
        g2 = hgg**2

        bmu = ((1.0_dp + g2) - &
          & ((1.0_dp - g2) / (1.0_dp - hgg + 2.0_dp * hgg * curand_uniform(ph%iseed)))**2) &
          & / (2.0_dp*hgg)
      else
        ! Isotropic scattering
        ph%cost = 2.0_dp * curand_uniform(ph%iseed) - 1.0_dp
        ph%nzp = ph%cost
        return
      end if

    case default
      print*, 'Invalid scattering phase function number: ', ph%iscat
      stop
    end select

    !! Change direction of packet given by sampled direction

    ! Avoid rare numerical issues with sampling bmu
    if (bmu > 1.0_dp) then
      bmu = 1.0_dp
    else if (bmu < -1.0_dp) then
      bmu = -1.0_dp
    end if

    !! Calculate change in direction in grid reference frame
    if (bmu >= 1.0_dp) then
      !! Packet directly forward scatters - no change in direction
      ! ph%cost = ph%cost
    else if (bmu <= -1.0_dp) then
      !! Packet directly backward scatters - negative sign for cost
      ph%cost = -ph%cost
    else
      !! Packet scatters according to sampled cosine and current direction
      ! Save current cosine direction of packet and calculate sine
      costp = ph%cost
      sintp = 1.0_dp - costp**2
      if (sintp <= 0.0_dp)then
        sintp = 0.0_dp
      else
        sintp = sqrt(sintp)
      endif

      ! Randomly decide if scattered in +/- quadrant direction and find new cosine direction
      sinbt = sqrt(1.0_dp - bmu**2)
      ri1 = twopi * curand_uniform(ph%iseed)
      if (ri1 > pi) then
        cosi3 = cos(twopi - ri1)
        ! Calculate new cosine direction
        ph%cost = costp * bmu + sintp * sinbt * cosi3
      else !(ri1 <= pi)
        cosi1 = cos(ri1)
        ph%cost = costp * bmu + sintp * sinbt * cosi1
      end if
    end if

    ! Give nzp the cost value
    ph%nzp = ph%cost

  end subroutine scatt

end module mc_scat_mod
