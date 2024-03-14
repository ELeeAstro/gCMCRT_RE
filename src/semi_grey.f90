module gCMCRT_semi_grey
  use mc_precision_mod, only : dp
  use mc_data_mod
  use mc_tauint_mod
  use mc_inc_mod
  use mc_emit_mod
  use mc_scat_mod
  use curand_device
  use cudafor
  implicit none


  integer, device :: Nph_d, Nph_irr_d, Nph_int_d
  real(dp), device :: Fint_d, Firr_d
  type(curandStateMRG32k3a), allocatable, dimension(:), device :: iseed

contains

  attributes(global) subroutine set_iseed(Nph)
    implicit none

    integer, intent(in) :: Nph
    integer(8) :: id, seed
    integer :: seq, offset

    id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (id > Nph) then
      return
    end if

    seed = id + id**2 + id/2
    seq = 0
    offset = 0
    call curand_init(seed, seq, offset, iseed(id))


  end subroutine set_iseed

  attributes(global) subroutine pp_semi_grey_irr_kernel(Nph)
    implicit none

    integer, intent(in) :: Nph
    type(pac) :: ph
    integer :: istat

    ph%id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (ph%id > Nph) then
      return
    end if
    ph%iseed = iseed(ph%id)

    ph%w = 1.0_dp
    ph%e0dt = Firr_d/real(Nph,dp)

    ph%b = 1

    ph%inc = 1

    call inc_stellar(ph)

    ph%flag = 0
    ph%nscat = 0
    ph%iscat = iscat_d

    do while (ph%flag == 0)

      ph%tau_p = -log(curand_uniform(ph%iseed))

      call tauint_1D_pp(ph)

      if (ph%flag == -2) then
        call emit_iso_surf(ph)
        ph%b = 2
        ph%flag = 0
        ph%inc = 0
        cycle
      end if
      if (ph%flag == 1) then
        exit
      end if
      if (curand_uniform(ph%iseed) < alb_d(ph%b,ph%zc)) then
        call scatt(ph)
        ph%nscat = ph%nscat + 1
      else
        call emit_iso(ph)
        ph%b = 2
        ph%inc = 0
      end if

    end do

    if (ph%flag == 1) then
      !! Add energy to OLR
      istat = atomicadd(OLR_d, ph%e0dt)
    end if

    ! Give back iseed to saved device array for next iteration with this ph%id
    iseed(ph%id) = ph%iseed

  end subroutine pp_semi_grey_irr_kernel
  
  attributes(global) subroutine pp_semi_grey_int_kernel(Nph)
    implicit none

    integer, intent(in) :: Nph
    type(pac) :: ph
    integer :: istat

    ph%id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (ph%id > Nph) then
      return
    end if
    ph%iseed = iseed(ph%id)

    ph%w = 1.0_dp
    ph%e0dt = Fint_d/real(Nph,dp)

    ph%b = 2

    call emit_iso_surf(ph)

    ph%flag = 0
    ph%nscat = 0
    ph%iscat = iscat_d

    ph%inc = 0

    do while (ph%flag == 0)

      ph%tau_p = -log(curand_uniform(ph%iseed))

      call tauint_1D_pp(ph)

      if (ph%flag == -2) then
        call emit_iso_surf(ph)
        ph%flag = 0
        cycle
      end if
      if (ph%flag == 1) then
        cycle
      end if
      if (curand_uniform(ph%iseed) < alb_d(ph%b,ph%zc)) then
        call scatt(ph)
        ph%nscat = ph%nscat + 1
      else
        call emit_iso(ph)
      end if

    end do

    if (ph%flag == 1) then
      !! Add energy to OLR
      istat = atomicadd(OLR_d, ph%e0dt)
    end if      

    ! Give back iseed to saved device array for next iteration with this ph%id
    iseed(ph%id) = ph%iseed

  end subroutine pp_semi_grey_int_kernel


end module gCMCRT_semi_grey

subroutine semi_grey()
  use mc_precision_mod, only : dp
  use mc_aux_mod
  use mc_data_mod
  use gCMCRT_semi_grey
  use cudafor
  implicit none

  integer :: u_nml
  integer :: Nph_irr, Nph_int, Nit, last_int_fac, nlay, iscat
  real(dp) :: Tirr, Tint, k_sw, k_lw, p_top, p_bot, T_IC, grav_const, Rd_gas, lw_a, lw_g, sw_a, sw_g, mu_z

  integer :: nlev, i, nb, n, u0, u3
  real(dp), allocatable, dimension(:) :: pe, dpe, pl, k_Ross, k_P, Tl, z, dze, Rd_air, grav, dTl
  real(dp), allocatable, dimension(:) :: rho, tau_R
  real(dp), allocatable, dimension(:,:) :: k_abs, k_sca, k_ext, alb, g, rhokap, tau
  real(dp) :: Fint, Firr

  integer :: istat, Nph
  type(dim3) :: blocks, threads

  namelist /semi_grey_nml/ Nph_irr, Nph_int, Nit, last_int_fac, nlay, iscat, &
    & Tirr, Tint, k_sw, k_lw, p_top, p_bot, T_IC, grav_const, Rd_gas, lw_a, lw_g, sw_a, sw_g, mu_z

  !! Read input variables from namelist
  open(newunit=u_nml, file='gCMCRT_RE.nml', status='old', action='read')
  read(u_nml, nml=semi_grey_nml)
  close(u_nml)

  nlev = nlay + 1

  Nph = max(Nph_irr, Nph_int)

  threads = dim3(128, 1, 1)
  blocks = dim3(ceiling(real(Nph,dp)/threads%x),1,1)
  allocate(iseed(Nph))
  Nph_d = Nph
  call set_iseed<<<blocks, threads>>>(Nph_d)


  !! Set up vertical grid

  p_bot = p_bot * 1e6_dp
  p_top = p_top * 1e6_dp

  ! Do log spacing of pressure - note, index 1 = bottom of atmosphere
  allocate(pe(nlev),pl(nlay),dpe(nlay))
  do i = 1, nlev
    pe(i) = 10.0_dp**(log10(p_bot)+real(i-1,dp)/real(nlev-1,dp)*log10(p_top/p_bot))
    !print*, i, pe(i)/1e6_dp
  end do

  do i = 1, nlay
    dpe(i) = pe(i+1) - pe(i)
    pl(i) = dpe(i) / log(pe(i+1)/pe(i))
  end do

  ! Constant grey opacity - semi-grey two bands
  nb = 2
  allocate(k_abs(nb,nlay), k_sca(nb,nlay), k_ext(nb,nlay), k_Ross(nlay), k_P(nlay))
  k_ext(1,:) = k_sw
  k_ext(2,:) = k_lw

  allocate(alb(nb,nlay),g(nb,nlay))
  alb(1,:) = sw_a
  alb(2,:) = lw_a
  g(1,:) = sw_g
  g(2,:) = lw_g

  ! Find the scattering and absorption opacity
  k_sca(:,:) = k_ext(:,:) * alb(:,:)
  k_abs(:,:) = k_ext(:,:) * (1.0_dp - alb(:,:))

  k_Ross(:) = k_abs(2,:)
  k_P(:) = k_Ross(:)

  ! Find initial temperature profile
  allocate(Tl(nlay), dTl(nlay))
  dTl(:) = 0.0_dp
  Tl(:) = mu_z*Tirr**4 + Tint**4
  Tl(:) = Tl(:)**(1.0_dp/4.0_dp)

  allocate(Rd_air(nlay),grav(nlay))
  Rd_air(:) = Rd_gas
  grav(:) = grav_const

  ! Find altitude grid
  allocate(z(nlev))
  call hypsometric(nlay,Rd_air,grav,Tl,pe,z)

  allocate(dze(nlay))
  do i = 1, nlay
    dze(i) = z(i+1) - z(i)
  end do

  !! Find grid opacity
  allocate(rho(nlay), rhokap(nb,nlay))
  allocate(tau(nb,nlev), tau_R(nlev))
  tau(:,nlev) = 0.0_dp
  tau_R(nlev) = 0.0_dp
  do i = nlay, 1, -1
    rho(i) = pl(i)/(Rd_air(i) * Tl(i))
    rhokap(:,i) = rho(i) * k_ext(:,i)
    tau(:,i) = tau(:,i+1) + rhokap(:,i) * dze(i)
    tau_R(i) = tau_R(i+1) + rho(i)*k_Ross(i)* dze(i)
    !print*, i, rho(i), rhokap(:,i), tau(:,i)
  end do

  !! Output IC
  open(newunit=u0,file='results/Tp_semi_grey_1D_pp.txt',action='readwrite')
  write(u0,*) Nit, nlay
  write(u0,*) 0
  do i = 1, nlay
    write(u0,*) i, pl(i)/1e6_dp, (tau_R(i)+tau_R(i+1))/2.0_dp, Tl(i)
  end do

  !! Begin RT calculations and RE calculation

  ! Allocate moments and absorption arrays
  allocate(Jdot(nlay),Hdot(nlay),Kdot(nlay))
  allocate(Jdot_d(nlay),Hdot_d(nlay),Kdot_d(nlay))
  allocate(Jdot_s(nlay),Hdot_s(nlay),Kdot_s(nlay))
  allocate(Jdot_s_d(nlay),Hdot_s_d(nlay),Kdot_s_d(nlay))
  allocate(Adot(nlay), Adot_d(nlay), Adot_s(nlay), Adot_s_d(nlay))

  !! Send constant grid variables to device and allocate needed device arrays
  nlay_d = nlay
  nlev_d = nlev
  nb_d = nb
  Nph_irr_d = Nph_irr
  Nph_int_d = Nph_int
  mu_z_d = mu_z
  iscat_d = iscat
  allocate(z_d(nlev),rhokap_d(nb,nlay),alb_d(nb,nlay),rho_d(nlay),k_abs_d(nb,nlay),g_d(nb,nlay))
  
  do n = 1, Nit

    print*, n, Nit

    ! Irradiation flux onto planet
    Firr = mu_z * sb_c * Tirr**4
    ! Internal flux of the planet
    Fint = sb_c * Tint**4

    ! Recalculate height grid
    call hypsometric(nlay,Rd_air,grav,Tl,pe,z)
    do i = 1, nlay
      dze(i) = z(i+1) - z(i)
    end do

    ! Calculate new optical depths and rhokap
    tau(:,nlev) = 0.0_dp
    tau_R(nlev) = 0.0_dp
    do i = nlay, 1, -1
      rho(i) = pl(i)/(Rd_air(i) * Tl(i))
      rhokap(:,i) = rho(i) * k_ext(:,i)
      tau(:,i) = tau(:,i+1) + rhokap(:,i) * dze(i)
      tau_R(i) = tau_R(i+1) + rho(i)*k_Ross(i)* dze(i)
      !print*, i, rho(i), rhokap(i), tau(i)
    end do

    ! Zero moment and absorption arrays for device
    Jdot_d(:) = 0.0_dp ; Hdot_d(:) = 0.0_dp ; Kdot_d(:) = 0.0_dp
    Jdot_s_d(:) = 0.0_dp; Hdot_s_d(:) = 0.0_dp;  Kdot_s_d(:) = 0.0_dp
    Adot_d(:) = 0.0_dp ; Adot_s_d(:) = 0.0_dp

    !! Send non-constant grid variables to device
    z_d(:) = z(:)
    rhokap_d(:,:) = rhokap(:,:)
    alb_d(:,:) = alb(:,:)
    g_d(:,:) = g(:,:)
    rho_d(:) = rho(:)
    k_abs_d(:,:) = k_abs(:,:)
    Firr_d = Firr
    Fint_d = Fint
    OLR_d = 0.0_dp


    print*, n,'Longwave'

    call pp_semi_grey_int_kernel<<<blocks, threads>>>(Nph_int_d)

    !istat = cudaDeviceSynchronize()

    print*, n,'Shortwave'

    call pp_semi_grey_irr_kernel<<<blocks, threads>>>(Nph_irr_d)

    istat = cudaDeviceSynchronize()


    ! Moment and absorption calculation
    !! Get values from device
    Jdot(:) = Jdot_d(:) ; Hdot(:) = Hdot_d(:) ; Kdot(:) = Kdot_d(:)
    Jdot_s(:) = Jdot_s_d(:) ; Hdot_s(:) = Hdot_s_d(:) ; Kdot_s(:) = Kdot_s_d(:)
    !! Calculate Moments
    Jdot(:) = Jdot(:)/dze(:)/fourpi ; Hdot(:) = Hdot(:)/dze(:)/fourpi ; Kdot(:) = Kdot(:)/dze(:)/fourpi
    Jdot_s(:) = Jdot_s(:)/dze(:)/fourpi ; Hdot_s(:) = Hdot_s(:)/dze(:)/fourpi ; Kdot_s(:) = Kdot_s(:)/dze(:)/fourpi

    Adot(:) = Adot_d(:) ;  Adot_s(:) = Adot_s_d(:)
    Adot(:) = Adot(:)/dze(:) ; Adot_s(:) = Adot_s(:)/dze(:)

    !! Do temperature adjustment with safety factor
    dTl(:) = (((Adot(:) + Adot_s(:))/(4.0_dp * sb_c * rho(:)*k_P(:)))**(1.0_dp/4.0_dp) - Tl(:)) * 0.8_dp
    Tl(:) = Tl(:) + dTl(:)

    write(u0,*) n
    do i = 1, nlay
      write(u0,*) i, pl(i)/1e6_dp, (tau_R(i)+tau_R(i+1))/2.0_dp, Tl(i)
    end do

    OLR = OLR_d
    print*, 'OLR: ', (OLR/sb_c)**(1.0_dp/4.0_dp), 'Tint: ', Tint, 'Tirr: ', Tirr, 'Tot: ', (mu_z*Tirr**4 + Tint**4)**(1.0_dp/4.0_dp)

  end do


  open(newunit=u3,file='results/estim_pp_semi_grey_1D_pp.txt',action='readwrite')
  do i = 1, nlay
    write(u3,*) i, (tau_R(i)+tau_R(i+1))/2.0_dp, &
    &  Jdot(i), Hdot(i), Kdot(i), Adot(i), Jdot_s(i), Hdot_s(i), Kdot_s(i), Adot_s(i)
  end do

  close(u0)
  close(u3)

end subroutine semi_grey

