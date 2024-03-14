program gCMCRT_RE
  implicit none

  integer :: u_nml
  character(len=100) :: exper

  namelist /main_nml/ exper

  !! Read input variables from namelist
  open(newunit=u_nml, file='gCMCRT_RE.nml', status='old', action='read')
  read(u_nml, nml=main_nml)
  close(u_nml)

  select case(trim(exper))

  case('grey')

    call grey()

  case('semi_grey')

    call semi_grey()

  case default
    print*, 'Invalid experiment string: ', trim(exper)
  end select



end program gCMCRT_RE
