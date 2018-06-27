
module inputoutput
  implicit none

  character(16) :: sysname
  integer :: nx1_m
  integer :: ny1_m
  integer :: nz1_m
  integer :: nx2_m
  integer :: ny2_m
  integer :: nz2_m
  integer :: nt
  real(8) :: hx_m
  real(8) :: hy_m
  real(8) :: hz_m
  real(8) :: dt
  real(8) :: ac_0
  real(8) :: epdir(1:3)
  real(8) :: t_pulse
  real(8) :: omega
  real(8) :: omega_l
  real(8) :: gamma_l
  real(8) :: chi_l0
  logical :: out_ac_bin
  logical :: out_ac_out
  integer :: ac_out_step
  real(8) :: angle
  logical :: inp_bc

contains



  subroutine read_input()
    implicit none

    namelist/input/ &
    & sysname, &
    & nx1_m, &
    & ny1_m, &
    & nz1_m, &
    & nx2_m, &
    & ny2_m, &
    & nz2_m, &
    & nt, &
    & hx_m, &
    & hy_m, &
    & hz_m, &
    & dt, &
    & ac_0, &
    & epdir, &
    & t_pulse, &
    & omega, &
    & omega_l, &
    & gamma_l, &
    & chi_l0, &
    & out_ac_bin, &
    & out_ac_out, &
    & ac_out_step, &
    & angle, &
    & inp_bc

    sysname="untitled"
    nx1_m=-500
    ny1_m=1
    nz1_m=1
    nx2_m=500
    ny2_m=1
    nz2_m=1
    nt=1000
    hx_m=250d0
    hy_m=250d0
    hz_m=250d0
    dt=1
    ac_0=1d0
    epdir=(/0d0,0d0,1d0/)
    t_pulse=440d0
    omega=0.057d0
    omega_l=1d1
    gamma_l=0d0
    chi_l0=1d0
    out_ac_bin=.false.
    out_ac_out=.true.
    ac_out_step=1000
    angle=0
    inp_bc=.true.

    read (*, nml=input)

    return
  end subroutine read_input



  subroutine var_dump(iunit)
    implicit none
    integer, intent(in) :: iunit
    
    write(iunit, '("#",4x,"sysname=",1(1X, a))') sysname
    write(iunit, '("#",4x,"nx1_m=",1(1X, i6))') nx1_m
    write(iunit, '("#",4x,"ny1_m=",1(1X, i6))') ny1_m
    write(iunit, '("#",4x,"nz1_m=",1(1X, i6))') nz1_m
    write(iunit, '("#",4x,"nx2_m=",1(1X, i6))') nx2_m
    write(iunit, '("#",4x,"ny2_m=",1(1X, i6))') ny2_m
    write(iunit, '("#",4x,"nz2_m=",1(1X, i6))') nz2_m
    write(iunit, '("#",4x,"nt=",1(1X, i6))') nt
    write(iunit, '("#",4x,"hx_m=",1(1X, es23.15e3))') hx_m
    write(iunit, '("#",4x,"hy_m=",1(1X, es23.15e3))') hy_m
    write(iunit, '("#",4x,"hz_m=",1(1X, es23.15e3))') hz_m
    write(iunit, '("#",4x,"dt=",1(1X, es23.15e3))') dt
    write(iunit, '("#",4x,"ac_0=",1(1X, es23.15e3))') ac_0
    write(iunit, '("#",4x,"epdir=",3(1X, es23.15e3))') epdir
    write(iunit, '("#",4x,"t_pulse=",1(1X, es23.15e3))') t_pulse
    write(iunit, '("#",4x,"omega=",1(1X, es23.15e3))') omega
    write(iunit, '("#",4x,"omega_l=",1(1X, es23.15e3))') omega_l
    write(iunit, '("#",4x,"gamma_l=",1(1X, es23.15e3))') gamma_l
    write(iunit, '("#",4x,"chi_l0=",1(1X, es23.15e3))') chi_l0
    write(iunit, '("#",4x,"out_ac_bin=",1(1X, l1))') out_ac_bin
    write(iunit, '("#",4x,"out_ac_out=",1(1X, l1))') out_ac_out
    write(iunit, '("#",4x,"ac_out_step=",1(1X, i6))') ac_out_step
    write(iunit, '("#",4x,"angle=",1(1X, es23.15e3))') angle
    write(iunit, '("#",4x,"inp_bc=",1(1X, l1))') inp_bc
    return
  end subroutine var_dump



end module inputoutput
