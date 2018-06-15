! 射法計算用テストモジュール
module fdtd1d1
  implicit none
  
  character(64) :: sysname
  
  ! system size
  ! n(x|y|z)(1|2)_m: lower and upper end of macroscopic region
  ! m(x|y|z)(1|2)_m: ... with taking account of overlapped region
  integer :: nx1_m, nx2_m, mx1_m, mx2_m
  integer :: ny1_m, ny2_m, my1_m, my2_m
  integer :: nz1_m, nz2_m, mz1_m, mz2_m
  
  real(8) :: hx_m, hy_m, hz_m
  real(8) :: dt
  
  real(8) :: e_ex, e_em
  
  integer :: iter, nt
  
  ! vector potential field
  real(8), allocatable :: ac_new_m(:, :, :, :)
  real(8), allocatable :: ac_cur_m(:, :, :, :)
  real(8), allocatable :: ac_old_m(:, :, :, :)
  
  ! matter current density
  real(8), allocatable :: jmat_new_m(:, :, :, :)
  real(8), allocatable :: jmat_cur_m(:, :, :, :)
  real(8), allocatable :: jmat_old_m(:, :, :, :)
  ! real(8), allocatable :: pmat_cur_m(:, :, :, :)
  
  real(8) :: ac_0
  real(8) :: epdir(3)
  real(8) :: t_pulse
  real(8) :: omega
  
  ! lorentz drude model
  real(8) :: omega_l, gamma_l, chi_l0
    
  real(8), parameter :: c_speed = 137.0
  real(8), parameter :: pi = 3.14159265
  
contains
  
  function sin2cos(t) result(f)
    implicit none
    real(8), intent(in) :: t
    real(8) :: f
    
    if ((0 < t) .and. (t < t_pulse)) then
      f = sin(pi * t / t_pulse) ** 2 * cos(omega * t)
    else
      f = 0d0
    end if
    return
  end function sin2cos
  
  
  
  subroutine init_fdtd()
    implicit none    
    integer :: ix_m, iy_m, iz_m
    real(8) :: x, t, f_cur, f_old
    
    namelist/input/ &
    & sysname, nx1_m, nx2_m, hx_m, nt, dt, ac_0, epdir, t_pulse, omega, omega_l, gamma_l, chi_l0
    
    sysname = "untitled"
    hx_m = 250; hy_m = 250; hz_m = 250;
    nt = 10000
    dt = 0.10
    omega = 1.55d0 / 13.6d0 / 2d0
    ac_0 = 1.00d0
    t_pulse = 400.0
    epdir(1:3) = (/0d0, 0d0, 1d0/)
    nx1_m = -1000
    nx2_m = +1000
    ny1_m = 1
    ny2_m = 1
    nz1_m = 1
    nz2_m = 1
    omega_l = 1d0
    gamma_l = 0.2
    chi_l0 = 1d0
    
    
    read (*, nml=input)
    
    mx1_m = nx1_m - 1; mx2_m = nx2_m + 1
    my1_m = ny1_m - 1; my2_m = ny2_m + 1
    mz1_m = nz1_m - 1; mz2_m = nz2_m + 1
    
    allocate(ac_old_m(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    allocate(ac_cur_m(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    allocate(ac_new_m(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    allocate(jmat_old_m(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    allocate(jmat_cur_m(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    allocate(jmat_new_m(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    ! allocate(pmat_cur_m(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
        
    ac_new_m = 0d0; ac_cur_m = 0d0; ac_old_m = 0d0
    jmat_new_m = 0d0; jmat_cur_m = 0d0; jmat_old_m = 0d0
    !pmat_cur_m = 0d0
    
    do ix_m = nx1_m, nx2_m
      x = ix_m * hx_m
      t = - x / c_speed
      f_cur = sin2cos(t)
      f_old = sin2cos(t - dt)
      
      ac_cur_m(1, ix_m, :, :) = epdir(1) * ac_0 * f_cur
      ac_cur_m(2, ix_m, :, :) = epdir(2) * ac_0 * f_cur
      ac_cur_m(3, ix_m, :, :) = epdir(3) * ac_0 * f_cur
      
      ac_old_m(1, ix_m, :, :) = epdir(1) * ac_0 * f_old
      ac_old_m(2, ix_m, :, :) = epdir(2) * ac_0 * f_old
      ac_old_m(3, ix_m, :, :) = epdir(3) * ac_0 * f_old
    end do
    iter = 0
    e_ex = 0d0; e_em = 0d0
    
    return
  end subroutine init_fdtd
  
  subroutine dt_evolve_ac()
    implicit none
    integer :: ix_m, iy_m, iz_m
    real(8) :: rlap(3)
    
    iy_m = ny1_m; iz_m = nz1_m
!$omp parallel do default(shared) private(ix_m, rlap)
    do ix_m = nx1_m, nx2_m
      rlap = ( &
      & + ac_cur_m(:, ix_m+1, iy_m, iz_m) &
      & + ac_cur_m(:, ix_m-1, iy_m, iz_m) &
      & - 2 * ac_cur_m(:, ix_m, iy_m, iz_m) &
      & ) / hx_m ** 2
      ac_new_m(:, ix_m, iy_m, iz_m) = ( &
      & + 2 * ac_cur_m(:, ix_m, iy_m, iz_m) & 
      & - ac_old_m(:, ix_m, iy_m, iz_m) &
      & + (4 * pi * dt**2) * (- jmat_cur_m(:, ix_m, iy_m, iz_m)) &
      & + (c_speed**2 * dt**2) * rlap(:) &
      & )
    end do
  !$omp end parallel do
  
  
    ! ac_old_m(:, :, :, :) = ac_cur_m(:, :, :, :)
    ! ac_cur_m(:, :, :, :) = ac_new_m(:, :, :, :)
    ! 
    
  end subroutine 
    
    
  subroutine current()
    implicit none
    integer :: ix_m, iy_m, iz_m
    real(8) :: e(3), de(3)
    
    real(8) :: f1, f2, f3, f4
    
    f1 = 1.0 / (1d0 + 0.5d0 * gamma_l * dt)
    f2 = (1d0 - 0.5d0 * gamma_l * dt) 
    f3 = 2d0 - (omega_l * dt) ** 2
    f4 = chi_l0 * omega_l ** 2
    
    iy_m = ny1_m; iz_m = nz1_m
    !$omp parallel do default(shared) private(ix_m)
    do ix_m = 1, nx2_m
      jmat_new_m(:, ix_m, iy_m, iz_m) = f1 * ( &
      & - f2 * jmat_old_m(:, ix_m, iy_m, iz_m) &
      & + f3 * jmat_cur_m(:, ix_m, iy_m, iz_m) &
      & + f4 * ac_old_m(:, ix_m, iy_m, iz_m) &
      & - f4 * ac_cur_m(:, ix_m, iy_m, iz_m) * 2 &
      & + f4 * ac_new_m(:, ix_m, iy_m, iz_m) &
      & ) 
    end do
    !$omp end parallel do
    
    return
  end subroutine current
  
  
  subroutine proceed_vars()
    implicit none
    ac_old_m = ac_cur_m; ac_cur_m = ac_new_m; ac_new_m = 0d0
    jmat_old_m = jmat_cur_m; jmat_cur_m = jmat_new_m; jmat_new_m = 0d0
    return
  end subroutine proceed_vars

  subroutine write_ac()
    implicit none
    integer :: ix_m, iy_m, iz_m
    character(64) :: file_ac_out
    
    write(file_ac_out, '(a, "_ac_", i6.6, ".data")') trim(sysname), iter
    write(*, '("# write_ac:", a)') trim(file_ac_out)
    write(*, '("# e_em:", e23.15e3)') e_em
    write(*, '("# e_ex:", e23.15e3)') e_ex
    write(*, '("# total:", e23.15e3)') e_ex + e_em
    
    open(unit=100, file=trim(file_ac_out))
    
    iy_m = ny1_m; iz_m = nz1_m
    do ix_m = nx1_m, nx2_m
      write(100, '(i6,3(1x,e23.15e3))') ix_m, ac_cur_m(:, ix_m, iy_m, iz_m)
    end do

    close(100)
    return
  end subroutine write_ac
  
  
  subroutine calc_elemag()
    implicit none
    integer :: ix_m, iy_m, iz_m
    real(8) :: elec(3), bmag(3)
    real(8) :: e_ex_dt
    

    e_em = 0d0
    e_ex_dt = 0d0
    
    iy_m = ny1_m; iz_m = nz1_m
!$omp parallel do default(shared) private(ix_m, elec, bmag) reduction(+:e_em,e_ex_dt)
    do ix_m = nx1_m, nx2_m
      elec = - ( &
      & + ac_new_m(:, ix_m, iy_m, iz_m) &
      & - ac_old_m(:, ix_m, iy_m, iz_m) &
      & ) * (0.5 / dt)
      bmag(1) = 0d0
      bmag(2) = - (0.5 * c_speed / hx_m) * ( &
      & ac_cur_m(3, ix_m+1, iy_m, iz_m) -  ac_cur_m(3, ix_m-1, iy_m, iz_m) &
      & )
      bmag(3) = + (0.5 * c_speed / hx_m) * ( &
      & ac_cur_m(2, ix_m+1, iy_m, iz_m) -  ac_cur_m(2, ix_m-1, iy_m, iz_m) &
      & )
      e_em = e_em + sum(elec**2 + bmag**2) * (hx_m * hy_m * hz_m / (8 * pi))
      e_ex_dt = e_ex_dt - sum(elec * jmat_cur_m(:, ix_m, iy_m, iz_m)) * (hx_m * hy_m * hz_m) * dt
    end do
!$omp end parallel do
    e_ex = e_ex + e_ex_dt
    
  end subroutine calc_elemag
  
  
  subroutine update(nstep)
    implicit none
    integer, intent(in) :: nstep
    integer :: i
    do i = 1, nstep
      call dt_evolve_ac()
      call calc_elemag()
      call current()
      call proceed_vars()
      iter = iter + 1
    end do
  end subroutine
  
  
  
end module fdtd1d1
    
