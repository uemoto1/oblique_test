! 射法計算用テストモジュール
module fdtd1d1
  implicit none
  
  character(64) :: sysname
  
  ! System Size
  ! n(x|y|z)(1|2)_m: lower and upper end of macroscopic region
  ! m(x|y|z)(1|2)_m: ... with taking account of overlapped region
  integer :: nx1_m, nx2_m, mx1_m, mx2_m
  integer :: ny1_m, ny2_m, my1_m, my2_m
  integer :: nz1_m, nz2_m, mz1_m, mz2_m
  
  real(8) :: hx_m, hy_m, hz_m
  real(8) :: dt
  
  real(8) :: E_ex, E_EM
  
  integer :: iter, nt
  
  ! Vector Potential Field
  real(8), allocatable :: Ac_new_m(:, :, :, :)
  real(8), allocatable :: Ac_cur_m(:, :, :, :)
  real(8), allocatable :: Ac_old_m(:, :, :, :)
  
  ! Matter Current Density
  real(8), allocatable :: Jmat_new_m(:, :, :, :)
  real(8), allocatable :: Jmat_cur_m(:, :, :, :)
  real(8), allocatable :: Jmat_old_m(:, :, :, :)
  ! real(8), allocatable :: Pmat_cur_m(:, :, :, :)
  
  real(8) :: Ac_0
  real(8) :: epdir(3)
  real(8) :: t_pulse
  real(8) :: omega
  
  ! Lorentz Drude Model
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
    & sysname, nx1_m, nx2_m, hx_m, nt, dt, Ac_0, epdir, t_pulse, omega, omega_l, gamma_l, chi_l0
    
    sysname = "untitled"
    hx_m = 250; hy_m = 250; hz_m = 250;
    nt = 10000
    dt = 0.10
    omega = 1.55d0 / 13.6d0 / 2d0
    Ac_0 = 1.00d0
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
    
    allocate(Ac_old_m(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    allocate(Ac_cur_m(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    allocate(Ac_new_m(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    allocate(Jmat_old_m(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    allocate(Jmat_cur_m(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    allocate(Jmat_new_m(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    ! allocate(Pmat_cur_m(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
        
    Ac_new_m = 0d0; Ac_cur_m = 0d0; Ac_old_m = 0d0
    Jmat_new_m = 0d0; Jmat_cur_m = 0d0; Jmat_old_m = 0d0
    !Pmat_cur_m = 0d0
    
    do ix_m = nx1_m, nx2_m
      x = ix_m * hx_m
      t = - x / c_speed
      f_cur = sin2cos(t)
      f_old = sin2cos(t - dt)
      
      Ac_cur_m(1, ix_m, :, :) = epdir(1) * Ac_0 * f_cur
      Ac_cur_m(2, ix_m, :, :) = epdir(2) * Ac_0 * f_cur
      Ac_cur_m(3, ix_m, :, :) = epdir(3) * Ac_0 * f_cur
      
      Ac_old_m(1, ix_m, :, :) = epdir(1) * Ac_0 * f_old
      Ac_old_m(2, ix_m, :, :) = epdir(2) * Ac_0 * f_old
      Ac_old_m(3, ix_m, :, :) = epdir(3) * Ac_0 * f_old
    end do
    iter = 0
    E_EX = 0d0; E_EM = 0d0
    
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
      & + Ac_cur_m(:, ix_m+1, iy_m, iz_m) &
      & + Ac_cur_m(:, ix_m-1, iy_m, iz_m) &
      & - 2 * Ac_cur_m(:, ix_m, iy_m, iz_m) &
      & ) / hx_m ** 2
      ac_new_m(:, ix_m, iy_m, iz_m) = ( &
      & + 2 * ac_cur_m(:, ix_m, iy_m, iz_m) & 
      & - ac_old_m(:, ix_m, iy_m, iz_m) &
      & + (4 * pi * dt**2) * (- Jmat_cur_m(:, ix_m, iy_m, iz_m)) &
      & + (c_speed**2 * dt**2) * rlap(:) &
      & )
    end do
  !$omp end parallel do
  
  
    ! Ac_old_m(:, :, :, :) = Ac_cur_m(:, :, :, :)
    ! Ac_cur_m(:, :, :, :) = Ac_new_m(:, :, :, :)
    ! 
    
  end subroutine 
    
    
  subroutine current()
    implicit none
    integer :: ix_m, iy_m, iz_m
    real(8) :: E(3), dE(3)
    
    real(8) :: f1, f2, f3, f4
    
    f1 = 1.0 / (1d0 + 0.5d0 * gamma_l * dt)
    f2 = (1d0 - 0.5d0 * gamma_l * dt) 
    f3 = 2d0 - (omega_l * dt) ** 2
    f4 = chi_l0 * omega_l ** 2
    
    iy_m = ny1_m; iz_m = nz1_m
    !$omp parallel do default(shared) private(ix_m)
    do ix_m = 1, nx2_m
      Jmat_new_m(:, ix_m, iy_m, iz_m) = f1 * ( &
      & - f2 * Jmat_old_m(:, ix_m, iy_m, iz_m) &
      & + f3 * Jmat_cur_m(:, ix_m, iy_m, iz_m) &
      & + f4 * Ac_old_m(:, ix_m, iy_m, iz_m) &
      & - f4 * Ac_cur_m(:, ix_m, iy_m, iz_m) * 2 &
      & + f4 * Ac_new_m(:, ix_m, iy_m, iz_m) &
      & ) 
    end do
    !$omp end parallel do
    
    return
  end subroutine current
  
  
  subroutine proceed_vars()
    implicit none
    Ac_old_m = Ac_cur_m; Ac_cur_m = Ac_new_m; Ac_new_m = 0d0
    Jmat_old_m = Jmat_cur_m; Jmat_cur_m = Jmat_new_m; Jmat_new_m = 0d0
    return
  end subroutine proceed_vars

  subroutine write_ac()
    implicit none
    integer :: ix_m, iy_m, iz_m
    character(64) :: file_ac_out
    
    write(file_ac_out, '(a, "_ac_", i6.6, ".data")') trim(sysname), iter
    write(*, '("# write_ac:", a)') trim(file_ac_out)
    write(*, '("# E_EM:", e23.15e3)') E_EM
    write(*, '("# E_EX:", e23.15e3)') E_EX
    write(*, '("# Total:", e23.15e3)') E_EX + E_EM
    
    open(unit=100, file=trim(file_ac_out))
    
    iy_m = ny1_m; iz_m = nz1_m
    do ix_m = nx1_m, nx2_m
      write(100, '(i6,3(1x,e23.15e3))') ix_m, Ac_cur_m(:, ix_m, iy_m, iz_m)
    end do

    close(100)
    return
  end subroutine write_ac
  
  
  subroutine calc_elemag()
    implicit none
    integer :: ix_m, iy_m, iz_m
    real(8) :: elec(3), bmag(3)
    
    
    E_EM = 0d0
    
    iy_m = ny1_m; iz_m = nz1_m
    do ix_m = nx1_m, nx2_m
      elec = - ( &
      & + ac_new_m(:, ix_m, iy_m, iz_m) &
      & - ac_old_m(:, ix_m, iy_m, iz_m) &
      & ) * (0.5 / dt)
      bmag(1) = 0d0
      bmag(2) = - (0.5 * c_speed / HX_m) * ( &
      & Ac_cur_m(3, ix_m+1, iy_m, iz_m) -  Ac_cur_m(3, ix_m-1, iy_m, iz_m) &
      & )
      bmag(3) = + (0.5 * c_speed / HX_m) * ( &
      & Ac_cur_m(2, ix_m+1, iy_m, iz_m) -  Ac_cur_m(2, ix_m-1, iy_m, iz_m) &
      & )
      E_em = E_em + sum(elec**2 + bmag**2) * (HX_m * HY_m * HZ_m / (8 * pi))
      E_ex = E_ex - sum(elec * Jmat_cur_m(:, ix_m, iy_m, iz_m)) * (HX_m * HY_m * HZ_m) * dt
    end do
    
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
    
