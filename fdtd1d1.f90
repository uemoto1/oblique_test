! 射法計算用テストモジュール
module fdtd1d1
  use inputoutput
  implicit none
  
  
  ! system size
  ! n(x|y|z)(1|2)_m: lower and upper end of macroscopic region
  ! m(x|y|z)(1|2)_m: ... with taking account of overlapped region
  integer :: mx1_m, mx2_m
  integer :: my1_m, my2_m
  integer :: mz1_m, mz2_m
  
  real(8) :: e_ex, e_em
  
  real(8) :: rinv_dx, rinv_dy, rinv_dz
  integer :: iter
  
  ! vector potential field
  real(8), allocatable :: ac_new_ms(:, :, :, :)
  real(8), allocatable :: ac_ms(:, :, :, :)
  real(8), allocatable :: ac_old_ms(:, :, :, :)
  
  ! matter current density
  real(8), allocatable :: jm_new_ms(:, :, :, :)
  real(8), allocatable :: jm_ms(:, :, :, :)
  real(8), allocatable :: jm_old_ms(:, :, :, :)
  ! real(8), allocatable :: pmat_cur_m(:, :, :, :)
  
  ! temporal region
  real(8), allocatable :: data(:, :, :, :, :)
    
  real(8), parameter :: c_light = 137.0
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
    
    call read_input()
    call var_dump()
    
    mx1_m = nx1_m - 1; mx2_m = nx2_m + 1
    my1_m = ny1_m - 1; my2_m = ny2_m + 1
    mz1_m = nz1_m - 1; mz2_m = nz2_m + 1
  
    rinv_dx = 1d0 / hx_m
    rinv_dy = 1d0 / hy_m
    rinv_dz = 1d0 / hz_m
    
    allocate(ac_old_ms(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    allocate(ac_ms(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    allocate(ac_new_ms(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    allocate(jm_old_ms(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    allocate(jm_ms(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    allocate(jm_new_ms(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    ! allocate(pmat_cur_m(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    
    allocate(data(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m, 0:nt))
        
    ac_new_ms = 0d0; ac_ms = 0d0; ac_old_ms = 0d0
    jm_new_ms = 0d0; jm_ms = 0d0; jm_old_ms = 0d0
    !pmat_cur_m = 0d0
    
    do ix_m = nx1_m, nx2_m
      x = ix_m * hx_m
      t = - x / c_light
      f_cur = sin2cos(t)
      f_old = sin2cos(t - dt)
      
      ac_new_ms(1, ix_m, :, :) = epdir(1) * ac_0 * f_cur
      ac_new_ms(2, ix_m, :, :) = epdir(2) * ac_0 * f_cur
      ac_new_ms(3, ix_m, :, :) = epdir(3) * ac_0 * f_cur
      
      ac_ms(1, ix_m, :, :) = epdir(1) * ac_0 * f_old
      ac_ms(2, ix_m, :, :) = epdir(2) * ac_0 * f_old
      ac_ms(3, ix_m, :, :) = epdir(3) * ac_0 * f_old
    end do
    iter = 0
    e_ex = 0d0; e_em = 0d0
    return
  end subroutine init_fdtd
  
  
  
  subroutine dt_evolve_ac()
    implicit none
    integer :: ix_m, iy_m, iz_m
    real(8) :: rr(3)
    
    do iz_m = nz1_m, nz2_m
      do iy_m = ny1_m, ny2_m
        do ix_m = nx1_m, nx2_m
          ! Calculate Rot Rot A
          rr(1) = - (rinv_dy**2) * Ac_ms(1, ix_m+0, iy_m-1, iz_m+0) &
                & - (rinv_dz**2) * Ac_ms(1, ix_m+0, iy_m+0, iz_m-1) &
                & + (2d0*(rinv_dy**2 + rinv_dz**2)) * Ac_ms(1, ix_m+0, iy_m+0, iz_m+0) &
                & - (rinv_dz**2) * Ac_ms(1, ix_m+0, iy_m+0, iz_m+1) &
                & - (rinv_dy**2) * Ac_ms(1, ix_m+0, iy_m+1, iz_m+0) &
                & + (rinv_dx*rinv_dy*0.25d0) * Ac_ms(2, ix_m-1, iy_m-1, iz_m+0) &
                & - (rinv_dx*rinv_dy*0.25d0) * Ac_ms(2, ix_m-1, iy_m+1, iz_m+0) &
                & - (rinv_dx*rinv_dy*0.25d0) * Ac_ms(2, ix_m+1, iy_m-1, iz_m+0) &
                & + (rinv_dx*rinv_dy*0.25d0) * Ac_ms(2, ix_m+1, iy_m+1, iz_m+0) &
                & + (rinv_dx*rinv_dz*0.25d0) * Ac_ms(3, ix_m-1, iy_m+0, iz_m-1) &
                & - (rinv_dx*rinv_dz*0.25d0) * Ac_ms(3, ix_m-1, iy_m+0, iz_m+1) &
                & - (rinv_dx*rinv_dz*0.25d0) * Ac_ms(3, ix_m+1, iy_m+0, iz_m-1) &
                & + (rinv_dx*rinv_dz*0.25d0) * Ac_ms(3, ix_m+1, iy_m+0, iz_m+1)
          rr(2) = + (rinv_dx*rinv_dy*0.25d0) * Ac_ms(1, ix_m-1, iy_m-1, iz_m+0) &
                & - (rinv_dx*rinv_dy*0.25d0) * Ac_ms(1, ix_m-1, iy_m+1, iz_m+0) &
                & - (rinv_dx*rinv_dy*0.25d0) * Ac_ms(1, ix_m+1, iy_m-1, iz_m+0) &
                & + (rinv_dx*rinv_dy*0.25d0) * Ac_ms(1, ix_m+1, iy_m+1, iz_m+0) &
                & - (rinv_dx**2) * Ac_ms(2, ix_m-1, iy_m+0, iz_m+0) &
                & - (rinv_dz**2) * Ac_ms(2, ix_m+0, iy_m+0, iz_m-1) &
                & + (2d0*(rinv_dx**2 + rinv_dz**2)) * Ac_ms(2, ix_m+0, iy_m+0, iz_m+0) &
                & - (rinv_dz**2) * Ac_ms(2, ix_m+0, iy_m+0, iz_m+1) &
                & - (rinv_dx**2) * Ac_ms(2, ix_m+1, iy_m+0, iz_m+0) &
                & + (rinv_dy*rinv_dz*0.25d0) * Ac_ms(3, ix_m+0, iy_m-1, iz_m-1) &
                & - (rinv_dy*rinv_dz*0.25d0) * Ac_ms(3, ix_m+0, iy_m-1, iz_m+1) &
                & - (rinv_dy*rinv_dz*0.25d0) * Ac_ms(3, ix_m+0, iy_m+1, iz_m-1) &
                & + (rinv_dy*rinv_dz*0.25d0) * Ac_ms(3, ix_m+0, iy_m+1, iz_m+1)
          rr(3) = + (rinv_dx*rinv_dz*0.25d0) * Ac_ms(1, ix_m-1, iy_m+0, iz_m-1) &
                & - (rinv_dx*rinv_dz*0.25d0) * Ac_ms(1, ix_m-1, iy_m+0, iz_m+1) &
                & - (rinv_dx*rinv_dz*0.25d0) * Ac_ms(1, ix_m+1, iy_m+0, iz_m-1) &
                & + (rinv_dx*rinv_dz*0.25d0) * Ac_ms(1, ix_m+1, iy_m+0, iz_m+1) &
                & + (rinv_dy*rinv_dz*0.25d0) * Ac_ms(2, ix_m+0, iy_m-1, iz_m-1) &
                & - (rinv_dy*rinv_dz*0.25d0) * Ac_ms(2, ix_m+0, iy_m-1, iz_m+1) &
                & - (rinv_dy*rinv_dz*0.25d0) * Ac_ms(2, ix_m+0, iy_m+1, iz_m-1) &
                & + (rinv_dy*rinv_dz*0.25d0) * Ac_ms(2, ix_m+0, iy_m+1, iz_m+1) &
                & - (rinv_dx**2) * Ac_ms(3, ix_m-1, iy_m+0, iz_m+0) &
                & - (rinv_dy**2) * Ac_ms(3, ix_m+0, iy_m-1, iz_m+0) &
                & + (2d0*(rinv_dx**2 + rinv_dy**2)) * Ac_ms(3, ix_m+0, iy_m+0, iz_m+0) &
                & - (rinv_dy**2) * Ac_ms(3, ix_m+0, iy_m+1, iz_m+0) &
                & - (rinv_dx**2) * Ac_ms(3, ix_m+1, iy_m+0, iz_m+0)
          Ac_new_ms(:,ix_m, iy_m, iz_m) = &
                & + (2 * Ac_ms(:,ix_m, iy_m, iz_m) - Ac_old_ms(:,ix_m, iy_m, iz_m) &
                & - Jm_ms(:,ix_m, iy_m, iz_m) * 4.0*pi*(dt**2) - rr(:)*(c_light*dt)**2 )
        end do
      end do
    end do

    do ix_m = mx1_m, mx2_m
      do iy_m = my1_m, my2_m
        Ac_new_ms(1:3, ix_m, iy_m, mz1_m) = Ac_new_ms(1:3, ix_m, iy_m, nz2_m)
        Ac_new_ms(1:3, ix_m, iy_m, mz2_m) = Ac_new_ms(1:3, ix_m, iy_m, nz1_m)
      end do
    end do
    
    do ix_m = mx1_m, mx2_m
      do iz_m = mz1_m, mz2_m
        Ac_new_ms(1:3, ix_m, my1_m, iz_m) = Ac_new_ms(1:3, ix_m, ny2_m, iz_m) 
        Ac_new_ms(1:3, ix_m, my2_m, iz_m) = Ac_new_ms(1:3, ix_m, ny1_m, iz_m) 
      end do
    end do
    
    do iy_m = my1_m, my2_m
      do iz_m = mz1_m, mz2_m
        Ac_new_ms(1:3, mx1_m, iy_m, iz_m) = Ac_new_ms(1:3, nx2_m, iy_m, iz_m)
        Ac_new_ms(1:3, mx2_m, iy_m, iz_m) = Ac_new_ms(1:3, nx1_m, iy_m, iz_m)
      end do
    end do
    return    
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
    
    do iz_m = nz1_m, nz2_m
      do iy_m = ny1_m, ny2_m
        do ix_m = 1, nx2_m
          jm_new_ms(:, ix_m, iy_m, iz_m) = f1 * ( &
          & - f2 * jm_old_ms(:, ix_m, iy_m, iz_m) &
          & + f3 * jm_ms(:, ix_m, iy_m, iz_m) &
          & + f4 * ac_old_ms(:, ix_m, iy_m, iz_m) &
          & - f4 * ac_ms(:, ix_m, iy_m, iz_m) * 2 &
          & + f4 * ac_new_ms(:, ix_m, iy_m, iz_m) &
          & )
        end do
      end do
    end do
    return
  end subroutine current
  
  
  
  subroutine proceed_vars()
    implicit none
    ac_old_ms = ac_ms; ac_ms = ac_new_ms; ac_new_ms = 0d0
    jm_old_ms = jm_ms; jm_ms = jm_new_ms; jm_new_ms = 0d0
    return
  end subroutine proceed_vars



  subroutine write_ac()
    implicit none
    integer :: ix_m, iy_m, iz_m
    character(64) :: file_ac_out
    
    write(file_ac_out, '(a, "_ac_", i6.6, ".data")') trim(sysname), iter
    print '(4x,"# export: ", a)', trim(file_ac_out)
    open(unit=100, file=trim(file_ac_out))
    do iz_m = nz1_m, nz2_m
      do iy_m = ny1_m, ny2_m
        do ix_m = nx1_m, nx2_m
          write(100, '(3(1x,i6),3(1x,e23.15e3))') ix_m, iy_m, iz_m, ac_ms(:, ix_m, iy_m, iz_m)
        end do
      end do
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
    
    do iz_m = nz1_m, nz2_m
      do iy_m = ny1_m, ny2_m
        do ix_m = nx1_m, nx2_m
          elec = - ( &
          & + ac_new_ms(:, ix_m, iy_m, iz_m) &
          & - ac_old_ms(:, ix_m, iy_m, iz_m) &
          & ) * (0.5 / dt)
          bmag(1) = 0d0
          bmag(2) = - (0.5 * c_light / hx_m) * ( &
          & ac_ms(3, ix_m+1, iy_m, iz_m) -  ac_ms(3, ix_m-1, iy_m, iz_m) &
          & )
          bmag(3) = + (0.5 * c_light / hx_m) * ( &
          & ac_ms(2, ix_m+1, iy_m, iz_m) -  ac_ms(2, ix_m-1, iy_m, iz_m) &
          & )
          e_em = e_em + sum(elec**2 + bmag**2) * (hx_m * hy_m * hz_m / (8 * pi))
          e_ex_dt = e_ex_dt - sum(elec * jm_ms(:, ix_m, iy_m, iz_m)) * (hx_m * hy_m * hz_m) * dt
        end do
      end do
    end do    
    e_ex = e_ex + e_ex_dt
    return
  end subroutine calc_elemag
  
  
  
  subroutine run_fdtd()
    implicit none
    do iter = 0, nt
      call proceed_vars()
      call dt_evolve_ac()
      call calc_elemag()
      call current()
      
      if (mod(iter, 1000) == 0) then
        print '("# iter=", i6)', iter
        print '(4x, "# E_ex=", es23.15e3)', e_ex
        print '(4x, "# E_em=", es23.15e3)', e_em
        print '(4x, "# E_tot=", es23.15e3)', e_ex + e_em
        call write_ac()
      end if
    end do
  end subroutine run_fdtd
  
  
  
end module fdtd1d1
    
