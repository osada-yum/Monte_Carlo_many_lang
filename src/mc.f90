module mc_m
  use, intrinsic :: iso_fortran_env
  implicit none
  integer(int64), parameter :: nx = 100, ny = 100, nall = nx * ny
  real(real64), parameter :: nall_inv = 1d0 / nall
  integer(int32), allocatable :: spins(:, :)
  real(real64), parameter :: kbt = 2.3d0, beta = 1 / kbt
  integer(int64), parameter :: dy(4) = [1_int64, 0_int64, -1_int64, 0_int64]
  integer(int64), parameter :: dx(4) = [0_int64, 1_int64, 0_int64, -1_int64]
contains
  impure subroutine init()
    allocate(spins(nx, ny))
  end subroutine init
  impure subroutine set_order_spin()
    spins(:, :) = 1_int32
  end subroutine set_order_spin
  impure subroutine metropolis()
    integer(int64) :: i, j
    do i = 1, ny
       do j = 1 + iand(i, b'1'), nx, 2
          call local_flip(j, i)
       end do
    end do
    do i = 1, ny
       do j = 2 - iand(i, b'1'), nx, 2
          call local_flip(j, i)
       end do
    end do
  end subroutine metropolis
  impure subroutine local_flip(x, y)
    integer(int64), intent(in) :: x, y
    integer(int32) :: d
    integer(int32) :: near_summ, delta_energy
    integer(int64) :: near_x, near_y
    real(real64) :: r
    near_summ = 0_int32
    do d = 1, size(dy)
       near_x = x + dx(d)
       if (near_x > nx) then
          near_x = 1
       else if (near_x < 1) then
          near_x = nx
       end if
       near_y = y + dy(d)
       if (near_y > ny) then
          near_y = 1
       else if (near_y < 1) then
          near_y = ny
       end if
       near_summ = near_summ + spins(near_x, near_y)
    end do
    delta_energy = 2 * spins(x, y) * near_summ
    if (delta_energy <= 0) then
       spins(x, y) = - spins(x, y)
    else
       call random_number(r)
       if (r < exp(-beta * delta_energy)) &
            & spins(x, y) = - spins(x, y)
    end if
  end subroutine local_flip
  pure real(real64) function calc_magne() result(res)
    res = sum(int(spins, int64)) * nall_inv
  end function calc_magne
end module mc_m

program spin_simulation
  use, intrinsic :: iso_fortran_env
  use mc_m
  implicit none
  integer(int32), parameter :: mcs = 1000, nsample = 100
  real(real64), allocatable :: magnes(:)
  integer(int32) :: i, j
  call init()
  allocate(magnes(mcs), source = 0d0)
  do j = 1, nsample
     call set_order_spin()
     write(error_unit, '(*(g0, 1x))') "sample", j
     do i = 1, mcs
        call metropolis()
        magnes(i) = magnes(i) + calc_magne()
     end do
  end do
  do i = 1, mcs
     magnes(i) = magnes(i) / nsample
     write(output_unit, '(*(g0 ,1x))') i, magnes(i)
  end do
  ! write(output_unit, '(*(g0, 1x))') spins(:, :)
end program spin_simulation
