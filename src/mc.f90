module mc_m
  use, intrinsic :: iso_fortran_env
  implicit none
  integer(int64), parameter :: nx = 100, ny = 100
  integer(int32), allocatable :: spins(:, :)
  real(real64), parameter :: kbt = 2.2d0, beta = 1 / kbt
  integer(int64), parameter :: dy(4) = [1_int64, 0_int64, -1_int64, 0_int64]
  integer(int64), parameter :: dx(4) = [0_int64, 1_int64, 0_int64, -1_int64]
contains
  impure subroutine init()
    allocate(spins(nx, ny))
    spins(:, :) = 1_int32
  end subroutine init
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
    integer(int64) :: d
    integer(int64) :: near_summ, near_x, near_y
    integer(int64) :: delta_energy
    real(real64) :: r
    near_summ = 0_int64
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
end module mc_m

program spin_simulation
  use, intrinsic :: iso_fortran_env
  use mc_m
  implicit none
  integer(int32), parameter :: mcs = 100000
  integer(int32) :: i
  call init()
  do i = 1, mcs
     call metropolis()
  end do
  write(output_unit, '(*(g0, 1x))') spins(:, :)
end program spin_simulation
