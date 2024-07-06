!> Module for deposit particle info on the grid
!>
module deposit
  implicit none
  private

  public :: deposit_charge_density

  contains
  subroutine deposit_charge_density(part, q, qm, nop, idimp, nxv)
    !> deposit charge density on the grid using linear interpolation
    implicit none
    integer, intent(in) :: nop, idimp, nxv
    !! nop: number of particles
    !! idimp: dimension of phase space = 2
    !! nxv: first dimension of charge array, must be >= nx+1

    real, intent(inout) :: part(idimp, nop), q(nxv)
    !! part: particle array, (1, n) is the position x of particle n
    !!       (2, n) is the velocity vx of particle n
    !! q: q(j) is the charge density at grid point j

    real, intent(in) :: qm
    !! qm: charge of particle

    real :: qdx
    !! dx: distance from particle to grid point 
    integer :: j, nn
    do j = 1, nop
      !> find the interpolation weights
      nn = part(1, j)
      qdx = qm * (part(1, j) - real(nn))
      nn = nn + 1
      !> deposit charge
      q(nn) = q(nn) + (qm - qdx)
      q(nn+1) = q(nn+1) + qdx
    end do
  end subroutine deposit_charge_density


  subroutine add_guard_cell(q, nx, nxe)
    implicit none
    integer, intent(in) :: nx, nxe
    real, intent(inout) :: q(nxe)

    q(1) = q(1) + q(nx + 1)
    q(nx+1) = 0
    
  end subroutine add_guard_cell

end module deposit