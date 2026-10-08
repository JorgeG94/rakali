!! Structured grid metadata
module rki_grid
   !! Pure metadata describing a structured Cartesian grid: physical and
   !! ghosted dimensions, ghost width, and cell sizes. No allocatables,
   !! no GPU mapping needed.
   use rki_constants, only: wp
   implicit none
   private

   public :: hgrid_t

   type :: hgrid_t
      integer  :: nx_total = 0
         !! Total cells in x including ghosts (nx_phys + 2*nghost)
      integer  :: ny_total = 0
         !! Total cells in y including ghosts
      integer  :: nx_phys = 0
         !! Physical (interior) cells in x
      integer  :: ny_phys = 0
         !! Physical (interior) cells in y
      integer  :: nghost = 2
         !! Ghost cell width on each side
      real(wp) :: dx = 0.0_wp
         !! Cell size in x (m)
      real(wp) :: dy = 0.0_wp
         !! Cell size in y (m)
   contains
      procedure, non_overridable :: init => grid_init
   end type hgrid_t

contains

   subroutine grid_init(this, nx_phys, ny_phys, nghost, dx, dy)
      class(hgrid_t), intent(inout) :: this
      integer, intent(in) :: nx_phys, ny_phys, nghost
      real(wp), intent(in) :: dx, dy

      this%nx_phys = nx_phys
      this%ny_phys = ny_phys
      this%nghost = nghost
      this%dx = dx
      this%dy = dy
      this%nx_total = nx_phys + 2*nghost
      this%ny_total = ny_phys + 2*nghost
   end subroutine grid_init

end module rki_grid
