!! This file contains the cute banner we will print when a simulation is run
module rki_banner
   use pic_logger, only: logger => global_logger
   use pic_knowledge, only: get_knowledge
   implicit none
   private
   public :: print_banner
contains

   subroutine print_banner(sim_type)
    !! Print a cute banner to the console when the simulation starts.
    !! `sim_type` is accepted (and ignored) so callers can pass the
    !! configured regime tag without branching.
      character(len=*), intent(in), optional :: sim_type
      character(len=54) :: subtitle   ! 54 = interior box width (the ═ border)

      subtitle = "       2D/3D Coastal Hydrodynamics Model"
      if (present(sim_type)) continue   ! single-regime build: no variants

      call logger%info("")
      call logger%info("  ╔══════════════════════════════════════════════════════╗")
      call logger%info("  ║                                                      ║")
      call logger%info("  ║    ██████╗  █████╗ ██╗  ██╗ █████╗ ██╗     ██╗       ║")
      call logger%info("  ║    ██╔══██╗██╔══██╗██║ ██╔╝██╔══██╗██║     ██║       ║")
      call logger%info("  ║    ██████╔╝███████║█████╔╝ ███████║██║     ██║       ║")
      call logger%info("  ║    ██╔══██╗██╔══██║██╔═██╗ ██╔══██║██║     ██║       ║")
      call logger%info("  ║    ██║  ██║██║  ██║██║  ██╗██║  ██║███████╗██║       ║")
      call logger%info("  ║    ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝╚═╝       ║")
      call logger%info("  ║                                                      ║")
      call logger%info("  ║        A GPU accelerated, Fortran native             ║")
      call logger%info("  ║"//subtitle//"║")
      call logger%info("  ║                                                      ║")
      call logger%info("  ║                                                      ║")
      call logger%info("  ╚══════════════════════════════════════════════════════╝")
      call logger%info("")
      call logger%info("  Contributors:")
      call logger%info("    Jorge Luis Galvez Vallejo")
      call logger%info("    Albert Sietze Thie")
      ! add your name here after you've contributed!
      call logger%info("")
      call get_knowledge()
   end subroutine print_banner

end module rki_banner
