!! Central Rakali namelist schema (strict validation + parameter dumps).
module rki_config_schema
   !! Re-export shim for Rakali's strict namelist schema. Re-exports
   !! `build_rakali_schema` (defined in `rki_config`, kept here to avoid a
   !! circular dependency) and owns the `nml_dirname` helper.
   use rki_config, only: build_rakali_schema
   implicit none
   private

   public :: build_rakali_schema
   public :: nml_dirname

contains

   pure function nml_dirname(path) result(dir)
      !! Directory part of `path` with a trailing '/', or '' if none.
      character(len=*), intent(in) :: path
      character(len=:), allocatable :: dir
      integer :: i
      dir = ""
      do i = len_trim(path), 1, -1
         if (path(i:i) == "/") then
            dir = path(1:i)
            return
         end if
      end do
   end function nml_dirname

end module rki_config_schema
