!! Rakali namelist knob-reference generator.
program nml_doc
   !! Host-only CLI tool: build a default `config_t`, register the full
   !! strict schema against it, and render every group + knob to Markdown
   !! via `render_markdown`.  Output path is arg 1 (default `nml_knobs.md`).
   !! The defaults captured are the pristine `config_t` field initialisers.
   use rki_config, only: config_t, build_rakali_schema
   use rki_nml_schema, only: nml_schema_t
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   implicit none

   type(config_t), target :: cfg
   type(nml_schema_t) :: schema
   character(len=512) :: out_path
   integer :: arg_len, io_stat, ngroups, nknobs, ig

   if (command_argument_count() >= 1) then
      call get_command_argument(1, out_path, arg_len, io_stat)
   else
      out_path = "nml_knobs.md"
   end if

   call build_rakali_schema(cfg, schema)
   call schema%render_markdown(trim(out_path))

   ngroups = size(schema%groups)
   nknobs = 0
   do ig = 1, ngroups
      if (allocated(schema%groups(ig)%keys)) nknobs = nknobs + size(schema%groups(ig)%keys)
   end do

   call logger%info("Wrote namelist knob reference to "//trim(out_path)// &
                    " ("//to_char(ngroups)//" groups, "//to_char(nknobs)//" knobs)")

end program nml_doc
