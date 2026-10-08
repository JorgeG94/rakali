!! Rakali — GPU-native shallow-water / multilayer coastal solver.
program rakali
   !! Thin entry-point wrapper.  Owns:
   !!   * command-line argument parsing + config load,
   !!   * output-directory bootstrap.
   !!
   !! All compute-side logic lives in `rki_driver%driver_run` so the
   !! same lifecycle is exercised from both `app/main.F90` and from
   !! unit tests (`tests/unit/test_driver_smoke.F90`).  This catches
   !! integration-level bugs (OpenACC partial-presence, init ordering,
   !! missing attaches) that per-kernel unit tests miss.
   use rki_config, only: config_t, read_config, validate_config
   use rki_nml_schema, only: nml_schema_t
   use rki_driver, only: driver_run, configure_log_level
   use pic_logger, only: logger => global_logger
   implicit none

   type(config_t), target :: cfg
   type(nml_schema_t) :: schema
   character(len=:), allocatable :: doc_dir
   character(len=256) :: input_file
   integer :: arg_len, io_stat
   logical :: input_exists

   ! Read input file from the command line.
   ! Running without a namelist used to fall back to a hard-coded default
   ! state, which then SIGBUS'd downstream once any output path was hit.
   ! Fail loudly here instead so the user gets a clear "no nml" message
   ! rather than a stale-default crash several seconds in.
   if (command_argument_count() < 1) then
      call logger%error("FATAL: no input file provided.")
      call logger%error("Usage: rakali <input_file.nml>")
      error stop 2
   end if

   call get_command_argument(1, input_file, arg_len, io_stat)
   inquire (file=trim(input_file), exist=input_exists)
   if (.not. input_exists) then
      call logger%error("FATAL: input file not found: "//trim(input_file))
      error stop 2
   end if
   ! read_config builds + strict-parses the schema internally (so any
   ! unknown-key / range / enum error error-stops here, with the full
   ! report); `schema` is returned for the parameter-doc dumps.
   call read_config(trim(input_file), cfg, schema)
   call validate_config(cfg)
   ! Apply log level from config
   call configure_log_level(cfg%log_level)
   call logger%info("Configuration read from: "//trim(input_file))

   ! Ensure output directory exists. Fortran's NetCDF/open calls don't create
   ! parent directories — mkdir -p is idempotent.
   if (len_trim(cfg%output_dir) > 0) then
      call execute_command_line('mkdir -p "'//trim(cfg%output_dir)//'"', &
                                wait=.true., exitstat=io_stat)
      if (io_stat /= 0) then
         call logger%warning("Could not create output directory: "// &
                             trim(cfg%output_dir))
      end if
   end if

   ! Emit the MOM6-style parameter-documentation dumps from the schema
   ! returned by read_config (already built + strict-parsed).
   if (len_trim(cfg%output_dir) > 0) then
      doc_dir = trim(cfg%output_dir)//"/"
   else
      doc_dir = ""
   end if
   call schema%write_doc_all(doc_dir//"rakali_parameter_doc.all")
   call schema%write_doc_short(doc_dir//"rakali_parameter_doc.short")

   call driver_run(cfg)

end program rakali
