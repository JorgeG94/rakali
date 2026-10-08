!! Closure -> SI seam kernels: k-epsilon / KPP interface viscosity fills for
!! the semi-implicit multilayer step (unstructured).
module rki_si_unstr_closures
   !! Q1 (SI vertical-mixing closures) for the unstructured backend.
   !! The SI column solve consumes ANY interface viscosity through
   !! `state%si%nu_iface`; closures are pure FILLS of that slot.  This module
   !! holds the flat explicit-shape fill kernels on the `(nz, ncells)` layout:
   !!
   !! * `si_unstr_fill_nu_iface_from_keps` — copy the k-epsilon `keps_nut`
   !!   (interface-shaped (nz+1, n), advanced by the UNCHANGED
   !!   `ml_keps_update_unstr` at SI cadence — the explicit kernel is already
   !!   backward-Euler/Patankar) into `nu_iface`, with the SI convention's
   !!   bed/surface zero rows and the vanished-layer decoupling gate.
   !! * `si_unstr_keps_copy_iface` — interface-array copy for the optional
   !!   `keps_si_picard > 1` save/restore of the pre-step k / eps.
   !! * `si_unstr_build_nu_interface_kpp` — PP81 interior viscosity with the
   !!   KPP (LMD94) shape-function overlay MAX'd in, written to `nu_iface` for
   !!   the IMPLICIT column tridiagonal.  NO explicit-flux CFL clamp — the
   !!   coastal explicit `ml_kpp_viscosity_unstr` clamp would strangle the
   !!   overlay at SI step sizes (dt = 120-480 s); the implicit solve is
   !!   unconditionally stable at any nu.
   !! * `si_unstr_kpp_diffusion_tracers` — S/T IMPLICIT (Thomas) vertical
   !!   diffusion with the PP81-cubic + KPP-overlay diffusivity, plus the
   !!   LMD94 non-local counter-gradient transport applied after the diffusion
   !!   (the coastal explicit order: vdiff-then-nonlocal), reading the
   !!   post-diffusion surface salinity for F_S.  Conserves column tracer
   !!   content (no-flux tridiagonal, conservative nonlocal divergence).
   !!
   !! Layer convention (load-bearing): BOTTOM-UP — k=1 bed, k=nz surface;
   !! interface m sits BELOW layer m (m=1 bed, m=nz+1 surface).  Cell-centric
   !! `do concurrent` over cells; per-column serial work in-body with
   !! NZ_STACK_MAX locals; explicit-shape dummies throughout (NVHPC rule).
#ifdef LFORTRAN_PASSING
   use rki_constants, only: wp, GRAVITY, H_VANISHED, DRY_TOLERANCE
#else
   use rki_constants, only: NZ_STACK_MAX, wp, GRAVITY, H_VANISHED, DRY_TOLERANCE
#endif
   use rki_ml_common_unstr, only: ML_MIN_DEPTH_UNSTR
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rki_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: si_unstr_fill_nu_iface_from_keps
   public :: si_unstr_keps_copy_iface
   public :: si_unstr_build_nu_interface_kpp
   public :: si_unstr_kpp_diffusion_tracers

   real(wp), parameter :: SHEAR2_EPS = 1.0e-20_wp
      !! Shear^2 floor in the interior gradient-Ri (matches PP81 + kpp kernels).

contains

   subroutine si_unstr_fill_nu_iface_from_keps(nz, nzp1, nct, nco, &
                                                    h_layer, keps_nut, nu_iface)
      !! Copy the k-epsilon interface viscosity into the SI column-solve slot.
      !! `keps_nut` and `nu_iface` share the interface convention (m=1 bed,
      !! m=nz+1 surface); interior interfaces copy 1:1; bed and surface rows
      !! are forced to 0 (no flux through bed/surface — wind and drag enter
      !! the column solve separately), and an interface bounded by a vanished
      !! layer decouples (nu = 0), mirroring `si_unstr_build_nu_interface`
      !! so the closure swap changes ONLY the viscosity values.
      integer, intent(in) :: nz
         !! Number of layers.
      integer, intent(in) :: nzp1
         !! nz + 1 (leading dim of interface arrays).
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: h_layer(nz, nct)
         !! Layer thicknesses [m], bottom-up (k=1 bed).
      real(wp), intent(in) :: keps_nut(nzp1, nct)
         !! k-epsilon interface eddy viscosity [m^2/s] (from ml_keps_update_unstr).
      real(wp), intent(out) :: nu_iface(nzp1, nct)
         !! SI interface viscosity [m^2/s], vdiff convention.

      integer :: c, k
      real(wp) :: hk, hkp1

      !$omp target teams distribute parallel do private(k, hk, hkp1)
      do c = 1, nco
         nu_iface(1, c) = 0.0_wp
         nu_iface(nz + 1, c) = 0.0_wp
         do k = 1, nz - 1
            hk = h_layer(k, c)
            hkp1 = h_layer(k + 1, c)
            if (hk <= H_VANISHED .or. hkp1 <= H_VANISHED) then
               nu_iface(k + 1, c) = 0.0_wp
            else
               nu_iface(k + 1, c) = keps_nut(k + 1, c)
            end if
         end do
      end do
      !$omp end target teams distribute parallel do
   end subroutine si_unstr_fill_nu_iface_from_keps

   subroutine si_unstr_keps_copy_iface(nzp1, nct, nco, src, dst)
      !! Interface-array copy (nzp1, nct) for the `keps_si_picard > 1`
      !! save/restore: the Picard loop re-runs `ml_keps_update_unstr` from
      !! the SAME pre-step k / eps with the refreshed nu_t production
      !! coefficient (see local_archive/prototypes/si_closures/DISCRETIZATION_SPEC.md
      !! section 3).
      integer, intent(in) :: nzp1
         !! Interface count (nz + 1).
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: src(nzp1, nct)
         !! Source interface array.
      real(wp), intent(out) :: dst(nzp1, nct)
         !! Destination interface array.

      integer :: c, m

      !$omp target teams distribute parallel do collapse(2)
      do m = 1, nzp1
      do c = 1, nco
         dst(m, c) = src(m, c)
      end do
      end do
      !$omp end target teams distribute parallel do
   end subroutine si_unstr_keps_copy_iface

   subroutine si_unstr_build_nu_interface_kpp(nz, nzp1, nct, nco, &
                                                   h_layer, hu_layer, hv_layer, &
                                                   rho_layer, rho_0, nu_max, &
                                                   nu_bg, ri_alpha, u_star, &
                                                   buoy_flux, bl_depth, nu_iface)
      !! PP81 interior interface viscosity with the KPP (LMD94) surface-
      !! boundary-layer overlay MAX'd in, for the SI implicit column solve:
      !!   nu = max( nu_max/(1+ri_alpha*Ri)^2 + nu_bg,
      !!             h_b * w_s * G(sigma) )        for interfaces inside h_b
      !! with G(sigma) = sigma(1-sigma)^2, sigma = d_iface/h_b, and the
      !! velocity scale w_s = sqrt(u_*^2 + w_*^2), w_*^3 = max(0,-B_0)*h_b.
      !! NO dt clamp: the consumer is the unconditionally-stable implicit
      !! tridiagonal.  Where the overlay is inactive (w_s = 0 / below the
      !! BL) the result is EXACTLY the PP81 `si_unstr_build_nu_interface` value.
      integer, intent(in) :: nz
         !! Number of layers.
      integer, intent(in) :: nzp1
         !! nz + 1 (leading dim of interface arrays).
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: h_layer(nz, nct)
         !! Layer thicknesses [m], bottom-up.
      real(wp), intent(in) :: hu_layer(nz, nct)
         !! Layer x-transports [m^2/s].
      real(wp), intent(in) :: hv_layer(nz, nct)
         !! Layer y-transports [m^2/s].
      real(wp), intent(in) :: rho_layer(nz, nct)
         !! Layer densities [kg/m^3].
      real(wp), intent(in) :: rho_0
         !! Reference density [kg/m^3].
      real(wp), intent(in) :: nu_max
         !! PP81 maximum turbulent viscosity [m^2/s].
      real(wp), intent(in) :: nu_bg
         !! Background viscosity floor [m^2/s].
      real(wp), intent(in) :: ri_alpha
         !! PP81 Richardson tuning coefficient.
      real(wp), intent(in) :: u_star(nct)
         !! Surface friction velocity [m/s] (from ml_kpp_boundary_layer_unstr).
      real(wp), intent(in) :: buoy_flux(nct)
         !! Surface buoyancy flux B_0 [m^2/s^3]; < 0 destabilizing.
      real(wp), intent(in) :: bl_depth(nct)
         !! Diagnosed boundary-layer depth [m].
      real(wp), intent(out) :: nu_iface(nzp1, nct)
         !! Interface viscosity [m^2/s], vdiff convention (bed/surface 0).

      integer :: c, k
      real(wp) :: h_col, hk, hkp1, h_avg, uk, ukp1, vk, vkp1
      real(wp) :: du, dv, drho, ri_num, denom_ri, nu_t
      real(wp) :: h_b, b_0, us, wstar3_val, w_star, w_s
      real(wp) :: d_below, d_iface, sig, g_shape, nu_kpp

      !$omp target teams distribute parallel do private(k, h_col, hk, hkp1, h_avg, uk, ukp1, vk, vkp1, du, dv, drho, ri_num, denom_ri, nu_t, h_b, b_0, us, wstar3_val, w_star, w_s, d_below, d_iface, sig, g_shape, nu_kpp)
      do c = 1, nco
         nu_iface(1, c) = 0.0_wp
         nu_iface(nz + 1, c) = 0.0_wp

         h_col = 0.0_wp
         do k = 1, nz
            h_col = h_col + h_layer(k, c)
         end do

         if (h_col > ML_MIN_DEPTH_UNSTR) then
            h_b = bl_depth(c)
            b_0 = buoy_flux(c)
            us = u_star(c)
            wstar3_val = max(0.0_wp, -b_0)*h_b
            w_star = wstar3_val**(1.0_wp/3.0_wp)
            w_s = sqrt(us*us + w_star*w_star)

            d_below = 0.0_wp
            do k = 1, nz - 1
               hk = h_layer(k, c)
               hkp1 = h_layer(k + 1, c)
               h_avg = 0.5_wp*(hk + hkp1)
               d_below = d_below + hk
               d_iface = h_col - d_below

               if (hk > H_VANISHED .and. hkp1 > H_VANISHED) then
                  uk = hu_layer(k, c)/hk
                  ukp1 = hu_layer(k + 1, c)/hkp1
                  vk = hv_layer(k, c)/hk
                  vkp1 = hv_layer(k + 1, c)/hkp1
                  du = ukp1 - uk
                  dv = vkp1 - vk
                  drho = rho_layer(k, c) - rho_layer(k + 1, c)
                  ri_num = GRAVITY*(drho/rho_0)*h_avg/(du*du + dv*dv + SHEAR2_EPS)
                  if (ri_num < 0.0_wp) ri_num = 0.0_wp
                  denom_ri = 1.0_wp + ri_alpha*ri_num
                  nu_t = nu_max/(denom_ri*denom_ri) + nu_bg
                  ! KPP overlay inside the boundary layer (max, LMD94).
                  if (w_s > 0.0_wp .and. h_b > 0.0_wp .and. d_iface < h_b) then
                     sig = d_iface/h_b
                     g_shape = sig*(1.0_wp - sig)*(1.0_wp - sig)
                     nu_kpp = h_b*w_s*g_shape
                     if (nu_kpp > nu_t) nu_t = nu_kpp
                  end if
                  nu_iface(k + 1, c) = nu_t
               else
                  nu_iface(k + 1, c) = 0.0_wp
               end if
            end do
         else
            do k = 1, nz - 1
               nu_iface(k + 1, c) = 0.0_wp
            end do
         end if
      end do
      !$omp end target teams distribute parallel do
   end subroutine si_unstr_build_nu_interface_kpp

   subroutine si_unstr_kpp_diffusion_tracers(nz, nzp1, nct, nco, h_layer, &
                                                  hu_layer, hv_layer, hS_layer, &
                                                  hT_layer, rho_layer, rho_0, &
                                                  nu_max, kappa_S_bg, kappa_T_bg, &
                                                  ri_alpha, cs_nonlocal, heat_flux, &
                                                  evap_minus_precip, cp, u_star, &
                                                  buoy_flux, bl_depth, dt)
      !! S/T vertical diffusion with the KPP overlay for the SI step:
      !! IMPLICIT (Thomas) column solve with the PP81-CUBIC interior
      !! diffusivity + the KPP shape-function overlay MAX'd in (same
      !! overlay as momentum; no dt clamp — implicit), then the LMD94
      !! non-local counter-gradient transport applied as a conservative
      !! interface divergence AFTER the diffusion (the coastal explicit
      !! order: vdiff-then-nonlocal), reading the post-diffusion surface
      !! salinity for F_S.  gamma = 0 at bed and surface => the divergence
      !! conserves the column integral exactly.  Conserves column tracer
      !! content (no-flux tridiagonal).
      integer, intent(in) :: nz
         !! Number of layers.
      integer, intent(in) :: nzp1
         !! nz + 1 (leading dim of interface arrays).
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: h_layer(nz, nct)
         !! Layer thicknesses [m], bottom-up.
      real(wp), intent(in) :: hu_layer(nz, nct)
         !! Layer x-transports [m^2/s] (shear for the PP81 Ri).
      real(wp), intent(in) :: hv_layer(nz, nct)
         !! Layer y-transports [m^2/s].
      real(wp), intent(inout) :: hS_layer(nz, nct)
         !! Salinity content [PSU m], updated in place.
      real(wp), intent(inout) :: hT_layer(nz, nct)
         !! Temperature content [degC m], updated in place.
      real(wp), intent(in) :: rho_layer(nz, nct)
         !! Layer densities [kg/m^3].
      real(wp), intent(in) :: rho_0
         !! Reference density [kg/m^3].
      real(wp), intent(in) :: nu_max
         !! PP81 maximum turbulent diffusivity scale [m^2/s].
      real(wp), intent(in) :: kappa_S_bg
         !! Background salinity diffusivity [m^2/s].
      real(wp), intent(in) :: kappa_T_bg
         !! Background temperature diffusivity [m^2/s].
      real(wp), intent(in) :: ri_alpha
         !! PP81 Richardson tuning coefficient.
      real(wp), intent(in) :: cs_nonlocal
         !! Non-local transport coefficient C_s (LMD94 eq 20).
      real(wp), intent(in) :: heat_flux(nct)
         !! Surface heat flux [W/m^2] (positive warms), per cell.
      real(wp), intent(in) :: evap_minus_precip(nct)
         !! Net evaporation [m/s], per cell.
      real(wp), intent(in) :: cp
         !! Specific heat of water [J/(kg K)].
      real(wp), intent(in) :: u_star(nct)
         !! Surface friction velocity [m/s].
      real(wp), intent(in) :: buoy_flux(nct)
         !! Surface buoyancy flux B_0 [m^2/s^3]; < 0 destabilizing.
      real(wp), intent(in) :: bl_depth(nct)
         !! Diagnosed boundary-layer depth [m].
      real(wp), intent(in) :: dt
         !! Time step [s].

      integer :: c, k
      real(wp) :: h_col, hk, hkp1, h_avg, uk, ukp1, vk, vkp1
      real(wp) :: du, dv, drho, ri_num, denom_ri, kap
      real(wp) :: h_b, b_0, us, wstar3_val, w_star, w_s
      real(wp) :: d_below, d_iface, sig, g_shape, nu_kpp
      real(wp) :: h_top, s_top, f_t_kin, f_s_kin, gamma_t, gamma_s
      real(wp) :: cd_val, cu_val, dz_dn, dz_up
      ! Per-site gfortran 15.1 gotcha: use separate locals for each write site
      ! inside do concurrent to avoid the "local reassigned across if/else"
      ! corruption bug.
      real(wp) :: kapS_k, kapT_k
      real(wp) :: kapS(NZ_STACK_MAX), kapT(NZ_STACK_MAX)
      real(wp) :: aa(NZ_STACK_MAX), bb(NZ_STACK_MAX), cc(NZ_STACK_MAX)
      real(wp) :: rS(NZ_STACK_MAX), rT(NZ_STACK_MAX), cpv(NZ_STACK_MAX)
      logical :: destabilizing

      if (nz < 2) return

      !$omp target teams distribute parallel do private(k, h_col, hk, hkp1, h_avg, uk, ukp1, vk, vkp1, du, dv, drho, ri_num, denom_ri, kap, h_b, b_0, us, wstar3_val, w_star, w_s, d_below, d_iface, sig, g_shape, nu_kpp, h_top, s_top, f_t_kin, f_s_kin, gamma_t, gamma_s, cd_val, cu_val, dz_dn, dz_up, kapS_k, kapT_k, kapS, kapT, aa, bb, cc, rS, rT, cpv, destabilizing)
      do c = 1, nco

         h_col = 0.0_wp
         do k = 1, nz
            h_col = h_col + h_layer(k, c)
         end do

         if (h_col > ML_MIN_DEPTH_UNSTR) then
            h_b = bl_depth(c)
            b_0 = buoy_flux(c)
            us = u_star(c)
            wstar3_val = max(0.0_wp, -b_0)*h_b
            w_star = wstar3_val**(1.0_wp/3.0_wp)
            w_s = sqrt(us*us + w_star*w_star)
            destabilizing = (b_0 < 0.0_wp)

            ! ---- interface diffusivities: PP81 cubic + KPP overlay ----
            ! kapS/kapT(m): interface m convention (1 = bed = 0 flux).
            kapS(1) = 0.0_wp
            kapT(1) = 0.0_wp
            kapS(nz + 1) = 0.0_wp
            kapT(nz + 1) = 0.0_wp
            d_below = 0.0_wp
            do k = 1, nz - 1
               hk = h_layer(k, c)
               hkp1 = h_layer(k + 1, c)
               h_avg = 0.5_wp*(hk + hkp1)
               d_below = d_below + hk
               d_iface = h_col - d_below
               if (hk > H_VANISHED .and. hkp1 > H_VANISHED) then
                  uk = hu_layer(k, c)/hk
                  ukp1 = hu_layer(k + 1, c)/hkp1
                  vk = hv_layer(k, c)/hk
                  vkp1 = hv_layer(k + 1, c)/hkp1
                  du = ukp1 - uk
                  dv = vkp1 - vk
                  drho = rho_layer(k, c) - rho_layer(k + 1, c)
                  ri_num = GRAVITY*(drho/rho_0)*h_avg/(du*du + dv*dv + SHEAR2_EPS)
                  if (ri_num < 0.0_wp) ri_num = 0.0_wp
                  denom_ri = 1.0_wp + ri_alpha*ri_num
                  kap = nu_max/(denom_ri*denom_ri*denom_ri)
                  nu_kpp = 0.0_wp
                  if (w_s > 0.0_wp .and. h_b > 0.0_wp .and. d_iface < h_b) then
                     sig = d_iface/h_b
                     g_shape = sig*(1.0_wp - sig)*(1.0_wp - sig)
                     nu_kpp = h_b*w_s*g_shape
                  end if
                  kapS_k = max(kap + kappa_S_bg, nu_kpp)
                  kapT_k = max(kap + kappa_T_bg, nu_kpp)
               else
                  kapS_k = 0.0_wp
                  kapT_k = 0.0_wp
               end if
               kapS(k + 1) = kapS_k
               kapT(k + 1) = kapT_k
            end do

            ! ---- implicit (Thomas) diffusion of S concentration ----
            do k = 1, nz
               hk = h_layer(k, c)
               if (hk > H_VANISHED) then
                  cd_val = 0.0_wp
                  if (k > 1) then
                     if (h_layer(k - 1, c) > H_VANISHED) then
                        dz_dn = 0.5_wp*(h_layer(k - 1, c) + hk)
                        cd_val = dt/(hk*dz_dn)
                     end if
                  end if
                  cu_val = 0.0_wp
                  if (k < nz) then
                     if (h_layer(k + 1, c) > H_VANISHED) then
                        dz_up = 0.5_wp*(hk + h_layer(k + 1, c))
                        cu_val = dt/(hk*dz_up)
                     end if
                  end if
                  aa(k) = -cd_val*kapS(k)
                  cc(k) = -cu_val*kapS(k + 1)
                  bb(k) = 1.0_wp + cd_val*kapS(k) + cu_val*kapS(k + 1)
                  rS(k) = hS_layer(k, c)/hk
               else
                  aa(k) = 0.0_wp
                  cc(k) = 0.0_wp
                  bb(k) = 1.0_wp
                  if (hk > 0.0_wp) then
                     rS(k) = hS_layer(k, c)/hk
                  else
                     rS(k) = 0.0_wp
                  end if
               end if
            end do
            call si_unstr_tridiag_noflux(aa, bb, cc, rS, cpv, nz)
            do k = 1, nz
               hS_layer(k, c) = rS(k)*h_layer(k, c)
            end do

            ! ---- implicit (Thomas) diffusion of T concentration ----
            do k = 1, nz
               hk = h_layer(k, c)
               if (hk > H_VANISHED) then
                  cd_val = 0.0_wp
                  if (k > 1) then
                     if (h_layer(k - 1, c) > H_VANISHED) then
                        dz_dn = 0.5_wp*(h_layer(k - 1, c) + hk)
                        cd_val = dt/(hk*dz_dn)
                     end if
                  end if
                  cu_val = 0.0_wp
                  if (k < nz) then
                     if (h_layer(k + 1, c) > H_VANISHED) then
                        dz_up = 0.5_wp*(hk + h_layer(k + 1, c))
                        cu_val = dt/(hk*dz_up)
                     end if
                  end if
                  aa(k) = -cd_val*kapT(k)
                  cc(k) = -cu_val*kapT(k + 1)
                  bb(k) = 1.0_wp + cd_val*kapT(k) + cu_val*kapT(k + 1)
                  rT(k) = hT_layer(k, c)/hk
               else
                  aa(k) = 0.0_wp
                  cc(k) = 0.0_wp
                  bb(k) = 1.0_wp
                  if (hk > 0.0_wp) then
                     rT(k) = hT_layer(k, c)/hk
                  else
                     rT(k) = 0.0_wp
                  end if
               end if
            end do
            call si_unstr_tridiag_noflux(aa, bb, cc, rT, cpv, nz)
            do k = 1, nz
               hT_layer(k, c) = rT(k)*h_layer(k, c)
            end do

            ! ---- non-local counter-gradient transport (post-diffusion) ----
            if (destabilizing .and. w_s > 0.0_wp .and. h_b > 0.0_wp) then
               h_top = h_layer(nz, c)
               if (h_top > DRY_TOLERANCE) then
                  s_top = hS_layer(nz, c)/h_top
               else
                  s_top = 0.0_wp
               end if
               f_t_kin = heat_flux(c)/(rho_0*cp)
               f_s_kin = evap_minus_precip(c)*s_top

               d_below = 0.0_wp
               do k = 1, nz - 1
                  hk = h_layer(k, c)
                  hkp1 = h_layer(k + 1, c)
                  d_below = d_below + hk
                  d_iface = h_col - d_below
                  if (hk > H_VANISHED .and. hkp1 > H_VANISHED .and. &
                      d_iface < h_b) then
                     sig = d_iface/h_b
                     g_shape = sig*(1.0_wp - sig)*(1.0_wp - sig)
                     gamma_t = cs_nonlocal*g_shape*f_t_kin
                     gamma_s = cs_nonlocal*g_shape*f_s_kin
                     hT_layer(k, c) = hT_layer(k, c) + dt*gamma_t
                     hT_layer(k + 1, c) = hT_layer(k + 1, c) - dt*gamma_t
                     hS_layer(k, c) = hS_layer(k, c) + dt*gamma_s
                     hS_layer(k + 1, c) = hS_layer(k + 1, c) - dt*gamma_s
                  end if
               end do
            end if
         end if
      end do
      !$omp end target teams distribute parallel do
   end subroutine si_unstr_kpp_diffusion_tracers

   pure subroutine si_unstr_tridiag_noflux(a, b, c, rhs, cpv, n)
      !$omp declare target
      !! In-place Thomas solve over rows 1..n with the boundary closure
      !! baked into a/b/c by the caller (the no-flux tridiagonal closure).
      !! `rhs` returns the solution.
      integer, intent(in) :: n
      real(wp), intent(in) :: a(n), b(n), c(n)
      real(wp), intent(inout) :: rhs(n), cpv(n)
      integer :: m
      real(wp) :: denom

      cpv(1) = c(1)/b(1)
      rhs(1) = rhs(1)/b(1)
      do m = 2, n
         denom = b(m) - a(m)*cpv(m - 1)
         cpv(m) = c(m)/denom
         rhs(m) = (rhs(m) - a(m)*rhs(m - 1))/denom
      end do
      do m = n - 1, 1, -1
         rhs(m) = rhs(m) - cpv(m)*rhs(m + 1)
      end do
   end subroutine si_unstr_tridiag_noflux

end module rki_si_unstr_closures
