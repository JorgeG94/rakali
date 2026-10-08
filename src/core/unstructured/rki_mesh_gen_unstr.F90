!! Native (Python-free) unstructured mesh generators.
module rki_mesh_gen_unstr
   !! Naive in-Fortran triangulators for idealized domains, so the
   !! unstructured path can be driven end-to-end without the Python +
   !! meshpy/Triangle toolchain (benchmarks, idealized test cases).  The
   !! output (vert_x, vert_y, tri_v) feeds straight into
   !! `mesh_unstr_init_from_triangles`.  These are NOT a substitute for
   !! meshpy on real coastlines -- no boundary conformity, no
   !! bathymetry-driven refinement.
   use rki_constants, only: wp
   implicit none
   private

   public :: gen_regular_triangulation
   public :: gen_delaunay_triangulation

contains

   subroutine gen_regular_triangulation(nx, ny, lx, ly, vert_x, vert_y, tri_v)
      !! Regular rectangular triangulation of [0,lx] x [0,ly]: an
      !! (nx+1) x (ny+1) vertex grid, each quad split into two CCW
      !! triangles (a lower (v00,v10,v01) + an upper (v10,v11,v01)).
      !! Yields (nx+1)(ny+1) vertices and 2*nx*ny triangles.  Vertices are
      !! numbered row-major (i fastest); the triangulation tiles the
      !! rectangle exactly with positively-oriented triangles.
      integer, intent(in) :: nx, ny
         !! Number of cells in x / y (>= 1).
      real(wp), intent(in) :: lx, ly
         !! Domain extents in x / y (m).
      real(wp), allocatable, intent(out) :: vert_x(:), vert_y(:)
      integer, allocatable, intent(out) :: tri_v(:, :)
         !! Triangle connectivity, shape (3, 2*nx*ny); 1-based vertex ids.

      integer :: nverts, ntri
      integer :: i, j, vid, tid
      integer :: v00, v10, v01, v11
      real(wp) :: dx, dy

      nverts = (nx + 1)*(ny + 1)
      ntri = 2*nx*ny
      allocate (vert_x(nverts), vert_y(nverts))
      allocate (tri_v(3, ntri))

      dx = lx/real(nx, wp)
      dy = ly/real(ny, wp)

      ! Vertices: row-major, i fastest.  vid(i,j) = j*(nx+1) + i + 1.
      do j = 0, ny
         do i = 0, nx
            vid = j*(nx + 1) + i + 1
            vert_x(vid) = real(i, wp)*dx
            vert_y(vid) = real(j, wp)*dy
         end do
      end do

      ! Each quad -> two CCW triangles.
      tid = 0
      do j = 0, ny - 1
         do i = 0, nx - 1
            v00 = j*(nx + 1) + i + 1               ! (i,   j)
            v10 = j*(nx + 1) + (i + 1) + 1         ! (i+1, j)
            v01 = (j + 1)*(nx + 1) + i + 1         ! (i,   j+1)
            v11 = (j + 1)*(nx + 1) + (i + 1) + 1   ! (i+1, j+1)

            tid = tid + 1
            tri_v(1, tid) = v00
            tri_v(2, tid) = v10
            tri_v(3, tid) = v01

            tid = tid + 1
            tri_v(1, tid) = v10
            tri_v(2, tid) = v11
            tri_v(3, tid) = v01
         end do
      end do
   end subroutine gen_regular_triangulation

   subroutine gen_delaunay_triangulation(nx, ny, lx, ly, vert_x, vert_y, tri_v)
      !! Acute-Delaunay variant of `gen_regular_triangulation` for the
      !! semi-implicit (TPFA) path.  The plain right-triangle split puts
      !! each circumcentre ON the shared hypotenuse midpoint — the two
      !! triangles of a quad get COINCIDENT circumcentres, `edge_delta = 0`,
      !! and the SI mesh-quality gate rightly rejects the mesh.  Here the
      !! interior vertices are jittered (deterministic hash — reproducible)
      !! and each quad is split along its locally-Delaunay diagonal
      !! (in-circle test), so circumcentres sit strictly inside the acute
      !! triangles and `edge_delta > 0` everywhere.  Boundary vertices stay
      !! ON the box so the domain is exactly rectangular for the wall BC.
      integer, intent(in) :: nx, ny
         !! Number of cells in x / y (>= 1).
      real(wp), intent(in) :: lx, ly
         !! Domain extents in x / y (m).
      real(wp), allocatable, intent(out) :: vert_x(:), vert_y(:)
      integer, allocatable, intent(out) :: tri_v(:, :)
         !! Triangle connectivity, shape (3, 2*nx*ny); 1-based vertex ids.

      real(wp), parameter :: JIT = 0.28_wp
      integer :: nverts, ntri, i, j, vid, tid
      integer :: v00, v10, v01, v11
      real(wp) :: dx, dy, jx, jy

      nverts = (nx + 1)*(ny + 1)
      ntri = 2*nx*ny
      allocate (vert_x(nverts), vert_y(nverts))
      allocate (tri_v(3, ntri))

      dx = lx/real(nx, wp)
      dy = ly/real(ny, wp)

      do j = 0, ny
         do i = 0, nx
            vid = j*(nx + 1) + i + 1
            jx = 0.0_wp
            jy = 0.0_wp
            ! Jitter interior vertices only (keep the boundary on the box).
            if (i > 0 .and. i < nx) jx = JIT*dx*(hash01(i, j, 1) - 0.5_wp)
            if (j > 0 .and. j < ny) jy = JIT*dy*(hash01(i, j, 2) - 0.5_wp)
            vert_x(vid) = real(i, wp)*dx + jx
            vert_y(vid) = real(j, wp)*dy + jy
         end do
      end do

      tid = 0
      do j = 0, ny - 1
         do i = 0, nx - 1
            v00 = j*(nx + 1) + i + 1               ! (i,   j)   bl
            v10 = j*(nx + 1) + (i + 1) + 1         ! (i+1, j)   br
            v01 = (j + 1)*(nx + 1) + i + 1         ! (i,   j+1) tl
            v11 = (j + 1)*(nx + 1) + (i + 1) + 1   ! (i+1, j+1) tr

            if (delaunay_diag_bl_tr(vert_x(v00), vert_y(v00), &
                                    vert_x(v10), vert_y(v10), &
                                    vert_x(v11), vert_y(v11), &
                                    vert_x(v01), vert_y(v01))) then
               tid = tid + 1
               tri_v(:, tid) = [v00, v10, v11]
               tid = tid + 1
               tri_v(:, tid) = [v00, v11, v01]
            else
               tid = tid + 1
               tri_v(:, tid) = [v00, v10, v01]
               tid = tid + 1
               tri_v(:, tid) = [v10, v11, v01]
            end if
         end do
      end do
   end subroutine gen_delaunay_triangulation

   pure function hash01(i, j, s) result(r)
      !! Deterministic pseudo-random value in [0,1) from a 3-integer key.
      integer, intent(in) :: i, j, s
      real(wp) :: r
      integer(kind=8) :: h
      h = int(i, 8)*73856093_8 + int(j, 8)*19349663_8 + int(s, 8)*83492791_8
      h = iand(h, 2147483647_8)
      r = real(modulo(h, 100003_8), wp)/100003.0_wp
   end function hash01

   pure function delaunay_diag_bl_tr(xbl, ybl, xbr, ybr, xtr, ytr, xtl, ytl) &
      result(use_bl_tr)
      !! In-circle test choosing the Delaunay diagonal of the quad
      !! (bl, br, tr, tl) CCW.  Returns .true. to split along bl--tr (the
      !! diagonal that keeps both triangles Delaunay), .false. for br--tl.
      real(wp), intent(in) :: xbl, ybl, xbr, ybr, xtr, ytr, xtl, ytl
      logical :: use_bl_tr
      real(wp) :: ax, ay, bx, by, cx, cy, incirc

      ax = xbl - xtl
      ay = ybl - ytl
      bx = xbr - xtl
      by = ybr - ytl
      cx = xtr - xtl
      cy = ytr - ytl
      incirc = (ax*ax + ay*ay)*(bx*cy - cx*by) &
               - (bx*bx + by*by)*(ax*cy - cx*ay) &
               + (cx*cx + cy*cy)*(ax*by - bx*ay)
      ! (bl,br,tr) CCW => incirc > 0 means tl is inside the circumcircle =>
      ! flip to br--tl; else keep bl--tr.
      use_bl_tr = (incirc <= 0.0_wp)
   end function delaunay_diag_bl_tr

end module rki_mesh_gen_unstr
