! hg-mesa-1108: extra history columns for the HG inner-boundary luminosity L_in(t).
! USE: in a star/work copy, replace "include 'standard_run_star_extras.inc'" in
! src/run_star_extras.f90 by the CONTENTS of $MESA_DIR/include/standard_run_star_extras.inc,
! delete its two routines how_many_extra_history_columns / data_for_extra_history_columns,
! and paste the three routines below in their place (inside the module, after "contains").
! The routines are not compiled on viper; fix any r15140 name mismatch (e.g. eps_nuc_categories
! index i3alf from chem_def) and say so in the NOTE.
! MESA convention: k = 1 is the surface; s%r(k), s%m(k), s%L(k) are at the OUTER face of cell k.

      real(dp) function face_interp(s, x, xq, y) result(yq)
         ! linear interpolation of the face quantity y at x = xq; x decreases with k
         type (star_info), pointer :: s
         real(dp), intent(in) :: x(:), y(:), xq
         integer :: k
         yq = -1d99
         if (xq >= x(1)) then
            yq = y(1); return
         end if
         do k = 1, s%nz - 1
            if (x(k) >= xq .and. x(k+1) < xq) then
               yq = y(k+1) + (y(k) - y(k+1))*(xq - x(k+1))/(x(k) - x(k+1))
               return
            end if
         end do
      end function face_interp

      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         how_many_extra_history_columns = 14
      end function how_many_extra_history_columns

      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ! enclosed masses of our boundaries in profile14 (profile14_crosscheck.txt)
         real(dp), parameter :: m_in3 = 2.627060d0, m_in8 = 2.719887d0   ! Msun
         integer :: k, k3a, khe, nz
         real(dp) :: e3a, x
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         nz = s%nz
         names(1) = 'L_r3';   vals(1) = face_interp(s, s%r(1:nz), 3d0*Rsun, s%L(1:nz))/Lsun
         names(2) = 'L_r8';   vals(2) = face_interp(s, s%r(1:nz), 8d0*Rsun, s%L(1:nz))/Lsun
         names(3) = 'L_m3';   vals(3) = face_interp(s, s%m(1:nz), m_in3*Msun, s%L(1:nz))/Lsun
         names(4) = 'r_m3';   vals(4) = face_interp(s, s%m(1:nz), m_in3*Msun, s%r(1:nz))/Rsun
         names(5) = 'L_m8';   vals(5) = face_interp(s, s%m(1:nz), m_in8*Msun, s%L(1:nz))/Lsun
         names(6) = 'r_m8';   vals(6) = face_interp(s, s%m(1:nz), m_in8*Msun, s%r(1:nz))/Rsun
         ! He-burning shell = cell of maximum triple-alpha energy generation
         k3a = 1; e3a = -1d99
         do k = 1, nz
            if (s%eps_nuc_categories(i3alf, k) > e3a) then
               e3a = s%eps_nuc_categories(i3alf, k); k3a = k
            end if
         end do
         names(7) = 'r_He3a';  vals(7) = s%r(k3a)/Rsun
         names(8) = 'm_He3a';  vals(8) = s%m(k3a)/Msun
         names(9) = 'L_He3a';  vals(9) = s%L(k3a)/Lsun
         names(10) = 'log_eps3a_max'; vals(10) = safe_log10(e3a)
         ! base of the He envelope = outermost face (from the surface down) where X(he4) < 0.5
         khe = nz
         do k = 1, nz
            x = s%xa(s%net_iso(ihe4), k)
            if (x < 0.5d0) then
               khe = k; exit
            end if
         end do
         names(11) = 'r_Hebase'; vals(11) = s%r(khe)/Rsun
         names(12) = 'm_Hebase'; vals(12) = s%m(khe)/Msun
         names(13) = 'L_Hebase'; vals(13) = s%L(khe)/Lsun
         ! L_nu integrated over the star (Lsun), for the energy budget of the contraction
         names(14) = 'L_neu_tot'; vals(14) = s%power_neutrinos
      end subroutine data_for_extra_history_columns
