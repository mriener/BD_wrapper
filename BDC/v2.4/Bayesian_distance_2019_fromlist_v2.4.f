      program Bayesian_distance_2019

c     Generates a distance probability density function for a source
c     given its (ell,bee,Vlsr) values.  It uses spiral arm models,
c     kinematic distances, nearby sources with parallax measurement,
c     and Galactic latitude inference.

c     This is for batch processing of survey data given (l,b,Vlsr), but
c     without parallax/p.m. information.  The program has the code to
c     use proper motion information, but it is hardwired to ignore it
c     for this application.

c     Input files:
c        galaxy_data_Univ.inp      ! Galaxy model parameters for kinematic distances, and
c                                  !    reliable longitude range for spiral arm distances

c        probability_controls.inp  ! Weights for each of the 5 probability density components

c        parallax_data.inp         ! Parallax source data

c        sources_info.inp          ! Source name, RA/longitude, Dec/latitude, Vlsr, +/-, P(far)

c                                  ! Individual arm (l,b,v, R, azimuth, D) traces
c        Out_lbvRBD.2019           ! Outer arm
c        Per_lbvRBD.2019           ! Perseus arm
c        Loc_lbvRBD.2019           ! Local Arm
c        CrN_lbvRBD.2019           ! Carina arm (near portion)
c        CrF_lbvRBD.2019           ! Carina arm (far portion)
c        SgN_lbvRBD.2019           ! Sagittarius arm (near portion)
c        SgF_lbvRBD.2019           ! Sagittarius arm (far portion)
c        CtN_lbvRBD.2019           ! Centaurus-Crux arm (near portion)
c        CtF_lbvRBD.2019           ! Centaurus-Crux arm (far portion)
c        ScN_lbvRBD.2019           ! Scutum arm (near portion)
c        ScF_lbvRBD.2019           ! Scutum arm (far portion)
c        OSC_lbvRBD.2019           ! Q1 and Q2 "Outer Scutum-Crux" arm
c        N1N_lbvRBD.2019           ! Norma arm Q1 near
c        N1F_lbvRBD.2019           ! Norma arm Q1 far
c        N4N_lbvRBD.2019           ! Norma arm Q4 near
c        N4F_lbvRBD.2019           ! Norma arm Q4 far
c        3kN_lbvRBD.2019           ! near 3-kpc arm
c        3kF_lbvRBD.2019           ! far 3-kpc arm
c        LoS_lbvRDB.2019           ! Local Spur
c        AqS_lbvRBD.2019           ! Aquila Spur
c        CnN_lbvRBD.2019           ! near Connecting arm
c        CnX_lbvRBD.2019           ! extension (bifurcation) of Connecting arm
c        AqR_lbvRBD.2019           ! Aquila Ridge


c     Output files:
c        summary.prt                            ! one line summary per source

c        Following are not made unless "lu_out = 7" specified below
c        sourcename.prt                         ! summary print out for one source
c        sourcename_final_distance_pdf.dat      ! full distance pdf values
c        sourcename_arm_latitude_pdf.dat        ! spiral-arm * latitude pdf
c        sourcename_kinematic_distance_pdf.dat  ! kinematic distance pdf
c        sourcename_parallaxes_pdf.dat          ! parallax-source association pdf
c        sourcename_pm_ell_distance_pdf.dat     ! proper motion (ell) pdf
c        sourcename_pm_bee_distance_pdf.dat     ! proper motion (ell) pdf
c        sourcename_arm_ranges.dat              ! distance ranges for possible arm associations


c     Printout:
c        one line per input source, including distance estimate for the density function
c            component with the greatest integrated probability


      implicit real*8 (a-h,o-z)

      character*48  control_file,sources_file,srcprt_file,summary_file
      character*48  data_file, out_file
      character*32  citation(1000), cite, pdf_name
      character*14  src, ref_src(1000), stripped_src
      character*12  near_far, ref, ref_arm(1000)
      character*13  Gname
      character*1   asign, bsign, M, arm, s1, nf, c1

      real*8        ref_ell(1000), ref_bee(1000)
      real*8        ref_vlsr(1000), ref_vunc(1000)
      real*8        ref_par(1000), ref_punc(1000)

      real*8        dist_prior(1001), prob_dist(1001)
      real*8        par_bins(1001), dist_bins(1001)
      real*8        prob_arm(1001), prob_lat(1001), prob_armlat(1001)
      real*8        prob_Dpm_ell(1001),prob_Dpm_bee(1001),prob_Dk(1001)
      real*8        prob_MWmodel(1001)

      real*8        peaks(25), peaks_low(25), peaks_high(25)
      real*8        peak_prob(25), peaks_width(25), prob_int(25)
      real*8        peak_dist(25), peak_dunc(25), peak_int(25)

      logical       use_peak(25)

      real*8        params(77), new_params(77), param_sigmas(77)
      character*16  parnames(77)
      integer       paramids(77)

c     Arm information arrays
      real*8        arm_probabilities(29)
      real*8        arm_ell(29,300), arm_bee(29,300), arm_vel(29,300)
      real*8        arm_Rgc(29,300), arm_beta(29,300), arm_dist(29,300)
      real*8        seg_ell_min(29,2), seg_ell_max(29,2)
      real*8        d_store(29), b_store(29)
      integer       iarm_entries(29)
      character*12  a_store(29)

      character*48  arm_file, galaxy_file

      character*12  arm_assigned, arm_indicated, unknown, p2_arm
      character*12  a_name, arm_name(29), arm_max, arm_max_bin(1001)

      character*20  extra_info
      character*2   questionable, q2

      logical       writeout, found_arm, accept, odd_source

c     Set some maximum dimension values, switches, and logical unit numbers...
      max_num_parallaxes = 1000  ! max number of parallaxes entered
      max_num_lbvs   = 300       ! max number of (l,b,v,R,beta,D) values to define an arm segment
      max_num_params = 77        ! max number of Gaussian fit parameters used for final PDF
      max_num_peaks  = 25        ! max number of peaks in final PDF
c     max_num_arms   = 29        ! max number of arm segments (hardwired in subroutine)

      lu_noprint   =-1
      lu_print     = 6

      lu_out       = 7           ! Make multiple PDF files for each source
c     lu_out       =-1           ! Don't make multiple PDF files for each source

      lu_data      = 8
      lu_control   = 9
      lu_sources   =10
      lu_srcprt    =12
      lu_summary   =99

      pi = 4.d0 * atan(1.d0)
      twopi = 2.d0 * pi
      factor = (1.d0/twopi)**1.5d0
      deg_to_rad = pi/180.d0
      hr_to_rad  = pi/12.d0

cc      write (lu_print,1000)
 1000 format(' Bayesian distance estimator: 2019: September 24, 2019')

c     =============================================================
c     Important numerical parameters
      num_bins = 1001         ! Dimension limit for par_prob array
      bin_size = 0.025d0      ! kpc

c     To calculate prob of being in an arm based on (l,b,v)-traces,
c     add Virial velocity uncertainty to other uncertainties (measured, peculiar)
      sig_vel     = 5.d0      ! km/s

c     Maximum difference in longitude to bother calculating arm prob
      ell_dif_max = 5.d0      ! deg

c     To calculate probability based on arm width vs Radius
c     Arm width formula  W = 0.17 + (R_kpc-3.5)*0.036 kpc for R>3.5kpc
c     (which comes from W = 0.336 + (R_kpc-8.15)*0.036; see Fig. 2 of Reid+2019)
      width_min  = 0.17d0     ! kpc  (Gaussian sigma)
      width_Rref = 3.5d0      ! kpc
      width_slope= 0.036d0    ! kpc/kpc

c     To calculate prob of being off the (warped) plane
c     Arm z-width formula  W_z = 0.03 + (R-7.0)*0.036 kpc for R>7kpc
      sigz    = 0.03d0               ! kpc      (1-sig z-height for Rgc<Rsigz)
      Rsigz   = 7.0d0                ! kpc
      sigzdot = 0.036d0              ! kpc/kpc  (slope of z-height with Rgc)

c     To calcuate prob of matching a source with one with measured parallax
      sig_GMC = 0.05d0        ! kpc (1-sig radius difference for sources in a GMC)

c     Smooth prob(d) based on arm assignment (which is "lumpy")
      smooth_kpc = 0.5d0      ! kpc (smooth dist_prob over this length)

c     =============================================================
c               Read in 4 logical program controls...
      call get_controls ( lu_control, lu_print,
     +                    P_max_SA, P_max_KD, P_max_GL, P_max_PS,
     +                    P_max_PM )

c     =============================================================
c               Read in Galactic/Solar parameters...
      galaxy_file = 'galaxy_data_Univ.inp'
      call galaxy_parameters_Univ( lu_data, galaxy_file, lu_print,
     +           Ro, a1, a2, a3, Uo, Vo, Wo, Us, Vs, Ws,
     +           glong_min, glong_max )

c     =============================================================
c               Read in spiral arm segment data...
      call get_arm_segments ( lu_data, lu_print, max_num_lbvs,
     +                        arm_name, seg_ell_min, seg_ell_max,
     +                        arm_ell, arm_bee, arm_vel,
     +                        arm_Rgc, arm_beta, arm_dist,
     +                        iarm_entries, num_arms )


c     =============================================================
c               Read in trig parallax results...
      call get_parallaxes ( lu_control, lu_print,
     +                      max_num_parallaxes, num_parallaxes,
     +                      ref_src, ref_arm, ref_ell,
     +                      ref_bee, ref_vlsr, ref_vunc,
     +                      ref_par, ref_punc )

c     =============================================================
c               Read in target sources for distance PDF...
c     Get source information from one file (one line per source)...
      sources_file = 'sources_info.inp'
      call open_ascii_file ( lu_sources, sources_file, lu_print )

      ieof = 0
      do while ( ieof .ge. 0 )

c      Check for "commented-out" source line (ie, starting with a "!")
       read (lu_sources,*,iostat=ieof) s1

       if ( s1.ne.'!' .and. ieof.ge.0 ) then

         backspace (unit=lu_sources)
         read (lu_sources,*) src, coord1, coord2, v_lsr, v_lsr_unc,
     +                       far_prob

c     Open summary output file...
c        Strip blanks out of source name in output file naming
      call strip_blanks_14 ( src, stripped_src, nch_src )
      write (summary_file,1180) stripped_src(1:nch_src)
 1180 format(a,'_summary.prt')
      open ( unit=lu_summary, file=summary_file )
c     Document summary output
      write (lu_summary,1100)
 1100 format('! Parallax-based distance estimator: Version 2',/
     +       '!',29x,'--------- 1st Peak --------      ',
     +       '-------- 2nd Peak ---------       ',/
     +       '! Long.   Lat.   Vlsr  +/-    ',
     +       'Dist.  +/-  Integrated  Arm      ',
     +       'Dist.  +/-  Integrated  Arm',/
     +       '! (deg)  (deg)  (km/s)        ',
     +       '(kpc)       Probability          ',
     +       '(kpc)       Probability          ')
cc      write (lu_print,1100)

c        Set source characteristics that are not entered in this version
c        (ie, assume proper motions are not measured)
         pm_x      = 0.d0
         pm_x_unc  = 0.d0
         pm_y      = 0.d0
         pm_y_unc  = 0.d0

c        Decide if input coordinates are (RA,Dec) or (ell,bee),
c        since input RA,Dec are in hhmmss and dd''"" formats.
         if ( coord1.gt.360.d0 .and. abs(coord2).gt.90.d0 ) then
c              Almost certainly (RA,Dec), so convert to (ell,bee)
               call hmsrad (coord1, ra_radians)
               call dmsrad (coord2, dec_radians)
               ra_hr   = ra_radians / hr_to_rad        ! hours
               dec_deg = dec_radians / deg_to_rad      ! deg
               call radec_to_galactic (ra_hr, dec_deg,
     +                                 ell, bee )
            else
c              Almost certainly (ell,bee) in degrees
               ell = coord1
               bee = coord2
c              Calculate (RA,Dec) from (ell,bee)
               call ellbee_to_radec ( ell, bee,
     +              ra_hr, dec_deg, ra_hhmmss, dec_ddmmss)
         endif

c        Make sure Prob(Far) ranges from 0 -> 1
         if ( far_prob .lt. 0.d0 ) far_prob = 0.d0
         if ( far_prob .gt. 1.d0 ) far_prob = 1.d0

c        Check if source in desired Galactic longitude range...
         call longitude_range ( ell, glong_min, glong_max,
     +                          accept )

c        Strip blanks out of source name in output file naming
         call strip_blanks_14 ( src,
     +                          stripped_src, nch_src )

         write (srcprt_file,1130) stripped_src(1:nch_src)
 1130    format(a,'.prt')
         open ( unit=lu_srcprt, file=srcprt_file )

c        Calculate (pm_Glong,pm_Glat) from (pm_x,pm_y)
c        (NB: increase shifts to avoid roundoff errors)
         x_motion = pm_x * 1000.d0                            ! uas in 1 yr
         y_motion = pm_y * 1000.d0                            ! uas in 1 yr
         call xy_to_galactic ( ra_hr, dec_deg, x_motion, y_motion,
     +                         ell_motion, bee_motion )
         pm_ell = ell_motion / 1000.d0                        ! mas/yr
         pm_bee = bee_motion / 1000.d0                        ! mas/yr

c        Calculate uncertainties in these Galactic motion components
         x_motion_unc = pm_x_unc * 1000.d0                    ! uas in 1 yr
         y_motion_unc = pm_y_unc * 1000.d0                    ! uas in 1 yr
         call xy_to_galactic ( ra_hr, dec_deg,
     +                         x_motion_unc, y_motion_unc,
     +                         ell_motion_unc, bee_motion_unc )
         pm_ell_unc = ell_motion_unc / 1000.d0                ! mas/yr
         pm_bee_unc = bee_motion_unc / 1000.d0                ! mas/yr

         write (lu_srcprt,1140) src, ell, bee, far_prob,
     +          v_lsr,v_lsr_unc, pm_ell,pm_ell_unc, pm_bee,pm_bee_unc
 1140    format(//'================================================',
     +            '================================',
     +      /'Source          Long.    Lat.    Pfar    Vlsr   +/- ',
     +       '  pm_ell  +/-   pm_bee  +/-',
     +      /a14,2f8.3,f7.2,2x,f7.2,f6.2,2x,2f6.2,2x,2f6.2,/1x)

c        ----------------------------------------------------------------
c        Pre-calculate non-distance based probability terms:
c              P(arm) and near/far kinematic distance values

c        Calculate P(arm) to be used later in P(d)=P(d|arm)*P(arm)
c        Use (ell,z,vel) to assign arm segment probabilities
         call assign_arm_probabilities ( ell, bee, v_lsr, v_lsr_unc,
     +            ell_dif_max,
     +            num_arms, arm_ell, arm_bee, arm_vel, iarm_entries,
     +            arm_Rgc, arm_beta, arm_dist,
     +            width_min, width_Rref, width_slope,
     +            sigz, Rsigz, sigzdot,
     +            seg_ell_min, seg_ell_max,
     +            sig_ell, sig_vel,
     +            arm_probabilities )

c        Print out absolute probabilities.
         write (lu_srcprt,1147)
 1147    format(' Arm  Prob(arm|l,b,v)')
         do n_a = 1, num_arms
            write (lu_srcprt,1148) arm_name(n_a),
     +                             arm_probabilities(n_a)
 1148       format(1x,a3,f13.5)
         enddo

c        Get "standard" kinematic distance(s) for informational use only.
c        (Later Prob_KD will be directly calculated from velocity differences,
c        avoiding complications associated with allowing for velocity uncertainties.)
         farnear_flag = 0.d0                              ! near distance flag
         call calc_Dk_Univ ( lu_out,
     +               Ro, a1, a2, a3, Uo, Vo, Wo, Us, Vs, Ws,
     +               ell, bee, farnear_flag, v_lsr,
     +               v_lsr_rev, Dk_near )

         write (lu_srcprt,1149)
 1149    format(' ')
         if ( Dk_near .gt. 0.d0 ) then
            write (lu_srcprt,1150) Dk_near
 1150       format(' Kinematic distance(s):',f7.2)
         endif

         farnear_flag = 1.d0                              ! far distance flag
         call calc_Dk_Univ ( lu_out,
     +               Ro, a1, a2, a3, Uo, Vo, Wo, Us, Vs, Ws,
     +               ell, bee, farnear_flag, v_lsr,
     +               v_lsr_rev, Dk_far )

         if ( Dk_far .gt. Dk_near ) then
            write (lu_srcprt,1150) Dk_far
         endif

c        Check for pathalogical values...
         if ( Dk_near.le.0.d0 .and. Dk_far.le.0.d0 ) then
            write (lu_srcprt,1160) src, v_lsr, ell
 1160       format(/' Source ',a12,': unlikely V(LSR),longitude pair (',
     +              2f8.2,').',
     +             /' Resetting kinematic distance to 5 kpc',
     +              ' with large uncertainty.')
            Dk_near       = 5.d0                              ! kpc
            Dk_far        = 0.d0                              ! kpc
         endif

c        ----------------------------------------------------------------
c           Start distance-based probability calculations
c        ----------------------------------------------------------------

c        Define distance bins and get warping model
c        Zero distance probability array and define bins
c        (eg, 0.025 kpc steps from 0.025 to 25.025 kpc with 1001 bins)
         do n_ps = 1, num_bins
            dist_prior(n_ps) = 0.d0
            dist_bins(n_ps)  = n_ps * bin_size                ! kpc
c           Generate 1/d parallax bins
            par_bins(n_ps)   = 1.d0 / dist_bins(n_ps)         ! mas
         enddo

c        Get warping model values for distance bins along ray from Sun through source...
         call get_warping ( num_bins, Ro, num_arms, iarm_entries,
     +                      arm_name, arm_probabilities,
     +                      arm_ell, arm_bee, arm_vel,
     +                      arm_Rgc, arm_beta, arm_dist,
     +                      ell, bee, dist_bins,
     +                      n_bees, d_store, b_store, a_store)

c        =========================================================
c        First calculate parallax-source based Prob_PS(d).
c        This requires examining hundreds of parallax sources and
c        calculating dist_prob(n) for each potential matching source.
c        Doing it now is much more efficient than later and going through
c        each distance bin and then searching for parallax-source match.

c        Find "nearby" in (l,b,v) parallax sources
         write (lu_srcprt,1170)
 1170    format(/' Priors from other source parallaxes:',
     +          /' Source         Vlsr      Parallax   +/-',
     +           '      Arm        Weight')
         n_pf = 0
         sum_weight = 0.d0

c        Check that we have some parallax sources and that we want to use them
         if ( num_parallaxes.ge.1 .and. P_max_PS.gt.0.d0 ) then
            do n_p = 1, num_parallaxes

c              Temporary variables for this parallax source...
               dist_par = 1.d0 / ref_par(n_p)
               ell_par  = ref_ell(n_p)
               bee_par  = ref_bee(n_p)
               vlsr_par = ref_vlsr(n_p)

c              Calculate differences between target and parallax source
               dell = abs( ell - ell_par )                       ! deg
               dbee = abs( bee - bee_par )                       ! deg
               dvel = abs( v_lsr - vlsr_par )                    ! km/s

c              Convert (longitude,latitude) to (rho,zee) in kpc units,
c              using the parallax-source's measured parallax distance
               drho = dell*deg_to_rad * dist_par                 ! kpc
               dzee = dbee*deg_to_rad * dist_par                 ! kpc

c              Both the target and parallax sources have significant V
c              uncertainty, the V-difference uncertainty will be larger
               sig_vel_dif = sqrt( sig_vel**2 + v_lsr_unc**2 )   ! km/s

c              Don't bother with this parallax source unless close in all parameters
               if ( drho.le.3.d0*sig_GMC     .and.
     +              dzee.le.3.d0*sig_GMC     .and.
     +              dvel.le.3.d0*sig_vel_dif       ) then

c                 Have a nearby parallax source...
                  n_pf = n_pf + 1

c                 Calculate weight
                  ssq = drho**2/(2.0d0*sig_GMC**2) +
     +                  dzee**2/(2.0d0*sig_GMC**2) +
     +                  dvel**2/(2.0d0*sig_vel_dif**2)
                  weight     = exp(-ssq)               ! unity for perfect match
                  sum_weight = sum_weight + weight

                  pref = ref_par(n_p)
                  eref = ref_punc(n_p)
                  vref = ref_vlsr(n_p)

                  write (lu_srcprt,1200) ref_src(n_p), vref,
     +                                pref, eref, ref_arm(n_p), weight
 1200             format(1x,a12,f7.1,5x,2f8.3,5x,a3,5x,f10.6)

c                 Add to prior parallax probability, converted to distance probability
                  do n_ps = 1, num_bins

c                    Calculate weighted parallax probability for this prior
                     p_bin = par_bins(n_ps)
                     call gauss_prob ( pref, eref, p_bin, prob )
                     p_prob = weight * prob

c                    Convert to distance probability (P(d)=P(pi)*pi^2)
c                    and sum in this distance bin...
                     d_prob = p_prob * p_bin**2
                     dist_prior(n_ps) = dist_prior(n_ps) + d_prob

                  enddo

               endif            ! prior source close to target source

            enddo            ! prior source loop

         endif          ! no parallaxes check

c        Do not want to assume entirely that the parallax sources give the distance
c        to the target.  So, add in a background "flat" pdf for all distances...
c        Use sum_weight, since we want total probability (not density)
c        Set maximum total parallax probability to P_max_PS
c        Normalize parallax association pdf
         call condition_pdf(num_bins, bin_size, P_max_PS, dist_prior)
         pdf_name = 'parallaxes_pdf'
         if ( lu_out .ge.7 ) then
            call output_pdf ( lu_out, pdf_name, src, num_bins,
     +                        dist_bins, dist_prior )
         endif

c        =============================================================
c        Now work towards a full Prob(d) by including
c        Prob_SA, Prob_KD, Prob_GL, and Prob_pm_ell...

         write (lu_srcprt,1900)
 1900    format(/'  Dist     P_PS    P_GalMod   P_posteriori',
     +           '   Arm',
     +          /'  (kpc)')

c        Run through distance bins and calculate various PDFs...
         do n_ps = 1, num_bins

c           Bin parallax value
            d_bin = dist_bins(n_ps)                               ! kpc

c           Calculate probabilities for spiral arm model, based on the
c           minimum separation of the target source (for a trial distance)
c           from the center of a given spiral arm
            call arm_model_prob ( Ro, num_arms, iarm_entries,
     +                            ell_dif_max,
     +                            arm_name, arm_probabilities,
     +                            arm_ell, arm_bee, arm_vel,
     +                            arm_Rgc, arm_beta, arm_dist,
     +                            width_min, width_Rref, width_slope,
     +                            ell, bee, d_bin,
     +                            arm_max, p_arm )

c           If user specified not to use Prob_SA or if the target's
c           Galactic longitude is out of range, reset p_arm to zero
            arm_max_bin(n_ps) = arm_max
            if ( P_max_SA.eq.0.d0 .or. .not.accept ) p_arm = 0.d0
            prob_arm(n_ps) = p_arm

c           Kinematic distance pdf...
c           Added measurement and Virial uncertainty in quadrature
            sig_vel_infl = sqrt( sig_vel**2 + v_lsr_unc**2 )      ! km/s
            call Dk_prob_density ( a1, a2, a3, Ro, To,
     +           Uo, Vo, Wo, Us, Vs, Ws,
     +           ell, bee, d_bin, v_lsr,
     +           Dk_near, Dk_far, far_prob, sig_vel_infl,
     +           p_Dk )
            if ( P_max_KD .eq. 0.d0 ) p_Dk = 0.d0
            prob_Dk(n_ps) = p_Dk

c           Galactic latitude pdf...
            call lat_prob_warped ( ell, bee, d_bin, Ro,
     +                             Rsigz, sigz, sigzdot,
     +                             n_bees, d_store, b_store,
     +                             p_lat, arg )
            if ( P_max_GL .eq. 0.d0 ) p_lat = 0.d0
            prob_lat(n_ps) = p_lat

c           Proper motion pdfs...
            prob_Dpm_ell(n_ps) = 0.d0
            prob_Dpm_bee(n_ps) = 0.d0

c           Caculate only if a non-zero PM component is entered
            if ( pm_x.ne.0.d0 .or. pm_y.ne.0.d0 ) then

c              Galactic longitude proper motion pdf...
               call pm_ell_prob_density ( a1, a2, a3, Ro, To,
     +                 Uo, Vo, Wo, Us, Vs, Ws,
     +                 ell, bee, d_bin, pm_ell, pm_ell_unc, sig_vel,
     +                 p_Dpm_ell )
               if ( P_max_PM .eq. 0.d0 ) p_Dpm_ell = 0.d0
               prob_Dpm_ell(n_ps) = p_Dpm_ell

c              Galactic latitude proper motion pdf...
               call pm_bee_prob_density ( Wo, ell, bee,
     +                                 pm_bee, pm_bee_unc,
     +                                 d_bin, sig_vel,
     +                                 p_Dpm_bee )
               if ( P_max_PM .eq. 0.d0 ) p_Dpm_bee = 0.d0
               prob_Dpm_bee(n_ps) = p_Dpm_bee

            endif

         enddo         ! distance bin loop

c        ------------------------------------------------------------
c        Normalize the component PDFs...

c        Arm model...
c        Smooth prob_arm to get rid of numerous secondary peaks caused by
c        using lbvRBD information, which has variations on small scales
         call smooth ( smooth_kpc, num_bins, dist_bins,
     +                 prob_arm )

c        Add a background probability density if Sum P(arm) < 1
         sum_arm_probs = 0.d0
         do i_a = 1, num_arms
            sum_arm_probs = sum_arm_probs + arm_probabilities(i_a)
         enddo

c        Final normalization, with background if needed
c        Use Sum{arm_probabilities) for P_max, but clip so that <P_max_SA
         P_max = min( P_max_SA, sum_arm_probs )
         call condition_pdf (num_bins, bin_size, P_max, prob_arm)

c        Latitude distance pdf...
         call condition_pdf (num_bins, bin_size, P_max_GL, prob_lat)

c        Combine spiral arm pdf and latitude pdf...
         do n_b = 1, num_bins
            prob_armlat(n_b) = prob_arm(n_b) * prob_lat(n_b)
         enddo
         call condition_pdf (num_bins, bin_size, P_max, prob_armlat)
         pdf_name = 'arm_latitude_pdf'
         if ( lu_out .ge.7 ) then
            call output_arm_pdf ( lu_out, pdf_name, src, num_bins,
     +                            dist_bins, prob_armlat, arm_max_bin )
         endif

c        Kinematic distance pdf...
         call condition_pdf (num_bins, bin_size, P_max_KD, prob_Dk)
         pdf_name = 'kinematic_distance_pdf'
         if ( lu_out .ge.7 ) then
            call output_pdf ( lu_out, pdf_name, src, num_bins,
     +                        dist_bins, prob_Dk )
         endif

c        Galactic longitude (ell) proper motion distance pdf...
         call condition_pdf (num_bins, bin_size, P_max_PM,
     +                       prob_Dpm_ell)
         pdf_name = 'pm_ell_distance_pdf'
         if ( lu_out .ge.7 ) then
            call output_pdf ( lu_out, pdf_name, src, num_bins,
     +                        dist_bins, prob_Dpm_ell )
         endif

c        Galactic latitude (bee) proper motion distance pdf...
         call condition_pdf (num_bins, bin_size, P_max_PM,
     +                       prob_Dpm_bee)
         pdf_name = 'pm_bee_distance_pdf'
         if ( lu_out .ge.7 ) then
            call output_pdf ( lu_out, pdf_name, src, num_bins,
     +                        dist_bins, prob_Dpm_bee )
         endif

c        Combine probabilities for arm+latitude, Dk, and PMs ...
c        to make a "Milky Way" model based pdf
         do n_ps = 1, num_bins
            prob_MWmodel(n_ps) = 0.d0
            if ( prob_armlat(n_ps)  .gt.1.d-09 .and.
     +           prob_Dk(n_ps)      .gt.1.d-09 .and.
     +           prob_Dpm_ell(n_ps) .gt.1.d-09 .and.
     +           prob_Dpm_bee(n_ps) .gt.1.d-09       ) then

               prob_MWmodel(n_ps) = prob_armlat(n_ps)*prob_Dk(n_ps)*
     +                         prob_Dpm_ell(n_ps)*prob_Dpm_bee(n_ps)
            endif
         enddo

c        Normalize this pdf
         P_max = 1.d0
         call condition_pdf(num_bins, bin_size, P_max, prob_MWmodel)

c        Combine the Milky Way and parallax source pdfs...
         do n_ps = 1, num_bins
            prob_dist(n_ps) = 0.d0
            if ( dist_prior(n_ps)  .gt.1.d-09 .and.
     +           prob_MWmodel(n_ps).gt.1.d-09       ) then
               p_dist          = dist_prior(n_ps)*prob_MWmodel(n_ps)
               prob_dist(n_ps) = p_dist      ! combined PDF for parallax bin
            endif
         enddo

c        Normalize
         P_max = 1.d0
         call condition_pdf (num_bins, bin_size, P_max, prob_dist)
         pdf_name = 'final_distance_pdf'
         if ( lu_out .ge.7 ) then
            call output_pdf ( lu_out, pdf_name, src, num_bins,
     +                        dist_bins, prob_dist )
         endif

c        Finished generating PDFs!
c        ---------------------------------------------------------------------
c        Now start to analyze the probability information...
c        First, print out PDF information
         do n_ps = 1, num_bins
            if ( prob_dist(n_ps) .gt. 0.5d-05 ) then
               write (lu_srcprt,2000) dist_bins(n_ps),
     +                   dist_prior(n_ps), prob_MWmodel(n_ps),
     +                   prob_dist(n_ps), arm_max_bin(n_ps)
 2000          format(f7.3,3f10.5,8x,a3)
            endif
         enddo

c        Find combined pdf peak(s)...to be used for initial values for fitting
c        individual distance components
         call find_probability_peaks(lu_srcprt, max_num_peaks,
     +                bin_size, num_bins, dist_bins, prob_dist,
     +                num_peaks, peak_prob,
     +                peaks, peaks_width)

         call edit_peaks ( num_peaks,
     +                     peak_prob, peaks, peaks_width,
     +                     use_peak )

         write (lu_srcprt,3000)
 3000    format(/'  Peak  Distance     +/-    Use?')
         if ( num_peaks .gt. 0 )
     +      write (lu_srcprt,3100) (n, peaks(n), peaks_width(n),
     +                              use_peak(n), n=1,num_peaks)
 3100       format(i5,2f10.2,5x,l1)

c        Initialize parameter info...
         do i = 1, max_num_params
            params(i)   = 0.d0
            paramids(i) = 0
         enddo

c        Include a flat "background" (baseline) when fitting
         params(1)   = 0.0d0                  ! baseline offset
         paramids(1) = 1
         params(2)   = 0.0d0                  ! baseline slope
         paramids(2) = 0

c        Fit multiple Gaussians to probability vs distance data,
c        provided there is at least one recognizable peak in the PDF
         if ( num_peaks .ge. 1 ) then

c           First, put initial guesses into parameters array
            index  = 2
            do j_p = 1, num_peaks

               index = index + 1
               params(index) = peak_prob(j_p)   ! amplitude (1/kpc)
               paramids(index) = 0
               if ( use_peak(j_p) ) paramids(index) = 1

               index = index + 1
               params(index) = peaks(j_p)       ! center (kpc)
               paramids(index) = 0
               if ( use_peak(j_p) ) paramids(index) = 1

               index = index + 1
c              convert from 1-sigma to FWHM...
               params(index) = 2.3548d0*peaks_width(j_p) ! FWHM (kpc)
               paramids(index) = 0
               if ( use_peak(j_p) ) paramids(index) = 1

            enddo
            num_params = index

            call encode_parnames ( num_peaks, parnames )

c           But, for sources with very odd (l,b,v) values, don't fit Gaussians;
c           instead use the initial guesses
            odd_source = .false.
            if ( Dk_near.gt.25.d0 ) then
c              Avoid fitting by setting via "solve-for" flags to zeros
               do n_p = 1, num_params
                  paramids(n_p) = 0
               enddo
               odd_source = .true.
            endif

            call fit_multiple_gaussians ( lu_srcprt,
     +            num_bins, num_params, parnames,
     +            params, paramids,
     +            dist_bins, prob_dist, arm_max_bin,
     +            peak_dist, peak_dunc, peak_int )

          else

cc            write (lu_print,3120) src
 3120       format(' No distance probability peaks found for ',a)

         endif                  ! num_peaks > 0 check

c        ===============================================================
c        Output to fort.3x depending on arm_assignment
c        Write separate files by arms for plotting

c        Find two greatest integrated probability peaks
         p_int_max = 0.d0
         n_max     = 0
         if ( num_peaks .ge. 1 ) then
c           Find greatest integrated probability peak...
            do n_p = 1, num_peaks
               if ( peak_int(n_p) .gt. p_int_max ) then
                  p_int_max = peak_int(n_p)
                  n_max     = n_p
               endif
            enddo
c           Find second greatest peak...
            p2_int_max = 0.d0
            n2_max     = 0
            if ( num_peaks .ge. 2 ) then
               do n_p = 1, num_peaks
                  if ( peak_int(n_p).gt.p2_int_max .and.
     +                 n_p.ne.n_max                     ) then
                     p2_int_max = peak_int(n_p)
                     n2_max     = n_p
                  endif
               enddo
            endif
         endif

c        If there is a maximum probability peak, output
         if ( n_max .gt. 0 ) then

            distance = peak_dist(n_max)
            call closest_arm ( num_bins, dist_bins, arm_max_bin,
     +                         distance,
     +                         arm_indicated )

            lu_arm = 30                                ! unassigned arm
            do i_a = 1, num_arms
               if (arm_indicated .eq. arm_name(i_a)) lu_arm=30+i_a
            enddo

c           Flag if we didn't fit Gaussians...
            questionable = ' '
            if ( sum_arm_probs .lt. 0.1d0 ) questionable = '? '
            if ( odd_source )               questionable = '??'

c           Write out results in several files...
            write (lu_arm,3200) ell, bee, v_lsr, v_lsr_unc,
     +          peak_dist(n_max),peak_dunc(n_max),peak_int(n_max),
     +          arm_indicated,questionable
 3200       format(2f7.2,f7.1,f5.1,f8.2,f7.2,f8.2,5x,a3,a2)

c           Get 2nd peak information...
            p2_dist    = 0.d0
            p2_dunc    = 0.d0
            p2_arm     = '...'
            q2         = '  '
            if ( n2_max .gt. 0 ) then
               p2_dist = peak_dist(n2_max)
               p2_dunc = peak_dunc(n2_max)
               d2 = peak_dist(n2_max)
               call closest_arm ( num_bins, dist_bins, arm_max_bin,
     +                            d2,
     +                            p2_arm )
               q2 = questionable
            endif

cc            write (lu_print,3210) ell, bee, v_lsr, v_lsr_unc,
cc     +          peak_dist(n_max),peak_dunc(n_max),peak_int(n_max),
cc     +          arm_indicated,questionable,
cc     +          p2_dist, p2_dunc, p2_int_max, p2_arm,q2
 3210       format(2f7.2,f7.1,f5.1,f8.2,f7.2,f8.2,5x,a3,a2,
     +                             f8.2,f7.2,f8.2,5x,a3,a2)

            write (lu_summary,3210) ell, bee, v_lsr, v_lsr_unc,
     +          peak_dist(n_max),peak_dunc(n_max),peak_int(n_max),
     +          arm_indicated,questionable,
     +          p2_dist, p2_dunc, p2_int_max, p2_arm,q2

         endif

         close ( unit=lu_srcprt )

       endif          ! "commented-out target" check

      enddo         ! target loop

      end

c====================================================================
      subroutine get_controls ( lu_control, lu_print,
     +               P_max_SA, P_max_KD, P_max_GL, P_max_PS,
     +               P_max_PM )

c     Read in maximum probabilities, used to add a background to the PDF

      implicit real*8 (a-h,o-z)

      character*48 control_file
      character*1  c1

      logical      bad_P_max

      control_file = 'probability_controls.inp'
cc      write (lu_print,1000) control_file
 1000 format(/' Reading in file "',a,'"')

      open ( unit=lu_control, file=control_file )

c     Flagged values (in case of an empty file)...
      P_max_SA     = -1.d0
      P_max_KD     = -1.d0
      P_max_GL     = -1.d0
      P_max_PS     = -1.d0
      P_max_PM     = -1.d0

      ieof    = 0
      do while ( ieof.ge.0 )

         read(lu_control,*,iostat=ieof) c1

c        Check for end of file, commented-line, or blank line...
         if ( ieof.ge.0 .and. c1.ne.'!'   ) then

            backspace ( unit=lu_control )
            read (lu_control,*) P_max_SA, P_max_KD, P_max_GL, P_max_PS,
     +                          P_max_PM

         endif

      enddo

      close(unit=lu_control)

cc      write (lu_print,2000) P_max_SA, P_max_KD, P_max_GL, P_max_PS,
cc     +                      P_max_PM
 2000 format( ' Max Prob_SA =',f5.2,
     +       /' Max Prob_KD =',f5.2,
     +       /' Max Prob_GL =',f5.2,
     +       /' Max Prob_PS =',f5.2,
     +       /' Max Prob_PM =',f5.2,/ )

c     Check for reasonable P_max values...
      bad_P_max = .false.
      if ( P_max_SA.lt.0.d0 .or. P_max_SA.gt.1.d0 ) bad_P_max = .true.
      if ( P_max_KD.lt.0.d0 .or. P_max_KD.gt.1.d0 ) bad_P_max = .true.
      if ( P_max_GL.lt.0.d0 .or. P_max_GL.gt.1.d0 ) bad_P_max = .true.
      if ( P_max_PS.lt.0.d0 .or. P_max_PS.gt.1.d0 ) bad_P_max = .true.
      if ( P_max_PM.lt.0.d0 .or. P_max_PM.gt.1.d0 ) bad_P_max = .true.

      if ( bad_P_max ) then
         stop ': P_max value(s) in "probability_controls.inp"'
      endif

      return
      end

c====================================================================
      subroutine longitude_range ( ell, glong_min, glong_max,
     +                             accept )

c     Assumes   0 < ell < 360 deg

c     Decides if source longitude (ell) is within range from
c     glong_min to glong_max.   There is a complication with
c     for glong_min (eg, if 350deg, need to treat as if -10 deg).

c     If in range, logical flag "accept" set to .true.

      implicit real*8 (a-h,o-z)

      logical  accept

      accept = .false.

      if ( (glong_min.ge.0.d0 .and. ell.ge.glong_min) .and.
     +      ell.le.glong_max                          ) accept=.true.

c     More complicated if glong_min < 0 deg
      if ( glong_min .lt. 0.d0 ) then

         if ( ell.ge.glong_min .and. ell.le.glong_max ) accept=.true.

         epm = ell - 360.d0
         if ( epm.ge.glong_min .and. epm.le.glong_max ) accept=.true.

      endif

      return
      end

c======================================================================

      subroutine output_pdf ( lu_plots, pdf_name, src, num_bins,
     +                        dist_bins, prob_arm )

      implicit real*8 (a-h,o-z)

      real*8        dist_bins(1001), prob_arm(1001)

      character*48  file_name
      character*14  src, stripped_src
      character*32  pdf_name, stripped_pdf

c     Strip blanks out of name in output file naming
      call strip_blanks_14 ( src,
     +                       stripped_src, nch_src )
      call strip_blanks_32 ( pdf_name,
     +                       stripped_pdf, nch_pdf )

      write (file_name,1000) stripped_src(1:nch_src),
     +                       stripped_pdf(1:nch_pdf)
 1000 format(a,'_',a,'.dat')

      open ( unit=lu_plots, file=file_name )

      write (lu_plots,2000) src, pdf_name
 2000 format('! Source: ',a14,'; ',a32,/
     +       '!   D(kpc) Probability')

      do n = 1, num_bins

         if ( prob_arm(n) .gt. 1.d-09 ) then

            write (lu_plots,3000) dist_bins(n), prob_arm(n)
 3000       format(f10.3,f10.6)

         endif

      enddo

      close ( unit=lu_plots )

      return
      end
c======================================================================

      subroutine output_arm_pdf( lu_plots, pdf_name, src, num_bins,
     +                           dist_bins, prob_arm, arm_max_bin )

c     Special output for spiral arm PDF
c     Also writes extra file with ranges of distances for arms

      implicit real*8 (a-h,o-z)

      real*8        dist_bins(1001), prob_arm(1001)

      character*48  file_name
      character*14  src, stripped_src
      character*32  pdf_name, stripped_pdf
      character*12  arm_max_bin(1001), arm_old

      real*8        range_start(50), range_stop(50)
      character*3   arm_in_range(50)

c     Strip blanks out of source name in output file naming
      call strip_blanks_14 ( src,
     +                       stripped_src, nch_src )
      call strip_blanks_32 ( pdf_name,
     +                       stripped_pdf, nch_pdf )

      write (file_name,1000) stripped_src(1:nch_src),
     +                       stripped_pdf(1:nch_pdf)
 1000 format(a,'_',a,'.dat')

      open ( unit=lu_plots, file=file_name )

      write (lu_plots,2000) src, pdf_name
 2000 format('! Source: ',a14,'; ',a32,/
     +       '!   D(kpc) Probability')

      arm_old = '   '
      n_range = 0
      dist_end = 20.d0
      do n = 1, num_bins

c        Cut-off arm regions for plotting at Prob=0.0001
         if ( prob_arm(n) .ge. 0.0001d0 ) then

            write(lu_plots,3000)dist_bins(n),prob_arm(n),arm_max_bin(n)
 3000       format(f10.3,f10.6,5x,a12)

c           Keep track of ranges for arms...
            if ( arm_max_bin(n).ne.arm_old .and.
     +           arm_max_bin(n).ne.'...'        ) then

c              Start new arm range...
               n_range = n_range + 1
               range_start(n_range) = dist_bins(n)
               arm_in_range(n_range)= arm_max_bin(n)
               arm_old              = arm_max_bin(n)

c              Finish previous arm range
               if ( n_range .gt. 1 ) then
                  n_previous = n_range - 1
                  range_stop(n_previous) = dist_end
               endif

            endif

            if ( arm_max_bin(n) .ne. '...' ) dist_end = dist_bins(n)

         endif

      enddo

c     Get last end of range...
      range_stop(n_range) = dist_end

      close ( unit=lu_plots )

c     -------------------------------------------------------
c     Now output new file with arm range information...
      write (file_name,4000) stripped_src(1:nch_src)
 4000 format(a,'_arm_ranges.dat')

      open ( unit=lu_plots, file=file_name )

      write (lu_plots,5000) src
 5000 format('! Arm distance ranges for ',a13,/
     +       '!  Start(kpc)  End(kpc)  Arm',/
     +       '  -10.0      -9.0      Dum')

      do n_r = 1, n_range

         if ( arm_in_range(n_r) .ne. '...') then

            write(lu_plots,6000) range_start(n_r),range_stop(n_r),
     +                           arm_in_range(n_r)
 6000       format(2f10.2,5x,a3)

          endif

      enddo

      close ( unit=lu_plots )

      return
      end

c======================================================================

      subroutine strip_blanks_14 ( in_name,
     +                             out_name, nchars )

      character*14  full_name, clean_name, in_name, out_name

      character*1   cin(14), cout(14)

      equivalence   (full_name,cin)
      equivalence   (clean_name,cout)

      full_name = in_name

      j = 0
      do i = 1, 14
         if ( cin(i) .ne. ' ') then
            j=j+1
            cout(j) = cin(i)
         endif
      enddo

      nchars   = j
      out_name = clean_name

      return
      end

c======================================================================

      subroutine strip_blanks_32 ( in_name,
     +                             out_name, nchars )

      character*32  full_name, clean_name, in_name, out_name

      character*1   cin(32), cout(32)

      equivalence   (full_name,cin)
      equivalence   (clean_name,cout)

      full_name = in_name

      j = 0
      do i = 1, 32
         if ( cin(i) .ne. ' ') then
            j=j+1
            cout(j) = cin(i)
         endif
      enddo

      nchars   = j
      out_name = clean_name

      return
      end

c======================================================================

      subroutine assign_arm_probabilities ( ell, bee, v_lsr, v_lsr_unc,
     +           ell_dif_max,
     +           num_arms, arm_ell, arm_bee, arm_vel, iarm_entries,
     +           arm_Rgc, arm_beta, arm_dist,
     +           width_min, width_Rref, width_slope,
     +           sigz, Rsigz, sigzdot,
     +           seg_ell_min, seg_ell_max,
     +           sig_ell, sig_vel,
     +           arm_probabilities )

c     Interpolates between tables of arm (l,b,v) values to find the
c     closest to the target values in order to assign a probability
c     that the target is in the arm: P(arm|l,b,v,ArmModels).

      implicit real*8 (a-h,o-z)

      real*8   arm_ell(29,300), arm_bee(29,300), arm_vel(29,300)
      real*8   arm_Rgc(29,300), arm_beta(29,300), arm_dist(29,300)
      real*8   seg_ell_min(29,2), seg_ell_max(29,2)

      integer  iarm_entries(29)

      real*8   arm_probabilities(29)

      pi     = 4.d0*atan(1.d0)
      twopi  = 2.d0*pi
      raddeg = 180.d0/pi
      degrad = pi/180.d0

c     Go through all arms...
      sum_prob_max = 0.d0
      do n_a = 1, num_arms

c        Go through table entries for this arm...
         prob_max = 0.d0
         i1 = iarm_entries(n_a) - 1
         do n_e = 1, i1

c           Adjacent table entries to be used for interpolation
            e0 = arm_ell(n_a,n_e)
            e1 = arm_ell(n_a,n_e+1)

c           Difference between longitude table entries
            diff_e  = (e1 - e0)

c           Resolve negative, positive longitude issues
            if ( diff_e.gt.180.d0 .and. ell.gt.180.d0 ) then
c              eg, e0=10, e1=350, ell=355
               e0 = e0 + 360.d0
            endif
            if ( diff_e.gt.180.d0 .and. ell.le.180.d0 ) then
c              eg, e0=10, e1=350, ell=5
               e1 = e1 - 360.d0
            endif
            if ( diff_e.lt.-180.d0 .and. ell.gt.180.d0 ) then
c              eg, e0=350, e1=10, ell=355
               e1 = e1 + 360.d0
            endif
            if ( diff_e.lt.-180.d0 .and. ell.le.180.d0 ) then
c              eg, e0=350, e1=10, ell=5
               e0 = e0 - 360.d0
            endif
c           Re-calculate difference
            diff_e  = (e1 - e0)
c           Check for duplicate entries to avoid a divide by zero
            if ( abs(diff_e) .gt. 0.001d0 ) then

c             Check that ell is near the [e0,e1] range,
c             now have arm information typically every 1 deg,
c             so limit to +/- "ell_dif_max" degrees
              if ( abs( ell-((e0+e1)/2.d0) ) .lt. ell_dif_max ) then

                denom = -99.9d0

                b0 = arm_bee(n_a,n_e)
                b1 = arm_bee(n_a,n_e+1)
                v0 = arm_vel(n_a,n_e)
                v1 = arm_vel(n_a,n_e+1)

                R0 = arm_Rgc(n_a,n_e)
                R1 = arm_Rgc(n_a,n_e+1)
                d0 = arm_dist(n_a,n_e)
                d1 = arm_dist(n_a,n_e+1)

c               Slopes of parameters vs e ...
                slope_b = (b1 - b0) / diff_e          ! latitude
                slope_v = (v1 - v0) / diff_e          ! velocity
                slope_R = (R1 - R0) / diff_e          ! Galoctocentric Radius
                slope_d = (d1 - d0) / diff_e          ! Heliocentric distance

c               Interpolate in steps of 0.1 times difference of longitudes
                n_i = abs(10.d0*diff_e) + 1
                if ( n_i .lt. 2 ) n_i = 2

                d_e = diff_e / (n_i - 1)
                do n = 1, n_i

                  e_step = (n - 1) * d_e
                  arm_ell_interp = e0 + e_step
                  arm_bee_interp = b0 + slope_b * e_step
                  arm_vel_interp = v0 + slope_v * e_step
                  arm_Rgc_interp = R0 + slope_R * e_step
                  arm_dist_interp= d0 + slope_d * e_step

c                 Only calc probability if within range of knowledge for arm
c                 (if arm not known for this ell, R and d should be zero)
                  if ( R0.gt.0.d0 .and. R1.gt.0.d0 .and.
     +                 d0.gt.0.d0 .and. d1.gt.0.d0       ) then

                     dell = abs( ell -   arm_ell_interp )
                     dbee = abs( bee -   arm_bee_interp )
                     dvel = abs( v_lsr - arm_vel_interp )

c                    Calculate sigma-z based on expected z-thickness of arm
c                    at the arm's distance from the Sun
                     if ( arm_Rgc_interp .lt. Rsigz ) then
                           sigma_z = sigz                             ! kpc
                        else
                           sigma_z = sigz +
     +                               sigzdot*(arm_Rgc_interp-Rsigz)   ! kpc
                     endif

c                    Convert this to sigma-b
                     sig_bee = ( sigma_z/arm_dist_interp ) * raddeg   ! deg

c                    Need similar calculation for sig_ell...
                     if ( arm_Rgc_interp .gt. width_Rref ) then
                           width = width_min +
     +                        (arm_Rgc_interp-width_Rref)*width_slope
                        else
                           width = width_min                        ! kpc
                     endif
c                    Convert width (kpc) to sigma_ell (deg)
                     sig_ell = (width/arm_dist_interp) * raddeg     ! deg

c                    Account for misalignment of arm with ray along source longitude
c                    by increasing "effective" sigma_ell by sec(misalignment).
c                    Estimate angle between ray and arm tangent...
                     t_d = 10.d0 ! arbitrary distance from Sun (0,0) to define slope of ray
                     call find_angle ( ell, t_d, e0, d0, e1, d1,
     +                                 angle_between )
                     cos_angle = abs( cos(angle_between) )
                     if ( cos_angle .gt. 0.1d0 ) then
                           sig_ell = sig_ell / cos_angle
                        else
c                          Allowing for curvature, cutoff at factor of 10...
                           sig_ell = sig_ell / 0.1d0
                     endif

c                    Increase sig_vel for uncertainty in the arm velocity,
c                    which jitters when following CO clouds in (l-v) plots.
                     sig_vel_inc = sqrt(2.d0) * sig_vel                  ! km/s
c                    And Vlsr (measurement) uncertainty in quadrature.
                     sig_vel_dif = sqrt( sig_vel_inc**2 + v_lsr_unc**2 )

c                    Calculate probability if position and velocity reasonably close
                     if ( dell.le.3.d0*sig_ell  .and.
     +                    dbee.le.3.d0*sig_bee  .and.
     +                    dvel.le.3.d0*sig_vel_dif     ) then

c                       Source possibly within this arm for this
c                       interpolated (l,b,v), calculate a probability.
                        pos_ssq = dell**2/(2.d0*sig_ell**2) +
     +                            dbee**2/(2.d0*sig_bee**2) +
     +                            dvel**2/(2.d0*sig_vel_dif**2)
                        prob    = exp(-pos_ssq)

                        if ( prob.gt.prob_max ) prob_max = prob

                     endif

                  endif         ! have arm info

                enddo        ! interpolation grid

              endif      ! ell close to [e0,e1] check

            endif   ! duplicate entries check (diff_e=0)

         enddo   ! table entries for one arm

         arm_probabilities(n_a) = prob_max
         sum_prob_max = sum_prob_max + prob_max

      enddo   ! different arms

c     If no clear arm segment association, that's OK...
c     will add a background to the final P_SA(d) later

      return
      end

c=====================================================================
      subroutine find_angle ( t_l,t_d, a_l_1,a_d_1, a_l_2,a_d_2,
     +                        angle_rad )

c     Given a "trial" longitude and distance (t_l,t_d) for a target from
c     the Sun at (0,0), and two points (s_l,s_d) for a spiral Arm segment
c     (assumed linear over small distances), find the minimum angle between
c     the two lines.

      implicit real*8 (a-h,o-z)

      pi      = 4.d0 * atan(1.d0)
      deg_rad = pi/180.d0

      t_dx = t_d*sin(t_l * deg_rad)
      t_dy = t_d*cos(t_l * deg_rad)

      a_dx = a_d_2*sin(a_l_2*deg_rad) - a_d_1*sin(a_l_1*deg_rad)
      a_dy = a_d_2*cos(a_l_2*deg_rad) - a_d_1*cos(a_l_1*deg_rad)

      if ( abs(t_dx) .lt. 1.d-9 ) then
         temp = t_dx
         if ( temp .lt. 0.d0 ) t_dx = -1.d-9
         if ( temp .ge. 0.d0 ) t_dx = +1.d-9
      endif

      if ( abs(a_dx) .lt. 1.d-9 ) then
         temp = a_dx
         if ( temp .lt. 0.d0 ) a_dx = -1.d-9
         if ( temp .ge. 0.d0 ) a_dx = +1.d-9
      endif

      angle_1 = atan(t_dy/t_dx) - atan(a_dy/a_dx)
      angle_1 = abs( angle_1 )                          ! rad

c     Now rotate system by 90 degrees and then take the
c     minimum of the absolute values of angles
      if ( abs(t_dy) .lt. 1.d-9 ) then
         temp = t_dy
         if ( temp .lt. 0.d0 ) t_dy = -1.d-9
         if ( temp .ge. 0.d0 ) t_dy = +1.d-9
      endif

      if ( abs(a_dy) .lt. 1.d-9 ) then
         temp = a_dy
         if ( temp .lt. 0.d0 ) a_dy = -1.d-9
         if ( temp .ge. 0.d0 ) a_dy = +1.d-9
      endif

      angle_2 = atan(t_dx/t_dy) - atan(a_dx/a_dy)
      angle_2 = abs( angle_2 )                          ! rad

      angle_rad = min( angle_1, angle_2 )               ! rad

      return
      end

c======================================================================

      subroutine closest_arm ( num_bins, dist_bins, arm_max_bin,
     +                         distance,
     +                         arm_indicated )

      implicit real*8 (a-h,o-z)

      real*8        dist_bins(1001)

      character*12  arm_indicated, arm_max_bin(1001)

      arm_indicated = '...'

      delta_min = 9.d9                        ! kpc
      do i = 1, num_bins

         delta_dist = abs( distance - dist_bins(i) )

         if ( delta_dist .lt. delta_min ) then

            delta_min = delta_dist

            arm_indicated = arm_max_bin(i)

         endif

      enddo

      return
      end

c======================================================================
      subroutine get_warping ( num_bins, Ro, num_arms,
     +                         iarm_entries,
     +                         arm_name, arm_probabilities,
     +                         arm_ell, arm_bee, arm_vel,
     +                         arm_Rgc, arm_beta, arm_dist,
     +                         ell, bee, dist_bins,
     +                         n_bees, d_store, b_store, a_store )


c     Calculates latitude of warping from spiral arm models, based on
c     the minimum separation of the target source (for a trial distance)
c     from the center of a given spiral arm

c     V2: Noted that this may fail if an arm wraps more than 360 deg
c     in azimuth, since two pairs of distance-latitude would be needed.
c     Could be fixed by picking largest distance point, since inner
c     Galaxy is not warped.

      implicit real*8 (a-h,o-z)

      real*8        dist_bins(1001), d_min(1001), b_min(1001)
      integer       i_a_min(1001)

      real*8        arm_probabilities(29)
      real*8        d_store(29), b_store(29)
      real*8        arm_d_min(29), arm_b_min(29)

      real*8        arm_ell(29,300), arm_bee(29,300), arm_vel(29,300)
      real*8        arm_Rgc(29,300), arm_beta(29,300), arm_dist(29,300)

      integer       iarm_entries(29)

      character*12  arm_name(29), arm_max, atemp
      character*12  a_store(29), arm_name_b_min


      do n_ps = 1, num_bins

         do n_a = 1, num_arms
            arm_d_min(n_a)=9.9d0
            arm_b_min(n_a)=0.d0
         enddo

         d_bin = dist_bins(n_ps)

         call arm_model_b ( Ro, num_arms, iarm_entries,
     +                      arm_name, arm_probabilities,
     +                      arm_ell, arm_bee, arm_vel,
     +                      arm_Rgc, arm_beta, arm_dist,
     +                      ell, bee, d_bin,
     +                      arm_d_min, arm_b_min )

c        Check each arm for minimum distance from bin point
         d_min(n_ps)   = 99.d9
         b_min(n_ps)   = 0.d0
         i_a_min(n_ps) = 0

         do i_a = 1, num_arms

            if ( arm_d_min(i_a) .lt. d_min(n_ps) ) then
c              Store values for closest arm
               d_min(n_ps)   = arm_d_min(i_a)
               b_min(n_ps)   = arm_b_min(i_a)
               i_a_min(n_ps) = i_a
            endif

         enddo            ! arms

      enddo           ! dist bin

c     Pick only 1 point from each arm among bins
      n_d = 0
      do i_a = 1, num_arms

         d_min_min = 9.d9
         b_min_min = 0.d0
         arm_name_b_min = 'No_Name'

         do n_p = 1, num_bins

            if ( i_a_min(n_p).eq.i_a .and.
     +           d_min(n_p)  .lt.d_min_min ) then

               arm_name_b_min = arm_name(i_a)
               d_min_min = d_min(n_p)
               b_min_min = b_min(n_p)
               d_bin_min = dist_bins(n_p)

            endif

         enddo

c        Accept this arm's latitude only if it "passes close by"
         if ( d_min_min .lt. 0.25d0 ) then

            n_d = n_d + 1
            d_store(n_d) = d_bin_min
            b_store(n_d) = b_min_min
            a_store(n_d) = arm_name_b_min

         endif

      enddo

c     Sort b-values by distance from Sun
      n_bees = n_d
      do n_d = 1, n_bees

         do m_d = n_d+1, n_bees

            if ( d_store(m_d) .lt. d_store(n_d ) ) then
               temp = d_store(n_d)
               d_store(n_d) = d_store(m_d)
               d_store(m_d) = temp
               temp = b_store(n_d)
               b_store(n_d) = b_store(m_d)
               b_store(m_d) = temp
               atemp = a_store(n_d)
               a_store(n_d) = a_store(m_d)
               a_store(m_d) = atemp
            endif

         enddo

      enddo

      return
      end


c======================================================================

      subroutine arm_model_b ( Ro, num_arms, iarm_entries,
     +                            arm_name, arm_probabilities,
     +                            arm_ell, arm_bee, arm_vel,
     +                            arm_Rgc, arm_beta, arm_dist,
     +                            ell, bee, dist,
     +                            arm_d_min, arm_b_min )

c     For each arm and given a trial source distance, calculate
c     the minimum separation point and the latitude at that point

      implicit real*8 (a-h,o-z)

      real*8   arm_probabilities(29)
      real*8   arm_b_min(29), arm_d_min(29)

      real*8   arm_ell(29,300), arm_bee(29,300), arm_vel(29,300)
      real*8   arm_Rgc(29,300), arm_beta(29,300), arm_dist(29,300)

      integer  iarm_entries(29)

      character*12  arm_name(29)

      do i_a = 1, num_arms

c           Calculate min distance to arm and latitude at that point
            i_arm = i_a
            call calc_arm_b( Ro, ell, bee, dist,
     +                       i_arm, iarm_entries,
     +                       arm_ell, arm_bee, arm_vel,
     +                       arm_Rgc, arm_beta, arm_dist,
     +                       Rsrc,
     +                       del_min, b_min )

c           Save b_min values
            arm_d_min(i_a) = del_min
            arm_b_min(i_a) = b_min

      enddo    ! all spiral arms

      return
      end

c===================================================================

      subroutine calc_arm_b( Ro, ell, bee, dist,
     +                       i_arm, iarm_entries,
     +                       arm_ell, arm_bee, arm_vel,
     +                       arm_Rgc, arm_beta, arm_dist,
     +                       Rsrc,
     +                       del_min_in, b_min )

c     Calculates minimum distance to an arm in the plane (del_min_in) and
c     the latitude (b_min) at its distance for arm number "i_arm"

      implicit real*8 ( a-h, o-z )

      real*8  arm_ell(29,300), arm_bee(29,300), arm_vel(29,300)
      real*8  arm_Rgc(29,300), arm_beta(29,300), arm_dist(29,300)

      integer iarm_entries(29)

      logical near_arm

      pi = 4.d0*atan(1.d0)
      deg_to_rad = pi/180.d0

c     Source values...
      ell_rad  = ell * deg_to_rad                      ! radians
      bee_rad  = bee * deg_to_rad                      ! radians
      x_galcen = dist * sin(ell_rad)                   ! kpc
      y_galcen = Ro - dist * cos(ell_rad)              ! kpc
      z_galcen = dist * sin(bee_rad)                   ! kpc
      Rsrc     = sqrt( x_galcen**2 + y_galcen**2 )     ! kpc
      beta_src = atan2( x_galcen, y_galcen )           ! radians
      beta_src_deg = beta_src / deg_to_rad             ! deg

c     Go through table of (l,b,V,R,beta,d) for arm "i_arm" to find the closest
c     azimuth (arm_beta) to the source_beta value.   Use azimuth (beta), not
c     longitude (ell), since it is much better behaved.

c     But, first check that ell is within the overall range for the arm segment.
      e_min = 999.d0
      e_max =-999.d0
      i_a_max = iarm_entries(i_arm)
      do i_a = 1, i_a_max
         e_val = arm_ell(i_arm,i_a)              ! deg
         if ( e_val .gt. e_max ) e_max = e_val
         if ( e_val .lt. e_min ) e_min = e_val
      enddo
      del_min_in = 999.d0
      b_min      = 0.d0
      if ( ell.lt.e_min .or. ell.gt.e_max ) return

c     If we get here, then target ell is within arm longitude range.
      del_beta_min = 999.d0
      i_a_min      = 0
      do i_a = 1, i_a_max

c        Check that arm entry has estimated values...
         if ( arm_Rgc(i_arm,i_a)  .gt. 0.d0 .and.
     +        arm_dist(i_arm,i_a) .gt. 0.d0      ) then

            bs = beta_src_deg                           ! deg
            ba = arm_beta(i_arm,i_a)                    ! deg

c           To allow beta values beyond -180 to +180, need 3 checks
            del_beta1 = abs( bs - ba )                  ! deg
            del_beta2 = abs( bs - ba - 360.d0 )         ! deg
            del_beta3 = abs( bs - ba + 360.d0 )         ! deg
            if ( del_beta1 .lt. del_beta_min ) then
               del_beta_min = del_beta1
               i_a_min = i_a
            endif
            if ( del_beta2 .lt. del_beta_min ) then
               del_beta_min = del_beta2
               i_a_min = i_a
            endif
            if ( del_beta3 .lt. del_beta_min ) then
               del_beta_min = del_beta3
               i_a_min = i_a
            endif

         endif

      enddo

c     Then get arm longitude of point with closest beta
      ell_arm_closest = arm_ell(i_arm,i_a_min)            ! deg
      del_ell_1       = abs( ell - ell_arm_closest )
c     Allow for possible 360 deg differences...
      del_ell_2       = abs( del_ell_1 - 360.d0 )
      del_ell_3       = abs( del_ell_1 + 360.d0 )

      del_min_in = 999.d0                                    ! kpc
      Rgc_min    = 999.d9
      b_min      = 0.d0

c     Longitudes must be close (<3 deg) to be considered for this arm
      near_arm = .false.
      if ( del_ell_1 .lt. 3.d0 ) then
         near_arm = .true.
         ambiguity = 0.d0
      endif
      if ( del_ell_2 .lt. 3.d0 ) then
         near_arm = .true.
         ambiguity = -360.d0
      endif
      if ( del_ell_3 .lt. 3.d0 ) then
         near_arm = .true.
         ambiguity = +360.d0
      endif

      if ( near_arm ) then

         i_a_0 = i_a_min - 3
         i_a_1 = i_a_min + 3
         if ( i_a_0 .lt. 1 ) i_a_0 = 1
         if ( i_a_1 .gt. i_a_max ) i_a_1 = i_a_max

c        Now go through in steps of 0.1 deg (actually 0.1 of
c        difference between longitude entries) and calculate
c        the minimum 2-D separation of source from arm center
         do i_a = i_a_0, i_a_1-1

          diff_e    = arm_ell(i_arm,i_a) - arm_ell(i_arm,i_a+1)
          diff_e    = diff_e + ambiguity

c         Check for duplicate entries (diff_e=0)
          if ( abs(diff_e) .gt. 0.001d0 ) then

            diff_b    = arm_bee(i_arm,i_a) - arm_bee(i_arm,i_a+1)
            diff_R    = arm_Rgc(i_arm,i_a) - arm_Rgc(i_arm,i_a+1)
            diff_beta = arm_beta(i_arm,i_a)- arm_beta(i_arm,i_a+1)
            diff_d    = arm_dist(i_arm,i_a)- arm_dist(i_arm,i_a+1)

c           Slopes...
            slope_Rgc  = diff_R   / diff_e                ! dR/dell
            slope_beta = diff_beta/ diff_e                ! dbeta/dell
            slope_b    = diff_b   / diff_e                ! dbee/dell
            slope_d    = diff_d   / diff_e                ! dd/dell

c           Interpolate in trial steps of 0.1 of entries
            do j_a = 1, 10

               step   = j_a * (0.1d0*diff_e)

               Rgc_t  = arm_Rgc(i_arm,i_a)  + slope_Rgc*step    ! kpc
               beta_t = arm_beta(i_arm,i_a) + slope_beta*step   ! deg
               beta_t_rad = beta_t * deg_to_rad                 ! rad

               dist_t = arm_dist(i_arm,i_a) + slope_d*step      ! kpc
               bee_t  = arm_bee(i_arm,i_a)  + slope_b*step      ! deg
               bee_t_rad  = bee_t * deg_to_rad                  ! rad

c              arm model Galactocentric cartesian coordinates...
               x_m    = Rgc_t  * sin( beta_t_rad )        ! kpc
               y_m    = Rgc_t  * cos( beta_t_rad )
               z_m    = dist_t * sin( bee_t_rad  )

c              difference source and arm model values...
               del_x = x_galcen - x_m                     ! kpc
               del_y = y_galcen - y_m
               del_z = z_galcen - z_m
               del_in   = sqrt( del_x**2 + del_y**2 )     ! kpc
               del_dist = sqrt( del_x**2 + del_y**2 + del_z**2 ) ! kpc

               if ( del_in .lt. del_min_in ) then
c                 Save best values...
                  del_min = del_dist                      ! kpc
                  del_min_in = del_in                     ! kpc
                  del_min_z  = del_z                      ! kpc
                  del_min = del_dist                      ! kpc
                  Rgc_min = Rgc_t                         ! kpc
                  b_min   = bee_t                         ! deg
               endif

            enddo

          endif    ! diff_e=0 check

         enddo

      endif  ! source within 3 degrees in longitude from arm

      return
      end

c======================================================================

      subroutine arm_model_prob ( Ro, num_arms, iarm_entries,
     +                            ell_dif_max,
     +                            arm_name, arm_probabilities,
     +                            arm_ell, arm_bee, arm_vel,
     +                            arm_Rgc, arm_beta, arm_dist,
     +                            width_min, width_Rref, width_slope,
     +                            ell, bee, dist,
     +                            arm_max, prob_max )

c     For each arm and given a trial source distance,
c     calculate the probability of it being in the arm via
c     P(d) = P(d|arm)*P(arm)

      implicit real*8 (a-h,o-z)

      real*8   arm_probabilities(29)
      real*8   arm_ell(29,300), arm_bee(29,300), arm_vel(29,300)
      real*8   arm_Rgc(29,300), arm_beta(29,300), arm_dist(29,300)

      integer  iarm_entries(29)

      character*12  arm_name(29), arm_max

      prob_max = 0.01d0        ! minimum (starting) value for Prob(d)
      cutoff   = 0.01d0        ! minimum value for P(arm)
      arm_max  = '...'

      do i_a = 1, num_arms

         if ( arm_probabilities(i_a) .gt. cutoff ) then

c           Calculate Prob(d|arm)
            i_arm = i_a
            call calc_arm_prob( Ro, ell, bee, dist,
     +                          i_arm, iarm_entries, ell_dif_max,
     +                          arm_ell, arm_bee, arm_vel,
     +                          arm_Rgc, arm_beta, arm_dist,
     +                          width_min, width_Rref, width_slope,
     +                          x_mod, y_mod, Rsrc,
     +                          del_min, prob_dist )

c           Multiply probabilities...
c           Prob(d) = Prob(d|arm) * Prob(arm)
            prob = prob_dist * arm_probabilities(i_a)

            if ( prob .gt. prob_max ) then
c              Store max probability and arm name
               prob_max = prob
               arm_max  = arm_name(i_a)
            endif

         endif     ! arm check

      enddo    ! all spiral arms

      return
      end

c===================================================================

      subroutine calc_arm_prob(Ro, ell, bee, dist,
     +                         i_arm, iarm_entries, ell_dif_max,
     +                         arm_ell, arm_bee, arm_vel,
     +                         arm_Rgc, arm_beta, arm_dist,
     +                         width_min, width_Rref, width_slope,
     +                         x_mod, y_mod, Rsrc,
     +                         del_min, prob)

c     Calculates probability density, Prob(d|arm), for arm number "i_arm"
c     at a given distance (dist)

      implicit real*8 ( a-h, o-z )

      real*8  arm_ell(29,300), arm_bee(29,300), arm_vel(29,300)
      real*8  arm_Rgc(29,300), arm_beta(29,300), arm_dist(29,300)

      integer iarm_entries(29)

      logical near_arm

      pi = 4.d0*atan(1.d0)
      deg_to_rad = pi/180.d0

c     Source values...
      ell_rad  = ell * deg_to_rad                      ! radians
      bee_rad  = bee * deg_to_rad                      ! radians
      x_galcen = dist * sin(ell_rad)                   ! kpc
      y_galcen = Ro - dist * cos(ell_rad)              ! kpc
      z_galcen = dist * sin(bee_rad)                   ! kpc
      Rsrc     = sqrt( x_galcen**2 + y_galcen**2 )     ! kpc
      beta_src = atan2( x_galcen, y_galcen )           ! radians
      beta_src_deg = beta_src / deg_to_rad             ! deg

c     Go through table of (l,b,V,R,beta,d) for arm "i_arm" to find
c     closest arm_beta to the source_beta value.   Use beta, not ell,
c     since it is much better behaved.
      i_a_max = iarm_entries(i_arm)
      del_beta_min = 999.d0
      i_a_min      = 0
      do i_a = 1, i_a_max

c        Check that arm entry has estimated values...
         if ( arm_Rgc(i_arm,i_a)  .gt. 0.d0 .and.
     +        arm_dist(i_arm,i_a) .gt. 0.d0      ) then

            bs = beta_src_deg                           ! deg
            ba = arm_beta(i_arm,i_a)                    ! deg

c           To allow beta values beyond -180 to +180, need 3 checks
            del_beta1 = abs( bs - ba )                  ! deg
            del_beta2 = abs( bs - ba - 360.d0 )         ! deg
            del_beta3 = abs( bs - ba + 360.d0 )         ! deg
            if ( del_beta1 .lt. del_beta_min ) then
               del_beta_min = del_beta1
               i_a_min = i_a
            endif
            if ( del_beta2 .lt. del_beta_min ) then
               del_beta_min = del_beta2
               i_a_min = i_a
            endif
            if ( del_beta3 .lt. del_beta_min ) then
               del_beta_min = del_beta3
               i_a_min = i_a
            endif

         endif

      enddo

c     Then get arm longitude of point with closest beta
      ell_arm_closest = arm_ell(i_arm,i_a_min)            ! deg
      del_ell_1       = abs( ell - ell_arm_closest )
c     Allow for possible 360 deg differences...
      del_ell_2       = abs( del_ell_1 - 360.d0 )
      del_ell_3       = abs( del_ell_1 + 360.d0 )

      del_min = 999.d0                                    ! kpc
      Rgc_min = 999.d9

c     Longitudes must be close to likely be in this arm
c     Use "ell_dif_max" degrees difference as dividing line.
      near_arm = .false.
      if ( del_ell_1 .lt. ell_dif_max ) then
         near_arm = .true.
         ambiguity = 0.d0
      endif
      if ( del_ell_2 .lt. ell_dif_max ) then
         near_arm = .true.
         ambiguity = -360.d0
      endif
      if ( del_ell_3 .lt. ell_dif_max ) then
         near_arm = .true.
         ambiguity = +360.d0
      endif

      i_ell_dif_max = int(ell_dif_max + 0.5)
      if ( near_arm ) then

         i_a_0 = i_a_min - i_ell_dif_max
         i_a_1 = i_a_min + i_ell_dif_max
         if ( i_a_0 .lt. 1 ) i_a_0 = 1
         if ( i_a_1 .gt. i_a_max ) i_a_1 = i_a_max

c        Now go through in steps of 0.1 deg (actually 0.1 of
c        difference between longitude entries) and calculate
c        the minimum 2-D separation of source from arm center

         i_a_1m1 = i_a_1-1
         do i_a = i_a_0, i_a_1m1

            diff_e    = arm_ell(i_arm,i_a) - arm_ell(i_arm,i_a+1)

c           Check for and skip a duplicate entry in lbvRBD data,
c           which would lead to diff_e=0 and then divide by zero...
            if ( abs(diff_e) .gt. 0.001d0 ) then

c              Don't rely on previously determined ambiguity
               call fix_ambiguity ( diff_e, ambiguity )

               diff_b    = arm_bee(i_arm,i_a) - arm_bee(i_arm,i_a+1)
               diff_R    = arm_Rgc(i_arm,i_a) - arm_Rgc(i_arm,i_a+1)
               diff_beta = arm_beta(i_arm,i_a)- arm_beta(i_arm,i_a+1)
               diff_d    = arm_dist(i_arm,i_a)- arm_dist(i_arm,i_a+1)

c              Slopes...
               slope_Rgc  = diff_R   / diff_e                       ! dR/dell
               slope_beta = diff_beta/ diff_e                       ! dbeta/dell
               slope_b    = diff_b   / diff_e                       ! dbee/dell
               slope_d    = diff_d   / diff_e                       ! dd/dell

c              Interpolate in trial steps of 10% of longitude entry difference
               do j_a = 1, 10

                  step   = j_a * (0.1d0*diff_e)

                  Rgc_t  = arm_Rgc(i_arm,i_a)  + slope_Rgc*step     ! kpc
                  beta_t = arm_beta(i_arm,i_a) + slope_beta*step    ! deg
                  beta_t_rad = beta_t * deg_to_rad                  ! rad

                  dist_t = arm_dist(i_arm,i_a) + slope_d*step       ! kpc
                  bee_t  = arm_bee(i_arm,i_a)  + slope_b*step       ! deg
                  bee_t_rad  = bee_t * deg_to_rad                   ! rad

c                 arm model Galactocentric cartesian coordinates...
                  cos_bee = cos( bee_t_rad )
                  x_m    = Rgc_t  * sin( beta_t_rad )*cos_bee       ! kpc
                  y_m    = Rgc_t  * cos( beta_t_rad )*cos_bee
                  z_m    = dist_t * sin( bee_t_rad  )

c                 difference source and arm model values...
                  del_x = x_galcen - x_m ! kpc
                  del_y = y_galcen - y_m
                  del_z = z_galcen - z_m
                  del_in   = sqrt( del_x**2 + del_y**2 ) ! kpc
                  del_dist = sqrt( del_x**2 + del_y**2 + del_z**2 ) ! kpc

                  if ( del_dist .lt. del_min ) then

c                    Save best values...
                     del_min_in = del_in                            ! kpc
                     del_min_z  = del_z                             ! kpc
                     Rgc_min    = Rgc_t                             ! kpc

                     del_min    = del_dist                          ! kpc

                  endif

               enddo

            endif     ! duplicate entry (diff_e=0) check

         enddo

      endif  ! source within "ell_dif_max" degrees in longitude from arm

c     Have minimum offset distance between source and center of arm,
c     now put a sign on this value with -'ve for Rsrc < Rmod
      if ( Rsrc .lt. Rgc_min ) del_min = - del_min          ! kpc

      prob = 0.d0
      if ( abs(del_min) .lt. 2.0d0 ) then

c           Now calculate probability based on arm width vs Radius
            if ( Rsrc .gt. width_Rref ) then
                  width = width_min +
     +                   (Rsrc-width_Rref)*width_slope      ! kpc
               else
                  width = width_min                         ! kpc
            endif
            zero  = 0.d0
            call gauss_prob (zero,width,del_min_in,prob_in)

c           Assume z-sigma is roughly one-third of the in-plane
c           of the spiral arm sigma
            width_z = width / 3.d0
            call gauss_prob (zero,width_z,del_min_z,prob_z)

            prob = prob_in * prob_z                     ! 1/kpc^2

      endif

      return
      end

c===================================================================
      subroutine fix_ambiguity ( ang_dif_deg, ambiguity )

c     Checks an angular difference (degrees), which should be near zero,
c     for +/-1 wrap ambiguity and returns the minimum difference

      implicit real*8 (a-h,o-z)

      ambiguity = 0.d0

      try1 = ang_dif_deg - 360.d0
      if ( abs(try1) .lt. 180.d0 ) ambiguity = -360.d0

      try2 = ang_dif_deg + 360.d0
      if ( abs(try2) .lt. 180.d0 ) ambiguity = +360.d0

      ang_dif_deg = ang_dif_deg + ambiguity

      return
      end

c===================================================================

      subroutine calc_Dk_Univ ( lu_out,
     +        Ro, a1, a2, a3, Uo, Vo, Wo, Us, Vs, Ws,
     +        gal_long, gal_lat, farnear, v_lsr,
     +        v_lsr_rev, Dk )

c     Calculate revised Vlsr by converting standard Vlsr back to
c     heliocentric, apply modern Solar Motion values (Uo,Vo,Wo), and
c     remove effects of average source non-ciruclar motion (Us,Vs,Ws).
c     Then calculate kinematic distance using the Persic "Universal"
c     rotation curve specified by Ro, a1, a2, and a3

      implicit real*8 (a-h,o-z)

      character*12   src

      pi      = 4.d0*atan(1.d0)
      deg_rad = pi/180.d0

      n_iter_max = 100

      gal_long_rad = gal_long * deg_rad       ! radians
      cos_l = cos( gal_long_rad )
      sin_l = sin( gal_long_rad )

      gal_lat_rad  = gal_lat * deg_rad        ! radians
      cos_b = cos( gal_lat_rad )
      sin_b = sin( gal_lat_rad )

c     -------------------------------------------------------------
c     Convert to true Heliocentric frame (v_rad)
c     Add back Sun's peculiar motion to radial velocity
c     Use old Standard Solar Motion (Sun moves at 20 km/s toward
c     18h, +30deg in 1900 coordinates), since this has been used to
c     define V_lsr:

      Uo_IAU = 10.27d0            ! km/s precessed to J2000
      Vo_IAU = 15.32d0
      Wo_IAU =  7.74d0

      v_helio = v_lsr - (Vo_IAU*sin_l + Uo_IAU*cos_l)*cos_b
     +                -  Wo_IAU*sin_b

c     -------------------------------------------------------------
c     Make "new" V(LSR) using best Solar Motion

      v_newlsr = v_helio + (Vo*sin_l + Uo*cos_l)*cos_b
     +                   +  Wo*sin_b

c     --------------------------------------------------------
c     Remove effects of common peculiar motions specified in
c     "galaxy_data_Univ.inp" file

c     In general, need to know distance to get source Galactocentric
c     radius (Rs) to evalute rotation curve.   So must iterate...
      n_iter =  0
      del_d  = 99.d0
      Dk     =  3.d0

      do while ( del_d.gt.0.01d0 .and. n_iter.lt.n_iter_max )

c        Save old value of kinematic distance
         Dk_old = Dk

c        Calculate "gamma" angle and projected Galactocentric radius
         d_proj = Dk * cos_b                     ! kpc in Gal Plane
         r_sq   = Ro**2 + d_proj**2 - 2.d0*Ro * d_proj * cos_l
         r_proj = sqrt( r_sq )                   ! kpc in Gal Plane

c        Calculate Galactocentric longitude (beta in paper)...
         sin_beta =   d_proj  * sin_l     / r_proj
         cos_beta = ( Ro - d_proj*cos_l ) / r_proj
         beta     = atan2( sin_beta, cos_beta )  ! radians
         beta_deg = beta / deg_rad               ! deg

c        Calculate Sun-Maser-GC angle...
         gamma = pi - gal_long_rad - beta        ! radians
         cos_gamma = cos( gamma )
         sin_gamma = sin( gamma )

         v_fixed = v_newlsr - (Vs*sin_gamma - Us*cos_gamma)*cos_b
     +                      -  Ws*sin_b          ! km/s

c        -----------------------------------------------------------------
c        Calculate a kinematic distance using best Ro, To and dTdr

         Rs = r_proj
         V_proj = v_fixed * cos_b
         call kinematic_distance_Univ ( v_proj, gal_long,
     +                                  Ro, a1, a2, a3,
     +                                  Rs,
     +                                  D_near, D_far )

         Dk = D_near
         if ( farnear .ne. 0.d0 ) Dk = D_far

c        Ignore "farnear" flag if one of the values is zero
         if ( D_near .le. 0.d0 .and. D_far  .gt. 0.d0 ) Dk = D_far
         if ( D_far  .le. 0.d0 .and. D_near .gt. 0.d0 ) Dk = D_near

         del_d = abs( Dk - Dk_old )
         n_iter = n_iter + 1

      enddo

      v_lsr_rev = v_fixed

      return
      end

c======================================================================

      subroutine Dk_prob_density ( a1, a2, a3, Ro, To,
     +           Uo, Vo, Wo, Us, Vs, Ws,
     +           ell, bee, dist, v_lsr,
     +           Dk_near, Dk_far, far_prob, sig_vel,
     +           prob_Dk )

c     Calculates probability density for a given trial distance.
c     This approach avoids actually calculating the kinematic distance(s),
c     which can be very complicated when a random velocity component is
c     considered for an uncertainty estimate.

c     Accounts for Prob(Far) entered by user.

c     V3 will incorporate a sig_vel that increases inward of 6 kpc from GC

      implicit real*8 (a-h,o-z)

      pi     = 4.d0*atan(1.d0)
      degrad = pi/180.d0

      sin_ell = sin(ell*degrad)
      cos_ell = cos(ell*degrad)
      cos_bee = cos(bee*degrad)

      lu_out = -1

c     Get circular rotation speed (Tr) at the Galactocentric radius (Rs)...
      d_proj = dist * cos_bee                             ! kpc in Gal Plane
      R_sq   = Ro**2 + d_proj**2 - 2.d0*Ro*d_proj*cos_ell
      Rs     = sqrt( R_sq )                               ! kpc in Gal Plane
      call Univ_RC_from_note ( Rs, a2, a3, Ro,
     +                         Tr)
c     Calculate To for this rotation curve...
      call Univ_RC_from_note ( Ro, a2, a3, Ro,
     +                         To)

c     Now get model observed circular velocity
      if ( Rs .gt. 0.01d0*Ro ) then
            V_mod = Tr*sin_ell*(Ro/Rs) - To*sin_ell       ! km/s
         else
            V_mod = 0.d0
      endif

c     Get revised V_lsr for the assumed Solar & average source parameters...
      call calc_revised_Vlsr ( lu_out,
     +                         Ro, Uo, Vo, Wo, Us, Vs, Ws,
     +                         ell, bee, v_lsr, dist,
     +                         v_lsr_rev )

c     Compare Vmodel with V_lsr_rev to calculate a distance PDF...
c     First calculate midpoint for two kinematic distances for weighting of
c     the near and far PDF components
      d_mid = 999.d0
      if ( Dk_near.gt.0.d0 .and. Dk_far.gt.Dk_near ) then
         d_mid = 0.5d0 * (Dk_near + Dk_far)
      endif

c     Increase expected velocity uncertainty interior to 6 kpc to
c     account for larger peculiar motions seen in the bar region.
c     Increase peculiar velocity term linearly between 6 and 4 kpc
c     and cap it at 25 km/s
      Rs_start = 6.d0
      if ( Rs .lt. Rs_start ) then
            vel_pec = (Rs_start - Rs) * 12.5d0           ! km/s
            if ( vel_pec .gt. 25.d0 ) vel_pec = 25.d0
            sig_vel_incr = sqrt( sig_vel**2 + vel_pec**2 )
         else
            sig_vel_incr = sig_vel
      endif

c     Increase velocity uncertainty for sources in Perseus arm known
c     to have Vlsr kinematically anomalies of ~20 km/s within a "box"
c     of 1.0 < Xgc < 3.2 kpc and 8.8 < Ygc < 10.0 kpc
      Xgc = dist * sin_ell * cos_bee                     ! kpc
      Ygc = Ro - dist * cos_ell * cos_bee
      vel_pec = 20.d0                                    ! km/s
      if ( Xgc.ge.1.0d0 .and. Xgc.le.3.2d0  .and.
     +     Ygc.ge.8.8d0 .and. Ygc.le.10.0d0      ) then
         sig_vel_incr = sqrt( sig_vel_incr**2 + vel_pec**2 )
      endif

      sig_vel_avg = sig_vel                              ! km/s
      call asym_conservative_prob( V_mod, sig_vel_incr, sig_vel_avg,
     +                             v_lsr_rev, prob )

c     Weight probability density by user-supplied near/far probability
c     if appropriate (in quadrants 1 and 4)
      if ( ell.lt.90.d0 .or. ell.gt.270.d0) then

            if ( dist .le. d_mid ) then
c                 PDF weighted for near distance
                  prob_Dk = prob * (1.d0 - far_prob)
               else
c                 PDF weighted for far distance
                  prob_Dk = prob * far_prob
            endif

         else

c           In Quadrants 2 and 3 there is no ambiguity
            prob_Dk = prob

      endif

      return
      end

c======================================================================

      subroutine pm_ell_prob_density ( a1, a2, a3, Ro, To,
     +           Uo, Vo, Wo, Us, Vs, Ws,
     +           ell, bee, dist, pm_ell, pm_ell_unc, sig_vel,
     +           p_Dpm_ell )


c     Calculates probability density based on proper motion in
c     Galactic longitude for a given trial distance (dist) toward source.

c     As for kinematic distances, this version avoids actually calculating
c     the "motion" distance(s), which can be very complicated when a
c     random velocity component is considered for an uncertainty estimate.

c     V3 will incorporate a sig_vel that increases inward of 6 kpc from GC

c     Currently NOT using the Us,Vs correction

      implicit real*8 (a-h,o-z)

      pi      = 4.d0*atan(1.d0)
      deg_rad = pi/180.d0

      ell_rad = ell*deg_rad
      sin_ell = sin(ell_rad)
      cos_ell = cos(ell_rad)

      cos_bee = cos(bee*deg_rad)

c     Get circular rotation speed (Tr) at the source's Galactocentric radius (Rs)...
      d_proj = dist * cos_bee                                        ! kpc in Gal Plane
      R_sq   = Ro**2 + d_proj**2 - 2.d0*Ro*d_proj*cos_ell
      Rs     = sqrt( R_sq )                                          ! kpc in Gal Plane

      call Univ_RC_from_note ( Rs, a2, a3, Ro,
     +                         Tr )
c     Get To from rotation curve
      call Univ_RC_from_note ( Ro, a2, a3, Ro,
     +                         To )

c     Now get model proper motion in Galactic longitude...
      if ( Rs .gt. 0.01d0*Ro ) then

c           Calculate Galactocentric azimuth (beta)
            cos_beta = (Ro - d_proj*cos_ell) / Rs
            sin_beta = d_proj*sin_ell / Rs
            beta_rad = atan2( sin_beta, cos_beta )                   ! radians
c           Calculate angle (phi) between pm_ell and circular rotation
            phi_rad = ell_rad + beta_rad                             ! radians
            cos_phi = cos( phi_rad )

            vel_ell = Tr*cos_phi - (To+Vo)*cos_ell + Uo*sin_ell      ! km/s
            pm_ell_mod_kmskpc  = vel_ell / d_proj                    ! km/s/kpc
            pm_ell_mod = pm_ell_mod_kmskpc / 4.74d0                  ! mas/yr

         else

            pm_ell_mod = 0.d0

      endif

c     Compare model pm_ell with measured value to calculate a
c     distance PDF, accounting for both measurement uncertainty
c     and virial motions...

c     Increase expected velocity uncertainty interior to 6 kpc to
c     account for larger peculiar motions seen in the bar region.
c     Increase peculiar velocity term linearly between 6 and 4 kpc
c     and cap it at 25 km/s
      Rs_start = 6.d0
      if ( Rs .lt. Rs_start ) then
            vel_pec = (Rs_start - Rs) * 12.5d0           ! km/s
            if ( vel_pec .gt. 25.d0 ) vel_pec = 25.d0
            sig_vel_incr = sqrt( sig_vel**2 + vel_pec**2 )
         else
            sig_vel_incr = sig_vel
      endif

c     Increase velocity uncertainty for sources in Perseus arm known
c     to have Vlsr kinematically anomalies of ~20 km/s within a "box"
c     of 1.0 < Xgc < 3.2 kpc and 8.8 < Ygc < 10.0 kpc
      Xgc = dist * sin_ell * cos_bee                     ! kpc
      Ygc = Ro - dist * cos_ell * cos_bee
      vel_pec = 20.d0                                    ! km/s
      if ( Xgc.ge.1.0d0 .and. Xgc.le.3.2d0  .and.
     +     Ygc.ge.8.8d0 .and. Ygc.le.10.0d0      ) then
         sig_vel_incr = sqrt( sig_vel_incr**2 + vel_pec**2 )
      endif

      sig_pm_incr = sig_vel_incr / d_proj / 4.74d0                 ! mas/yr
      sig_pm      = sqrt ( sig_pm_incr**2 + pm_ell_unc**2 )        ! mas/yr

      sig_pm_avg  = sqrt( sig_vel**2 * pm_ell_unc**2 )             ! km/s
      call asym_conservative_prob( pm_ell_mod, sig_pm, sig_pm_avg,
     +                             pm_ell, prob )

      p_Dpm_ell = prob


      return
      end

c======================================================================

      subroutine pm_bee_prob_density ( Wo, ell, bee,
     +                                 pm_bee, pm_bee_unc,
     +                                 dist, sig_vel,
     +                                 p_Dpm_bee )

c     Calculates probability density based on Galactic latitude motion

c     Will NOT incorporate a sig_vel that increases inward of 6 kpc
c     from GC, since one expects enhanced peculiar motions to be almost
c     entirely in the Galactic plane

      implicit real*8 (a-h,o-z)

      pi      = 4.d0*atan(1.d0)
      deg_rad = pi/180.d0

      d_proj = dist * cos(bee*deg_rad)                               ! kpc

      sig_pm_virial = sig_vel / d_proj / 4.74d0                      ! mas/yr
      sig_pm     = sqrt( sig_pm_virial**2 + pm_bee_unc**2 )          ! mas/yr
      sig_pm_avg = sig_pm

c     Expect to see reflex of Solar Motion in z-direction...
      pm_bee_mod = -Wo / d_proj / 4.74d0                             ! mas/yr
      call asym_conservative_prob( pm_bee_mod, sig_pm, sig_pm_avg,
     +                             pm_bee, prob )

      p_Dpm_bee = prob


      return
      end

c===================================================================

      subroutine calc_revised_Vlsr ( lu_out,
     +                          Ro, Uo, Vo, Wo, Us, Vs, Ws,
     +                          gal_long, gal_lat, v_lsr, dist,
     +                          v_lsr_rev )

c     Calculate revised Vlsr by
c     1) converting standard Vlsr back to Heliocentric,
c     2) appling modern Solar Motion values (Uo,Vo,Wo), and
c     3) removing effects of average source non-circular motion (<Us>,<Vs>,<Ws>).

      implicit real*8 (a-h,o-z)

      pi     = 4.d0 * atan(1.d0)
      deg_rad = pi/180.d0

      gal_long_rad = gal_long * deg_rad       ! radians
      cos_l = cos( gal_long_rad )
      sin_l = sin( gal_long_rad )

      gal_lat_rad  = gal_lat * deg_rad        ! radians
      cos_b = cos( gal_lat_rad )
      sin_b = sin( gal_lat_rad )

c     -------------------------------------------------------------
c     Convert to true Heliocentric frame (v_rad)
c     Add back Sun's peculiar motion to radial velocity
c     Use old Standard Solar Motion (Sun moves at 20 km/s toward
c     18h, +30deg in 1900 coordinates), since this has been used to
c     define V_lsr:

      Uo_IAU = 10.27d0            ! km/s precessed to J2000
      Vo_IAU = 15.32d0
      Wo_IAU =  7.74d0

      v_helio = v_lsr - (Vo_IAU*sin_l + Uo_IAU*cos_l)*cos_b
     +                -  Wo_IAU*sin_b

c     -------------------------------------------------------------
c     Make "new" V(LSR) using best Solar Motion

      v_newlsr = v_helio + (Vo*sin_l + Uo*cos_l)*cos_b
     +                   +  Wo*sin_b

c     --------------------------------------------------------
c     Remove effects of common peculiar motions specified in
c     "parameter_file.inp"

c     Calculate "gamma" angle and projected Galactocentric radius
      d_proj = dist * cos_b                     ! kpc in Gal Plane
      r_sq   = Ro**2 + d_proj**2 - 2.d0*Ro * d_proj * cos_l
      r_proj = sqrt( r_sq )                   ! kpc in Gal Plane

c     Calculate Galactocentric longitude (beta in paper)...
      sin_beta =   d_proj  * sin_l     / r_proj
      cos_beta = ( Ro - d_proj*cos_l ) / r_proj
      beta     = atan2( sin_beta, cos_beta )  ! radians
      beta_deg = beta / deg_rad               ! deg

c     Calculate Sun-Maser-GC angle...
      gamma = pi - gal_long_rad - beta        ! radians
      cos_gamma = cos( gamma )
      sin_gamma = sin( gamma )

      v_fixed = v_newlsr - (Vs*sin_gamma - Us*cos_gamma)*cos_b
     +                   -  Ws*sin_b          ! km/s

      v_lsr_rev = v_fixed

      return
      end


c     =======================================================

      subroutine kinematic_distance_Univ ( Vlsr, gal_long,
     +                                     Ro, a1, a2, a3,
     +                                     Rs,
     +                                     D_near, D_far )

c     Caluclate kinematic distance given Vlsr (projected in Gal plane)
c     and information required to construct the kinematic model.
c     Returns both near and far distances (when appropriate).
c     Uses Persic's "Universal" rotation curve formulation.

      implicit real*8 (a-h,o-z)

      pi      = 4.0d0*atan(1.0d0)
      deg_rad = pi/180.d0	     ! convert degrees to radians

      glongrad = gal_long*deg_rad

      cos_l = cos(glongrad)
      sin_l = sin(glongrad)

      Rosinl = Ro * sin_l
      Rocosl = Ro * cos_l

      call Univ_RC_from_note ( Rs, a2, a3, Ro,
     +                         Tr )
c     Get To from rotation curve
      call Univ_RC_from_note ( Ro, a2, a3, Ro,
     +                         To )

      Tosinl = To * sin_l
      Trsinl = Tr * sin_l

      rootterm = Rocosl**2 + ( Trsinl / ( Tosinl/Ro + Vlsr/Ro ) )**2
     +         - Ro**2
      if ( rootterm .lt. 0 ) rootterm = 0.d0

      if ( gal_long.ge.0.d0   .and. gal_long.lt.90.d0 ) then
         D_near = Rocosl - sqrt( rootterm )
         D_far  = Rocosl + sqrt( rootterm )
      endif

      if ( gal_long.ge.90.d0  .and. gal_long.le.270.d0 ) then
         D_near = Rocosl + sqrt( rootterm )
         D_far  = D_near
      endif

      if ( gal_long.gt.270.d0 .and. gal_long.lt.360.d0 ) then
         D_near = Rocosl - sqrt( rootterm )
         D_far  = Rocosl + sqrt( rootterm )
      endif

      return
      end

c============================================================================

      subroutine get_Vmax ( ell, Ro, a1, a2, a3,
     +                      V_max )

      implicit real*8 (a-h,o-z)

c     Calculates the maximum line of sight velocity for a given longitude
c     for a Persic "Universal" rotation curve.  This may differ from the
c     "tangent point" velocity.

      V_max = 0.d0

c     Only valid for quadrants 1 and 4 (and avoid edges)
      if ( (ell.gt.0.1d0   .and. ell.lt.89.9d0 ) .or.
     +     (ell.gt.270.1d0 .and. ell.lt.359.9d0)      ) then

         pi = 4.d0*atan(1.d0)
         deg_rad = pi/180.d0

         ell_rad = ell*deg_rad
         sin_l = sin(ell_rad)
         cos_l = cos(ell_rad)

c        Walk outward from Sun and find maximum V
c        (positive in Q1; negative in Q4)
         do i_d = 1, 100

            d = i_d*0.1d0                          ! kpc

            y = Ro - d*cos_l                       ! kpc
            x = d*sin_l
            Rs= sqrt( x**2 + y**2 )

            sin_beta = x / Rs
            cos_beta = y / Rs
            beta     = atan2( sin_beta, cos_beta ) ! rad
            ellbeta  = ell_rad + beta              ! rad

            call Univ_RC_from_note ( Rs, a2, a3, Ro,
     +                               Tr )
c           Get To from rotation curve
            call Univ_RC_from_note ( Ro, a2, a3, Ro,
     +                               To )

            V_los = Tr*sin(ellbeta) - To*sin_l     ! km/s

            if ( abs(V_los).gt.abs(V_max) ) V_max = V_los

         enddo         ! trial distance loop

      endif       ! quadrant check

      return
      end

c============================================================================
           subroutine Univ_RC_from_note ( r, a2, a3, Ro,
     +                                    Tr )

           implicit real*8 (a-h,o-z)

           real*8   lambda, log_lam

c          Disk plus halo parameterization of rotation curve...
c          see Persic, Salucci and Stel 1996 "Note Added in Proof"

c          NB: doesn't use "a1" parameter

           lambda = (a3/1.5d0)**5                          ! L/L*

           Ropt   = a2 * Ro     ! a2 = Ropt/Ro; (NB: Ropt=3.2Rd encloses 83% of light)
                                ! where Rd=scale length ~ 2.5 kpc if Ro ~ 8.1
                                ! P.C. van der Kruit & K.C. Freeman ARAA, 2011, 49, 301

           rho    = r/Ropt

c          Calculate Tr...
           log_lam = log10( lambda )

           term1 = 200.d0*lambda**0.41d0

           top   = 0.75d0*exp(-0.4d0*lambda)
           bot   = 0.47d0 + 2.25d0*lambda**0.4d0
           term2 = sqrt( 0.80d0 + 0.49d0*log_lam + (top/bot) )

           top   = 1.97d0*rho**1.22d0
           bot   = (rho**2 + 0.61d0)**1.43d0
           term3 = (0.72d0 + 0.44*log_lam) * (top/bot)

           top   = rho**2
           bot   = rho**2 + 2.25*lambda**0.4d0
           term4 = 1.6d0*exp(-0.4d0*lambda) * (top/bot)

           Tr    = (term1/term2) * sqrt(term3 + term 4)            ! km/s

c
           return
           end

c==============================================================================
c The following has been replaced by   subroutine Univ_RC_from_note
      subroutine Univ_RC ( a1, a2, a3, Ro, Rs,
     +                     Tr, To )

      implicit real*8 (a-h,o-z)

c     Disk plus halo parameterization of rotation curve...
c     see http://ned.ipac.caltech.edu/level5/March01/Battaner/node7.html
c     referenced as following  Persic and Salucci (1997)

c     Returns T(Rs) and T(Ro)

      beta   = 0.72d0      ! 0.72 + 0.44 ln(L/L*)
      V2_Ropt= a1**2       ! a1 = V(Ropt)
      Ropt   = a2 * Ro     ! a2 = Ropt/Ro; (NB: Ropt=3.2Rd encloses 83% of light)
                           ! where Rd=scale length ~ 4.5 kpc (but often ~3 kpc used)
                           ! P.C. van der Kruit & K.C. Freeman ARAA, 2011, 49, 301
      a      = a3          ! "a"     = 1.5(L/L*)^(1/5)

c     So, for L = L*, beta=0.72; a1~To; a2~3.2*(3.0to4.5)/8.3 = 1.2to1.7; a3~1.5

      do i_rho = 1, 2

c        Calculate T(Rs) and then T(Ro)...
         rho    = Rs/Ropt
         if ( i_rho .eq. 2 ) rho = Ro/Ropt

c        Exponential thin disk component...
         top = 1.97d0*rho**1.22
         bot = (rho**2 + 0.78d0**2)**1.43d0
         V2_disk = V2_Ropt * beta * top / bot

c        Halo component...
         top = rho**2
         bot = rho**2 + a**2
         V2_halo = V2_Ropt * (1.d0-beta) * (1.d0+a**2)*top/bot

c        Total...
         Theta = sqrt( V2_disk + V2_halo ) ! km/s

         if ( i_rho .eq. 1 ) then
               Tr = Theta
            else
               To = Theta
         endif

      enddo

      return
      end

c======================================================================

      subroutine lat_prob_warped ( ell, bee, dist, Ro,
     +                             Rsigz, sigz, sigzdot,
     +                             n_bees, d_store, b_store,
     +                             prob, arg )

c     Calculates latitude Prob(d|b,sigz,warping_info) that a source is at a
c     given distance (dist)

      implicit real*8 (a-h,o-z)

      real*8   d_store(29), b_store(29)

      pi = 3.141592654d0
      twopi = 2.d0 * pi
      degrad = pi / 180.d0

c     Calculate expected z-sigma for arm at Rgc (from dist and Ro)
c     using a constant for Rgc<Rsigz and sloping up for Rgc>Rsigz
      cosl = cos( ell*degrad )
      Rgc = sqrt( Ro**2 + dist**2 - 2.d0*Ro*dist*cosl ) ! kpc
      if ( Rgc .lt. Rsigz ) then
            sigma_z = sigz                              ! kpc
         else
            sigma_z = sigz + sigzdot*(Rgc-Rsigz)        ! kpc
      endif

c     First, get warping value..
      call interpolate_warping ( n_bees, d_store, b_store, dist,
     +                           b_warp )

      del_bee     = bee - b_warp                        ! deg
      del_bee_rad = del_bee * degrad                    ! rad
      del_z       = del_bee_rad * dist                  ! kpc

c     Calculate Prob(d|b) = Prob(b|d) * Prob(b)/Prob(d)
c     For flat priors on b and d drop the ratio term.
c     Next, Prob(b|d) = Prob(z|d) dz/db = Prob(z|d)*d
      prefact = sqrt(1.d0/twopi) / sigma_z              ! 1/kpc
      arg     = del_z**2 / (2.d0*sigma_z**2)
      prob    = prefact * exp(-arg) * dist              ! kpc/kpc

      return
      end

c=======================================================================

      subroutine interpolate_warping ( n_bees, d_store, b_store, dist,
     +                                 b_warp )

c     Linearly interpolates between (d,b) values.
c     Returns latitude warp at specified dist.
c     Complicated because of edges

      implicit real*8 (a-h,o-z)

      real*8   d_store(29), b_store(29)

c     Default: no warp
      b_warp = 0.d0

      if ( n_bees .ge. 1 ) then

         n_b = 0
         do while ( d_store(n_b+1).lt.dist .and.
     +              n_b .lt. n_bees-1           )
            n_b = n_b + 1
         enddo
c        n_b is first index of (d,b) arrays for which source dist > d
c        It could be zero if source dist < first arm d

         if ( n_b.eq.0 ) then
c           interpolate/extrapolate from zero and first value
            dif_d  = dist                                  ! kpc
            dif_a  = d_store(n_b+1)                        ! kpc
            dif_b  = b_store(n_b+1)                        ! deg
            b_warp = (dif_d/dif_a)*dif_b                   ! deg
         endif

         if ( n_b.ge.1 .and. n_bees.gt.n_b ) then
c           interpolate between pairs of values
            dif_d  = dist - d_store(n_b)                   ! kpc
            dif_a  = d_store(n_b+1) - d_store(n_b)         ! kpc
            dif_b  = b_store(n_b+1) - b_store(n_b)         ! deg
            b_warp = b_store(n_b) + (dif_d/dif_a)*dif_b    ! deg
         endif

         if ( n_b.eq.1 .and. n_bees.eq.1 ) then
c           extrapolate past the only value
            dif_d  = dist                                  ! kpc
            dif_a  = d_store(n_b)                          ! kpc
            dif_b  = b_store(n_b)                          ! deg
            b_warp = b_store(n_b) + (dif_d/dif_a)*dif_b    ! deg
         endif

         if ( n_b.eq.n_bees .and. n_bees.ge.2 ) then
c           extrapolate past last (of 2 or more) values
            dif_d  = dist - d_store(n_b-1)                 ! kpc
            dif_a  = d_store(n_b) - d_store(n_b-1)         ! kpc
            dif_b  = b_store(n_b) - b_store(n_b-1)         ! deg
            b_warp = b_store(n_b) + (dif_d/dif_a)*dif_b    ! deg
         endif

      endif

      return
      end

c======================================================================

      subroutine open_ascii_file ( lu_in, ascii_file, lu_print )

c     Opens ascii file and reads and prints comments (if first
c     character is "!").  Leaves file open, ready to read data.

      character*48       ascii_file
      character*80       comment
      character*1        c_1
      equivalence       (c_1, comment)

      logical            first_comment

c     Open input ascii file...
      first_comment = .true.
      if ( lu_print .gt. 0 ) then
cc         write (lu_print,1000) ascii_file
 1000    format(/' Openning input ascii file "',a32,'"')
      endif

      open (unit=lu_in, file=ascii_file, status='old')

c     Read all comment lines (but stop at data)
      if ( lu_print .gt. 0 ) then
cc         write (lu_print,1010)
 1010    format(' Comment lines follow:')
      endif

      c_1 = '!'
      do while ( c_1 .eq. '!' )
	 read  (lu_in,1020) comment
 1020    format(a80)

         if ( c_1.eq.'!') then
            if ( lu_print .gt. 0 ) then
cc               write (lu_print,1030) comment
 1030          format(1x,a80)
            endif
          else
             backspace (unit=lu_in)
         endif
      enddo

      if ( lu_print .gt. 0 ) then
cc         write (lu_print,1040)
 1040    format(' End of Comment lines',/1x)
      endif

      return
      end

c===================================================================

      subroutine open_sources_info ( lu_sources, lu_print,
     +                               sources_file,
     +                               num_sources )

      implicit real*8 (a-h,o-z)

      character*48  sources_file

      character*14  src

      call open_ascii_file ( lu_sources, sources_file, lu_print )

      ieof = 0
      n    = 0
      do while ( ieof .ge. 0 )

         read (lu_sources,*,iostat=ieof) src, ell, bee, vlsr

         if ( ieof .ge. 0 ) n = n + 1

      enddo

      num_sources = n

cc      write (lu_print,1000) sources_file, num_sources
 1000 format(/'Opened source information file: "',a48,
     +       /'Found',i5,' lines of source information')

      rewind (unit=lu_sources)

      return
      end

c     ========================================================

      subroutine hmsrad ( hms, hms_rad )

C     Converts hhmmss.sss format to radians

      implicit real*8 (a-h,o-z)

      pi = 4.d0 * atan(1.d0)

      xhms = abs(hms)

      ih = xhms/10000.d0 + 1.d-10
      im = dmod(xhms,10000.d0)/100.d0 + 1.d-10
      s  = dmod(xhms,100.d0)

      hms_rad = dfloat(ih) + dfloat(im)/60.d0 + s/3600.d0
      hms_rad = hms_rad * pi / 12.d0

c     Shouldn't have negative hms, but just in case...
      if ( hms .lt. 0.d0 )  hms_rad = hms_rad + 2.d0*pi

      return
      end

c     ========================================================

      subroutine dmsrad ( dms, dms_rad )

c     Converts dd''"".""" format to radians

      implicit real*8 (a-h,o-z)

      pi = 4.d0 * atan(1.d0)
      deg_rad = pi / 180.d0

      xdms = abs(dms)

      id = xdms/10000.d0
      im = dmod(xdms,10000.d0)/100.d0
      s  = dmod(xdms,100.d0)

      dms_rad = dfloat(id) + dfloat(im)/60.d0 + s/3600.d0
      dms_rad = dms_rad * deg_rad
      if ( dms .lt. 0.d0 )  dms_rad = -dms_rad

      return
      end

c=======================================================================
        subroutine ellbee_to_radec ( l_II, b_II,
     +             ra_hr, dec_deg, ra_hhmmss, dec_ddmmss)

c       Converts Galactic coordinates (l_II,b_II) in degrees
c       to (RA,Dec) in radians in J2000 coordinates.

c       Also gives RA in hhmmss and Dec in ddmmss format

        implicit real*8 (a-h,o-z)

        real*8    l_II

	pi 	= 4.d0*atan(1.d0)
	hr_rad	= pi/12.d0	  !	/* convert hours to radians	*/
	deg_rad = pi/180.d0	  !	/* convert degrees to radians	*/
        hr_deg  = 15.d0           !     /* convert hours to degrees     */

c       J2000 Galactic plane-pole data...
        g_pole_ra = 12.d0 + 51.d0/60.d0 + 26.2817d0/3600.d0 ! hrs
        g_alpha   = (g_pole_ra + 6.d0) * hr_deg ! deg
                      ! = 282.859507

        g_pole_dec= 27.d0 + 07.d0/60.d0 + 42.013d0/3600.d0 ! deg
        g_dec     = 90.d0 - g_pole_dec ! deg
                      ! = 62.8716631
        g_33      = 33.d0 - 0.068d0 ! Longitude of ascending node of gal plane
                      ! = 32.932
                      ! Determined g_33 by requiring that the
                      ! B1950 defined G.C., when precessed to J2000,
                      ! still gives l=0.0 and b=0.0.  (works to 0.0001 deg)

c       Now convert to RA, Dec...
        cos_b = cos( b_II*deg_rad )
        sin_b = sin( b_II*deg_rad )

        cos_l = cos( l_II*deg_rad )
        sin_l = sin( l_II*deg_rad )

        sin_l33 = sin( (l_II-g_33)*deg_rad )
        cos_l33 = cos( (l_II-g_33)*deg_rad )

        sin_62  = sin( g_dec*deg_rad )
        cos_62  = cos( g_dec*deg_rad )

        sin_d = cos_b*sin_l33*sin_62 + sin_b*cos_62
        dec   = asin( sin_d )
        cos_d = cos( dec )
        dec_deg = dec/deg_rad

        temp = cos_b * sin_l33 - sin_d * sin_62
        sin_a282 = temp / ( cos_d * cos_62 )
        cos_a282 = cos_b * cos_l33 / cos_d

        a282_rad = atan2( sin_a282, cos_a282 )         ! radians
        a282_deg = a282_rad / deg_rad                  ! degrees

        ra_deg   = a282_deg + g_alpha                  ! degrees
        if ( ra_deg .lt. 0.d0   ) ra_deg = ra_deg + 360.d0
        if ( ra_deg .gt. 360.d0 ) ra_deg = ra_deg - 360.d0

        ra_hr = ra_deg / 15.d0                         ! hours

        ra_rad = ra_deg * deg_rad
        call radians_to_hhmmss ( ra_rad, ra_hhmmss)

        dec_rad = dec_deg * deg_rad
        call radians_to_ddmmss ( dec_rad, dec_ddmmss)

c        write (6,1000)  l_II,b_II, ra_hhmmss,dec_ddmmss
c 1000   format(1x,2f12.8,5x,2f15.6)

        return
        end

c     ========================================================

      subroutine radec_to_galactic ( ra, dec, ell_II, bee_II )

c     inputs:    ra     (decimal hours)  in J2000  !!!
c                dec    (decimal degrees)
c     outputs:   ell_II (galactic longitude in degrees)
c                bee_II (galactic latitude in degrees)

      implicit real*8 (a-h,o-z)

      pi        = 4.0d0*atan(1.0d0)
      hr_rad    = pi/12.d0      ! convert hours to radians
      deg_rad   = pi/180.d0     ! convert degrees to radians
      hr_deg    = 15.d0         ! convert hours to degrees

c     Galactic plane-pole data for B1950 (not used)
c     g_alpha  = 282.25d0
c     g_dec    = 62.6d0
c     g_33     = 33.0d0
c     J2000 Galactic plane-pole data...
      g_pole_ra = 12.d0 + 51.d0/60.d0 + 26.2817d0/3600.d0     ! hrs
      g_alpha   = (g_pole_ra + 6.d0) * hr_deg                 ! deg
                       ! = 282.859507
      g_pole_dec= 27.d0 + 07.d0/60.d0 + 42.013d0/3600.d0      ! deg
      g_dec     = 90.d0 - g_pole_dec                          ! deg
                       ! = 62.8716631
      g_33      = 33.d0 - 0.068d0 ! Longitude of ascending node of gal plane
                       ! = 32.932

c     Useful quantities...
      cos_a     = cos( ra*hr_rad - g_alpha*deg_rad )
      sin_a     = sin( ra*hr_rad - g_alpha*deg_rad )
      cos_d     = cos( dec*deg_rad )
      sin_d     = sin( dec*deg_rad )
      cos_gd    = cos( g_dec*deg_rad )
      sin_gd    = sin( g_dec*deg_rad )

c     Now calculate central galactic coordinates given ra,dec...
      sin_bII   = sin_d*cos_gd - cos_d*sin_a*sin_gd
      bee       = asin( sin_bII )
      cos_bII   = cos( bee )

      cos_ell   = cos_d*cos_a / cos_bII
      sin_ell   = ( cos_d*sin_a*cos_gd + sin_d*sin_gd ) / cos_bII
      ell       = atan2( sin_ell,cos_ell )

      ell_II    = ell/deg_rad + g_33
      bee_II    = bee/deg_rad

      if ( ell_II .lt. 0.d0 ) ell_II = ell_II + 360.d0

      return
      end

c     ========================================================
      subroutine xy_to_galactic ( ra, dec, x_offset, y_offset,
     +                            l_II_offset, b_II_offset)

c     Converts (X,Y) position offsets to Galactic offsets.
c     Can be used to convert X,Y motions to Galactic motions.

c     Takes (ra,dec) in decimal hours,degrees,
c     converts to Galactic (ell,bee) in degrees.

c     Then it adds X-offset (ie, RA*cos(Dec) and Y-offset
c     in arcseconds to (ra,dec) and again converts to Galactic.

c     Finally, differences the two Galactic coordinates to get
c     a shift in Galactic coordinates.


      implicit real*8 (a-h,o-z)

      real*8    l_II, l_II_new, l_II_offset

      pi 	= 4.0d0*atan(1.0d0)
      hr_rad	= pi/12.d0	! convert hours to radians
      deg_rad   = pi/180.d0	! convert degrees to radians
      hr_deg    = 15.d0         ! convert hours to degrees


c     Galactic plane-pole data for B1950...
c      g_alpha	= 282.25d0	! RA(?) of Galactic pole (deg)
c      g_dec	= 62.6d0	! Dec(?) of Galactic pole (deg)
c      g_33	= 33.0d0	! l(II) - l(I) (?) (deg)
c     J2000 Galactic plane-pole data...
      g_pole_ra = 12.d0 + 51.d0/60.d0 + 26.2817d0/3600.d0     ! hrs
      g_alpha   = (g_pole_ra + 6.d0) * hr_deg                 ! deg
                       ! = 282.859507

      g_pole_dec= 27.d0 + 07.d0/60.d0 + 42.013d0/3600.d0      ! deg
      g_dec     = 90.d0 - g_pole_dec                          ! deg
                       ! = 62.8716631
      g_33      = 33.d0 - 0.068d0 ! Longitude of ascending node of gal plane
                       ! = 32.932
                       ! Determined g_33 by requiring that the
                       ! B1950 defined G.C., when precessed to J2000,
                       ! still gives l=0.0 and b=0.0.  (works to 0.0001 deg)

c     Useful quantities...
      cos_a	= cos( ra*hr_rad - g_alpha*deg_rad )
      sin_a	= sin( ra*hr_rad - g_alpha*deg_rad )
      cos_d	= cos( dec*deg_rad )
      sin_d	= sin( dec*deg_rad )
      cos_gd	= cos( g_dec*deg_rad )
      sin_gd	= sin( g_dec*deg_rad )

c     Now calculate central galactic coordinates given ra,dec...
      sin_bII	= sin_d*cos_gd - cos_d*sin_a*sin_gd
      bee	= asin( sin_bII )
      cos_bII	= cos( bee )

      cos_ell	= cos_d*cos_a / cos_bII
      sin_ell	= ( cos_d*sin_a*cos_gd + sin_d*sin_gd ) / cos_bII
      ell	= atan2( sin_ell,cos_ell )

      l_II	= ell/deg_rad + g_33
      b_II	= bee/deg_rad

c      write (6,1000) ra,dec, l_II,b_II
c 1000 format('XY_TO_GALACTIC: ra,dec',2f10.4,'; ell,bee',2f10.2)

c     Now add in offsets (eg, in arcseconds)

      ra_new = ra + x_offset/3600.d0/15.d0/cos_d
      dec_new= dec+ y_offset/3600.d0

c     Useful quantities...
      cos_a	= cos( ra_new*hr_rad - g_alpha*deg_rad )
      sin_a	= sin( ra_new*hr_rad - g_alpha*deg_rad )
      cos_d	= cos( dec_new*deg_rad )
      sin_d	= sin( dec_new*deg_rad )

c     Now calculate new galactic coordinates given ra,dec...
      sin_bII	= sin_d*cos_gd - cos_d*sin_a*sin_gd
      bee	= asin( sin_bII )
      cos_bII	= cos( bee )

      cos_ell	= cos_d*cos_a / cos_bII
      sin_ell	= ( cos_d*sin_a*cos_gd + sin_d*sin_gd ) / cos_bII
      ell	= atan2( sin_ell,cos_ell )

      l_II_new= ell/deg_rad + g_33
      b_II_new= bee/deg_rad

      l_II_offset = (l_II_new - l_II)*3600.
      b_II_offset = (b_II_new - b_II)*3600.

      return
      end

c     ========================================================

      subroutine standard_vlsr_to_helio ( ra_rad, dec_rad, V_proj )

c     Enter (RA,Dec) in radians

c     Gives V_proj = V(Heliocentric) - V(LSR)  in km/s
c         ie, add V_proj to V(LSR) to get V(Heliocentric)

c     Uses Standard Solar Motion (old values)

      implicit real*8 (a-h,o-z)

      pi = 4.d0 * atan( 1.d0 )

c     LSR defined by removing peculiar Solar Motion of
c     20.0 km/s toward 18.0 hours, +30.0 degrees (RA,Dec)
c     (Defined in 1900.0 system, so precess to J2000)

      Vsun     = 20.d0                  ! km/s
c     RA_Vsun  = 18.d0                  ! hours      1900.0
c     Dec_Vsun = 30.d0                  ! degrees
      RA_Vsun  = 18.0640d0              ! hours      J2000.0
      Dec_Vsun = 30.0024d0              ! degrees

      cos_RA = cos( RA_Vsun * pi / 12.d0 )
      sin_RA = sin( RA_Vsun * pi / 12.d0 )

      cos_Dec= cos( Dec_Vsun* pi /180.d0 )
      sin_Dec= sin( Dec_Vsun* pi /180.d0 )

c     In equatorial Cartesian frame the Solar Motion is

      Xo = Vsun * cos_RA * cos_Dec       ! toward RA=0h in equat plane
      Yo = Vsun * sin_RA * cos_Dec       ! toward RA=6h in equat plane
      Zo = Vsun * sin_Dec                ! toward Dec=90 (North)

c     Projection of Solar Motion toward a source at RA,Dec

      cos_alpha = cos( ra_rad )
      sin_alpha = sin( ra_rad )

      cos_delta = cos( dec_rad )
      sin_delta = sin( dec_rad )

      V_proj = -Xo*cos_alpha*cos_delta - Yo*sin_alpha*cos_delta -
     +          Zo*sin_delta

      return
      end

c======================================================================

      subroutine galactic_name ( ell, bee, Gname )

      implicit real*8 (a-h,o-z)

      character*13 Gname
      character*1  bsign

c     Generate full galactic source name (truncating last digit per IAU
c     naming conventions!)

      i_ell       = ell
      decimal_ell = ell - i_ell
      i_d_ell     = 100.d0*decimal_ell

      bsign       = '+'
      if ( bee .lt. 0.d0 ) bsign = '-'
      abs_bee     = abs(bee)
      i_bee       = abs_bee
      decimal_bee = abs_bee - i_bee
      i_d_bee     = 100.d0*decimal_bee

      write (Gname,1000) i_ell,i_d_ell, bsign,i_bee,i_d_bee
 1000 format('G',i3.3,'.',i2.2,a1,i2.2,'.',i2.2)


      return
      end

c======================================================================

      subroutine encode_parnames ( num_peaks, parnames )

      character*16 parnames(77)

      parnames(1) = 'Constant        '
      parnames(2) = 'Slope           '

      index = 2
      do i = 1, num_peaks

         index = index + 1
         write (parnames(index),1000) i
 1000    format('Peak #',i1,'  (1/kpc)')

         index = index + 1
         write (parnames(index),1010) i
 1010    format('Distance #',i1,'(kpc)')

         index = index + 1
         write (parnames(index),1020) i
 1020    format('FWHM #    ',i1,'(kpc)')

      enddo

      return
      end

c======================================================================

      subroutine fit_multiple_gaussians (lu_print,
     +               num_data, num_params, parnames,
     +               params, paramids,
     +               vel, data, arm_max_bin,
     +               peak_dist, peak_dunc, peak_int )

      implicit real*8 (a-h,o-z)

      real*8        vel(1001), data(1001)
      real*8        model(1001), resids(1001), res_err(1001)
      real*8        partls(1001,77)

      character*12  arm_max_bin(1001), arm_indicated

      real*8        params(77), new_params(77), param_sigmas(77)
      real*8        params_original(77)
      character*16  parnames(77)
      integer       paramids(77)

      real*8        peak_dist(25), peak_dunc(25), peak_int(25)

      logical       print_cor, print_resids, converged


      max_num_data = 1001

      itermax = 50
      idebug  =  0

      print_resids = .false.
      converged    = .false.

c     ------------------------------------------------------------
c                       Setup for Fitting

c     Save original parameter values in case fitting fails...
      do n_p = 1, num_params
         params_original(n_p) = params(n_p)
      enddo

c     Get initial data residuals...
      call calc_residuals ( num_params, params,
     +                      num_data, vel, data,
     +                      model, resids )

c     Count number of parameters solved-for
      num_solved_for = 0
      do n = 1, num_params
         if (paramids(n).gt.0) num_solved_for=num_solved_for+1
      enddo

c     Assign errors for data; will be adjusted later by sigma_pdf
      do n = 1, num_data
         res_err(n) = 1.d0/num_data
      enddo

C     Calculate initial "sigma-squared"...
      sigsq     = 9999.d0
      call calc_sigsq (num_data, num_solved_for,
     +                 resids, res_err,
     +                 sigsq_old, sigsq, sigma_pdf)


      num_deg_freedom = num_data - num_solved_for
      if ( idebug .gt. 0 )
     +   write (6,1300) sigma_pdf, num_deg_freedom
1300     format(/' Pre-fit sigma_pdf =',f10.3,' for',i4,
     +           ' degrees of freedom.')

C     ------------------------------------------------------------
C                      Iterate Least-Squares Fitting

      print_cor = .false.
      iterpass  = 1
      itest2    = 0
      do while ( iterpass.lt.itermax .and. .not.converged .and.
     +           itest2.eq.0 )

         iterpass = iterpass + 1

         call calc_partials ( num_params, params, paramids,
     +                        num_data, vel,
     +                        partls )

         call scale_errors( num_data, sigma_pdf, res_err )

         call least_squares_fit( idebug,print_cor,max_num_data,
     +                 num_params,num_data,
     +                 parnames,params,paramids,resids,
     +                 partls,res_err,new_params,param_sigmas,
     +                 itest2 )

c        Check that matrix inversion worked...
         if ( itest2 .eq. 0 ) then
            call update_params ( num_params, params,
     +                           new_params, converged )

            call calc_residuals ( num_params, params,
     +                            num_data, vel, data,
     +                            model, resids )

            call calc_sigsq ( num_data, num_solved_for,
     +                        resids, res_err,
     +                        sigsq_old, sigsq, sigma_pdf)
         endif

      enddo

c     If matrix inversion failed, restore original parameter values...
      if ( itest2 .gt. 0 ) then
         do n_p = 1, num_params
            params(n_p) = params_original(n_p)
c           But, don't allow width param >0.5 kpc as this can dominate probability
            if ( mod(n_p-2,3).eq.0 .and. params(n_p).gt.0.5d0 ) then
               params(n_p) = 0.5d0
            endif
         enddo
      endif

c     Integrals of probabilities for each Gaussian...
      pi = 4.d0*atan(1.d0)
      factor = sqrt( pi / (4.d0*log(2.d0) ) )

      num_gaussians = (num_params - 2)/3
      peak_int_sum = 0.d0
      do n_g = 1, num_gaussians

         n_a = 3 + (n_g-1)*3              ! amp index
         n_c = n_a + 1                    ! center (distance) index
         n_w = n_a + 2                    ! width index

         peak_int(n_g) = factor * params(n_a) * params(n_w)
         peak_int_sum = peak_int_sum + peak_int(n_g)

      enddo

c     Normalize integrated probabilities...
cc      write (lu_print,2540)
 2540 format(/28x,'Dist(kpc)    +/-   Probability   Arm')
      do n_g = 1, num_gaussians
         n_a = 3 + (n_g-1)*3              ! amp index
         n_c = n_a + 1                    ! center (distance) index
         n_w = n_a + 2                    ! width index

         peak_dist(n_g)= params(n_c)      ! distance estimate (kpc)
         peak_dunc(n_g)= params(n_w)/2.3548d0  ! converted FWHM to 1-sigma
         peak_int(n_g) = peak_int(n_g) /  peak_int_sum


         distance = peak_dist(n_g)
         num_bins = max_num_data
         call closest_arm ( num_bins, vel, arm_max_bin,
     +                      distance,
     +                      arm_indicated )


cc         write (lu_print,2550) n_g, peak_dist(n_g), peak_dunc(n_g),
cc     +                         peak_int(n_g), arm_indicated
 2550    format(' Probability component',i2,': ',3f10.4,5x,a3)
      enddo


C     ------------------------------------------------------------
C                        Other outputs
C     Printout residuals?
      if ( print_resids ) then
cc         write (lu_print,3000)
3000     format (/'         Data       Model    Residual      ',
     +            ' Error     Res/Err'/)
         do n = 1, num_data

            wtres = resids(n) / res_err(n)
cc            write (lu_print,3100) data(n), model(n), resids(n),
cc     +                                     res_err(n), wtres
3100        format (1x,5f12.4)
         enddo
      endif

      return
      end

c======================================================================

      subroutine scale_errors( num_data, sigma_pdf, res_err )

      implicit real*8 (a-h,o-z)

      real*8   res_err(1001)

      do n = 1, num_data

         res_err(n) = res_err(n) * sigma_pdf

      enddo

      return
      end

c======================================================================

      subroutine smooth ( smooth_kpc, num_bins, dist_bins,
     +                    prob_dist )

c     Smoothes "prob_dist" and returns smoothed values in same array
c     Will use boxcar smoothing with width = smooth_kpc

      implicit real*8 (a-h,o-z)

      real*8    dist_bins(1001), prob_dist(1001), temp(1001)

      do i = 1, num_bins
         temp(i) = prob_dist(i)
      enddo

      bin_width = dist_bins(2) - dist_bins(1)               ! kpc
      n_smooth = smooth_kpc / bin_width + 0.5d0

c     Want odd number of channels to center the boxcar on the bin
      if ( mod(n_smooth,2) .eq. 0 ) n_smooth = n_smooth + 1
      n_half   = n_smooth / 2

c     Don't smooth edge channels (need to keep them for some nearby sources)
      n0 = n_half + 1
      n1 = num_bins - n_half
      do n = 1, num_bins

         if ( n.ge.n0 .and. n.le.n1 ) then

               n_start = n - n_half - 1
               sum = 0.d0
               do nn = 1, n_smooth
                  nnn = n_start + nn
                  sum = sum + temp(nnn)
               enddo

               prob_dist(n) = sum / n_smooth

            else

               prob_dist(n) = temp(n)

         endif

      enddo

      return
      end

c======================================================================

      subroutine find_probability_peaks(lu_print, max_num_peaks,
     +                bin_size, num_bins, dist_bins, prob_dist,
     +                num_peaks, peak_prob,
     +                peaks, peaks_width)

c     Goes through a PDF, finds up to max_num_peaks=10 (local) peaks, and
c     provides the peak bin value and +/-34% integrated probability bins from
c     the peak.   Peaks must be greater than a value currently set to
c     ~1/(num_bins*bin_size) or 4%

c     Of course for overlapping PDF this is problematic and will
c     be handled by fitting multiple Gaussians to the full PDF later.

      implicit real*8 (a-h,o-z)

      real*8    dist_bins(1001), prob_dist(1001)
      real*8    peaks(25), peaks_low(25), peaks_high(25)
      real*8    peak_prob(25), peaks_width(25)

      integer   npeaks(25)

      pmin   = 1.01d0/(num_bins*bin_size)   ! set just above pure flat PDF
      pratio = 0.0d0                        ! max drop off next to peak

      onesig_down = exp(-0.5d0)

c     Zero outputs
      num_peaks = 0
      do i_p = 1, max_num_peaks
         npeaks(i_p)     = 0
         peak_prob(i_p)  = 0.d0
         peaks(i_p)      = 0.d0
         peaks_low(i_p)  = 0.d0
         peaks_high(i_p) = 0.d0
      enddo

c     --------------------------------------------------
c     First scan for the locations of multiple peaks...
      n_p = 0
      n0  = 2
      n1 = num_bins - 1
      do n = n0, n1
         pcen  = prob_dist(n)
         plow  = prob_dist(n-1)
         phigh = prob_dist(n+1)
         if ( pcen .ge. plow        .and.
     +        pcen .ge. phigh       .and.
     +        pcen .gt. pmin        .and.
     +        plow .gt. pratio*pcen .and.
     +        phigh.gt. pratio*pcen      ) then
c           Have a significant local peak...
            if ( n_p .lt. max_num_peaks ) then
               n_p = n_p + 1
               peak_prob(n_p) = prob_dist(n)
               peaks(n_p)     = dist_bins(n)
               npeaks(n_p)    = n
            endif
         endif
      enddo
      num_peaks = n_p

c     -------------------------------------------------
c     For each peak, estimate 1-sigma widths

      do n_p = 1, num_peaks
c        Start at peak and increase bin index until
c        reaching 0.61 of peak probability (approx 1 sigma)
         onesig_amp = onesig_down * peak_prob(n_p)
         n_b  = npeaks(n_p)
         n_onesig = n_b
         do while ( prob_dist(n_b).gt.onesig_amp .and.
     +              n_b.le.num_bins )
            n_b = n_b + 1
            n_onesig = n_b
         enddo
         peaks_low(n_p) = dist_bins(n_onesig)
      enddo

      do n_p = 1, num_peaks
c        Start at peak and decrease bin index until
c        reaching 0.61 of peak probability (approx 1 sigma)
         onesig_amp = onesig_down * peak_prob(n_p)
         n_b  = npeaks(n_p)
         n_onesig = n_b
         do while ( prob_dist(n_b).gt.onesig_amp .and.
     +              n_b.ge.1 )
            n_b = n_b - 1
            n_onesig = n_b
         enddo
         peaks_high(n_p) = dist_bins(n_onesig)
      enddo

c     Now estimate the 1sigma width as the minimum
c     of the low and high estimates (because blended lines
c     can lead to erratic estimates)
      do n_p = 1, num_peaks
         sig_low = abs( peaks(n_p) - peaks_low(n_p) )
         sig_high= abs( peaks_high(n_p) - peaks(n_p))
         peaks_width(n_p) = min( sig_low, sig_high )
      enddo

      return
      end

c==================================================================

      subroutine edit_peaks ( num_peaks,
     +                        peak_prob, peaks, peaks_width,
     +                        use_peak )

c     Discards unimportant peaks in Prob(dist) distribution
c     by setting "use_peak" logical flag array.

c     Two parameters hardwired here:
c        minimum accepted peak is >3% of maximum peak
c        minimum separation between peaks is 0.5 kpc

c     Also will combine peaks if they overlap significantly.

      implicit real*8 (a-h,o-z)

      real*8        peaks(25)
      real*8        peak_prob(25), peaks_width(25)

      logical       use_peak(25)

      if ( num_peaks .eq. 1 ) use_peak(1) = .true.

c     Skip if only 1 peak found...
      if ( num_peaks .ge. 2 ) then

c        Find maximum peak probability; needed for setting "peak_min"
         peak_max = 0.d0
         do n = 1, num_peaks
            if ( peak_prob(n).gt.peak_max ) peak_max = peak_prob(n)
         enddo

c        Set minimum percentage of maximum peak to bother fitting...
         peak_min = 0.01d0*peak_max        ! minimum peak to do L-S fitting

c        Set minimum absolute separation
         separation_min = 0.25d0           ! minimum spacing in between peaks (kpc)

c        Set minimum separation to 25% of width
         sep_width_fract = 0.25d0

c        Check that peak prob and separation are greater than minima...
         do n = 1, num_peaks
            use_peak(n) = .false.
            if ( peak_prob(n) .gt. peak_min ) then
               if ( n .gt. 1 ) then
                     width_max = max(peaks_width(n-1),peaks_width(n))
                     sep_width = sep_width_fract * width_max
                     del = abs( peaks(n) - peaks(n-1) )
                     if ( del.gt.separation_min .and.
     +                    del.gt.sep_width          ) use_peak(n)=.true.
                  else
                     use_peak(n)=.true.
               endif
            endif
         enddo

c        Remove peaks not used...
         n_used = 0
         do n = 1, num_peaks
            if ( use_peak(n) ) then
c              Re-order used ones...
               n_used = n_used + 1
               peaks(n_used) = peaks(n)
               peak_prob(n_used) = peak_prob(n)
               peaks_width(n_used) = peaks_width(n)
c              Reset flags...
               use_peak(n)      = .false.
               use_peak(n_used) = .true.
            endif
         enddo

         num_peaks = n_used

      endif        ! num_peaks >0 check

      return
      end

c==================================================================

	subroutine calc_residuals ( num_params, params,
     +                              num_data, vel, data,
     +                              model, resids )

C	Calculates models and resduals for entire data set...

	implicit real*8 ( a-h,o-z )

	real*8	params(77)

	real*8	vel(1001), data(1001), model(1001), resids(1001)

        num_peaks = ( num_params - 2 ) / 3

        do n = 1, num_data

           data_vel = vel(n)
           call calc_model ( num_peaks, params, data_vel,
     +                       flux_density )

           model(n)  = flux_density
           resids(n) = data(n) - model(n)

	enddo

	return
	end
c=========================================================

      subroutine calc_model ( num_peaks, params, data_vel,
     +                        flux )

c     This routine was written for modeling spectra, hence
c     the unusual parameter names.

c     Calculates a model "flux density" given a "velocity"
c     and the model parameters for "num_peaks" Gaussians

	implicit real*8 ( a-h, o-z )

        real*8    params(77)

c       Baseline contribution...
c       (NB: baseline pivots about v=0)
        flux = params(1) + params(2) * data_vel

        index = 2
        do n = 1, num_peaks

           index     = index + 1
           peak_flux = params(index)

           index     = index + 1
           center    = params(index)

           index     = index + 1
           fwhm      = params(index)

           if ( peak_flux .gt. 0.d0 ) then

              call gaussian_profile ( peak_flux, center,
     +                                fwhm, data_vel,
     +                                g_flux )

              flux = flux + g_flux

           endif

        enddo

	return
	end
c=========================================================

      subroutine gaussian_profile(peak_flux,center,fwhm,data_vel,
     +                            flux)

c     Calculates "flux density" (flux) for a single Gaussian line
c     profile at a specified velocity (data_vel)

c     Note: flux set to 0 if exponential term < 10^-9

	implicit real*8 (a-h, o-z)

        dv_sq    = (data_vel - center)**2
        two_sigma_sq = 2.d0 * fwhm**2 / ( 8.d0 * log(2.d0) )
        arg  = dv_sq / two_sigma_sq
        if ( arg .lt. 20.7d0 ) then
              flux = peak_flux * exp( -arg )
           else
              flux = 0.d0
        endif

	return
	end
c=========================================================

      subroutine calc_partials ( num_params, params, paramids,
     +                 num_data, vel,
     +                 partls )

c     Numerically calculates partial derivatives of the model
c     for all parameters

      implicit real*8 ( a-h,o-z )

      real*8		params(77), params_wiggled(77)
      real*8		partls(1001,77), vel(1001)

      integer		paramids(77)

      num_peaks = ( num_params - 2 ) / 3

c     Set secondary parameter array and zero partials array...
      do i_p = 1, num_params

	params_wiggled(i_p) = params(i_p)

	do i_d = 1, num_data
	   partls(i_d,i_p) = 0.d0
	enddo

      enddo

c     Calculate numerical partials...
      do i_d = 1, num_data

         data_vel = vel(i_d)

         call calc_model ( num_peaks, params, data_vel, flux )

         do i_p = 1, num_params

            if ( paramids(i_p) .eq. 1 ) then

c              Change parameters by 1 part in 10^4; but never
c	       by less than an absolute value of 10^-8...

               del_param = abs( params(i_p) ) * 1.d-04
               if ( del_param .lt. 1.d-08 ) del_param = 1.d-08

               params_wiggled(i_p) = params(i_p) + del_param

               call calc_model(num_peaks,params_wiggled,data_vel,
     +                         flux_wiggled )

               partls(i_d,i_p)  = (flux_wiggled - flux)/del_param

               params_wiggled(i_p) = params(i_p)

            endif

         enddo

      enddo

      return
      end

c=========================================================

      subroutine update_params (num_params, params,
     +                          new_params, converged )

c     Takes parameter values before (params) and after (new_params)
c     least squares fit; adjusts the changes by the "gain";
c     and replaces the adjusted new parameters in the "params" array.

      implicit real*8 (a-h, o-z)

      real*8	params(77), new_params(77)

      logical   converged

      gain = 0.1d0                      ! hardwired gain parameter
      change_max  = 0.2d0               ! max fractional change per iteration
      convergence = 0.01d0              ! convergence criterion (max changes)

      converged = .true.
      do j = 1, num_params

         delta_param = new_params(j) - params(j)
         fractional_change = abs( delta_param/params(j) )
         if ( fractional_change .gt. change_max) then
            delta_param = delta_param*change_max/fractional_change
         endif

         trial_param = params(j) + gain * delta_param
         params(j) = trial_param

         if ( fractional_change .gt. convergence ) then
            converged=.false.
         endif

      enddo

      return
      end

c=========================================================

      subroutine calc_sigsq (num_data, num_solved_for,
     +			resids, res_err,
     +			sigsq_old, sigsq, sigma_pdf)

C	Calculate new "sigsq"...

	implicit real*8	(a-h,o-z)

	real*8		resids(1001), res_err(1001)

	sigsq_old = sigsq
	sigsq     = 0.d0

	do i_d = 1, num_data
           sigsq = sigsq + ( resids(i_d)/res_err(i_d) )**2
	enddo

	num_deg_freedom = num_data - num_solved_for

	if ( num_deg_freedom .gt. 1 ) then
              sigsq_pdf = sigsq / float(num_deg_freedom)
           else
              sigsq_pdf = 0.d0
        endif

	sigma_pdf = sqrt( sigsq_pdf )

	return
	end


c===================================================================

      subroutine galaxy_parameters_Univ( lu_control,
     +           file_name, lu_print,
     +           Ro, a1, a2, a3, Uo, Vo, Wo, Us, Vs, Ws,
     +           glong_min, glong_max )

c     Opens file and reads parameters needed for revised kinematic distance
c     using the Persic "Universal" rotation curve formulation

      implicit real*8 (a-h,o-z)

      character*48    file_name

      call open_ascii_file ( lu_control, file_name, lu_print )

c     Galactic structure parameters
      read (lu_control,*) Ro          ! Ro [kpc]
      read (lu_control,*) a1          ! To [km/s]; not used in 2 parameter version
      read (lu_control,*) a2, a3      ! Univ rotation curve parameters near 0.9, 1.5

c     Solar Motion parameters
      read (lu_control,*) Uo       ! Uo toward G.C. [km/s]
      read (lu_control,*) Vo       ! Vo toward Gal. rotation [km/s]
      read (lu_control,*) Wo       ! Wo toward North Gal.Pole [km/s]

c     Source peculiar motion in its local Galactocentric frame
c     (G.C.is as viewed by source, not Sun)
      read (lu_control,*) Us       ! Us toward G.C. [km/s]
      read (lu_control,*) Vs       ! Vs toward Gal. rotation [km/s]
      read (lu_control,*) Ws       ! Ws toward North Gal.Pole [km/s]

c     Minimum/maximum Galactic longitude to consider (where models are valid)
      read (lu_control,*) glong_min
      read (lu_control,*) glong_max

      close(unit=lu_control)

cc      write (lu_print,1000)
 1000 format(' Galaxy model:')
cc      write (lu_print,1050) Ro, a1, a2, a3
 1050 format('  Ro, a1, a2, a3:                   ',f10.2,3f10.3)
cc      write (lu_print,1100) Uo, Vo, Wo
 1100 format('  Solar Motion (Uo,Vo,Wo):          ',3f10.2)
cc      write (lu_print,1200) Us, Vs, Ws
 1200 format('  Source Peculiar Motion (Us,Vs,Ws):',3f10.2)

      return
      end

c==============================================================================

      subroutine get_arm_segments ( lu_data, lu_print, max_num_lbvs,
     +                              arm_name, seg_ell_min, seg_ell_max,
     +                              arm_ell, arm_bee, arm_vel,
     +                              arm_Rgc, arm_beta, arm_dist,
     +                              iarm_entries, num_arms )

      implicit real*8 (a-h,o-z)

c     Arm information arrays
      real*8        arm_ell(29,300), arm_bee(29,300), arm_vel(29,300)
      real*8        arm_Rgc(29,300), arm_beta(29,300), arm_dist(29,300)
      real*8        seg_ell_min(29,2), seg_ell_max(29,2)

      integer       iarm_entries(29)

      character*48  arm_file
      character*12  a_name, arm_name(29)
      character*1   c1

c     Spiral arm segment names; used to point to (lbvRBD) files

      num_arms    = 23            ! NB: current dimension limit is 29

      arm_name(1) = 'Out'         ! Outer arm

      arm_name(2) = 'Per'         ! Perseus arm

      arm_name(3) = 'Loc'         ! Local arm

      arm_name(4) = 'CrN'         ! Carina near side of tangent point
      arm_name(5) = 'CrF'         ! Carina far side
      arm_name(6) = 'SgN'         ! Sagittarius near
      arm_name(7) = 'SgF'         ! Sagittarius far

      arm_name(8) = 'CtN'         ! Centaurus (Crux) near side of tangent point
      arm_name(9) = 'CtF'         ! Centaurus (Crux) far side
      arm_name(10)= 'ScN'         ! Scutum near
      arm_name(11)= 'ScF'         ! Scutum far
      arm_name(12)= 'OSC'         ! "Outer Scutum-Centaurus" arm

      arm_name(13)= 'N1N'         ! Norma arm Q1 near
      arm_name(14)= 'N1F'         ! Norma arm Q1 far
      arm_name(15)= 'N4N'         ! Norma arm Q4 near
      arm_name(16)= 'N4F'         ! Norma arm Q4 far

      arm_name(17)= '3kN'         ! 3-kpc near side of GC
      arm_name(18)= '3kF'         ! 3-kpc far side

      arm_name(19)= 'AqS'         ! Aquila Spur

      arm_name(20)= 'CnN'         ! Connecting arm near side (of GC)
      arm_name(21)= 'CnX'         ! Connecting arm near side extension

      arm_name(22)= 'AqR'         ! Aquila Rift

      arm_name(23)= 'LoS'         ! Local spur


c     Open files (for each arm) used for plotting and document with comment line
      lu = 30
      a_name = 'Unk'
      write (30,1010) lu, a_name
      do n_a = 1, num_arms
         lu = 30 + n_a
         write (lu,1010) lu, arm_name(n_a)
 1010    format('! File fort.',i2,' for spiral arm segment "',a3,'"',
     +         /'! Long.   Lat.   Vlsr  +/-    Dist.  +/-  ',
     +          'Integrated  Arm',
     +         /'! (deg)  (deg)  (km/s)        (kpc)       ',
     +          'Probability')

      enddo

c    ----------------------------------------------------------
c               Read in spiral arm (l,b,v,R,beta,d) information
cc      write (lu_print,1014)
 1014 format(/' Getting Spiral arm segment data:',/' ')
      do n_a = 1, num_arms

         write (arm_file,1015) arm_name(n_a)
 1015    format(a3,'_lbvRBD.2019')

         lu_noprint = -1
         call open_ascii_file(lu_data, arm_file, lu_noprint)

c        Will keep track of max and min segment longitudes,
c        which is complicated by crossing from 359 to 001 degrees.
c        Need to allow for 2 ranges.
         seg_ell_min(n_a,1) = 9.d9
         seg_ell_max(n_a,1) =-9.d9
         seg_ell_min(n_a,2) = 9.d9
         seg_ell_max(n_a,2) =-9.d9
         i_r = 1

         ieof = 0
         n_v  = 0
         do while ( ieof.ge.0 .and. n_v.lt.max_num_lbvs )

            read (lu_data,1016,iostat=ieof) c1
 1016       format(a1)

            if ( c1.ne.'!' .and. ieof.ge.0 ) then

               backspace (unit=lu_data)
               read (lu_data,*) a_ell,a_bee,a_vel,a_R,a_beta,a_d

               n_v = n_v + 1
               arm_ell(n_a,n_v) = a_ell        ! longitude [deg]
               arm_bee(n_a,n_v) = a_bee        ! latitude  [deg]
               arm_vel(n_a,n_v) = a_vel        ! Vlsr      [km/s]
               arm_Rgc(n_a,n_v) = a_R          ! Galocentric Radius [kpc]
               arm_beta(n_a,n_v)= a_beta       !      "     Azimuth [deg]
               arm_dist(n_a,n_v)= a_d          ! Distance from Sun  [kpc]

               if (n_v.gt.1 .and. abs(a_ell-a_ell_old).gt.180.d0) i_r=2
               if (a_ell.lt.seg_ell_min(n_a,i_r))
     +                                     seg_ell_min(n_a,i_r) = a_ell
               if (a_ell.gt.seg_ell_max(n_a,i_r))
     +                                     seg_ell_max(n_a,i_r) = a_ell
               a_ell_old = a_ell

            endif

         enddo

         iarm_entries(n_a) = n_v

cc         write (lu_print,1017) arm_file, n_v
 1017    format(/'  Read file ',a20,' containing',i4,
     +           ' (l,b,v,R,beta,D) values.')
cc         write (lu_print,1018) seg_ell_min(n_a,1), seg_ell_max(n_a,1)
 1018    format(7x,'Spanning longitude range:',f7.1,' to',f7.1,' deg.')
         if ( seg_ell_max(n_a,2) .gt. -900.d0 ) then
cc           write (lu_print,1019) seg_ell_min(n_a,2), seg_ell_max(n_a,2)
 1019      format(12x,'and longitude range:',f7.1,' to',f7.1,' deg.')
         endif

         close (unit=lu_data)

         call enforce_ccw_order ( n_a, n_v, arm_ell, arm_bee, arm_vel,
     +                            arm_Rgc, arm_beta, arm_dist )

      enddo

      return
      end

c====================================================================
      subroutine enforce_ccw_order(n_a, n_v, arm_ell,arm_bee,arm_vel,
     +                             arm_Rgc, arm_beta, arm_dist)

c     Check that spiral arm azimuths are in counter-clockwise order;
c     if not, then reverse array entries.

      implicit real*8 (a-h,o-z)

c     Arm information arrays
      real*8  arm_ell(29,300), arm_bee(29,300), arm_vel(29,300)
      real*8  arm_Rgc(29,300),arm_beta(29,300),arm_dist(29,300)

c     Temporary arrays
      real*8  ell(29,300), bee(29,300), vel(29,300)
      real*8  Rgc(29,300),beta(29,300),dist(29,300)

      logical reverse1, reverse2

c     Check both beginning and ending longitudes for desired CCW order
      reverse1 = .false.
      reverse2 = .false.
      if ( arm_beta(n_a,2).gt.arm_beta(n_a,1) )       reverse1 = .true.
      if ( arm_beta(n_a,n_v).gt.arm_beta(n_a,n_v-1) ) reverse2 = .true.

      if ( reverse1 .neqv. reverse2 ) then
        print*,'STOP: subroutine enforce_ccw_order found mixed orders'
        STOP
      endif

      if ( reverse1 ) then

cc         write (6,1000)
 1000    format('  Reversing array order to enforce CCW azimuths')

c        Copy to temporary arrays in reversed order...
         do n = 1, n_v
            n_rev = n_v + 1 - n
            ell(n_a,n_rev) = arm_ell(n_a,n)
            bee(n_a,n_rev) = arm_bee(n_a,n)
            vel(n_a,n_rev) = arm_vel(n_a,n)
            Rgc(n_a,n_rev) = arm_Rgc(n_a,n)
            beta(n_a,n_rev)= arm_beta(n_a,n)
            dist(n_a,n_rev)= arm_dist(n_a,n)
         enddo

c        Transfer back to original arrays
         do n = 1, n_v
            arm_ell(n_a,n) = ell(n_a,n)
            arm_bee(n_a,n) = bee(n_a,n)
            arm_vel(n_a,n) = vel(n_a,n)
            arm_Rgc(n_a,n) = Rgc(n_a,n)
            arm_beta(n_a,n)= beta(n_a,n)
            arm_dist(n_a,n)= dist(n_a,n)
         enddo

      endif

      return
      end

c======================================================================

      subroutine get_parallaxes ( lu_control, lu_print,
     +                            max_num_parallaxes, num_parallaxes,
     +                            ref_src, ref_arm, ref_ell,
     +                            ref_bee, ref_vlsr, ref_vunc,
     +                            ref_par, ref_punc )

c     Reads in trigonometric parallax data...

      implicit real*8 (a-h,o-z)

      character*48  data_file
      character*14  src, ref_src(1000)
      character*12  arm, ref_arm(1000)
      character*1   c1

      real*8        ref_ell(1000), ref_bee(1000)
      real*8        ref_vlsr(1000), ref_vunc(1000)
      real*8        ref_par(1000), ref_punc(1000)

c     Read all parallax data and store essential info in arrays...
      data_file = 'parallax_data.inp'
      call open_ascii_file ( lu_control, data_file, lu_print )

      n_p  = 0
      ieof = 0
      do while ( ieof .ge. 0 )

         read (lu_control,1000,iostat=ieof) c1
 1000    format(a1)

         if ( ieof.ge.0 .and. c1.ne.'!' .and.
     +        n_p.lt.max_num_parallaxes       ) then

c           Input data for a source ...
            backspace (unit=lu_control)
            read (lu_control,*) src, arm, ell, bee,
     +                          v_lsr, err_v, par, err_par

c           Only use if fractional parallax uncertainty <20%
            fract_par_unc = err_par / par
            if ( fract_par_unc .lt. 0.20d0 ) then

               n_p = n_p + 1

               ref_src(n_p) = src
               ref_arm(n_p) = arm
               ref_ell(n_p) = ell
               ref_bee(n_p) = bee
               ref_vlsr(n_p)= v_lsr
               ref_vunc(n_p)= err_v
               ref_par(n_p) = par
               ref_punc(n_p)= err_par

            endif

         endif

      enddo

      close (unit=lu_control)

      if ( n_p .ge. max_num_parallaxes ) then
         write (6,2000) max_num_parallaxes
 2000    format(/' **********************************************',
     +          /' parallax_data.inp file contains more than',i5,
     +           ' entries (beyond dimension limit).  STOP')
         stop
      endif

      num_parallaxes = n_p
cc      write (lu_print,1100) num_parallaxes
 1100 format('  Read in',i4,' sources with measured parallaxes')

      return
      end

c======================================================================

      subroutine gauss_prob ( c_d, c_s, d, prob )

c     Calculates a normalized Gaussian probability (prob) density for
c     value d, offset from the center c_d, given a sigma c_s

c     Note: returns zero prob if exponential term is <10^-9

      implicit real*8 (a-h,o-z)

      pi    = 4.d0*atan(1.d0)
      twopi = 2.d0 * pi

      prefactor = sqrt(1.d0/twopi) / abs(c_s)

      arg = (d - c_d)**2 / (2.d0*c_s**2)

      if ( arg .lt. 20.7d0 ) then
            prob = prefactor * exp(-arg)
         else
c           NB: exp(-20.7) = 10^-9; if this term is <10^-9
c           call it zero to avoid underflows later
            prob = 0.d0
      endif

      return
      end

c======================================================================

      subroutine asym_gauss_prob ( c_d, c_s, c_s_avg, d, prob )

c     Calculates an approximately normalized Gaussian probability (prob)
c     density for value d, offset from the center c_d, given a sigma c_s

c     This version "smoothes" the top of the probability when there are
c     asymmetric errors by using the average uncertainty in the prefactor.

c     Note: returns zero prob if exponential term is <10^-9

      implicit real*8 (a-h,o-z)

      pi    = 4.d0*atan(1.d0)
      twopi = 2.d0 * pi

      prefactor = sqrt(1.d0/twopi) / abs(c_s_avg)

      arg = (d - c_d)**2 / (2.d0*c_s**2)

      if ( arg .lt. 20.7d0 ) then
            prob = prefactor * exp(-arg)
         else
c           NB: exp(-20.7) = 10^-9; if this term is <10^-9
c           call it zero to avoid underflows later
            prob = 0.d0
      endif

      return
      end

c======================================================================

      subroutine asym_conservative_prob( c_d, c_s, c_s_avg, d, prob )

c     Calculates a non-Gaussian ("conservative or error tolerant") probability
c     (prob) density for values d, offset from the center c_d, and sigma c_s

c     This version "smoothes" the top of the probability when there are
c     asymmetric errors by using the average uncertainty in the prefactor.

c     Note: exponential term is cutoff so prob=0 for values less than 10^-9

      implicit real*8 (a-h,o-z)

      pi    = 4.d0*atan(1.d0)
      twopi = 2.d0 * pi

      prefactor = sqrt(1.d0/twopi) / abs(c_s_avg)

      resid = (d - c_d) / c_s
      rsq   = resid**2

      exp_term = 0.d0
      if ( rsq .lt. 41.4d0 ) exp_term = exp(-rsq/2.d0)

      if ( rsq .lt. 1.d09 ) then
            prob = prefactor * ( 1.d0 - exp_term ) / rsq
         else
            prob = 0.d0
      endif

      return
      end

c=========================================================

      subroutine condition_pdf(num_bins, bin_size, P_max, pdf)

c     "Conditions" a "pdf" with "num_bins" elements spaced by "bin_size"
c     A flat background can be added, controlled by "P_max" {0 -> 1}:
c          P_max = 1.0 means add no background
c                = 0.5 means add a background with 50% of the area
c                = 0.0 means make a totally flat PDF (tossing the input values)
c          NB, if Sum(pdf)=0; a flat background will automatically be used
c     A normalized pdf is replaced, such that
c          Sum{ pdf(n)*bin_size } = 1


      implicit real*8 (a-h,o-z)

      real*8   pdf(1001)

      if ( P_max.lt.0.d0 .or. P_max.gt.1.d0 ) then
         write(6,1000) P_max
 1000    format(/'Subroutine condition_pdf: invalid P_max=',f10.3)
         stop
      endif

c     Calculate sum of PDF
      pdf_sum = 0.d0
      do n = 1, num_bins
         pdf_sum = pdf_sum + pdf(n)
      enddo
c     Multiply by bin_size
      pdf_sum = pdf_sum * bin_size

      scale = 1.d0
      if ( pdf_sum .gt. 1.d-99 ) scale = P_max / pdf_sum

      background = (1.d0 - P_max) / (bin_size*num_bins)

c     Normalize with background included...
      do n =  1, num_bins
         if ( pdf_sum .gt. 1.d-99 ) then
               pdf(n) = scale*pdf(n) + background
            else
               pdf(n) = 1.d0 / (bin_size*num_bins)
         endif
      enddo

      return
      end

c       =========================================================
        subroutine radians_to_hhmmss ( ra_rad, ra_hhmmss)

c       Input :  ra_rad       in radians
c       Output:  ra_hhmmss    as hhmmss.sssss

        implicit real*8 (a-h, o-z)

        pi        = 4.0d0 * atan(1.0d0)
        rad_to_hr = 12.d0/ pi

        ra_hr = ra_rad * rad_to_hr

        ihr  = ra_hr
        imin = (ra_hr - ihr*1.d0) * 60.d0
        sec  = (ra_hr - ihr*1.0d0 - imin/60.0d0) * 3600.d0

        ra_hhmmss  = ihr*10000.d0 + imin*100.d0 + sec

        return
        end


c       =========================================================
        subroutine radians_to_ddmmss ( dec_rad, dec_ddmmss)

c       Input :  dec_rad       in radians
c       Output:  dec_ddmmss    as hhmmss.sssss

        implicit real*8 (a-h, o-z)

        pi         = 4.0d0 * atan(1.0d0)
        rad_to_deg = 180.d0/ pi

        if ( dec_rad .ge. 0.d0 ) then
             dec_sign =  1.d0
          else
             dec_sign = -1.d0
        endif

        dec_deg = abs( dec_rad ) * rad_to_deg

        ideg = dec_deg
        imin = (dec_deg - ideg*1.d0) * 60.d0
        sec  = (dec_deg - ideg*1.0d0 - imin/60.0d0) * 3600.d0

        dec_ddmmss  = ideg*10000.d0 + imin*100.d0 + sec
        dec_ddmmss  = dec_ddmmss * dec_sign

        return
        end

c=========================================================

      SUBROUTINE LEAST_SQUARES_FIT ( IDEBUG,PRINT_COR,MAXDAT,NUMXES,
     +				N,VNAME,X,ID,S,C,E,XHAT2,ERROR,itest2 )

c     The following code was written before you were born,
c     not using structured code!

C     SUBROUTINE FOR LINEAR LEAST SQUARE FITTING...
C     IT REQUIRES INPUT OF THE FOLLOWING:
C        IDEBUG    = DEBUGGING PARAMETER
C                    0 = SUPRESS PRINT OUT BEYOND SOLUTION AND CORRELATIONS
C                    1 = PRINT OUT ALL DATA AND MATRICES
C	 PRINT_COR = LOGICAL FLAG: PRINT CORRELATION MATRIX
C	 MAXDAT    = MAX NUM DATA POINTS.  DIMENSION LIMITS IN MAIN PROG.
C        NUMXES    = NUMBER OF PARAMETERS IN MODEL (WHETHER SOLVED FOR OR NOT)
C        N         = NUMBER OF EQUATIONS
C        VNAME     = ALPHANUMERIC 'NAME' FOR THE PARAMETERS
C        X         = 1-D ARRAY OF INITIAL PARAMETER VALUES
C        ID        = 1-D ARRAY OF 'SOLVE FOR' CODES
C                  0 = HOLD PARAMETER CONSTANT
C                  1 = SOLVE FOR PARAMETER
C        S         = 1-D ARRAY OF DATA (E.G. VLB DIFFERENCED PHASES)
C        C         = 2-D ARRAY OF PARTIAL DERIVATIVES OF MODEL
C        E         = 1-D ARRAY OF ERRORS ON DATA POINTS
C
C     THE OUTPUT INCLUDES:
C        XHAT2     = 1-D ARRAY OF FINAL LEAST SQUARE VALUES FOR ALL PARAMS
C        THE SOLUTION AND CORRELATIONS ARE AUTOMATICALLY PRINTED OUT
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(MAXDAT,1)
      DIMENSION S(1),E(1),X(1),ID(1),XHAT2(1)

      character*16  VNAME(1)

      REAL*8        B(77,77),BSTORE(77,77),XIDENT(77,77),
     +              COR(77,77),ERROR(77),STDEV(77),PR(77),
     +              XHAT(77),SNEW(1001),PRODCT(1001)

      DIMENSION LL(77),MM(77)
      DATA      MSIZE /77/

C
      LOGICAL PRINT_COR

      character*8   QUEST(2)
      DATA QUEST/'YES     ','NO      '/

 9999 FORMAT (1X,77(1PD13.5))
 9998 FORMAT (1X,20F6.2)
 9417 FORMAT (/1X)

      IF (IDEBUG.ge.0) GO TO 502
C
      WRITE (6,5522)
 5522 FORMAT (////50X,' C MATRIX')
      DO 500 I=1,N
      WRITE (6,9999) (C(I,J),J=1,NUMXES)
  500 CONTINUE
      WRITE (6,5523)
 5523 FORMAT (////40X,'        S(I)           E(I)')
      DO 501 I=1,N
      WRITE (6,9999) S(I),E(I)
  501 CONTINUE
  502 CONTINUE
      M=0
      DO 4 I=1,NUMXES
    4 M=M+ID(I)
C        M IS THE NUMBER OF 'SOLVE FOR' PARAMETERS
      IF ( M.GT.N )  GO TO 998
C
      JNEW=0
      DO 12  J=1,NUMXES
	      IF (ID(J).EQ.0) GO TO 12
	      JNEW=JNEW+1
	      DO 11 I=1,N
   11		 C(I,JNEW)=C(I,J)
   12 CONTINUE
C
C        WEIGHT EQUATIONS BY DIVIDING EACH BY THE ERROR-E(I)-(SEE SECT 2
C        OF TEXT)
      DO 13 I=1,N
	      SNEW(I)=S(I)/E(I)
	      DO 13 J=1,M
   13		 C(I,J)=C(I,J)/E(I)
      IF (IDEBUG.ge.0) GO TO 21
      WRITE (6,2006)
 2006 FORMAT ('1'/////50X,' CSTAR MATRIX')
      DO 14 I=1,N
   14 WRITE (6,9999) (C(I,J),J=1,M)
      WRITE (6,2009)
 2009 FORMAT (/////50X,'     SSTAR')
      DO 15 I=1,N
   15 WRITE (6,9999) SNEW(I)
   21 CONTINUE
      ITEST2=0
      BMIN10=0.D0
      BMAX10=0.D0
      DO 22 I=1,M
      DO 22 J=1,M
   22 	B(I,J)=0.D0
      DO 24 I=1,M
      DO 24 J=1,M
	      DO 23 L=1,N
   23		 B(I,J)=B(I,J) + C(L,I)*C(L,J)
      IF (B(I,J).EQ.0.D0) GO TO 24
	      B1=DLOG10( DABS( B(I,J) ) )
	      IF (BMIN10.GT.B1) BMIN10=B1
	      IF (BMAX10.LT.B1) BMAX10=B1
   24 BSTORE(I,J)=B(I,J)
      IF (IDEBUG.ge.0) GO TO 761
      WRITE (6,2010)
 2010 FORMAT ('1'/////50X,' C*TRANSPOSE C')
      DO 25 I=1,M
   25 WRITE (6,9999) (B(I,J),J=1,M)
  761 CONTINUE
C
C        THE SUBROUTINE 'MINV' IS A DOUBLE PRECISION MATRIX INVERSION
C        IT INVERTS THE MATRIX 'BB' AND STORES IT AS 'BB' (I.E. THE ORIGIN
C        MATIRX IS DESTROYED
C        *******************************
C
C
      IF (DABS(BMAX10-BMIN10).LT.24.D0) GO TO 962
c      WRITE (6,9622) BMIN10,BMAX10
c 9622 FORMAT (//1X,'********************     BMIN10 = ',F6.1,'
c     1 BMAX10 = ',F6.1,'     ********************',/////1X)
  962 BADJST=10.D0**( (BMAX10+BMIN10)/2.D0 )
      DO 963 I=1,M
      DO 963 J=1,M
  963	 B(I,J)=B(I,J)/BADJST
      CALL MINV(B,M,DETERM,LL,MM,MSIZE)
      DO 964 I=1,M
      DO 964 J=1,M
  964	 B(I,J)=B(I,J)/BADJST
C
C
C        *******************************
C
      IF (IDEBUG.ge.0) GO TO 226
      WRITE (6,2011) DETERM
 2011 FORMAT ('1'////'  THE DETERMINANT IS',1PD13.5)
      WRITE (6,2022)
 2022  FORMAT (////45X,' (C*TRANSPOSE C*) INVERSE')
      DO 26 I=1,M
   26 WRITE (6,9999) (B(I,J),J=1,M)
  226 CONTINUE
      DO 27 I=1,M
      DO 27 J=1,M
   27 	XIDENT(I,J)=0.D0
      DO 30 I=1,M
      DO 30 J=1,M
	      DO 28 L=1,M
   28		 XIDENT(I,J)=XIDENT(I,J) + B(I,L)*BSTORE(L,J)
C        XIDENT = B INVERSE  B     WHICH SHOULD BE THE IDENTITY MATRIX I
C        INVERSION ROUTINE WORKED.  ITEST2 CHECKS THAT XIDENT IS THE IDE
C        MATRIX
	      IF (I.EQ.J) GO TO 29
		      IF (DABS(XIDENT(I,J)).GT.1.D-06) ITEST2=1
		      GO TO 30
   29	      CONTINUE
	      IF (DABS(XIDENT(I,I)-1.D0).GT.1.D-06) ITEST2=-1
   30 CONTINUE
      IF (ITEST2.NE.0) GO TO 50
	      DO 33 I=1,M
		      XHAT(I)=0.D0
		      DO 32 J=1,N
		      PRODCT(J)=0.D0
			      DO 31 L=1,M
   31				 PRODCT(J)=PRODCT(J) + B(I,L)*C(J,L)
   32		      XHAT(I)=XHAT(I)+PRODCT(J)*SNEW(J)
   33	      CONTINUE
C        XHAT'S ARE (IF ADDED TO THE X'S) THE UNCONSTRAINED LEAST SQUARE
C        SOLUTION OF THE 'SOLVE FOR' PARAMETERS ONLY.
	      IN=0
C        XHAT2(J) = THE LEAST SQUARES SOLTIONS (INCLUDING NON SOLVE FOR
C        PARAMETERS
	      DO 43 J=1,NUMXES
		      IF (ID(J).EQ.0) GO TO 42
		      IN=IN+1
		      XHAT2(J)=XHAT(IN) + X(J)
		      GO TO 43
   42 			XHAT2(J)=X(J)
   43	      CONTINUE
C
   50 CONTINUE
      DO 51 I=1,M
      DO 51 J=1,M
   51	 COR(I,J)=-BSTORE(I,J)/DSQRT(BSTORE(I,I)*BSTORE(J,J))
C
C
      IF (ITEST2.NE.0) GO TO 999
      INOTE=0
      DO 53 I=1,M
	      IF (B(I,I).GT.0.D0) GO TO 53
	      B(I,I)=-B(I,I)
	      INOTE = INOTE + 1
   53 CONTINUE
      IF (INOTE.EQ.0) GO TO 71
	WRITE (6,2071) INOTE
 2071	FORMAT (/////' ***** THERE ARE ',I2,' NEGATIVE VARIANCES'/'
     1             YOUR SOLUTION IS UNPHYSICAL')
   71 CONTINUE
      DO 54 J=1,M
   54	 STDEV(J)=DSQRT(B(J,J))
C        REARRANGE CALCULATED 1 SIGMA ERRORS TO APPLY TO THE SOLVED FOR
C        PARAMETERS
      IN=0
      DO 59 J=1,NUMXES
	      IF (ID(J).EQ.0) GO TO 58
		      IN=IN+1
		      ERROR(J)=STDEV(IN)
		      GO TO 59
   58	      ERROR(J)=0.D0
   59 CONTINUE

      IF ( PRINT_COR )  THEN
C        OUTPUT FOR THE UNCONSTRAINED LEAST SQUARES SOLUTION
         WRITE (6,2040)
C
 2040    FORMAT (/1X,'    PARAMETER        ORIGINAL VALUE   ',
     1        'LEAST SQUARE VALUE  1 SIGMA ERRORS   SOLVED FOR?')
C
C
          DO 44 J=1,NUMXES
	      L=1
	      IF (ID(J).EQ.0) L=2
  44	  WRITE (6,2041) VNAME(J), X(J), XHAT2(J), ERROR(J), QUEST(L)
 2041	  FORMAT (2X,A16,5X,3(F13.6,5X),4X,A8)
C
C
C	        CALCULATE THE CORRELATION COEFFICIENTS (COR)...
	      DO 55 I=1,M
	      DO 55 J=1,M
   55		 COR(I,J)=B(I,J)/DSQRT(B(I,I)*B(J,J))
C
	      WRITE (6,2056)
 2056	      FORMAT(/10X,' THE CORRELATION COEFFICIENTS')
	      DO 57 I=1,M
   57		      WRITE (6,9998) (COR(I,J),J=1,M)

C  	      THE MULTIPLE CORRELATION COEFFICIENTS (PR) ARE CALCULATED
	      DO 60 I=1,M
   60		 PR(I)= 1.D0 - (1.d0/(BSTORE(I,I)*B(I,I)))
	      WRITE (6,2060)
 2060         FORMAT (///10X,'THE MULTIPLE CORRELATION COEFFICIENTS'/)

C			not sure if the true PR is the sqrt(PR)??
              I = 0
              DO 61 J=1,NUMXES
              	IF ( ID(J).EQ.0 ) GO TO 61
              		I = I + 1
                        WRITE (6,2061) VNAME(J),PR(I)
 2061                   FORMAT (10X,A16,2X,F10.5)
   61         CONTINUE

  648 END IF
C
C
      GO TO 989
  998	WRITE (6,2094)  N,M
 2094	FORMAT (////' LEAST_SQUARES_FIT: # DATA POINTS (',I4,
     +              ') < # PARAMETERS (',I2,').  STOP' )
	STOP
C
  999 	continue

c        WRITE (6,2095)  ITEST2
c 2095 	FORMAT (////'  ***** ITEST2 =',I2,' ***** '
c     +             /'       MATRIX INVERSION FAILED.')
c	STOP
C
  989 RETURN
      END

c=========================================================

      SUBROUTINE MINV(A,N,D,L,M,MSIZE)

      IMPLICIT REAL*8 (A-H,O-Z)

C     IBM SCIENTIFIC SUBROUTINE PACAKAGE PAGE 118
C     ******************************************************************
C
C     PURPOSE  INVERT A MATRIX
C
C     DESCRIPTION OF PARAMETERS
C        A  INPUT MATRIX, DESTROYED IN COMPUTATION AND REPLACED BY INVER
C        N  ORDER OF MATRIX A
C        D  RESULTANT DETERMINANT
C        L  WORK VECTOR OF LENGTH N
C        M  WORK VECTOR OF LENGTH N
C     MSIZE ORDER OF TWO DIMENSIONAL MATRIX A IN MAIN PROGRAM
C
C     METHOD
C        THE STANDARD GAUSS-JORDAN METHOD IS USED. THE DETERMINANT IS A
C        CALCULATED. A DETERMINANT OF ZERO INDICATES THAT THE MATRIX IS
C        SINGULAR.
C
C     ******************************************************************
C
      DIMENSION A(1),L(1),M(1)
C
C     STORE MATRIX A IN VECTOR FORM (COLUMN BY COLUMN)
C     SEE SSP PROGRAM ARRAY  P.98
C
      IF(MSIZE.EQ.N) GO TO 2
      NI=MSIZE-N
      IJ=0
      NM=0
      DO 230 K=1,N
      DO 225 LL=1,N
      IJ=IJ+1
      NM=NM+1
  225 A(IJ)=A(NM)
  230 NM=NM+NI
    2 CONTINUE
C
C     SEARCH FOR LARGEST ELEMENT
C
      D=1.0D0
      NK=-N
      DO 80 K=1,N
      NK=NK+N
      L(K)=K
      M(K)=K
      KK=NK+K
      BIGA=A(KK)
      DO 20 J=K,N
      IZ=N*(J-1)
      DO 20 I=K,N
      IJ=IZ+I
   10 IF(ABS(BIGA)-ABS(A(IJ))) 15,20,20
   15 BIGA=A(IJ)
      L(K)=I
      M(K)=J
   20 CONTINUE
C
C     INTERCHANGE ROWS
C
      J=L(K)
      IF(J-K) 35,35,25
   25 KI=K-N
      DO 30 I=1,N
      KI=KI+N
      HOLD=-A(KI)
      JI=KI-K+J
      A(KI)=A(JI)
   30 A(JI)=HOLD
C
C     INTERCHANGE COLUMNS
C
   35 I=M(K)
      IF(I-K) 45,45,38
   38 JP=N*(I-1)
      DO 40 J=1,N
      JK=NK+J
      JI=JP+J
      HOLD=-A(JK)
      A(JK)=A(JI)
   40 A(JI)=HOLD
C
C     DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS CONTAINED
C     BIGA)
C
   45 IF(BIGA) 48,46,48
   46 D=0.0D0
      RETURN
   48 DO 55 I=1,N
      IF(I-K) 50,55,50
   50 IK=NK+I
      A(IK)=A(IK)/(-BIGA)
   55 CONTINUE
C
C     REDUCE MATRIX
C
      DO 65 I=1,N
      IK=NK+I
      HOLD=A(IK)
      IJ=I-N
      DO 65 J=1,N
      IJ=IJ+N
      IF(I-K) 60,65,60
   60 IF(J-K) 62,65,62
   62 KJ=IJ-I+K
      A(IJ)=HOLD*A(KJ)+A(IJ)
   65 CONTINUE
C
C     DIVIDE ROW BY PIVOT
C
      KJ=K-N
      DO 75 J=1,N
      KJ=KJ+N
      IF(J-K) 70,75,70
   70 A(KJ)=A(KJ)/BIGA
   75 CONTINUE
C
C     PRODUCT OF PIVOTS
C
      D=D*BIGA
C
C     REPLACE PIVOT BY RECIPROCAL
C
      A(KK)=1.0D0/BIGA
   80 CONTINUE
C
C     FINAL ROW AND COLUMN INTERCHANGE
C
      K=N
  100 K=K-1
      IF(K) 150,150,105
  105 I=L(K)
      IF(I-K) 120,120,108
  108 JQ=N*(K-1)
      JR=N*(I-1)
      DO 110 J=1,N
      JK=JQ+J
      HOLD=A(JK)
      JI=JR+J
      A(JK)=-A(JI)
  110 A(JI)=HOLD
  120 J=M(K)
      IF(J-K) 100,100,125
  125 KI=K-N
      DO 130 I=1,N
      KI=KI+N
      HOLD=A(KI)
      JI=KI-K+J
      A(KI)=-A(JI)
  130 A(JI)=HOLD
      GO TO 100
  150 CONTINUE
C
C     PUT MATRIX BACK INTO SQUARE FORM
C
      IF(MSIZE.EQ.N) GO TO 4
      IJ=N*N+1
      NM=N*MSIZE+1
      DO 210 K=1,N
      NM=NM-NI
      DO 210 LL=1,N
      IJ=IJ-1
      NM=NM-1
  210 A(NM)=A(IJ)
    4 CONTINUE
      RETURN
      END
