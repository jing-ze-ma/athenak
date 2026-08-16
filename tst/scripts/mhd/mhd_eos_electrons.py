# Regression test for the EOS electron fraction and ohmic_resistivity = eos
#
# The feature under test is the free-electron fraction the general EOS carries -- Saha
# equilibrium over hydrogen, helium and, when eos_metal_ionization is on, the alkalis and
# iron -- together with the ohmic diffusivity built from it,
#
#     eta = 230 sqrt(T)/x_e + 5.2e11 lnL/T^1.5
#
# It is tested in two halves, because they fail independently and the first does NOT imply
# the second. That is not a hypothetical: the tabulated x_e surface was once wired into a
# run whose startup banner reported the right electron fraction while the resistivity used
# a different number entirely, and only a fluid observable could tell.
#
#   MODEL half -- the startup banner reports mu and n_e/n_tot at the table's grid centre,
#   evaluated from the composition model itself. Zero-cycle runs pin the physics: metals
#   dominate x_e by nine orders of magnitude at 2000 K, [M/H] scales it, condensation
#   removes the donors once the gas is cold enough to rain them out and does nothing when
#   it is not, and none of it moves mu by more than a whisker.
#
#   SURFACE half -- the tabulated log10(x_e) patch is read by nothing except the
#   resistivity, so it is measured through the fluid. An Alfven wave damps as
#   exp(-eta k^2 t) (the rate is eta k^2, not 2 eta k^2, because only the magnetic half of
#   the wave energy is dissipated), so the transverse magnetic energy in the history file
#   is a reading of eta. Two things are asserted from it: that the decay agrees with an
#   ohmic_resistivity = constant run at the eta the banner's x_e predicts, which pins the
#   ABSOLUTE electron fraction and its units, and that changing [M/H] changes the decay by
#   exactly the factor the banner says it changes x_e by, which pins the propagation of
#   composition into the table without depending on the damping law at all.
#
# The input file places the table's grid centre on the wave's own (rho,T) so that the
# banner value and the value the wave sees are the same number; see the header there.

# Modules
import logging
import math
import re
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name

_input = 'tests/mhd_eos_electrons.athinput'

# The wave's temperature, which the input file has made the table's grid centre. Used only
# to turn the banner's x_e into an expected eta.
_twave = 10.0**3.297756

# A temperature range whose centre is 700 K, cold enough that KCl and Na2S have condensed
# out. 0.5*(2.545 + 3.145) = 2.845 and 10^2.845 = 700 K.
_cold_logt = ['mhd/eos_logt_min=2.545', 'mhd/eos_logt_max=3.145']

# Tolerances. The two ratio tests compare one measurement against another from the same
# run and are held tight; the absolute test carries the interpolation error of the x_e
# surface at the grid centre and the wave's own truncation error, so it is held to 10%.
_rtol_ratio = 0.03
_rtol_abs = 0.10
# The ideal-gas comparison carries, on top of the above, the small mismatch between the
# temperature Rgas implies and the one the table gives this background.
_rtol_ideal = 0.10

# Results are handed from run() to analyze() through these.
_banner = {}
_decay = {}


def _parse_banner(out):
    """mu and n_e/n_tot at the grid centre, from the general EOS startup report."""
    m = re.search(r'mu at the grid centre\s*=\s*(\S+),\s*n_e/n_tot\s*=\s*(\S+)', out)
    if m is None:
        raise athena.AthenaError('Could not find the EOS report in the output')
    return float(m.group(1)), float(m.group(2))


def _banner_run(label, arguments):
    out = athena.run_output(_input, ['job/basename=eosel_' + label,
                                     'time/nlim=0'] + arguments)
    mu, xe = _parse_banner(out)
    _banner[label] = (mu, xe)
    logger.info('banner {0}: mu = {1:g}, n_e/n_tot = {2:g}'.format(label, mu, xe))
    return xe


def _fluid_run(label, arguments):
    athena.run(_input, ['job/basename=eosel_' + label] + arguments)
    data = athena_read.hst('build/src/eosel_{0}.mhd.hst'.format(label))
    # by0 = bz0 = 0, so the whole of 3-ME is the wave and none of it is the background
    return data['3-ME'][-1]


# Run AthenaK
def run(**kwargs):
    logger.debug('Runnning test ' + __name__)

    # ------------------------------------------------------------- the model, at 2000 K
    _banner_run('nometal', ['mhd/eos_metal_ionization=false'])
    xe0 = _banner_run('mh00', ['mhd/eos_metal_mh=0.0'])
    _banner_run('mh05', ['mhd/eos_metal_mh=0.5'])
    _banner_run('hotcond', ['mhd/eos_metal_condensation=true'])

    # -------------------------------------------------------------- the model, at 700 K
    _banner_run('cold', _cold_logt)
    _banner_run('coldcond', _cold_logt + ['mhd/eos_metal_condensation=true'])

    # ------------------------------------------------------------------------- the fluid
    # eta the banner's electron fraction implies at the wave's temperature: the
    # electron-neutral term plus Spitzer with the Coulomb logarithm folded in at lnL = 20,
    # exactly as ResistivityEOS() forms it
    eta = 230.0*math.sqrt(_twave)/xe0 + 5.2e11*20.0/(_twave*math.sqrt(_twave))
    logger.info('expected eta at the wave = {0:g} cm^2/s'.format(eta))

    # the same wave with no resistivity at all, to measure how much of the decay is the
    # scheme's own truncation error rather than eta
    _decay['ideal'] = _fluid_run('ideal', ['mhd/ohmic_resistivity=constant',
                                           'mhd/eta_ohm_const=0.0'])
    _decay['const'] = _fluid_run('const', ['mhd/ohmic_resistivity=constant',
                                           'mhd/eta_ohm_const={0!r}'.format(eta)])
    _decay['mh00'] = _fluid_run('fmh00', ['mhd/eos_metal_mh=0.0'])
    _decay['mh05'] = _fluid_run('fmh05', ['mhd/eos_metal_mh=0.5'])

    # ------------------------------------------------- the same x_e on an IDEAL gas
    # An ideal gas can have the tabulated electron fraction without the tabulated
    # thermodynamics: the table is built for x_e alone. The input's Rgas is chosen so the
    # ideal gas sits at the same 1985 K the table gives this background, so with the same
    # x_e it must produce the same eta and damp the wave the same way. That is what makes
    # this a test of the wiring rather than of the ideal gas.
    _decay['idealgas'] = _fluid_run('fideal', ['mhd/eos=ideal'])

    # ------------------------------------------------------------------- the refusals
    # An isothermal EOS has no composition and no table is built for it, so asking for the
    # EOS electron fraction has to be fatal rather than silently resistivity-free.
    _decay['isothermal_aborts'] = True
    try:
        athena.run(_input, ['job/basename=eosel_abort', 'time/nlim=0',
                            'mhd/eos=isothermal'])
        _decay['isothermal_aborts'] = False
    except athena.AthenaError:
        pass


# Analyze outputs
def analyze():
    logger.debug('Analyzing test ' + __name__)
    status = True

    mu_off, xe_off = _banner['nometal']
    mu_0, xe_0 = _banner['mh00']
    mu_05, xe_05 = _banner['mh05']
    mu_hc, xe_hc = _banner['hotcond']
    mu_c, xe_c = _banner['cold']
    mu_cc, xe_cc = _banner['coldcond']

    # ---- the metals are what ionizes at 2000 K -----------------------------------------
    # hydrogen's 13.6 eV contributes essentially nothing here, so switching the alkalis
    # and iron off has to cost many orders of magnitude, not a factor of a few
    if not xe_0 > 1.0e6*xe_off:
        logger.warning('eos_metal_ionization barely changes the electron fraction '
                       '({0:g} with, {1:g} without); the metal donors are not reaching '
                       'the composition model'.format(xe_0, xe_off))
        status = False

    # ---- [M/H] scales x_e, and it is trace Saha so the scaling is sublinear ------------
    ratio = xe_05/xe_0
    logger.info('[M/H] 0 -> 0.5 scales n_e/n_tot by {0:.4f}'.format(ratio))
    if not 1.5 < ratio < 2.2:
        logger.warning('[M/H] = 0.5 scales the electron fraction by {0:g}; trace Saha '
                       'over a 3.16x abundance change should give close to sqrt(3.16)'
                       .format(ratio))
        status = False

    # ---- condensation is a cliff, not a taper ------------------------------------------
    # K and Na condense last and dominate n_e, so above their condensation temperature
    # turning condensation on must do nothing at all
    if not abs(xe_hc/xe_0 - 1.0) < 1.0e-3:
        logger.warning('eos_metal_condensation changes the electron fraction by {0:g} at '
                       '{1:g} K, where nothing should have condensed yet'
                       .format(xe_hc/xe_0 - 1.0, _twave))
        status = False
    # and below it must bite
    logger.info('condensation at 700 K scales n_e/n_tot by {0:.4f}'.format(xe_cc/xe_c))
    if not xe_cc < 0.8*xe_c:
        logger.warning('eos_metal_condensation leaves the electron fraction at {0:g} of '
                       'its uncondensed value at 700 K; the alkalis should be raining out'
                       .format(xe_cc/xe_c))
        status = False

    # ---- none of it is allowed to move the mean molecular weight -----------------------
    # the donors were already inside the inert metal lump and give one heavy particle each
    # either way; all they add is their electrons
    for name, mu in (('nometal', mu_off), ('mh05', mu_05), ('hotcond', mu_hc)):
        if not abs(mu/mu_0 - 1.0) < 1.0e-4:
            logger.warning('mu moves from {0:g} to {1:g} between the mh00 and {2} cases; '
                           'the metal donors should not change it'
                           .format(mu_0, mu, name))
            status = False
    if not abs(mu_cc/mu_c - 1.0) < 1.0e-4:
        logger.warning('mu moves from {0:g} to {1:g} when condensation is switched on at '
                       '700 K'.format(mu_c, mu_cc))
        status = False

    # ---- the fluid sees the same electron fraction the banner reports ------------------
    # log decay of the transverse magnetic energy, with the resistivity-free run divided
    # out so that what is left is the ohmic part alone
    me_ideal = _decay['ideal']
    damp = {}
    for key in ('const', 'mh00', 'mh05'):
        damp[key] = math.log(me_ideal/_decay[key])
        logger.info('{0}: 3-ME {1:g} vs {2:g} without resistivity, ohmic log decay {3:g}'
                    .format(key, _decay[key], me_ideal, damp[key]))

    if not damp['const'] > 0.1:
        logger.warning('a constant resistivity at the expected eta damps the wave by '
                       'only {0:g} in the log; either eta_ohm_const is not being applied '
                       'or the test state is wrong'.format(damp['const']))
        return False

    # ABSOLUTE: eos must reproduce the constant run at the eta its own x_e implies
    err = damp['mh00']/damp['const'] - 1.0
    logger.info('ohmic_resistivity = eos damps {0:+.2%} relative to constant at the same '
                'eta'.format(err))
    if not abs(err) < _rtol_abs:
        logger.warning('ohmic_resistivity = eos damps the wave {0:+.2%} differently from '
                       'ohmic_resistivity = constant at the eta the banner predicts; the '
                       'tabulated x_e is not the composition model\'s x_e'.format(err))
        status = False

    # RELATIVE: and the metallicity has to propagate through the table by the same factor
    fluid_ratio = damp['mh00']/damp['mh05']
    err = fluid_ratio/ratio - 1.0
    logger.info('[M/H] 0 -> 0.5 changes the damping by {0:.4f}, the banner x_e by {1:.4f}'
                .format(fluid_ratio, ratio))
    if not abs(err) < _rtol_ratio:
        logger.warning('raising [M/H] to 0.5 changes the ohmic damping by {0:g} but the '
                       'electron fraction by {1:g}; the metallicity is not reaching the '
                       'tabulated surface'.format(fluid_ratio, ratio))
        status = False

    # ---- an ideal gas reading the same table must damp the same way --------------------
    damp_ideal = math.log(me_ideal/_decay['idealgas'])
    err = damp_ideal/damp['mh00'] - 1.0
    logger.info('ideal gas + tabulated x_e: ohmic log decay {0:g} against {1:g} for the '
                'general EOS, {2:+.2%}'.format(damp_ideal, damp['mh00'], err))
    if not abs(err) < _rtol_ideal:
        logger.warning('an ideal gas at the same temperature reading the same x_e table '
                       'damps the wave {0:+.2%} differently from the general EOS; the '
                       'ideal-gas path is not reaching the same electron fraction'
                       .format(err))
        status = False

    # ---- and the refusal -------------------------------------------------------------
    if not _decay['isothermal_aborts']:
        logger.warning('ohmic_resistivity = eos with eos = isothermal did not abort; no '
                       'table is built for it, so there is no electron fraction to give')
        status = False

    return status
