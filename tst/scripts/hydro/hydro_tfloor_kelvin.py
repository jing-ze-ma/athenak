# Regression test for <block>/tfloor_kelvin, the temperature floor stated in kelvin
#
# `tfloor` is a floor on the CODE temperature, and that is not one quantity: an ideal gas
# floors p/d, whose kelvin scale carries the fixed <units>/mu, while a general EOS runs at
# mu_ref = 1 with composition inside the EOS, so its scale is the same one with that mu
# divided back out. A value copied from an ideal-gas input into a general-EOS one is
# therefore wrong by a factor of mu, silently. tfloor_kelvin removes the trap by letting
# both inputs carry the same physical number.
#
# What is checked:
#   * the conversion, in BOTH modes, against the two scales computed here independently of
#     the code -- the ideal one through <units>/mu and the general one without it. The
#     input sets mu = 2.32408 precisely so the two answers differ and the mu factor is
#     visible rather than cancelling.
#   * that the converted floor is APPLIED and not merely reported: the input's background
#     sits at ~1985 K under both EOSs, so a 2500 K floor has to change the solution and a
#     100 K floor has to leave it alone.
#   * that the parameter is refused where `tfloor` does not mean this quantity at all
#     (isothermal EOS).
#
# Not covered here, because the harness cannot express it: setting both tfloor and
# tfloor_kelvin is a fatal error, but an input may declare only one of the two, and the
# command line can only override parameters the input already declares.

# Modules
import logging
import re
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name

_input = 'tests/tfloor_kelvin.athinput'

# The input's <units>: cgs throughout, so one code temperature is m_u/k_B kelvin for the
# general EOS and mu times that for an ideal gas.
_mu = 2.32408
_k_per_code_general = 1.66053906660e-24/1.380649e-16
_k_per_code_ideal = _mu*_k_per_code_general

_tk = 2500.0        # above the ~1985 K background: the floor bites
_tk_low = 100.0     # far below it: the floor does nothing

# The converted floor is read back from the startup report, which prints at the default
# six significant digits, so that -- not the arithmetic -- is what sets the tolerance.
_rtol = 2.0e-5

_res = {}


def _run(label, arguments):
    """Run, returning (reported code tfloor or None, L1 error of the wave)."""
    out = athena.run_output(_input, ['job/basename=tfk_' + label] + arguments)
    m = re.search(r'tfloor_kelvin\s*=\s*\S+\s*K\s*->\s*tfloor\s*=\s*(\S+)', out)
    tfloor = float(m.group(1)) if m else None
    data = athena_read.error_dat('build/src/tfk_{0}-errs.dat'.format(label))
    return tfloor, data[-1][4]


# Run AthenaK
def run(**kwargs):
    logger.debug('Runnning test ' + __name__)
    _res['general'] = _run('general', ['hydro/eos=general'])
    _res['ideal'] = _run('ideal', ['hydro/eos=ideal'])
    _res['general_low'] = _run('general_low',
                               ['hydro/eos=general',
                                'hydro/tfloor_kelvin={0!r}'.format(_tk_low)])
    _res['ideal_low'] = _run('ideal_low',
                             ['hydro/eos=ideal',
                              'hydro/tfloor_kelvin={0!r}'.format(_tk_low)])

    # an isothermal EOS has no temperature floor for this to set, so asking for one is
    # fatal rather than quietly ignored
    _res['isothermal_aborts'] = True
    try:
        athena.run(_input, ['job/basename=tfk_iso', 'hydro/eos=isothermal',
                            'time/nlim=0'])
        _res['isothermal_aborts'] = False
    except athena.AthenaError:
        pass


# Analyze outputs
def analyze():
    logger.debug('Analyzing test ' + __name__)
    status = True

    expect = {'general': _tk/_k_per_code_general, 'ideal': _tk/_k_per_code_ideal}
    for mode in ('general', 'ideal'):
        got = _res[mode][0]
        if got is None:
            logger.warning('%s: the run reported no tfloor_kelvin conversion at all; '
                           'the parameter was ignored', mode)
            return False
        err = got/expect[mode] - 1.0
        logger.info('%s: tfloor_kelvin = %g K -> %g (expected %g, %+.2e)',
                    mode, _tk, got, expect[mode], err)
        if not abs(err) < _rtol:
            logger.warning('%s: tfloor_kelvin = %g K became %g in code temperature, but '
                           'this units block makes it %g', mode, _tk, got, expect[mode])
            status = False

    # the ideal and general scales must differ by exactly <units>/mu -- the whole reason
    # the parameter exists
    ratio = _res['general'][0]/_res['ideal'][0]
    logger.info('general/ideal code tfloor = %.6f, <units>/mu = %.6f', ratio, _mu)
    if not abs(ratio/_mu - 1.0) < _rtol:
        logger.warning('the same kelvin floor gave code values differing by %g, but the '
                       'two temperature scales differ by <units>/mu = %g', ratio, _mu)
        status = False

    # and it has to be applied, not just printed: 2500 K is above the ~1985 K background
    # and 100 K is below it, so the two must not give the same answer
    for mode in ('general', 'ideal'):
        hi, lo = _res[mode][1], _res[mode + '_low'][1]
        logger.info('%s: L1 with the floor at %g K = %g, at %g K = %g',
                    mode, _tk, hi, _tk_low, lo)
        if not lo > 0.0:
            logger.warning('%s: zero L1 error, the test is vacuous', mode)
            return False
        if abs(hi/lo - 1.0) < 1.0e-10:
            logger.warning('%s: a 2500 K floor and a 100 K floor give the same answer on '
                           'a ~1985 K background, so the converted floor is reported but '
                           'never applied', mode)
            status = False

    if not _res['isothermal_aborts']:
        logger.warning('tfloor_kelvin with an isothermal EOS did not abort, although an '
                       'isothermal C2P has no temperature floor to set')
        status = False

    return status
