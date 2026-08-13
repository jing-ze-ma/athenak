# Regression test for the TABULATED general EOS in non-relativistic hydro
#
# hydro_general_eos.py pins the general EOS *interface* down by requiring it to reproduce
# the ideal-gas path; that only works because the general EOS defaults to evaluating a
# gamma law (<hydro>/general_eos = gamma). This test covers the other mode,
# <hydro>/general_eos = table, where the EOS is the tabulated H2 + Saha + radiation one
# and there is no ideal-gas run to compare against.
#
# What it checks is convergence. The background state sits in the hydrogen
# partial-ionization zone, where Gamma_1 is about 1.30 rather than 5/3, and a sound wave
# is an eigenmode of the linearized equations only if the sound speed, the eigenvector
# used to build the initial condition, and the pressure the Riemann solver reconstructs
# all come from the SAME thermodynamics. If any of them fell back to the gamma law the
# initial state would contain an admixture of the other modes, and the L1 error after one
# wave period would stall at the amplitude rather than converging away.
#
# So: run the same wave at three resolutions and require second-order convergence. A
# convergence test alone would still pass if `general_eos = table` silently did nothing,
# so the same wave is also run in gamma mode and the two error levels are required to
# differ -- the tabulated EOS carries dissociation and ionization energy, so its total
# energy, and with it the error scale, is about twice the gamma-law one.

# Modules
import logging
import math
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name

_res = [64, 128, 256]
_input = 'tests/linear_wave_hydro_geneos_table.athinput'
# second order in space and time; the margin allows for the leading correction at these
# resolutions, but is nowhere near enough to let a first-order result through
_min_rate = 1.8
# how far the gamma-law run must sit from the tabulated one for the table to count as
# actually in use (the observed ratio is about 2)
_min_contrast = 1.25


# Run AthenaK
def run(**kwargs):
    logger.debug('Runnning test ' + __name__)
    for n in _res:
        arguments = ['job/basename=hydro_gen_eos_table',
                     'mesh/nx1={0}'.format(n),
                     'meshblock/nx1={0}'.format(n),
                     'hydro/reconstruct=plm',
                     'hydro/rsolver=hllc']
        athena.run(_input, arguments)

    # the same wave under the gamma-law mode of the same general EOS, for contrast
    athena.run(_input, ['job/basename=hydro_gen_eos_gammamode',
                        'mesh/nx1={0}'.format(_res[0]),
                        'meshblock/nx1={0}'.format(_res[0]),
                        'hydro/reconstruct=plm',
                        'hydro/rsolver=hllc',
                        'hydro/general_eos=gamma'])


# Analyze outputs
def analyze():
    logger.debug('Analyzing test ' + __name__)
    data = athena_read.error_dat('build/src/hydro_gen_eos_table-errs.dat')
    gamma_data = athena_read.error_dat('build/src/hydro_gen_eos_gammamode-errs.dat')
    analyze_status = True

    if data.shape[0] != len(_res):
        logger.warning("Expected {0} rows in the error file, found {1}"
                       .format(len(_res), data.shape[0]))
        return False

    # column 4 is the RMS L1 error over all conserved variables
    errs = [row[4] for row in data]
    if not errs[0] > 0.0:
        logger.warning("Zero L1 error at the coarsest resolution, test is vacuous")
        return False

    for i in range(len(_res) - 1):
        rate = math.log(errs[i]/errs[i+1], 2.0)
        logger.info("nx1 {0} -> {1}: L1 {2:g} -> {3:g}, order {4:.2f}"
                    .format(_res[i], _res[i+1], errs[i], errs[i+1], rate))
        if rate < _min_rate:
            logger.warning("Tabulated EOS linear wave converges at order {0:.2f} between "
                           "nx1 = {1} and {2}, expected at least {3}"
                           .format(rate, _res[i], _res[i+1], _min_rate))
            analyze_status = False

    err_gamma = gamma_data[-1][4]
    contrast = max(errs[0]/err_gamma, err_gamma/errs[0])
    logger.info("table vs gamma-law L1 at nx1 = {0}: {1:g} vs {2:g}"
                .format(_res[0], errs[0], err_gamma))
    if not contrast > _min_contrast:
        logger.warning("general_eos = table gives the same answer as general_eos = gamma "
                       "({0:g} vs {1:g}); the table does not appear to be in use"
                       .format(errs[0], err_gamma))
        analyze_status = False

    return analyze_status
