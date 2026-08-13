# Regression test for the TABULATED general EOS in non-relativistic MHD
#
# Companion to hydro/hydro_general_eos_table.py; see that file for why the tabulated EOS
# is checked by convergence rather than against the ideal-gas path.
#
# The extra thing this covers is the MHD eigensystem. In the primitive variables the ratio
# of specific heats enters the fast/slow/Alfven eigenvectors only through the sound speed,
# so the general-EOS fast wave is the ideal-gas one with Gamma_1 substituted for gamma. If
# that substitution were wrong the initial state would not be an eigenmode and the error
# would stall instead of converging.

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
_input = 'tests/linear_wave_mhd_geneos_table.athinput'
_min_rate = 1.8
_min_contrast = 1.25


# Run AthenaK
def run(**kwargs):
    logger.debug('Runnning test ' + __name__)
    for n in _res:
        arguments = ['job/basename=mhd_gen_eos_table',
                     'mesh/nx1={0}'.format(n),
                     'meshblock/nx1={0}'.format(n),
                     'mhd/reconstruct=plm',
                     'mhd/rsolver=hlld']
        athena.run(_input, arguments)

    # the same wave under the gamma-law mode of the same general EOS, for contrast
    athena.run(_input, ['job/basename=mhd_gen_eos_gammamode',
                        'mesh/nx1={0}'.format(_res[0]),
                        'meshblock/nx1={0}'.format(_res[0]),
                        'mhd/reconstruct=plm',
                        'mhd/rsolver=hlld',
                        'mhd/general_eos=gamma'])


# Analyze outputs
def analyze():
    logger.debug('Analyzing test ' + __name__)
    data = athena_read.error_dat('build/src/mhd_gen_eos_table-errs.dat')
    gamma_data = athena_read.error_dat('build/src/mhd_gen_eos_gammamode-errs.dat')
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
            logger.warning("Tabulated EOS fast wave converges at order {0:.2f} between "
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
