# Regression test for the general EOS interface in non-relativistic MHD
#
# Companion to hydro/hydro_general_eos.py. The general EOS path (<mhd>/eos = general)
# evaluates a gamma law for now, so it MUST reproduce the ideal-gas path run for run: the
# same linear wave problem is run twice, changing nothing but the EOS, and the L1 errors
# written to mhd_gen_eos-errs.dat are required to match.
#
# As in hydro, hlle is exempt from equality -- its ideal-gas wave speeds use a Roe average
# that has no general-EOS analogue, so the general path falls back to a Davis/Einfeldt
# bound and the two differ at truncation level by design.
#
# NOTE: these equality checks stop being meaningful once the analytic partial-ionization
# EOS replaces the gamma-law placeholder, and must be revisited then.

# Modules
import logging
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name

_recon = ['plm', 'wenoz']
_flux_exact = ['llf', 'hlld']   # must reproduce the ideal path
_flux_close = ['hlle']          # differs by construction, only required to be close
_flux = _flux_exact + _flux_close
_eos = ['ideal', 'general']
_grid = ['uniform', 'smr']
_input = 'tests/linear_wave_mhd_geneos.athinput'


# Run AthenaK
def run(**kwargs):
    logger.debug('Runnning test ' + __name__)
    for gv in _grid:
        if gv == 'smr':
            grid_args = ['mesh_refinement/refinement=static',
                         'meshblock/nx1=8', 'meshblock/nx2=8', 'meshblock/nx3=8']
        else:
            grid_args = ['mesh_refinement/refinement=none',
                         'meshblock/nx1=32', 'meshblock/nx2=16', 'meshblock/nx3=16']
        for rv in _recon:
            for fv in _flux:
                for ev in _eos:
                    arguments = grid_args + ['job/basename=mhd_gen_eos',
                                             'time/tlim=1.0',
                                             'time/nlim=200',
                                             'mesh/nghost=4',
                                             'mesh/nx1=32',
                                             'mesh/nx2=16',
                                             'mesh/nx3=16',
                                             'mhd/reconstruct=' + rv,
                                             'mhd/rsolver=' + fv,
                                             'mhd/eos=' + ev,
                                             'problem/wave_flag=0',
                                             'problem/amp=1.0e-4']
                    athena.run(_input, arguments)


# Analyze outputs
def analyze():
    logger.debug('Analyzing test ' + __name__)
    data = athena_read.error_dat('build/src/mhd_gen_eos-errs.dat')
    data = data.reshape([len(_grid), len(_recon), len(_flux), len(_eos),
                         data.shape[-1]])
    analyze_status = True

    for gi, gv in enumerate(_grid):
        for ri, rv in enumerate(_recon):
            for fi, fv in enumerate(_flux):
                l1_ideal = data[gi][ri][fi][_eos.index('ideal')][4]
                l1_gen = data[gi][ri][fi][_eos.index('general')][4]
                if not l1_ideal > 0.0:
                    logger.warning("Zero L1 error for {0}+{1}+{2}, test is vacuous"
                                   .format(gv, rv, fv))
                    analyze_status = False
                if fv in _flux_exact:
                    if l1_gen != l1_ideal:
                        logger.warning("General EOS does not reproduce ideal EOS for "
                                       "{0}+{1}+{2}: ideal {3:.17g} general {4:.17g}"
                                       .format(gv, rv, fv, l1_ideal, l1_gen))
                        analyze_status = False
                else:
                    if not 0.5 < l1_gen/l1_ideal < 2.0:
                        logger.warning("General EOS hlle error too far from ideal for "
                                       "{0}+{1}: ideal {2:g} general {3:g}"
                                       .format(gv, rv, l1_ideal, l1_gen))
                        analyze_status = False

    return analyze_status
