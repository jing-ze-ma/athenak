# Regression test for the general EOS interface in non-relativistic hydro
#
# The general EOS path (<hydro>/eos = general) evaluates a gamma law for now, so it MUST
# reproduce the ideal-gas path (<hydro>/eos = ideal) run for run. That is what this test
# pins down: the same linear wave problem is run twice, changing nothing but the EOS, and
# the L1 errors the executable writes to hydro_gen_eos-errs.dat are required to match.
#
# Two solvers are deliberately excluded from the "must match" set:
#   * hlle  -- its ideal-gas wave speeds use a Roe average, which has no general-EOS
#              analogue without Vinokur-Montagne averaging, so the general path falls back
#              to a Davis/Einfeldt bound. The two differ at truncation level BY DESIGN, so
#              they are checked to be close rather than equal.
#   * roe   -- refused outright under a general EOS, because the Roe average is only
#              defined for an ideal gas. The refusal itself is tested.
#
# NOTE: when the analytic (partial-ionization) EOS replaces the gamma-law placeholder,
# the equality checks below stop being meaningful and this test has to be revisited --
# either by giving the general EOS a gamma-law option to select here, or by replacing
# these checks with convergence tests against known thermodynamics.

# Modules
import logging
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name

_recon = ['plm', 'wenoz']
_flux_exact = ['llf', 'hllc']   # must reproduce the ideal path
_flux_close = ['hlle']          # differs by construction, only required to be close
_flux = _flux_exact + _flux_close
_eos = ['ideal', 'general']
_grid = ['uniform', 'smr']
_input = 'tests/linear_wave_hydro_geneos.athinput'

# roe under a general EOS must abort rather than silently produce a wrong answer
_roe_refused = None


# Run AthenaK
def run(**kwargs):
    logger.debug('Runnning test ' + __name__)
    global _roe_refused
    for gv in _grid:
        # SMR needs small MeshBlocks so the refined region is resolvable
        if gv == 'smr':
            grid_args = ['mesh_refinement/refinement=static',
                         'meshblock/nx1=8', 'meshblock/nx2=8', 'meshblock/nx3=8']
        else:
            grid_args = ['mesh_refinement/refinement=none',
                         'meshblock/nx1=32', 'meshblock/nx2=16', 'meshblock/nx3=16']
        for rv in _recon:
            for fv in _flux:
                for ev in _eos:
                    arguments = grid_args + ['job/basename=hydro_gen_eos',
                                             'time/tlim=1.0',
                                             'time/nlim=200',
                                             'mesh/nghost=4',
                                             'mesh/nx1=32',
                                             'mesh/nx2=16',
                                             'mesh/nx3=16',
                                             'hydro/reconstruct=' + rv,
                                             'hydro/rsolver=' + fv,
                                             'hydro/eos=' + ev,
                                             'problem/wave_flag=0',
                                             'problem/amp=1.0e-4']
                    athena.run(_input, arguments)

    # The Roe solver has no general-EOS form and must exit rather than approximate one.
    try:
        athena.run(_input, ['job/basename=hydro_gen_eos_roe',
                            'time/tlim=1.0', 'time/nlim=1',
                            'mesh/nghost=4',
                            'mesh/nx1=32', 'mesh/nx2=16', 'mesh/nx3=16',
                            'mesh_refinement/refinement=none',
                            'meshblock/nx1=32', 'meshblock/nx2=16',
                            'meshblock/nx3=16',
                            'hydro/rsolver=roe', 'hydro/eos=general'])
        _roe_refused = False
    except athena.AthenaError:
        _roe_refused = True


# Analyze outputs
def analyze():
    logger.debug('Analyzing test ' + __name__)
    data = athena_read.error_dat('build/src/hydro_gen_eos-errs.dat')
    data = data.reshape([len(_grid), len(_recon), len(_flux), len(_eos),
                         data.shape[-1]])
    analyze_status = True

    for gi, gv in enumerate(_grid):
        for ri, rv in enumerate(_recon):
            for fi, fv in enumerate(_flux):
                l1_ideal = data[gi][ri][fi][_eos.index('ideal')][4]
                l1_gen = data[gi][ri][fi][_eos.index('general')][4]
                # a wave that did not evolve would make the comparison vacuous
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
                    # hlle: different wave-speed estimate, so only require the errors to
                    # stay within a factor of two of each other
                    if not 0.5 < l1_gen/l1_ideal < 2.0:
                        logger.warning("General EOS hlle error too far from ideal for "
                                       "{0}+{1}: ideal {2:g} general {3:g}"
                                       .format(gv, rv, l1_ideal, l1_gen))
                        analyze_status = False

    if not _roe_refused:
        logger.warning("The Roe solver did not abort under eos = general; its Roe "
                       "average is only defined for an ideal gas and must be refused")
        analyze_status = False

    return analyze_status
