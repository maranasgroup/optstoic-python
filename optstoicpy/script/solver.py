import pulp
from optstoicpy.script.utils import create_logger

def load_pulp_solver(
    solver_names=['SCIP_CMD', 'GUROBI', 'GUROBI_CMD', 'CPLEX_CMD', 'GLPK_CMD'],
    logger=None):
    """Load a pulp solver based on what is available.

    Args:
        solver_name (`list` of `str`, optional): A list of solver names in the order of
            loading preferences.
        logger (None, optional): A logging.Logger object

    Returns:
        TYPE: Description
    """
    if logger is None:
        logger = create_logger('optstoic.load_pulp_solver')

    if isinstance(solver_names, str):
        solver_names = [solver_names]
    elif not isinstance(solver_names, list):
        raise Exception("Argument solver_names must be a list!")

    GUROBI_CMD_OPTIONS = [('Threads', 2), ('TimeLimit', 1800),
                          ('MIPGapAbs', 1e-6), ('MIPGap', 1e-6), ('CliqueCuts', 2)]

    CPLEX_CMD_OPTIONS = ['mip tolerances mipgap 1e-6','mip tolerances absmipgap 1e-6']

    GLPK_CMD_OPTIONS  = ['--clique', '--pcost', '--gomory', '--mipgap', '1e-6']

    pulp_solvers = {
        'SCIP_CMD': pulp.solvers.SCIP_CMD(
            keepFiles=0,
            mip=True,
            msg=True),
        'GUROBI': pulp.solvers.GUROBI(
            mip=True,
            msg=True,
            timeLimit=1800,
            MIPGapAbs=1e-6),
        'GUROBI_CMD': pulp.solvers.GUROBI_CMD(
            path=None,
            keepFiles=0,
            mip=1,
            msg=1,
            options=GUROBI_CMD_OPTIONS),
        'CPLEX_CMD': pulp.solvers.CPLEX_CMD(
            path=None,
            keepFiles=0,
            mip=1,
            msg=1,
            options=CPLEX_CMD_OPTIONS,
            timelimit=1800),
        'GLPK_CMD': pulp.solvers.GLPK_CMD(
            msg=1,
            mip=1,
            options=GLPK_CMD_OPTIONS)
        }

    # Load solvers in the order of preferences
    pulp_solver = None

    for solver_name in solver_names:
        ps = pulp_solvers.get(solver_name, None)

        if ps is None:
            raise Exception("The solver %s is not defined in the function!"%solver_name)

        if ps.available():
            pulp_solver = ps
            logger.warning("Pulp solver set to %s."%solver_name)

            if hasattr(pulp_solver,'tmpDir'):
                pulp_solver.tmpDir = './'

            if solver_name == 'GLPK_CMD':
                logger.warning("GLPK takes a significantly longer time to solve "
                                "OptStoic. Please be patient.")
            break

    if pulp_solver is None:
        logger.warning("No solver is available!")
        return None

    return pulp_solver