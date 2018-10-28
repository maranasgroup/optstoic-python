import copy
import json
import pulp
from nose.tools import (
    assert_equal)
from sympy import Matrix, S, nsimplify
from optstoicpy.core.database import (
    load_custom_reactions_to_be_excluded, 
    load_base_reaction_db
    )
from optstoicpy.script.utils import create_logger
from optstoicpy.script.solver import load_pulp_solver

def blocked_reactions_analysis(
    database, 
    pulp_solver,
    user_defined_export_rxns_Sji,
    specific_bounds,
    custom_flux_constraints,
    excluded_reactions=None, 
    target_reactions_list=None,
    logger=None):
    """
    Perform flux variability analysis on the database,
    based on the overall reaction equation of optstoic.
    If a reaction cannot carry flux (i.e., -eps <= v(j) <= eps, where eps = 1e-8), 
    then the reaction is considered as a blocked reaction.
    The blocked reactions are then eliminated from the database S matrix.
    Next, the internal loops (excluding cofactors) are identified. 
    Then, optStoic analysis can be performed for pathway prospecting.
    
    max/min v(j)
    
    subject to:
            sum(j, S(i,j) * v(j)) = 0, for all i
            custom_flux_constraints
    
        Note: The glycolysis study was done using the GAMS version of this code.
        This is written in attempt to port find_blocked_reactions.gms from GAMS to Python,
        as a part of effort to generalize optstoic analysis.
    
    Args:
        database (:obj:`BaseReactionDatabase`): The default reaction database
            without blocked reactions/loops.
        pulp_solver (TYPE): The solver for PuLP.
        user_defined_export_rxns_Sji (dict): The list of export reactions that
            need to be added to the model for metabolite exchange (i.e., any metabolite
            that participate in the design equation)
        specific_bounds (dict): LB and UB for exchange reactions which defined the
            overall design equations. E.g. {'Ex_glc': {'LB': -1, 'UB':-1}}
        custom_flux_constraints (TYPE): The custom constraints that need to be 
            added to the model formulation.
        excluded_reactions (None, optional): The list of reactions that are manually
            selected to be excluded from optstoic solution.
        target_reactions_list (None, optional): If provided, the blocked reaction analysis is performed
            only on a subset of the reaction provided. If None, the blocked reaction analysis
            will be performed on all reactions in the database. The excluded_reactions set
            can be subtracted(e.g., set(database.reactions) - excluded_reactions), since
            they are blocked reactions.
        logger (:obj:`logging.logger`, optional): The logging instance
    
    Returns:
        TYPE: Description
    
    Raises:
        ValueError: Description
    """
    if logger is None:
        logger = create_logger(name="Blocked reaction analysis")

    logger.warning("This process may take a long time to run. It is recommended to be run in a batch script.")
    
    M = 1000
    EPS = 1e-8

    # Initialize variables
    v = pulp.LpVariable.dicts("v", database.reactions,
                              lowBound=-M, upBound=M, cat='Continuous')

    for j in database.reactions:
        if database.rxntype[j] == 0:
            # Forward irreversible
            v[j].lowBound = 0
            v[j].upBound = M

        elif database.rxntype[j] == 1:
            # Reversible
            v[j].lowBound = -M
            v[j].upBound = M

        elif database.rxntype[j] == 2:
            # Reverse irreversible
            v[j].lowBound = -M
            v[j].upBound = 0

        elif database.rxntype[j] == 4:
            v[j].lowBound = 0
            v[j].upBound = 0

        else:
            raise ValueError("Reaction type for reaction %s is unknown."%j)

    if excluded_reactions is not None:
        for j in excluded_reactions:
            v[j].lowBound = 0
            v[j].upBound = 0

    # Fix stoichiometry of source/sink metabolites
    for j, bounds in specific_bounds.iteritems():
        v[j].lowBound = bounds['LB']
        v[j].upBound = bounds['UB']

    FVA_res = {}
    blocked_reactions = []    
    lp_prob = None

    if target_reactions_list is None:
        target_reactions_list = database.reactions
    num_rxn = len(target_reactions_list)

    for ind, j1 in enumerate(target_reactions_list):
        logger.debug("%s/%s"%(ind,num_rxn))
        FVA_res[j1] = {}

        for obj in ['min', 'max']:

            # Variables (make a copy)
            vt = copy.deepcopy(v)
            del lp_prob

            # Objective function
            if obj == 'min':
                lp_prob = pulp.LpProblem("FVA%s"%obj, pulp.LpMinimize)
                lp_prob += vt[j1], "FVA_min"
            elif obj == 'max':
                lp_prob = pulp.LpProblem("FVA%s"%obj, pulp.LpMaximize)
                lp_prob += vt[j1], "FVA_max"

            # Constraints
            # Mass_balance
            for i in database.metabolites:
                # If metabolites not involve in any reactions
                if i not in database.S:
                    continue
                label = "mass_balance_%s" % i
                dot_S_v = pulp.lpSum([database.S[i][j] * vt[j]
                                      for j in database.S[i].keys()])
                condition = dot_S_v == 0
                lp_prob += condition, label

            if custom_flux_constraints is not None:
                logger.info("Adding custom constraints...")

                for group in custom_flux_constraints:
                    lp_prob += pulp.lpSum(vt[rxn] for rxn in group['reactions']) <= group['UB'], "%s_UB"%group['constraint_name']
                    lp_prob += pulp.lpSum(vt[rxn] for rxn in group['reactions']) >= group['LB'], "%s_LB"%group['constraint_name']

            lp_prob.solve(solver=pulp_solver)

            FVA_res[j1][obj] = pulp.value(lp_prob.objective)

        if (FVA_res[j1]['max'] < EPS) and (FVA_res[j1]['min'] > -EPS):
            blocked_reactions.append(j1)

        json.dump(FVA_res,
              open("temp_FVA_result.json", 'w+'),
              sort_keys=True,
              indent=4)

    return blocked_reactions, FVA_res


def internal_loop_analysis():
    """
    M = Matrix([[16, 2, 3,13],
    [5,11,10, 8],
    [9, 7, 6,12],
    [4,14,15, 1]])

    print(nsimplify(M, rational=True).nullspace())
    """
    pass


def test_blocked_reactions_analysis():

    logger = create_logger(name="Test blocked_reactions_analysis")

    user_defined_export_rxns_Sji = {
        'EX_glc': {'C00031': -1.0},
        'EX_nad': {'C00003': -1.0},
        'EX_adp': {'C00008': -1.0},
        'EX_phosphate': {'C00009': -1.0},
        'EX_pyruvate': {'C00022': -1.0},
        'EX_nadh': {'C00004': -1.0},
        'EX_atp': {'C00002': -1.0},
        'EX_h2o': {'C00001': -1.0},
        'EX_hplus': {'C00080': -1.0},
        'EX_nadp': {'C00006': -1.0},
        'EX_nadph': {'C00005': -1.0}
        }

    custom_flux_constraints = [
        {'constraint_name': 'nadphcons1',
         'reactions': ['EX_nadph', 'EX_nadh'],
         'UB': 2,
         'LB': 2},
        {'constraint_name': 'nadphcons2',
        'reactions': ['EX_nadp', 'EX_nad'],
        'UB': -2,
        'LB': -2},
        {'constraint_name': 'nadphcons3',
        'reactions': ['EX_nadh', 'EX_nad'],
        'UB': 0,
        'LB': 0},
        {'constraint_name': 'nadphcons4',
        'reactions': ['EX_nadph', 'EX_nadp'],
        'UB': 0,
        'LB': 0}]

    specific_bounds = {'EX_glc': {'LB': -1, 'UB': -1},
                    'EX_pyruvate': {'LB': 2, 'UB': 2},
                    'EX_nad': {'LB': -2, 'UB': 0},
                    'EX_nadh': {'LB': 0, 'UB': 2},
                    'EX_nadp': {'LB': -2, 'UB': 0},
                    'EX_nadph': {'LB': 0, 'UB': 2},
                    'EX_adp': {'LB': -5, 'UB': -1}, # range 1-5 ATP
                    'EX_phosphate': {'LB': -5, 'UB': -1},
                    'EX_atp': {'LB': 1, 'UB': 5}, # range 1-5 ATP
                    'EX_h2o': {'LB': 1, 'UB': 5},
                    'EX_hplus': {'LB': -10, 'UB': 10}} #pulp/gurobi has issue with "h+"

    pulp_solver = load_pulp_solver(
        solver_names=['SCIP_CMD', 'GUROBI', 'GUROBI_CMD', 'CPLEX_CMD', 'GLPK_CMD'],
        logger=logger)

    exclude_reactions = load_custom_reactions_to_be_excluded()

    db = load_base_reaction_db(
        user_defined_export_rxns_Sji=user_defined_export_rxns_Sji
        )

    blocked_reactions_list, FVA_res = blocked_reactions_analysis(
        database=db,
        pulp_solver=pulp_solver,
        user_defined_export_rxns_Sji=user_defined_export_rxns_Sji,
        specific_bounds=specific_bounds,
        custom_flux_constraints=custom_flux_constraints,        
        excluded_reactions=exclude_reactions,
        target_reactions_list=['R01266', 'R07882', 'R00658', 'R01059'])

    assert_equal(set(blocked_reactions_list), set(['R01266', 'R07882']))

    return blocked_reactions_list, FVA_res


if __name__ == '__main__':
    blocked_reactions_list, FVA_res = test_find_blocked_reactions()
