"""
pip install nose

nosetests optstoic_tests.py
"""
from nose.tools import (
    assert_equal,
    assert_in,
    assert_not_equal,
    nottest,
    set_trace
    )
from optstoicpy.core.database import (
    load_db_v3,
    Database)
from optstoicpy.script.utils import create_logger
from optstoicpy.script.solver import load_pulp_solver
import optstoicpy.script.optstoic_glycolysis as osgly
import optstoicpy.script.optstoic as optstoic

class TestOptStoic:
    def setUp(self):
        self.DB = self.test_load_database()
        self.pulp_solver = self.test_pulp_solver_loading()

    def test_pulp_solver_loading(self):
        pulp_solver = load_pulp_solver(
            solver_names=['SCIP_CMD', 'GUROBI_CMD', 'GLPK_CMD'])
        assert_not_equal(pulp_solver, None)
        return pulp_solver

    def test_load_database(self):
        DB = load_db_v3()
        assert_equal(len(DB.metabolites), 5969)
        assert_equal(len(DB.reactions), 7175)

        return DB

    @nottest
    def test_optstoic_setup(self):
        model = osgly.OptStoic(database=self.DB,
                        objective='MinFlux',
                        nATP=1,
                        zlb=10,
                        max_iteration=1,
                        pulp_solver=self.pulp_solver,
                        result_filepath=None,
                        M=1000)

        lp_prob, pathways = model.solve(
            outputfile='test_optstoic.txt')

    @nottest
    def test_general_optstoic_setup(self):

        logger = create_logger(name='Test generalized optstoic ')

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
                        'EX_adp': {'LB': -1, 'UB': -1},
                        'EX_phosphate': {'LB': -1, 'UB': -1},
                        'EX_atp': {'LB': 1, 'UB': 1},
                        'EX_h2o': {'LB': 1, 'UB': 1},
                        'EX_hplus': {'LB': -10, 'UB': 10}} #pulp/gurobi has issue with "h+"

        model = optstoic.OptStoic(database=self.DB,
                        objective='MinFlux',
                        zlb=None,
                        specific_bounds=specific_bounds,
                        custom_flux_constraints=custom_flux_constraints,
                        add_loopless_constraints=False,
                        max_iteration=1,
                        pulp_solver=self.pulp_solver,
                        result_filepath='./result/',
                        M=1000,
                        logger=logger)

        lp_prob, pathways = model.solve(
            outputfile='test_optstoic_general.txt')

        assert_equal(pathways[1]['modelstat'], 'Optimal')
