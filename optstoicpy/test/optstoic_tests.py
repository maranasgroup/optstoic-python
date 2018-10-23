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
    def setup(self):
        self.logger = create_logger(name='Test generalized optstoic')
        self.DB = self.test_load_database()
        self.pulp_solver = self.test_pulp_solver_loading()        

    def test_pulp_solver_loading(self):
        self.logger.info("Test loading PuLP solver(s)")

        pulp_solver = load_pulp_solver(
            solver_names=['SCIP_CMD', 'GUROBI', 'GUROBI_CMD', 'CPLEX_CMD', 'GLPK_CMD'])
        assert_not_equal(pulp_solver, None)
        return pulp_solver

    def test_load_database(self):
        self.logger.info("Test loading Database")

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

        DB = load_db_v3(
            reduce_model_size=True,
            user_defined_export_rxns_Sji=user_defined_export_rxns_Sji)
        DB.validate()
        # assert_equal(len(DB.metabolites), 5969)
        # assert_equal(len(DB.reactions), 7175)
        return DB

    @nottest
    def test_optstoic_setup(self):
        model = osgly.OptStoicGlycolysis(
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
                        logger=self.logger)

        lp_prob, pathways = model.solve(
            outputfile='test_optstoic_general.txt')

        assert_equal(pathways[1]['modelstat'], 'Optimal')
