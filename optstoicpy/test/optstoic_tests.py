"""
pip install nose

nosetests optstoic_tests.py
"""
from nose.tools import (
    assert_equal,
    assert_in,
    assert_not_equal
    )
from optstoicpy.core.database import (
    load_db_v3,
    Database)
from optstoicpy.script.optstoic_glycolysis import (
    load_pulp_solver,
    OptStoic,
    DATA_DIR)

class TestOptStoic:
    def setUp(self):
        self.DB = self.test_load_database()
        self.pulp_solver = self.test_pulp_solver_loading()

    def test_pulp_solver_loading(self):
        pulp_solver = load_pulp_solver(
            solver_names=['GLPK_CMD'])
        assert_not_equal(pulp_solver, None)

    def test_load_database(self):
        DB = load_db_v3()
        assert_equal(len(DB.metabolites), 5969)
        assert_equal(len(DB.reactions), 7175)

        return DB

    def test_optstoic_setup(self):
        test = OptStoic(database=self.DB,
                        objective='MinFlux',
                        nATP=1,
                        zlb=10,
                        max_iteration=1,
                        pulp_solver=self.pulp_solver,
                        result_filepath=None,
                        M=1000)

