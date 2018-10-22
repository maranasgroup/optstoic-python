#/usr/bin/python
"""
Loopless OptStoic program to identify glycolytic pathway
(glucose to pyruvate) for n ATP production.
It read input files that are used for GAMS.
Currently, it has been tested with SCIP, GLPK, Gurobi and CPLEX solvers.

Author: Chiam Yu

To-do list:
2. generalize the input for designing other pathways
3. Fix MinRxn
4. Add feature to run MinRxn for a range of zlb

Known issues:
1. MinRxn takes a very long time to solve (gap = 60% after 20 mins)
2. Pulp is having issue with gurobi700

Tip:
Change the number of Thread for gurobi if needed.

"""
import os
import time
import sys
import copy
import random
import string  # to generate random hex code
import pulp
# import cPickle as pickle
import gams_parser
import json
#import pdb
from optstoicpy.core import database
from optstoicpy.script.utils import create_logger
from optstoicpy.script.solver import load_pulp_solver
from gurobi_command_line_solver import *

# Global variables/solver options
EPS = 1e-5
GUROBI_OPTIONS = 'Threads=2 TimeLimit=1800 MIPGapAbs=1e-6 MIPGap=1e-6 CliqueCuts=2'


class OptStoic(object):
    """An OptStoic problem Class to identify glycolytic pathways
        using either minFlux or minRxn algorithm.
    """

    def __init__(self,
                 database,
                 objective='MinFlux',
                 nATP=1,
                 zlb=8,
                 max_iteration=2,
                 pulp_solver=None,
                 result_filepath=None,
                 M=1000,
                 logger=None):
        """
        Args:
            database (TYPE): An optStoic Database object (equivalent to GSM model)
            objective (str, optional): The mode for optStoic pathway prospecting.
                Options available are: ['MinFlux', 'MinRxn']
            nATP (int, optional): The number of ATP
            zlb (int, optional): The lowerbound on objective value z
            max_iteration (int, optional): The default maximum number of iteration
            pulp_solver (None, optional): A pulp.solvers object (load any of the user-defined solver)
            result_filepath (str, optional): Filepath for result
            M (int, optional): The maximum flux bound (default 1000)
            logger (:obj:`logging.Logger`, optional): A logging.Logger object
        """
        if logger is None:
            self.logger = create_logger(name='optstoic.OptStoic')
        else:
            self.logger = logger

        self.objective = objective
        self.nATP = nATP
        self.zlb = zlb
        self.M = M

        # When nATP is not integer, change variables v, vf and vb to continuous variables
        if float(nATP).is_integer():
            self._varCat = 'Integer'
        else:
            self._varCat = 'Continuous'

        self.max_iteration = max_iteration

        if result_filepath is None:
            result_filepath = './result'
        self.result_filepath = result_filepath

        if not os.path.exists(result_filepath):
            os.makedirs(result_filepath)

        self.database = database
        self.pathways = {}
        self.iteration = 1
        self.lp_prob = None
        self.pulp_solver = pulp_solver
        self.lp_prob_fname = "OptStoic_{0}".format(self.generate_random_string(6))

    @staticmethod
    def generate_random_string(N):
        """LP file is appended with random string when using command line mode
            to prevent overwrite/read issues.
        """
        return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(N))

    def change_objective(self, new_objective):
        if new_objective not in ['MinFlux', 'MinRxn']:
            raise ValueError("The objective for OptStoic is "
                             "not correctly defined. "
                             "Please use either 'MinFlux' or 'MinRxn'.")
        self.objective = new_objective

    def change_nATP(self, new_nATP):
        self.nATP = new_nATP

    def create_minflux_problem(self):
        """
        Create minflux/minRxn LP problem (a pulp LpProblem object)
        minimize z = sum(j\j_exchange, vf(j) + vb(j))
        subject to:
            z = zlb
            v(j) = vf(j) - vb(j)
            sum(j, S(i,j) * v(j)) = 0, for all i
            vf(j) >= yf(j) * eps, for all j
            vf(j) <= yf(j) * M, for all j
            vb(j) >= yb(j) * eps, for all j
            vb(j) <= yb(j) * M, for all j
            yf(j) + yb(j) <= 1, for all j

            for reaction in jloop (not in jblock)
            sum(j, Nint(l,j) * G(j)) = 0, for all l
            G(j) >= -M * a(j) + (1 - a(j))
            G(j) <= -a(j) + M * (1 - a(j))
            v(j) >= -M * (1 - a(j))
            v(j) <= M * a(j)

            v('EX_nadph') + v('EX_nadh') = 2;
            v('EX_nadp') + v('EX_nad') = -2;
            v('EX_nadh') + v('EX_nad') = 0;
            v('EX_nadph') + v('EX_nadp') = 0;
        """
        self.logger.info("Formulating problem...")
        # Scalar
        M = self.M

        # Initialize variables
        v = pulp.LpVariable.dicts("v", self.database.reactions,
                                  lowBound=-M, upBound=M, cat=self._varCat)
        vf = pulp.LpVariable.dicts("vf", self.database.reactions,
                                   lowBound=0, upBound=M, cat=self._varCat)
        vb = pulp.LpVariable.dicts("vb", self.database.reactions,
                                   lowBound=0, upBound=M, cat=self._varCat)
        yf = pulp.LpVariable.dicts("yf", self.database.reactions,
                                   lowBound=0, upBound=1, cat='Binary')
        yb = pulp.LpVariable.dicts("yb", self.database.reactions,
                                   lowBound=0, upBound=1, cat='Binary')
        a = pulp.LpVariable.dicts("a", self.database.reactions,
                                  lowBound=0, upBound=1, cat='Binary')
        G = pulp.LpVariable.dicts("G", self.database.reactions,
                                  lowBound=-M, upBound=M, cat='Continuous')

        for j in self.database.reactions:

            if j in self.database.all_excluded_reactions:
                v[j].lowBound = 0
                v[j].upBound = 0
                vf[j].lowBound = 0
                vf[j].upBound = 0
                yf[j].lowBound = 0
                yf[j].upBound = 0
                yb[j].lowBound = 0
                yb[j].upBound = 0

                continue

            if self.database.rxntype[j] == 0:
                # Forward irreversible
                v[j].lowBound = 0
                v[j].upBound = M
                yb[j].upBound = 0
                vb[j].upBound = 0

            elif self.database.rxntype[j] == 1:
                # Reversible
                v[j].lowBound = -M
                v[j].upBound = M

            elif self.database.rxntype[j] == 2:
                # Reverse irreversible
                v[j].lowBound = -M
                v[j].upBound = 0
                vf[j].upBound = 0
                yf[j].upBound = 0

            elif self.database.rxntype[j] == 4:
                v[j].lowBound = 0
                v[j].upBound = 0

        # Fix stoichiometry of source/sink metabolites
        # Allow user to change this for generalization
        v['EX_glc'].lowBound = -1
        v['EX_glc'].upBound = -1
        v['EX_pyruvate'].lowBound = 2
        v['EX_pyruvate'].upBound = 2
        v['EX_nad'].lowBound = -2
        v['EX_nad'].upBound = 0
        v['EX_nadh'].lowBound = 0
        v['EX_nadh'].upBound = 2
        v['EX_nadp'].lowBound = -2
        v['EX_nadp'].upBound = 0
        v['EX_nadph'].lowBound = 0
        v['EX_nadph'].upBound = 2
        v['EX_adp'].lowBound = -self.nATP
        v['EX_adp'].upBound = -self.nATP
        v['EX_phosphate'].lowBound = -self.nATP
        v['EX_phosphate'].upBound = -self.nATP
        v['EX_atp'].lowBound = self.nATP
        v['EX_atp'].upBound = self.nATP
        v['EX_h2o'].lowBound = self.nATP
        v['EX_h2o'].upBound = self.nATP
        v['EX_hplus'].lowBound = -10
        v['EX_hplus'].upBound = 10

        LB = {}
        UB = {}
        for j in self.database.reactions:
            LB[j] = v[j].lowBound
            UB[j] = v[j].upBound

        lp_prob = pulp.LpProblem("OptStoic", pulp.LpMinimize)

        # Min-Rxn objective
        if self.objective == 'MinRxn':
            condition = pulp.lpSum([yf[j] + yb[j]
                                    for j in self.database.reactions
                                    if self.database.rxntype[j] != 4])
            lp_prob += condition, "MinRxn"

        # Min-Flux objective
        elif self.objective == 'MinFlux':
            condition = pulp.lpSum([vf[j] + vb[j]
                                    for j in self.database.reactions
                                    if self.database.rxntype[j] != 4])
            lp_prob += condition, "MinFlux"
            #fix lower bound
            lp_prob += condition == self.zlb, 'zLowerBound'

        # Constraints
        # Mass_balance
        for i in self.database.metabolites:
            # If metabolites not involve in any reactions
            if i not in self.database.S:
                continue
            label = "mass_balance_%s" % i
            dot_S_v = pulp.lpSum([self.database.S[i][j] * v[j]
                                  for j in self.database.S[i].keys()])
            condition = dot_S_v == 0
            lp_prob += condition, label

        # if self.objective == 'MinRxn':
        # for j in self.database.reactions:
        #     lp_prob += v[j] >= y[j]*LB[j], "cons1_%s"%j
        #     lp_prob += v[j] >= y[j]*LB[j], "cons1_%s"%j
        #     lp_prob += v[j] <= y[j]*UB[j], "cons2_%s"%j

        if self.objective == 'MinFlux':
            for j in self.database.reactions:
                lp_prob += (v[j] == vf[j] - vb[j]), "flux_%s" % j

                # These constraints ensure that when yf=0 and yb=0 ,
                # no flux goes through the reaction
                lp_prob += vf[j] >= yf[j] * 0.5, "cons1_%s" % j
                lp_prob += vf[j] <= yf[j] * M, "cons2_%s" % j
                lp_prob += vb[j] >= yb[j] * 0.5, "cons3_%s" % j
                lp_prob += vb[j] <= yb[j] * M, "cons4_%s" % j
                # Ensure that either yf or yb can be 1, not both
                lp_prob += yf[j] + yb[j] <= 1, 'cons5_%s' % j

        loop_rxn = list(set(self.database.internal_rxns) - set(self.database.blocked_rxns))

        # Loopless contraints
        for l in self.database.loops:
            label = "loopless_cons_%s" % l
            dot_N_G = pulp.lpSum([self.database.Ninternal[l][j] * G[j]
                                  for j in self.database.Ninternal[l].keys()])
            condition = dot_N_G == 0
            lp_prob += condition, label

        for j in loop_rxn:
            lp_prob += G[j] >= -M * a[j] + (1 - a[j]), "llcons1_%s" % j
            lp_prob += G[j] <= -a[j] + M * (1 - a[j]), "llcons2_%s" % j
            lp_prob += v[j] >= -M * (1 - a[j]), "llcons3_%s" % j
            lp_prob += v[j] <= M * a[j], "llcons4_%s" % j

        #lp_prob += v['R00756'] == 1, 'enforce_pfk'

        #Fix nad(p)h production and consumption
        lp_prob += v['EX_nadph'] + v['EX_nadh'] == 2, 'nadphcons1'
        lp_prob += v['EX_nadp'] + v['EX_nad'] == -2, 'nadphcons2'
        lp_prob += v['EX_nadh'] + v['EX_nad'] == 0, 'nadphcons3'
        lp_prob += v['EX_nadph'] + v['EX_nadp'] == 0, 'nadphcons4'

        return lp_prob, v, vf, vb, yf, yb, a, G

    def solve(self, exclude_existing_solution=False, outputfile="OptStoic_pulp_result.txt", max_iteration=None):
        """
        Solve OptStoic problem using pulp.solvers interface

        Keyword Arguments:
            outputfile: name of outpufile
            max_iteration: Externally specified maximum number of pathway to be found using OpStoic.
                            If not specified, it will set to the internal max iterations.
        """
        if self.objective not in ['MinFlux', 'MinRxn']:
            raise ValueError("The objective for OptStoic is not correctly defined. Please use either 'MinFlux' or 'MinRxn'.")

        if max_iteration is None:
            max_iteration = self.max_iteration

        self.logger.info("Finding multiple pathways using Optstoic %s...", self.objective)
        lp_prob, v, vf, vb, yf, yb, a, G = self.create_minflux_problem()

        # Create integer cut for existing pathways
        if exclude_existing_solution and bool(self.pathways):
            self.iteration = max(self.pathways.keys()) + 1
            if self.iteration > max_iteration:
                raise ValueError('Max iteration is less than current '
                                 'iteration. Increase max_iteration '
                                 'before solving!')

            for ind, entry in self.pathways.iteritems():
                rxnlist = list(set(entry['reaction_id']) -
                               set(self.database.user_defined_export_rxns))
                condition = pulp.lpSum(
                    [(1 - yf[j] - yb[j]) for j in rxnlist]) >= 1
                lp_prob += condition, "IntegerCut_%d" % ind


        self.logger.info("Solving problem...")
        if self.iteration == 1:
            result_output = open(os.path.join(self.result_filepath, outputfile),"w+")
        else:
            result_output = open(os.path.join(self.result_filepath, outputfile),"a+")

        while True and self.iteration <= max_iteration:
            self.logger.info("Iteration %s", self.iteration)
            #lp_prob.writeLP("OptStoic.lp", mip=1)  #optional
            e1 = time.time()
            lp_prob.solve(self.pulp_solver)
            e2 = time.time()
            self.logger.info("This iteration solved in %.3f seconds.", (e2-e1))

            # The solution is printed if it was deemed "optimal
            if pulp.LpStatus[lp_prob.status] == "Optimal":
                self.logger.info("Writing result to output file...")
                result_output.write("\nIteration no.: %d\n" %self.iteration)
                result_output.write("\nModelstat: %s\n" %pulp.LpStatus[lp_prob.status])

                res = {}
                res['reaction_id'] = []
                res['flux'] = []
                res['iteration'] = self.iteration
                res['time'] = (e2-e1)
                res['modelstat'] = "Optimal"

                for j in self.database.reactions:
                    if v[j].varValue is not None:
                        if v[j].varValue > EPS or v[j].varValue < -EPS:
                            res['reaction_id'].append(j)
                            res['flux'].append(v[j].varValue)
                            result_output.write("%s %.8f\n" %(v[j].name, v[j].varValue))

                result_output.write("%s = %.8f\n" % (self.objective, pulp.value(lp_prob.objective)))
                result_output.write("----------------------------------\n\n")

                integer_cut_reactions = list(set(res['reaction_id']) - set(self.database.user_defined_export_rxns))

                self.pathways[self.iteration] = res
                json.dump(self.pathways, open(os.path.join(self.result_filepath,'temp_pathways.json'),'w+'), sort_keys=True, indent=4)

                # Integer cut constraint is added so that
                # the same solution cannot be returned again
                condition = pulp.lpSum([(1 - yf[j] - yb[j])
                                        for j in integer_cut_reactions]) >= 1
                lp_prob += condition, "IntegerCut_%d" % self.iteration
                self.iteration += 1

            # If a new optimal solution cannot be found, end the program
            else:
                break
        result_output.close()

        self.lp_prob = lp_prob

        return self.lp_prob, self.pathways

    def add_existing_pathways(self, user_defined_pathways):
        """
        Add list of existing solutions (pathways) to be
        excluded from being identified.

        Keyword Arguments:
        user_defined_pathways -- pathways output from solve_gurobi_cl()
                                 or solve() or pathway in dictionary format
                                 e.g. {1: 'reaction_id': []}
        """
        if (isinstance(user_defined_pathways, dict) and
           ('reaction_id' in user_defined_pathways.values()[0])):
            self.pathways = copy.deepcopy(user_defined_pathways)
        else:
            raise ValueError("user_defined_pathways must be a "
                             "pathways dictionary "
                             "{1: 'reaction_id': ['R00001', 'R00002']}")

    def reset_pathways(self):
        """
        Reset self.pathways to empty dictionary
        """
        self.pathways = {}

    def solve_gurobi_cl(self, exclude_existing_solution=False,
                        outputfile="OptStoic_pulp_result_gcl.txt",
                        max_iteration=None, cleanup=True,
                        gurobi_options=GUROBI_OPTIONS):
        """
        Solve OptStoic problem using Gurobi command line (gurobi_cl)
        when pulp.solvers.GUROBI_CMD failed.
        Require the module "gurobi_command_line_solver.py".

        Keyword Arguments:
        exclude_existing_solution -- If true and if self.pathway is not None,
                                    exclude the pathways from being identified.
        outputfile -- name of outpufile
        max_iteration -- Externally specified maximum number of pathway
                         to be found using OpStoic. If not specified,
                         it will set to the internal max iterations.
        """
        if self.objective not in ['MinFlux', 'MinRxn']:
            raise ValueError("The objective for OptStoic is not correctly "
                             "defined. Please use either 'MinFlux' or "
                             "'MinRxn'.")

        if max_iteration is None:
            max_iteration = self.max_iteration

        t1 = time.time()

        self.logger.info("Finding multiple pathways using"
                     " Optstoic %s and Gurobi CL...", self.objective)
        lp_prob, v, vf, vb, yf, yb, a, G = self.create_minflux_problem()

        # Create integer cut for existing pathways
        if exclude_existing_solution and bool(self.pathways):
            self.iteration = max(self.pathways.keys()) + 1
            if self.iteration > max_iteration:
                raise ValueError('Max iteration is less than current '
                                 'iteration. Increase max_iteration '
                                 'before solving!')

            for ind, entry in self.pathways.iteritems():
                rxnlist = list(set(entry['reaction_id']) -
                               set(self.database.user_defined_export_rxns))
                condition = pulp.lpSum(
                    [(1 - yf[j] - yb[j]) for j in rxnlist]) >= 1
                lp_prob += condition, "IntegerCut_%d" % ind

        # Solve problem
        self.logger.info("Solving problem...")

        if self.iteration == 1:
            result_output = open(os.path.join(
                self.result_filepath, outputfile), "w+")
        else:
            result_output = open(os.path.join(
                self.result_filepath, outputfile), "a+")

        while True and self.iteration <= max_iteration:
            self.logger.info("Iteration %s", self.iteration)
            lp_prob.writeLP(self.lp_prob_fname + ".lp", mip=1)
            e1 = time.time()
            lp_status, solver_message = solve_with_gurobi_cl_debug(
                self.lp_prob_fname, options=gurobi_options)
            e2 = time.time()
            self.logger.info("This iteration solved in %.3f seconds.", (e2 - e1))

            # The solution is printed if it was deemed "optimal
            if lp_status in ["Optimal", "Time_limit"]:
                objective_function, varValue = parse_gurobi_sol(self.lp_prob_fname)

                res = {}
                res['reaction_id'] = []
                res['flux'] = []
                res['iteration'] = self.iteration
                res['time'] = (e2 - e1)
                res['modelstat'] = lp_status
                res['solvestat'] = solver_message

                result_output.write("\nIteration no.: %d\n" %self.iteration)
                result_output.write("\nModelstat: %s\n" %lp_status)

                for j in self.database.reactions:
                    if 'v_'+j in varValue:
                        v = varValue['v_'+j]
                        if v > EPS or v < -EPS :
                            res['reaction_id'].append(j)
                            res['flux'].append(v)
                            result_output.write("%s %.8f\n" %(j, v))

                result_output.write("%s = %.8f\n" %(self.objective, objective_function))
                result_output.write("----------------------------------\n\n")

                integer_cut_reactions = list(set(res['reaction_id']) - set(self.database.user_defined_export_rxns))

                self.pathways[self.iteration] = res
                # Keep a copy of pathways in case program terminate midway
                json.dump(self.pathways, open(os.path.join(
                    self.result_filepath, 'temp_pathways.json'), 'w+'),
                    sort_keys=True, indent=4)

                # Integer cut constraint is added so that
                # the same solution cannot be returned again
                condition = pulp.lpSum([(1 - yf[j] - yb[j])
                                        for j in integer_cut_reactions]) >= 1
                lp_prob += condition, "IntegerCut_%d" % self.iteration
                self.iteration += 1

            # If a new optimal solution cannot be found, end the program
            else:
                break

        result_output.close()
        # Clean up directory
        if cleanup:
            self.logger.debug("Cleaning up directory...")
            os.remove("./" + self.lp_prob_fname + ".lp")
            os.remove("./" + self.lp_prob_fname + ".sol")
            os.remove("./gurobi.log")

        self.lp_prob = lp_prob

        return self.lp_prob, self.pathways

    def __repr__(self):
        return "<OptStoic(nATP='%s', objective='%s')>" % (self.nATP, self.objective)



if __name__ == '__main__':

    logger = create_logger(name='script.optstoic.main')

    db3 = database.load_db_v3()

    #logger.debug('Testing optstoic output filepath: %s', res_dir)

    pulp_solver = load_pulp_solver(
        solver_names=['GLPK_CMD', 'GUROBI', 'GUROBI_CMD', 'CPLEX_CMD'],
        logger=logger)

    test = OptStoic(database=db3,
                    objective='MinFlux',
                    nATP=1,
                    zlb=10,
                    max_iteration=1,
                    pulp_solver=pulp_solver,
                    result_filepath='./result/',
                    M=1000,
                    logger=logger)

    if sys.platform == 'cygwin':
        lp_prob, pathways = test.solve_gurobi_cl(outputfile='test_optstoic_cyg.txt', cleanup=False)
        #test.max_iteration = test.max_iteration + 2
        #lp_prob, pathways = test.solve_gurobi_cl(outputfile='test_optstoic_cyg.txt', exclude_existing_solution=True, cleanup=False)
    else:
        lp_prob, pathways = test.solve(outputfile='test_optstoic.txt')
        #test.max_iteration = test.max_iteration + 1
        #lp_prob, pathways = test.solve(outputfile='test_optstoic.txt', exclude_existing_solution=True)
