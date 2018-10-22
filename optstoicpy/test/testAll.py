import sys, os, logging
current_dir =  os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.normpath(os.path.join(current_dir,'../'))
res_dir = os.path.join(current_dir,'test_results/')
sys.path.append(parent_dir)
from optstoicpy.script.optstoic_glycolysis import *
from optstoicpy.core import database
from optstoicpy.core.reaction import Reaction
from optstoicpy.core.pathway import Pathway, generate_kegg_model
from optstoicpy.core import drawpathway

logging.basicConfig(
        format='%(asctime)s %(levelname)s %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG
    )

if __name__ =="__main__":

    if not os.path.exists(res_dir):
         os.makedirs(res_dir)

    #load database version 3
    db3 = database.load_db_v3()

    #create OptStoic object
    test = OptStoic(db3, objective='MinFlux', nATP=1,
                    zlb=9,
                    max_iteration=2,
                    pulp_solver=pulp_solver,
                    data_filepath=data_dir,
                    result_filepath=res_dir,
                    M=1000)

    if sys.platform == 'cygwin':
        lp_prob, pathways = test.solve_gurobi_cl(outputfile='test_optstoic.txt')
        test.max_iteration = test.max_iteration + 2
        lp_prob, pathways = test.solve_gurobi_cl(outputfile='test_optstoic.txt', exclude_existing_solution=True)
        print "Running OptStoic in cygwin using gurobi_cl: Pass!"
    else:
        lp_prob, pathways = test.solve(outputfile='test_optstoic.txt')
        test.max_iteration = test.max_iteration + 2
        lp_prob, pathways = test.solve(outputfile='test_optstoic.txt', exclude_existing_solution=True)
        print "Running OptStoic using pulp solvers: Pass!"

    # Creating kegg model and drawing pathways.
    f = open(os.path.join(res_dir, 'test_KeggModel.txt'), 'w+')
    pathway_objects = []

    for ind, res in sorted(test.pathways.iteritems()):
        p = Pathway(id=ind, name='OptStoic', reaction_ids=res['reaction_id'], fluxes=res['flux'])
        p.rearrange_reaction_order()
        pathway_objects.append(p)
        generate_kegg_model(p, filehandle=f)
        graph_title = "{0}_{1}ATP_P{2}".format(p.name, p.nATP, p.id)
        drawpathway.draw_pathway(p, os.path.join(res_dir+'/pathway_{0:03d}'.format(p.id)),
                    imageFormat='png', graphTitle=graph_title)
    f.close()
    print "Generate kegg_model and draw pathway: Pass!"