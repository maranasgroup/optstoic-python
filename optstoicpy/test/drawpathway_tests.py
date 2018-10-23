"""
pip install nose

nosetests drawpathway_tests.py
"""
import os
from nose.tools import (
    assert_equal,
    assert_in,
    assert_not_equal,
    nottest
    )
from optstoicpy.core.pathway import Pathway
from optstoicpy.core.drawpathway import (
    draw_pathway)


class PathwayTest:
    def setup(self):
        self.logger = create_logger(name='Test core.Pathway')
        self.pathway_fixture = {'flux': [-1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 1.0,
                                1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 2.0, 1.0,
                                1.0, 1.0, -1.0, 1.0],
                                   'iteration': 1,
                                   'reaction_id': ['R00200', 'R00300', 'R00658', 'R01059',
                                                   'R01063', 'R01512', 'R01518', 'R01519',
                                                   'R01538', 'R08570', 'EX_glc', 'EX_nad',
                                                   'EX_adp', 'EX_phosphate', 'EX_pyruvate',
                                                   'EX_nadh', 'EX_atp', 'EX_h2o', 'EX_nadp',
                                                   'EX_nadph']}

        self.p1 = Pathway(id=1, 
                     name='OptStoic',
                     reaction_ids=self.pathway_fixture['reaction_id'],
                     fluxes=self.pathway_fixture['flux'])

    @nottest
    def test_rearrange_pathway(self):
        self.logger.info("Test rearranging reaction order")
        p1.rearrange_reaction_order()

    def test_kegg_model_generation(self):
        
        self.logger.info("Creating 'res' folder in the current directory if not exist...")
        # outputFilepath = 'res'
        # outputFilename = 'OptStoic'
        # try:
        #     os.makedirs(outputFilepath)
        # except OSError:
        #     if not os.path.isdir(outputFilepath):
        #         raise Exception

        self.logger.info("Test create KEGG model file")

        filename = "./test_kegg_model_generation.txt"
        f = open(filename, 'a+')
        kegg_model_text = generate_kegg_model(self.p1, filehandle=f)
        print kegg_model_text
        assert_in('R01512', kegg_model_text)
        assert_in('R01512', kegg_model_text)
        f.close()

        if os.path.exists(filename):
            os.remove(filename)


class TestDrawPathway:
    def draw_pathway(self):
        # Sample input
        reaction_list = [
         'R00200',
         'R00217',
         'R00299',
         'R00346',
         'R00658',
         'R00764',
         'R00771',
         'R01015',
         'R01061',
         'R01068',
         'R01512',
         'R01514',
         'R01748'
         ]

        fluxes = [-1, 1, 1, -1, 2, 1, 1, -1, 2, 1, -2, -2, -2]

        test_pathway = Pathway(name='EMP_pathway',
                            reaction_ids=reaction_list,
                            fluxes=fluxes)

        # Create png image
        draw_pathway(test_pathway,
            imageFileName='test_EMP_pathway',
            imageFormat='png',
            graphTitle=test_pathway.name,
            cleanup=True,
            darkBackgroundMode=False)

        fname = 'test_EMP_pathway.png'

        assert_equal(os.path.exists(fname), True)

        if os.path.exists(fname):
            os.remove(fname)