"""
pip install nose

nosetests drawpathway_tests.py
"""
import os
from nose.tools import (
    assert_equal,
    assert_in,
    assert_not_equal
    )
from optstoicpy.core.pathway import Pathway
from optstoicpy.core.drawpathway import (
    draw_pathway)



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



