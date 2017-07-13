import sys
sys.path.append('../')
from optstoicpy.core.reaction import *
from optstoicpy.core.drawpathway import *
from optstoicpy.core.pathway import *
import pandas as pd


#Read the "test" worksheet in the excel file "any_reaction_list.xlsx"
df = pd.read_excel("any_reaction_list.xlsx", sheetname="test")

#This create a "Pathway" instance
test_pathway = Pathway(name='EMP',
                    reaction_ids=df.reaction_id.tolist(),
                    fluxes=df.flux.tolist())

#Create png image
draw_pathway(test_pathway, imageFileName='test', imageFormat='png', graphTitle=test_pathway.name,
    cleanup=True, darkBackgroundMode=False)

#Create svg image
draw_pathway(test_pathway, imageFileName='test', imageFormat='svg', graphTitle=test_pathway.name,
    cleanup=True, darkBackgroundMode=False)

#-----------------------------------------------------------------------
###alternately, if you do not want to read from an excel file
###You need to have two list (reaction IDs from KEGG and fluxes value)
###See the excel file in optstoicpy/data/optstoic_db_v3 for list of reactions and equations

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

test_pathway2 = Pathway(name='EMP_test',
                    reaction_ids=reaction_list,
                    fluxes=fluxes)

#Create png image
draw_pathway(test_pathway2, imageFileName='EMP_test', imageFormat='png', graphTitle=test_pathway.name,
    cleanup=True, darkBackgroundMode=False)
