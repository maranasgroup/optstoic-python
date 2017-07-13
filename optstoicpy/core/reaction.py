# Author: Chiam Yu
# Inspired by the Reaction class of cobrapy package
# Reaction class
from .config import rxnSji


class Reaction(object):
    """Reaction class"""

    def __init__(self, rid=None, flux=1, metabolites={},
                 equation='', reversible=True):
        self.rid = rid
        self.flux = flux
        self.metabolites = metabolites
        self.equation = equation
        self.reversible = reversible

    def get_rid(self):
        return self.rid

    def get_metabolites(self):
        return self.metabolites

    def get_equation(self):
        return self.equation

    def autoset_metabolites(self):
        if len(self.metabolites) > 0:
            print "WARNING: Metabolites exists!"
        else:
            print "Retrieving metabolites from default database"
            self.metabolites = rxnSji[self.rid]
        return self.metabolites

    def flux(self):
        return self.flux

    @property
    def reactants(self):
        if self.flux > 0:
            return [k for k, v in self.metabolites.items() if v < 0]
        else:
            return [k for k, v in self.metabolites.items() if v > 0]

    @property
    def products(self):
        if self.flux > 0:
            return [k for k, v in self.metabolites.items() if v > 0]
        else:
            return [k for k, v in self.metabolites.items() if v < 0]

    def __str__(self):
        return "Reaction('%s')" % self.rid

    def __repr__(self):
        return "Reaction('%s')" % self.rid

    def set_equation(self):
        """Write equation in the direction of the flux"""
        if len(self.equation) != 0:
            print "WARNING: Equation exists!"
        else:
            if len(self.metabolites) == 0:
                print "Metabolites are not available!"
                print "Auto-updating metabolites..."
                self.autoset_metabolites()

            temp_list = []
            for cpd in sorted(self.reactants):
                coeff = abs(self.metabolites[cpd])
                if coeff == 1:
                    temp_list.append(cpd)
                else:
                    temp_list.append('%1.0f %s' % (coeff, cpd))

            eqStr = ' + '.join(temp_list)
            eqStr += ' <=> '

            temp_list = []
            for cpd in sorted(self.products):
                coeff = abs(self.metabolites[cpd])
                if coeff == 1:
                    temp_list.append(cpd)
                else:
                    temp_list.append('%1.0f %s' % (coeff, cpd))
            eqStr += ' + '.join(temp_list)
            self.equation = eqStr
        return self.equation

    @classmethod
    def create_Reaction_list_from_dict(cls, dataDict, excludeExchangeRxn=True):
        """
        Make a list of Reaction object from dataDict,
        excluding exchange reaction

        E.g.
        dataDict = {'reaction_id': ['R00001', 'R00002'], 'flux': [-1, 1]}
        output = [Reaction('R00001'), Reaction('R00002')]

        Keyword arguments:
        dataDict -- dictionary with reaction_id and flux
        excludeExchangeRxn -- Exclude all exchange reactions in the list
                              (default true)
        """
        RxnObjList = []
        for i in range(len(dataDict['reaction_id'])):
            if excludeExchangeRxn:
                if 'EX_' in dataDict['reaction_id'][i]:
                    continue
            tempRxn = cls(dataDict['reaction_id'][i], dataDict['flux'][i])
            # Get the metabolites dictionary {'C00001': -1, ...} for each
            # reaction
            tempRxn.metabolites = rxnSji[tempRxn.rid]
            RxnObjList.append(tempRxn)
        return RxnObjList