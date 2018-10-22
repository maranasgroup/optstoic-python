import os
import json
import copy
import pandas as pd
from optstoicpy.script import gams_parser
from optstoicpy.script.utils import create_logger

current_dir = os.path.dirname(os.path.abspath(__file__))

data_dir = os.path.normpath(
    os.path.join(current_dir, '../data/', 'optstoic_db_v3')
)

REACTION_TYPE = {0: 'Forward irreversible', 1: 'Reversible',
                 2: 'Reverse irreverisble', 4: 'Export reaction'}


class Database(object):
    """optstoic Database class: loading database from GAMS input file"""

    def __init__(
        self,
        version='',
        data_filepath='../data/',
        dbdict={},
        logger=None):
        """Summary

        Args:
            version (str, optional): Description
            data_filepath (str, optional): Description
            dbdict (dict, optional): Description
            logger (None, optional): Description
        """
        if logger is None:
            self.logger = create_logger('core.Database')
        else:
            self.logger = logger

        self.version=version
        self.data_filepath = data_filepath
        self.dbdict = dbdict

        #initalize
        self.reactions = []
        self.metabolites = []
        self.S = {}
        self.Sji = {}
        self.rxntype = []
        self.loops = []
        self.Ninternal = {}
        self.all_excluded_reactions = []
        self.excluded_reactions = dbdict.get('excluded_reactions_list') or []
        self.user_defined_export_rxns = []
        self.blocked_rxns = []

    def load(self):

        # Load S matrix dictionary (S(i,j)) from json
        if 'S' not in self.dbdict:
            self.Sji = json.load(open(os.path.join(self.data_filepath,
                                                   self.dbdict['Sji']), 'r+'))

            self.S = self.transpose_S(self.Sji)
        else:
            try:
                self.S = json.load(open(os.path.join(self.data_filepath,
                                                     self.dbdict['S']), 'r+'))

            except:
                # TODO add function to differentiate json/txt input
                self.S = gams_parser.convert_parameter_table_to_dict(
                    os.path.join(self.data_filepath,
                                 '20160616_optstoic_Sij.txt')
                )
            self.Sji = self.transpose_S(self.S)

        # Load reactions
        self.logger.debug('Reading reaction file...')
        self.reactions = gams_parser.convert_set_to_list(
            os.path.join(self.data_filepath, self.dbdict['reaction'])
        )

        self.internal_rxns = copy.deepcopy(self.reactions)

        self.logger.debug('Reading metabolite file...')
        self.metabolites = gams_parser.convert_set_to_list(
            os.path.join(self.data_filepath, self.dbdict['metabolite'])
        )

        self.logger.debug('Reading blocked reactions file...')
        if 'blocked_rxns' in self.dbdict:
            self.blocked_rxns = gams_parser.convert_set_to_list(
                os.path.join(self.data_filepath, self.dbdict['blocked_rxns'])
            )

        self.all_excluded_reactions = list(
            set(self.excluded_reactions + self.blocked_rxns)
        )

        self.logger.debug('Reading reaction type file...')
        self.rxntype = gams_parser.convert_parameter_list_to_dict(
            os.path.join(self.data_filepath, self.dbdict['reactiontype']),
            datadict=None
        )

        self.logger.debug('Reading loop file...')
        self.loops = gams_parser.convert_set_to_list(
            os.path.join(self.data_filepath, self.dbdict['loops']))

        self.logger.debug('Reading Nint(loop, j) file...')
        self.Ninternal = gams_parser.convert_parameter_table_to_dict(
            os.path.join(self.data_filepath, self.dbdict['Nint']))

    @staticmethod
    def transpose_S(Sji):
        """Tranpose Sji into Sij and also Sij to Sji dictionary."""
        # Update to pandas 0.19 (using sparse dataframe)
        df_Sji = pd.DataFrame(Sji).T
        Sij = dict(
            (k, v.dropna().to_dict()) for k, v in pd.compat.iteritems(df_Sji)
        )
        return Sij

    @staticmethod
    def to_json(Sdict, filepath):
        with open(filepath, 'w+') as fp:
            json.dump(Sdict, fp, sort_keys=True, indent=4)

    def to_mat_file():
        """write to matlab file"""
        pass

    def get_reaction_type(self, rid, verbose=True):
        try:
            self.rxntype[rid]
        except:
            self.logger.warning("Reaction %s not in database!" % rid)
            return None
        else:
            if verbose:
                print "Reaction: {0} is ({1}) {2}".format(
                    rid, self.rxntype[rid], REACTION_TYPE[self.rxntype[rid]]
                )
            return self.rxntype[rid]

    def set_reaction_type(self, rid, rxntype):
        try:
            t0 = self.rxntype[rid]
            self.rxntype[r] = rxntype
        except KeyError:
            self.logger.warning('Reaction %s not in database!' % r)
        else:
            self.logger.info('Reaction %s has been updated from %s to %s.'
                         % (r, REACTION_TYPE[t0], REACTION_TYPE[rtype])
                         )

    def extend_S_from_file(self, filename='Sij_extension_for_glycolysis.txt'):
        self.S = gams_parser.convert_parameter_table_to_dict(
            os.path.join(self.data_filepath, filename), Sdict=self.S)

    def update_S(self, extension_dict):
        temp_rxn = []
        for met, entries in extension_dict.iteritems():
            if met not in self.S:
                self.S[met] = {}
                self.metabolites.append(met)
            for rxn, coeff in entries.iteritems():
                self.S[met][rxn] = float(coeff)
                if rxn not in self.reactions:
                    self.reactions.append(rxn)
                    temp_rxn.append(rxn)
        return self.S, temp_rxn

    def set_database_export_reaction(self, export_reactions_Sij_dict):
        _, temp_rxn = self.update_S(export_reactions_Sij_dict)
        if len(self.user_defined_export_rxns) != 0:
            self.logger.warning("Warning: The current list of export reactions\
                will be replaced! %s" % str(self.user_defined_export_rxns))
        self.user_defined_export_rxns = temp_rxn

        for rxn in self.user_defined_export_rxns:
            self.set_reaction_type(rxn, 4)

    def update_rxntype(self, new_reaction_type_dict):
        for (r, rtype) in new_reaction_type_dict.iteritems():
            self.set_reaction_type(r, rtype)
        return self.rxntype

    # def create_excluded_reactions_list(self, use_default=True, user_defined_reactions=None):
    #     """
    #     Create a list of reactions (from the reaction database) that user
    #     wants to exclude from OptStoic.

    #     Args:
    #     use_default: If True, use the default list, else user need to supply own list of reactions
    #     user_defined_reactions: list of reactions to be excluded (default is empty list)
    #                             load default list if no input list is given
    #     """
    #     if not use_default:
    #         assert isinstance(user_defined_reactions, list), 'user_defined_reactions should be a list!'
    #         self.all_excluded_reactions.extend(user_defined_reactions)

    #     else:
    #         self.all_excluded_reactions = self.database.all_excluded_reactions
    #         return self.all_excluded_reactions


    def __repr__(self):
        return "OptStoic Database(Version='%s')" % self.version


def load_db_v3(
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
    ):
    """Load OptStoic database v3

    Returns:
        TYPE: Description
    """
    NTP_involving_rxns = gams_parser.convert_set_to_list(
        os.path.join(data_dir, 'NTP_and_AMP_reactions.txt')
    )
    cofactor_only_rxns = gams_parser.convert_set_to_list(
        os.path.join(data_dir, 'cofactor_only_reactions.txt')
    )
    cofactor_only_rxns.append('R10092')

    methylglyoxal_rxns = [
        'R00203',
        'R00205',
        'R01016',
        'R02260',
        'R02527',
        'R02528',
        'R02529',
        'R02530',
        'R02531',
        'R07183',
        'R09796',
        'R10049',
        'R10050'
    ]
    atp_cycle_rxns = [
        'R00085',
        'R00086',
        'R00087',
        'R00088',
        'R00122',
        'R00123',
        'R00125',
        'R00127',
        # These are a part of a cycle
        'R01150',
        'R05226',
        'R07302',
        'R07651',
        # part of NADH/NADPH cycle
        'R00112'
    ]

    other_undesirable_rxns = [
        # Bicarbonate and pyrrole cycle
        'R09794',
        'R09795',
        # Undesirable glucose uptake loop
        'R00305',
        'R00874',
        'R07359',
        'R00837',
        'R09749',
        'R03075',
        'R02985',
        'R02558',
        'R01555',
        'R02727',
        'R00010',
        'R02778',
        'R08946',
        'R00306'
    ]

    excluded_reactions = (NTP_involving_rxns + methylglyoxal_rxns +
                          cofactor_only_rxns + atp_cycle_rxns +
                          other_undesirable_rxns)

    all_excluded_reactions = list(set(excluded_reactions))

    dbdict = {
        'Sji': 'optstoic_v3_Sji_dict.json',
        'reaction': 'optstoic_v3_reactions.txt',
        'metabolite': 'optstoic_v3_metabolites.txt',
        'reactiontype': 'optstoic_v3_reactiontype.txt',
        'loops': 'optstoic_v3_loops_nocofactor.txt',
        'Nint': 'optstoic_v3_null_sij_nocofactor.txt',
        'blocked_rxns': 'optstoic_v3_blocked_reactions_0to5ATP.txt',
        'excluded_reactions_list': all_excluded_reactions
    }

    DB = Database(version='3', data_filepath=data_dir, dbdict=dbdict)
    DB.load()

    # Update reaction type
    # #Method 1

    # DB.extend_S_from_file('Sij_extension_for_glycolysis.txt')
    # DB.rxntype = gams_parser.convert_parameter_list_to_dict(
    #     os.path.join(data_filepath, "rxntype_extension_for_glycolysis.txt"),
    #     datadict=DB.rxntype
    # )

    # #Method 2
    # Update reaction type  = 0
    irreversible_fwd_rxns = gams_parser.convert_set_to_list(os.path.join(
        data_dir, 'optstoic_v3_ATP_irreversible_forward_rxns.txt')
    )

    new_reaction_type_dict = dict(zip(
        irreversible_fwd_rxns, [0] * len(irreversible_fwd_rxns))
    )
    # Update reaction type  =  2
    irreversible_bwd_rxns = gams_parser.convert_set_to_list(os.path.join(
        data_dir, 'optstoic_v3_ATP_irreversible_backward_rxns.txt')
    )

    new_reaction_type_dict.update(dict(
        zip(irreversible_bwd_rxns, [2] * len(irreversible_bwd_rxns)))
    )

    DB.update_rxntype(new_reaction_type_dict)

    # user_defined_export_rxns = ['EX_glc', 'EX_nad', 'EX_adp',
    #                             'EX_phosphate', 'EX_pyruvate', 'EX_nadh',
    #                             'EX_atp', 'EX_h2o', 'EX_hplus', 'EX_nadp',
    #                             'EX_nadph']

    # Add a list of export reactions and the metabolites
    if user_defined_export_rxns_Sji is not None:
        user_defined_export_rxns_Sij = Database.transpose_S(
            user_defined_export_rxns_Sji
        )

        DB.set_database_export_reaction(user_defined_export_rxns_Sij)

    return DB


if __name__ == "__main__":
    DB = load_db_v3()
