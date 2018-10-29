'''
SYNBIOCHEM (c) University of Manchester. 2018

PathwayGenie is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=broad-except
# pylint: disable=invalid-name
# pylint: disable=no-member
import json
import sys
from rdkit import Chem, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import pandas as pd


def analyse(filename_in, filename_out):
    '''Analyse predictions.'''
    with open(filename_in) as json_data:
        data = json.load(json_data)

        for entry in data['reaction']:
            entry['similarity'] = _get_similarity(entry['reactant_smile'],
                                                  entry['product_smile'])

    # Write updated json:
    with open(filename_out + '.json', 'w') as outfile:
        json.dump(data, outfile)

    # Write Dataframe:
    df = pd.DataFrame(data['reaction'])
    df.sort_values(by=['similarity'], ascending=False, inplace=True)
    df.to_csv(filename_out + '.csv', encoding='utf-8-sig', index=False)


def _get_similarity(smiles1, smiles2):
    '''Get similarity.'''
    try:
        return DataStructs.FingerprintSimilarity(_get_fingerprint(smiles1),
                                                 _get_fingerprint(smiles2))
    except Exception:
        return float('NaN')


def _get_fingerprint(smiles):
    '''Get fingerprint.'''
    mol = Chem.MolFromSmiles(smiles)

    if mol.GetNumAtoms(onlyExplicit=False) == 1:
        raise ValueError(smiles)

    return FingerprintMols.FingerprintMol(mol)


def main(args):
    '''main method.'''
    analyse(*args)


if __name__ == '__main__':
    main(sys.argv[1:])
