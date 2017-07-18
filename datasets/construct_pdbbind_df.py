"""
Contains methods for generating a pdbbind dataset mapping
  complexes (protein + ligand) to experimental binding measurement.
"""
from __future__ import print_function
import pickle
import os
import pandas as pd
from rdkit import Chem
from glob import glob
import re
from sklearn.externals import joblib
from deepchem.utils import rdkit_util


def extract_labels(pdbbind_label_file):
  """Extract labels from pdbbind label file."""
  assert os.path.isfile(pdbbind_label_file)
  labels = {}
  with open(pdbbind_label_file) as f:
    content = f.readlines()
    for line in content:
      if line[0] == "#":
        continue
      line = line.split()
      # lines in the label file have format
      # PDB-code Resolution Release-Year -logKd Kd reference ligand-name
      #print line[0], line[3]
      labels[line[0]] = line[3]
  return labels

def construct_df(datekind, pdb_stem_directory, pdbbind_label_file, pdbbind_df_joblib, pdb_exception_files):
  """
  Takes as input a stem directory containing subdirectories with ligand
    and protein pdb/mol2 files, a pdbbind_label_file containing binding
    assay data for the co-crystallized ligand in each pdb file,
    and a pdbbind_df_pkl to which will be saved a pandas DataFrame
    where each row contains a pdb_id, smiles string, unique complex id,
    ligand pdb as a list of strings per line in file, protein pdb as a list
    of strings per line in file, ligand mol2 as a list of strings per line in
    mol2 file, and a "label" containing the experimental measurement.
  """
  labels = extract_labels(pdbbind_label_file)
  df_rows = []
  os.chdir(pdb_stem_directory)
  pdb_directories = [pdb.replace('/', '') for pdb in glob('*/')]    # Only select the directory and convert it to pdb_id

  exception_sets = [pdb_no_subdir, ligand_error, protein_error] = [[] for i in xrange(3)]
  num_complete = num_incomplete = 0
  for pdb_id in labels.iterkeys():
    print('==> Start to prcess {}'.format(pdb_id))
    if pdb_id in pdb_directories:
      print("About to extract ligand and protein input files")
      protein_pdb_withouth = pdb_id + '_protein.pdb'
      ligand_pdb = None
      protein_pdb = None
      for f in os.listdir(pdb_id):
        if re.search("_ligand_hyd.pdb$", f):
          ligand_pdb = f
        elif re.search("_protein_hyd.pdb$", f):
          protein_pdb = f
        elif re.search("_ligand.mol2$", f):
          ligand_mol2 = f

      if not ligand_pdb:
        print("Generate ligand with hydrogen file.")
        ligand_mol = Chem.MolFromMol2File(os.path.join(pdb_stem_directory, pdb_id, ligand_mol2))
        try:
          ligand_mol_withh = Chem.AddHs(ligand_mol)
          Chem.MolToPDBFile(ligand_mol_withh, os.path.join(pdb_stem_directory, pdb_id, pdb_id + '_ligand_hyd.pdb'))
          ligand_pdb = [i for i in os.listdir(pdb_id) if i.endswith('_ligand_hyd.pdb')][0]
        except:
          num_incomplete += 1
          ligand_error.append(pdb_id)
          print('%s the ligand with hydrogen cannot be created.' % pdb_id)
          continue

      if not protein_pdb:
        print("Generate protein with hydrogen file.")
        try:
          protein_pdb_withh = rdkit_util.load_molecule(os.path.join(pdb_stem_directory, pdb_id, protein_pdb_withouth))[1]
          writer = Chem.PDBWriter(os.path.join(pdb_stem_directory, pdb_id, pdb_id + '_protein_hyd.pdb'))
          writer.write(protein_pdb_withh)
          writer.close()
          protein_pdb = [i for i in os.listdir(pdb_id) if i.endswith('_protein_hyd.pdb')][0]
        except:
          num_incomplete += 1
          protein_error.append(pdb_id)
          print('%s the protein with hydrogen cannot be created.' % pdb_id)
          continue

      print('ligand_pdb: {}, protein_pdb: {}, ligand_mol2: {}'.format(ligand_pdb, protein_pdb, ligand_mol2))
      if not ligand_pdb or not protein_pdb:
        raise ValueError("Required files not present for %s" % pdb_id)

      ligand_pdb_path = os.path.join(pdb_id, ligand_pdb)
      protein_pdb_path = os.path.join(pdb_id, protein_pdb)
      ligand_mol2_path = os.path.join(pdb_id, ligand_mol2)

      with open(protein_pdb_path, "rb") as f:
        protein_pdb_lines = f.readlines()

      with open(ligand_pdb_path, "rb") as f:
        ligand_pdb_lines = f.readlines()

      try:
        with open(ligand_mol2_path, "rb") as f:
          ligand_mol2_lines = f.readlines()
      except:
        ligand_mol2_lines = []

      print("About to compute ligand smiles string.")
      ligand_mol = Chem.MolFromPDBFile(ligand_pdb_path)
      if ligand_mol is None:
        raise ValueError('The ligand_mol is None, can not convert to smiles.')
        # continue
      smiles = Chem.MolToSmiles(ligand_mol, canonical=True)

      complex_id = "%s%s" % (pdb_id, smiles)
      label = labels[pdb_id]
      df_rows.append([pdb_id, smiles, complex_id, protein_pdb_lines,
                      ligand_pdb_lines, ligand_mol2_lines, label])

      print('==> Complete num: %d, pdb_code: %s' % (num_complete, pdb_id))
      num_complete += 1
    else:
      print('==> Incomplete num: %d, no sub-directory for pdb_code: %s' % (num_incomplete, pdb_id))
      num_incomplete += 1
      pdb_no_subdir.append(pdb_id)

  # for i, pdb_dir in enumerate(pdb_directories):
  #   print("About to extract ligand and protein input files")
  #   pdb_id = os.path.basename(pdb_dir)  # pdb_dir and pdb_id are same
  #   # print(pdb_dir, pdb_id)
  #
  #   protein_pdb_withouth = pdb_id + '_protein.pdb'
  #   ligand_pdb = None
  #   protein_pdb = None
  #   for f in os.listdir(pdb_dir):
  #     if re.search("_ligand_hyd.pdb$", f):
  #       ligand_pdb = f
  #     elif re.search("_protein_hyd.pdb$", f):
  #       protein_pdb = f
  #     elif re.search("_ligand.mol2$", f):
  #       ligand_mol2 = f
  #
  #   print("Extracted Input Files:")
  #   if not ligand_pdb:
  #     ligand_mol = Chem.MolFromMol2File(os.path.join(pdb_stem_directory, pdb_id, ligand_mol2))
  #     ligand_mol_withh = Chem.AddHs(ligand_mol)
  #     ligand_pdb = Chem.MolToPDBFile(ligand_mol_withh, os.path.join(pdb_stem_directory, pdb_id , pdb_id + '_ligand_hyd.pdb'))
  #
  #   if not protein_pdb:
  #     protein_pdb_withh = rdkit_util.load_molecule(os.path.join(pdb_stem_directory, pdb_id, protein_pdb_withouth))[1]
  #     writer = Chem.PDBWriter(os.path.join(pdb_stem_directory, pdb_id, pdb_id + '_protein_hyd.pdb'))
  #     protein_pdb = writer.write(protein_pdb_withh)
  #     writer.close()
  #
  #   if not ligand_pdb or not protein_pdb:
  #       raise ValueError("Required files not present for %s" % pdb_dir)
  #   print('ligand_pdb: {}, protein_pdb: {}, ligand_mol2: {}'.format(ligand_pdb, protein_pdb, ligand_mol2))
  #
  #   ligand_pdb_path = os.path.join(pdb_dir, ligand_pdb)
  #   protein_pdb_path = os.path.join(pdb_dir, protein_pdb)
  #   ligand_mol2_path = os.path.join(pdb_dir, ligand_mol2)
  #
  #   with open(protein_pdb_path, "rb") as f:
  #     protein_pdb_lines = f.readlines()
  #
  #   with open(ligand_pdb_path, "rb") as f:
  #     ligand_pdb_lines = f.readlines()
  #
  #   try:
  #     with open(ligand_mol2_path, "rb") as f:
  #       ligand_mol2_lines = f.readlines()
  #   except:
  #     ligand_mol2_lines = []
  #
  #   print("About to compute ligand smiles string.")
  #   ligand_mol = Chem.MolFromPDBFile(ligand_pdb_path)
  #   if ligand_mol is None:
  #     continue
  #     # ligand_mol = Chem.MolFromMol2File(ligand_mol2_path)
  #   smiles = Chem.MolToSmiles(ligand_mol)
  #   complex_id = "%s%s" % (pdb_id, smiles)
  #   label = labels[pdb_id]
  #   df_rows.append([pdb_id, smiles, complex_id, protein_pdb_lines,
  #                   ligand_pdb_lines, ligand_mol2_lines, label])
  #
  #   print('==> Complete num: %d, pdb_code: %s' % (i, pdb_id))

  pdbbind_df = pd.DataFrame(df_rows, columns=('pdb_id', 'smiles', 'complex_id',
                                              'protein_pdb', 'ligand_pdb',
                                              'ligand_mol2', 'label'))

  joblib.dump(pdbbind_df, pdbbind_df_joblib)
  print('For dataset: {}, complete sample : {}, incomplete sample: {}'.format(data_kind, num_complete, num_incomplete))

  print('Store the exception data to disk.')
  for i, x in enumerate(exception_sets):
    with open(pdb_exception_files[i], 'w') as f:
      for item in x:
        f.write('%s\n' % item)

  print('The task: %s is done.\n******' % data_kind)

if __name__ == '__main__':
  pdb_data_path = '/home/yanjun/develop/deepchem/examples/pdbbind/v2015/'
  data_kinds = ('core', 'refined', 'general')
  pdb_label_files = ('INDEX_core_data.2013', 'INDEX_refined_data.2015', 'INDEX_general_PL_data.2015')
  kind_id = 1
  data_kind = data_kinds[kind_id]
  pdb_label_file = pdb_label_files[kind_id]
  save_file = pdb_data_path + str(data_kind)
  exception_file_suffix = ['_no_subdir.txt', '_ligand_exception.txt', '_protein_exception.txt']
  pdb_exception_files = map(lambda x: pdb_data_path + str(data_kind) + str(x), exception_file_suffix)

  pdb_label_file = os.path.join(pdb_data_path, pdb_label_file)
  construct_df(data_kind, pdb_data_path, pdb_label_file, save_file, pdb_exception_files)

  # 1mh5: https://sourceforge.net/p/rdkit/mailman/message/29471913/