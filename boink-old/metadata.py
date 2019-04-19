import datetime
import os
import pandas as pd
from jinja2 import Environment, PackageLoader, select_autoescape

BOINK_DIR = os.path.abspath(os.path.join(__file__, os.pardir))
DATA_DIR =  os.path.join(BOINK_DIR, 'data')
TEMPLATE_DIR = os.path.join(BOINK_DIR, 'templates')
CUR_TIME = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M")

__version__ = '0.1'

def get_template_env():
    return Environment(loader=PackageLoader('boink', 'templates'))

def get_typegen_env():
    return Environment(loader=PackageLoader('boink', 'typegen'))

def parse_mmetsp_metadata(sraruninfo_path,
                          assembly_dir='/mnt/research/ged/data/mmetsp/mmetsp_assemblies_trinity2.2.0_zenodo/',
                          assembly_suffix='.trinity_out_2.2.0.Trinity.fasta'):

    sra_df = pd.read_csv(sraruninfo_path)
    sra_df['ReadsDir'] = sra_df.ScientificName.str.replace(' ', '_').apply(lambda name: os.path.join(DATA_DIR, name))
    sra_df['ReadsDir'] = sra_df.apply(lambda row: os.path.join(row.ReadsDir, row.Run), axis=1)
    sra_df['AssemblyPath'] = sra_df.apply(lambda row: os.path.join(assembly_dir, row.SampleName + assembly_suffix), axis=1)

    return sra_df

