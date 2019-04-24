import os
from shutil import rmtree
import pandas as pd
import jinja2
from doit.task import clean_targets
from boink.metadata import parse_mmetsp_metadata, DATA_DIR, BOINK_DIR, TEMPLATE_DIR, get_template_env

METADATA = parse_mmetsp_metadata(os.path.join(DATA_DIR, 'MMETSP_SraRunInfo_subset.csv'))
SUBMIT = True
templates = get_template_env()

def clean_folder(target):
    '''Function for doit task's `clean` parameter to remove a folder.

    Args:
        target (str): The folder to remove.
    '''

    try:
        rmtree(target)
    except OSError:
        pass


def write_partition_script(partition_dir, script_path, samples, 
                           time='04:00:00', ppn=1, mem='16gb',
                           account='ged', email='camille.scott.w@gmail.com',
                           pairing='single'):

    template = templates.get_template('pbs-partition.tpl')
    with open(script_path, 'wb') as fp:
        print('\tWriting', script_path)
        pbs = template.render(samples=samples, stats_dir=partition_dir, time=time,
                              ppn=ppn, mem=mem, account=account, 
                              pairing=pairing)
        fp.write(pbs.encode())

def task_partition_mmetsp_assemblies():
    for _, row in METADATA.iterrows():
        if not os.path.exists(row.AssemblyPath):
            continue
        script_path = os.path.join(row.ReadsDir, 'partition-assembly.sh')
        partition_dir = os.path.join(row.ReadsDir, 'partitioned-assembly')
        partition_done_path = os.path.join(partition_dir, 'global-stats.csv')
        actions = [(write_partition_script, [partition_dir, script_path, [row.AssemblyPath]])]
        if SUBMIT:
            actions.append('qsub -V {0}'.format(script_path))

        yield {'name': row.AssemblyPath,
                'actions': actions,
                'file_dep': [row.AssemblyPath],
                'targets': [script_path,
                            partition_done_path],
                'clean': [(clean_folder, [partition_dir]),
                          clean_targets]}


def task_partition_mmetsp_reads():
    for _, row in METADATA.iterrows():
        for K in [19, 23, 27, 31]:
            script_path = os.path.join(row.ReadsDir, 'partition-reads-{0}.sh'.format(K))
            partition_dir = os.path.join(row.ReadsDir, 'partitioned-reads-{0}'.format(K))
            partition_done_path = os.path.join(partition_dir, 'global-stats.csv')
            samples = [os.path.join(row.ReadsDir, '{0}_{1}.fastq'.format(row.Run, n)) for n in [1,2]]
            
            actions = [(write_partition_script, 
                        [partition_dir, script_path, samples],
                        {'time': '36:00:00', 'pairing': 'split'})]

            if SUBMIT:
                actions.append('qsub -V {0}'.format(script_path))

            yield {'name': '{0}-{1}'.format(row.SampleName, K),
                    'actions': actions,
                    'file_dep': samples,
                    'targets': [script_path,
                                partition_done_path],
                    'clean': [(clean_folder, [partition_dir]),
                              clean_targets]}
