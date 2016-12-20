#!/usr/bin/env python
"""Menace - Metagenomic Estimation of Relative Cell Periods
After verifying the command line arguments and the directory
structure, this module generates the necessary scripts for running the
whole pipeline (for all samples). Adapted from rnaseq_pipeline by
Matthias Nilsson <mattiasn@chalmers.se>.

@author Daniel Hermansson hedani@chalmers.se
"""
import argparse
import configparser
import codecs
import os
import sys
from sys import exit
from subprocess import Popen, PIPE, STDOUT
from ftplib import FTP
import re
from math import ceil
from platform import system
from jinja2 import Environment, FileSystemLoader
from shutil import copy

from lib.Community import local_conf

CWD = os.getcwd()
CODE_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = os.path.join(os.path.dirname(CODE_DIR), "templates")
DEFAULT_CONFIG = os.path.join(CWD, "project.conf")

def execute(command):
    process = Popen(command, shell=True, stdout=PIPE, stderr=STDOUT)

    # Poll process for new output until finished
    while True:
        nextline = process.stdout.readline()
        if nextline == '' and process.poll() is not None:
            break
        sys.stdout.write(nextline)
        sys.stdout.flush()

    output = process.communicate()[0]
    exitCode = process.returncode

    if (exitCode == 0):
        return output
    else:
        raise ProcessException(command, exitCode, output)

def read_config(args):
    """Read project config."""
    config = configparser.ConfigParser(allow_no_value=True, inline_comment_prefixes=(';',))
    config.optionxform = str
    config.readfp(codecs.open(args.config_file, "r", "utf8"))
    #config.read(args.config_file)

    if type(args.email) is str:
        config['Other']['Email'] = args.email
    #elif 'Email' in config['Other'].keys():
    elif not 'Email' in list(config['Other'].keys()):
        raise NameError('No email specified. Edit email in the project config file or pass it as an argument with "-e".')

    return config

def get_parser():
    """Get command line argument parser."""
    parser = argparse.ArgumentParser(prog='menace',
                                     description=("Perform tasks (individually or grouped) to extract relative cell periods from a metagenomic cohort. Before using the software, please edit the necessary information in the supplied configuration file after running 'menace init'."))

    parser.add_argument('-e', '--email', dest='email',
                        type=str, required=False,
                        help="Contact email address (can also be set in the config.)")

    parser.add_argument('-c', '--config', dest='config_file', type=str,
                        default=DEFAULT_CONFIG,
                        help="Config file for project specific settings (default: './project.conf').")

    subparsers = parser.add_subparsers(dest='subparser_name')

   
    #parser_t.add_argument("-l", action='store_true',help='Test a local config.')
    #parser_t.add_argument("-c", action='store_true',help='Test a minimal cluster config.')    

    parser_feSeq = subparsers.add_parser('fetch-data', help='download metagenomic sequences (using FTP)')
    parser_feSeq.add_argument("-t", dest='feSeq_threads', help='The number of simultaneous wget processes to use.',default=1)

    parser_fe = subparsers.add_parser('fetch-references', help='download references from NCBI')
    parser_fe.add_argument("-s", dest='fe_srch_file', help='A file with search strings to download.',default='./searchStrings')

    parser_b = subparsers.add_parser('build-index', help='build a mapping index')
    parser_b.add_argument("-t", dest='b_threads', help='The number of cpu threads to use.',default=1)
    #parser_b.add_argument("--opt4", action='store_true')

    #parser_m = subparsers.add_parser('make', help='generate an sbatch jobscript')
    #parser_m.add_argument("-e", action='store_true')

    parser_s = subparsers.add_parser('run', help='submit job to slurm or run locally')
    #parser_s.add_argument("-opt9", action='store_true')
    #parser_s.add_argument("--opt10", action='store_true')

    parser_f = subparsers.add_parser('full', help='perform all the above steps')
    parser_f.add_argument("-t1", dest='feSeq_threads', help='The number of simultaneous wget processes to use.',default=1)
    parser_f.add_argument("-t2", dest='b_threads', help='The number of cpu threads to use.',default=1)
    parser_f.add_argument("-s1", dest='fe_srch_file', help='A file with search strings to download.',default='./searchStrings')

    #parser_f.add_argument("-opt5", action='store_true')
    #parser_f.add_argument("--opt6", action='store_true')

    parser_c = subparsers.add_parser('collect', help='collect output cluster data')
    parser_c.add_argument("-o", dest='c_out_path', help='Path to store output.',default="./")
    parser_c.add_argument("--minOrics", dest='min_orics', help='The minimum number of fitted origins needed to make an end estimate of a reference origin.',default=1)
    #parser_c.add_argument("--opt8", action='store_true')

    parser_i = subparsers.add_parser('init', help='init an empty project in current directory')
    parser_i.add_argument("-l", action='store_true',help='Generate a minimal local config.')
    parser_i.add_argument("-c", action='store_true',help='Generate a minimal cluster config.')

    parser_t = subparsers.add_parser('test', help='generate a test project and run it on a small example data set in current directory')

    #parser_n = subparsers.add_parser('notebook', help='open a jupyter notebook of the project in the current directory with data import and analysis examples (requires jupyter)')

    return parser

    # parser.add_argument('-m', '--mode', dest='mode', type=str,
    #                     choices=['single', 'paired'], default='paired',
    #                     help="Sequencing read mode")

def parse_args(args):
    """Parse command line arguments"""
    args = get_parser().parse_args(args)
    if(args.subparser_name=='full'):
        args.fe_srch_file='./searchStrings'
        args.b_threads=1
        args.feSeq_threads=1
    return args


# def format_args(args):
#     """Format command line arguments (adding defaults)."""
#     if hasattr(args, 'contact_email'):
#         contact_email = args.contact_email
#     else
#         contact_email = args.contact_email
    
#     config_file = args.config_file

#     return {
#         'contact_email': contact_email,
#         'config_file': config_file
#     }


# def generate_scripts(project, config):
#     """Generate all pipeline scripts."""
#     env = Environment(loader=FileSystemLoader("templates"))
#     generate_job_scripts(project, config, env)
#     generate_de_script(project, config, env)
#     generate_pipeline_script(project, config, env)


# def generate_job_scripts(project, config, env):
#     """Generate the scripts needed to run a single sample."""
#     template = env.get_template("run_sample.sh")
#     with open("run_sample.sh", "w") as f:
#         output = template.render(
#             project=config['ProjectID'],
#             job_name=config['JobName'],
#             job_nodes=config['JobNodes'],
#             estimated_time=config['EstimatedTime'],
#             contact_email=project.email,
#         )
#         f.write(output)

#     template = env.get_template("run_job.py")
#     with open("run_job.py", "w") as f:
#         output = template.render(
#             code_dir=CODE_DIR,
#             genome_dir=project.genome_dir,
#             fasta_file=project.reference_genome,
#             annotation_file=project.annotations,
#             read_mode=project.mode,
#             adapter=project.adapter,
#             sample_ids=[sample.sample_id for sample in project.samples],
#             samples=[sample.tojson() for sample in project.samples],
#             provider=project.provider.get_provider_name(),
#             project_dir=project.project_dir
#         )1
#         f.write(output)

def get_current_dir_subdirs(f):
    ret = []
    f.dir("",ret.append)
    ret = [x.split()[-1] for x in ret if x.startswith("d")]
    return ret 

def compile_config(args,config):    
    #try:
    #    start_ind1, data_prefix=get_data_prefix(config)
    #except e:
    #    print("The data accession prefix (eg. ERR525) could not be retrieved. Enter manually? (y/n)")

    if(config['Other']['DataPrefix']=="" and not config['Other']['FtpURL']==""):
        start_ind1, data_prefix = get_data_prefix(config)
    else:
        start_ind1 = 0
        data_prefix = config['Other']['DataPrefix']

    email=args.email if args.email else config['Other']['Email']
    if (config['Other']['StartInd']==""):
        start_ind=start_ind1
    else:
        start_ind=config['Other']['StartInd']
    
    job_range = ceil(float(config['Other']['NrSamples'])/float(config['Other']['SamplesPerNode']))-1
    conf={
        'project': config['Project']['ProjectID'],
        'cluster': config['Project']['Cluster'],
        'job_name': config['Project']['JobName'],
        'job_nodes': config['Project']['JobNodes'],
        'cpu_cores': config['Project']['CpuCores'],
        'estimated_time': config['Project']['EstimatedTime'],

        'node_path': config['Directories']['Node'],
        'ref_path': config['Directories']['References'],
        'data_path': config['Directories']['Data'],
        'output_path': config['Directories']['Output'],
        'doric_path': config['Directories']['DoriC'],

        'mapper': config['Other']['Mapper'],
        'ref_name': config['Other']['RefName'],
        'nr_samples': config['Other']['NrSamples'],
        'samples_per_node': config['Other']['SamplesPerNode'],
        'email': email,
        'data_prefix': data_prefix,
        'start_ind': start_ind,
        'job_range': job_range,

        'ftp_url': config['Other']['FtpURL'],
        'data_url': config['Other']['DataURL']
    }
    return conf

        # 'fe_srch_file': args.fe_srch_file,
        # 'b_threads': args.b_threads,
        # 'c_out_path': args.c_out_path,
        # 'feSeq_threads': feSeq_threads

def save_config(config):
    render_conf('project.conf',config,CWD,'project.conf')

def generate_jobscript(config):
    template = 'jobscript'

    if (config['cluster']==''):
        template = 'jobscript_local'

    render_conf(template,config,CWD,'jobscript')

def render_conf(template,config,path,name):
    path=os.path.join(path,name)
    env = Environment(loader=FileSystemLoader(os.path.join(CODE_DIR,'templates')))
    template = env.get_template(template)
    
    with open(path, 'w') as f:
            output = template.render(
                project=config['project'],
                cluster=config['cluster'],
                job_name=config['job_name'],
                job_nodes=config['job_nodes'],
                cpu_cores=config['cpu_cores'],
                estimated_time=config['estimated_time'],
                
                node_path=config['node_path'],
                ref_path=config['ref_path'],
                data_path=config['data_path'],
                output_path=config['output_path'],
                doric_path=config['doric_path'],

                mapper=config['mapper'],
                ref_name=config['ref_name'],
                nr_samples=config['nr_samples'],
                samples_per_node=min(int(config['samples_per_node']),int(config['nr_samples'])),
                email=config['email'],
                data_prefix=config['data_prefix'],
                start_ind=config['start_ind'],
                job_range=config['job_range']
            )
            f.write(output)
            os.chmod(path,0o777)

# def generate_pipeline_script(project, config, env):
#     """Generate script for running complete pipeline."""
#     template = env.get_template("run_pipeline.sh")
#     with open("run_pipeline.sh", "w") as f:
#         output = template.render(
#             trimming_cmd=generate_sbatch_command(project),
#             de_script="run_de_and_generate_report.sh",
#             project_id=project.provider.get_project_id()
#         )
#         f.write(output)
#     permissions = os.stat('run_pipeline.sh')
#     os.chmod("run_pipeline.sh", (permissions.st_mode |
#                                  stat.S_IXUSR |
#                                  stat.S_IXGRP |
#                                  stat.S_IXOTH))

def generate_bt2_build_command(args,config):
    """Generate command for building a bowtie2 index."""
    #for file in os.listdir(os.path.join(config['Directories']['References'],'Fasta')):
        #if file.endswith(".fasta"):
        #    print(file)
    d=os.path.join(config['ref_path'],'Fasta');
    files=[f for f in os.listdir(d)]
    referenceList=",".join(files)
    #os.environ['reflist'] = referenceList
    #proc = Popen("for f in " + d + "/*.fasta" + '; do cat "$f" >> tmpfasta', shell=True)
    #print(referenceList)

    cmd = "(cd " + d + " && bowtie2-build --large-index " + "--threads " + str(args.b_threads) + " " + referenceList + " " + os.path.join(config['ref_path'],"Index",config['ref_name']) + ")"
    return cmd

def generate_gem_build_command(args,config):
    """Generate command for building a gem index."""
    #for file in os.listdir(os.path.join(config['Directories']['References'],'Fasta')):
        #if file.endswith(".fasta"):
        #    print(file)
    d=os.path.join(config['ref_path'],'Fasta');
    #files=[f for f in os.listdir(d)]
    #referenceList=",".join(files)
    #os.environ['reflist'] = referenceList
    #proc = Popen("for f in " + d + "/*.fasta" + '; do cat "$f" >> tmpfasta', shell=True)
    #print(referenceList)
    if (system()=='Linux'):
        cmd = "(cd " + d + " && gem-indexer " + "-T " + str(args.b_threads) + " -i " + os.path.join(d,"multi.fasta") + " -o " + os.path.join(config['ref_path'],"Index",config['ref_name']) + ")"
    elif (system()=='Darwin'): # for Mac OS X assume older (2010) binaries
        cmd = "(cd " + d + " && gem-do-index " + " -i " + os.path.join(config['ref_path'],"Fasta","multi.fasta") + " -o " + os.path.join(config['ref_path'],"Index",config['ref_name']) + ")"
    return cmd

def get_data_prefix(config):
    f = FTP(config['Other']['FtpURL'])
    f.login()
    f.cwd(config['Other']['DataURL'])
    dirs=get_current_dir_subdirs(f)
    data_prefix=os.path.commonprefix(dirs)
    startind=re.match(data_prefix+"(.+)",dirs[0])
    startind=startind.group(1)
    return (startind,data_prefix)

def generate_fetch_seq_command(args,config):
    nr_digits=len(config['start_ind'])
    end_ind=int(config['start_ind'])+int(config['nr_samples'])-1
    url=os.path.join(config['data_url'],config['data_prefix'])
    url=config['ftp_url']+url
    
    if(not re.match("ftp://",url)):
        url="ftp://"+url

    cmd=" | parallel -j " + str(args.feSeq_threads) + " wget -c -r -nH -np -nd -R index.html* -P " + config['data_path'] + " " + url + "{}/"
    #cmd="$(echo {"+ str(config['start_ind']) + ".." + str(end_ind) + "})"+cmd
    cmd='seq -f %0' + str(nr_digits) + 'g ' + str(config['start_ind']) + " " + str(end_ind) + cmd
    #print(cmd)
    return cmd

def generate_fetch_ref_command(args,config):
    """Generate command for downloading reference fastas and headers."""
    #dir_path = os.path.dirname(os.path.realpath(__file__))
    cmd = "fetch_seq -e {email} -t True -d {ref_path} -s " + args.fe_srch_file
    return cmd.format(**config)

def generate_sbatch_command(config):
    """Generate command for scheduling all sample runs."""
    cmd = "sbatch --array=0-{0} jobscript"
    return cmd.format(config['job_range'])

def generate_local_command(config):
    """Generate command for performing a local run."""
    cmd = os.path.join(CWD,'jobscript '+ CODE_DIR)
    return cmd

def generate_collect_command(config,args):
    koremLoc=os.path.join(CODE_DIR,'extra/accLoc.csv')
    cmd="python -W ignore " + CODE_DIR + "/bin/PTRMatrix.py {output_path} {ref_path} " + str(args.min_orics) + " {output_path} {doric_path} " + koremLoc
    return cmd.format(**config)

def generate_notebook_command(args,config):
    copy(os.path.join(CODE_DIR,'extra','menace_run.ipynb'),CWD)
    cmd="jupyter notebook menace_run.ipynb &"
    return cmd.format(**config)

def run_test_command(args):
    conf = local_conf(CWD,'bowtie2','','1')
    if args.email:
        conf['email']=args.email
    make_dirs(conf)
    copy(os.path.join(CODE_DIR,'test','comm00_1.fastq'),os.path.join(CWD,'Data'))
    copy(os.path.join(CODE_DIR,'test','comm00_2.fastq'),os.path.join(CWD,'Data'))
    copy(os.path.join(CODE_DIR,'test','searchStrings'),CWD)
    save_config(conf)
    return conf

def run_init_command(args):
    config=local_conf(CWD,'bowtie2','','1')
    if args.email:
        config['email']=args.email
    copy(os.path.join(CODE_DIR,'test', 'searchStrings'),CWD)
    make_dirs(config)
    save_config(config)
    return config

def make_dirs(conf):
    key_arr=['ref_path','data_path','output_path','doric_path']
    for k in key_arr:
        if not os.path.exists(conf[k]):
            os.makedirs(conf[k])

    ref_sub=['Fasta','Headers','Index']
    rp=conf['ref_path']
    for d in ref_sub:
        if not os.path.exists(os.path.join(rp,d)):
            os.makedirs(os.path.join(rp,d))


# def print_instructions(config):
#     """Print instructions for executing pipeline."""
#     print("The necessary scripts have now been generated. To execute the "
#           "pipeline, run the following command:")
#     print("")
#     print("\t./run_pipeline.sh")
#     print("")
#     print("The results will be saved in "
#           "$SNIC_NOBACKUP/alignment_{0}/.".format(23))
#     print("")
#     print("Once you have started the pipeline, you can run 'deactivate' to "
#           "exit the virtual environment.")
#     print("")

def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def main2(args,config):
    if(args.subparser_name=='init'):
        run_init_command(args)
        exit()
    
    if(args.subparser_name=='notebook'):
        if (not exists(join(CWD,'project.conf'))):
            run_init_command(args)
        process = Popen(generate_notebook_command(args,config), shell=True)
        exit()

    if(args.subparser_name=='full' or args.subparser_name=='fetch-data'):
        ans='run'
        if (config['DataURL']=='' or config['FtpURL']==''):
            ans='skip'
            ans=raw_input('DataURL or FtpURL empty in project.conf, skip fetch-data or end? [skip]/end: ')
        #if ans=='skip':
            #process = Popen(generate_fetch_seq_command(args,config), shell=True)
        #    break
            #process.wait()
        if ans=='end':
            exit()

        if not ans=='skip':
            execute(generate_fetch_seq_command(args,config))

    if(args.subparser_name=='test'):
        config=run_test_command(args)
        
        #streamdata = process.communicate()[0]
        #if(process.returncode==0):
        args.subparser_name='full'
        args.fe_srch_file='searchStrings'
        args.b_threads='1'
        #args.

    if(args.subparser_name=='full' or args.subparser_name=='fetch-references'):
        #print(os.path.join(CWD,'searchStrings'))
        #ans='run'
        #if (is_non_zero_file(os.path.join(CWD,'searchStrings'))):
        #    ans='skip'
        #    ans=raw_input('searchStrings empty or non existent, skip fetch-references or end? [skip]/end: ')
        #if ans=='end':
        #    exit()

        #if not ans=='skip':
        process = Popen(generate_fetch_ref_command(args,config), shell=True)
        process.wait()

    if(args.subparser_name=='full' or args.subparser_name=='build-index'):
        process = Popen(CODE_DIR+"/bin/changeTID.sh " + config['ref_path'], shell=True)
        process.wait()
        
        tmp_dir=os.path.join(config['ref_path'],"Index")
        
        if (not os.path.isdir(tmp_dir)):
            os.makedirs(tmp_dir)
        
        if(  os.listdir(tmp_dir)!=[] and args.subparser_name=='full'):
            pass
        else:
            if(config['mapper']=='bowtie2'):
                process = Popen(generate_bt2_build_command(args,config), shell=True)
                process.wait()
            elif(config['mapper']=='gem'):
                process = Popen(CODE_DIR+"/extra/concatenateFasta.sh " + os.path.join(config['ref_path'],'Fasta'), shell=True)
                process.wait()

                process = Popen(generate_gem_build_command(args,config), shell=True)
                process.wait()

    #if((args.subparser_name=='full' and not config.cluster=='') or args.subparser_name=='make'):
    #    generate_jobscript(config)

    if(args.subparser_name=='full' or args.subparser_name=='run'):
        generate_jobscript(config)

        if (config['cluster']==''):
            process = Popen(generate_local_command(config), shell=True)
            process.wait()
        else:
            process = Popen(generate_sbatch_command(config), shell=True)
            process.wait()

    if(args.subparser_name=='collect'):
        process = Popen(generate_collect_command(config,args), shell=True)
        process.wait()

    # prov = provider.get_provider(args['provider_name'], args['project_dir'])
    # p = project.Project(args, prov)
    # p.verify_directory_structure()
    # p.read_samples()
    # p.validate_de_input()
    # generate_scripts(p, config)
    # print_instructions(p)

def main():
    args = parse_args(sys.argv[1:])

    if(not (args.subparser_name=='init' or args.subparser_name=='test')):
        config = read_config(args)
        config = compile_config(args,config)
    else:
        config=[]
    main2(args,config)

if __name__ == "__main__":
    main()
