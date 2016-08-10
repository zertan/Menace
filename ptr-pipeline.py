#!/usr/bin/env python3
"""Generate scripts for running PTR pipeline on C3SE cluster.
After verifying the command line arguments and the directory
structure, this module generates the necessary scripts for running the
whole pipeline (for all samples). Adapted from rnaseq_pipeline by
Matthias Nilsson <mattiasn@chalmers.se>.
"""

import argparse
import configparser
import os
import stat
import sys
import subprocess
import ftplib
import re
from math import ceil

from jinja2 import Environment, FileSystemLoader

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = os.path.join(os.path.dirname(CODE_DIR), "templates")
DEFAULT_CONFIG = os.path.join(CODE_DIR, "project.conf")


def read_config(args):
    """Read project config."""
    config = configparser.ConfigParser(allow_no_value=True, inline_comment_prefixes=(';',))
    config.read(args.config_file)

    if type(args.email) is str:
        config['Other']['Email'] = args.email
    #elif 'Email' in config['Other'].keys():
    elif not 'Email' in config['Other'].keys():
        raise NameError('No email specified. Edit email in the project config file or pass it as an argument with "-e".')

    return config

def get_parser():
    """Get command line argument parser."""
    parser = argparse.ArgumentParser(prog='ptr-pipeline',
                                     description=("Perform tasks (individually or grouped) to extract PTR values of a metagenomic cohort. Before using the script, please edit the necessary information in the supplied configuration file."))

    parser.add_argument('-e', '--email', dest='email',
                        type=str, required=False,
                        help="Contact email address (can also be set in the config.)")

    parser.add_argument('-c', '--config', dest='config_file', type=str,
                        default=DEFAULT_CONFIG,
						help="Config file for project specific settings (default: './project.conf').")

    subparsers = parser.add_subparsers(dest='subparser_name')

    parser_feSeq = subparsers.add_parser('fetch-data', help='download metagenomic sequences (using FTP)')
    parser_feSeq.add_argument("-t", dest='feSeq_threads', help='The number of simultaneous wget processes to use.',default=1)
    

    parser_fe = subparsers.add_parser('fetch-references', help='download references from NCBI')
    parser_fe.add_argument("-s", dest='fe_srch_file', help='A file with search strings to download.',default='./searchStrings')

    parser_b = subparsers.add_parser('build-index', help='build a bowtie2 index')
    parser_b.add_argument("-t", dest='b_threads', help='The number of cpu threads to use.',default=1)
    #parser_b.add_argument("--opt4", action='store_true')

    parser_m = subparsers.add_parser('make', help='generate an sbatch jobscript')
    #parser_m.add_argument("-e", action='store_true')

    parser_s = subparsers.add_parser('submit', help='submit job to slurm')
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
    #parser_c.add_argument("--opt8", action='store_true')

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
    #	start_ind1, data_prefix=get_data_prefix(config)
    #except e:
    #	print("The data accession prefix (eg. ERR525) could not be retrieved. Enter manually? (y/n)")

    if(config['Other']['DataPrefix']==""):
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

def generate_jobscript(config):
    env = Environment(loader=FileSystemLoader("templates"))
    template = env.get_template("jobscript")

    with open("jobscript", "w") as f:
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
    #proc = subprocess.Popen("for f in " + d + "/*.fasta" + '; do cat "$f" >> tmpfasta', shell=True)
	#print(referenceList)

    cmd = "(cd " + d + " && bowtie2-build --large-index " + "-t " + str(args.b_threads) + " " + referenceList + " " + os.path.join(config['ref_path'],"Index",config['ref_name']) + ")"
    return cmd

def generate_gem_build_command(args,config):
    """Generate command for building a bowtie2 index."""
    #for file in os.listdir(os.path.join(config['Directories']['References'],'Fasta')):
        #if file.endswith(".fasta"):
        #    print(file)
    #d=os.path.join(config['ref_path'],'Fasta');
    #files=[f for f in os.listdir(d)]
    #referenceList=",".join(files)
    #os.environ['reflist'] = referenceList
    #proc = subprocess.Popen("for f in " + d + "/*.fasta" + '; do cat "$f" >> tmpfasta', shell=True)
    #print(referenceList)

    cmd = "(cd " + d + " && gem-indexer " + "-T " + str(args.b_threads) + " -i " + os.path.join(config['ref_path'],"Fasta","multi.fasta") + " -o " + os.path.join(config['ref_path'],"Index",config['ref_name']) + ")"
    return cmd

def get_data_prefix(config):
    f = ftplib.FTP(config['Other']['FtpURL'])
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
    cmd = "bin/fetchSeq.py -e {email} -t True -d {ref_path} -s " + args.fe_srch_file
    return cmd.format(**config)

def generate_sbatch_command(config):
    """Generate command for scheduling all sample runs."""
    cmd = "sbatch --array=0-{0} jobscript"
    return cmd.format(config['job_range'])

def generate_collect_command(config):
	koremLoc=os.path.join(CODE_DIR,"extra/accLoc.csv") 
	cmd="bin/PTRMatrix.py {output_path} {ref_path} {doric_path} " + koremLoc
	return cmd.format(**config)

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

def main():
	args = parse_args(sys.argv[1:])
	config = read_config(args)
	config = compile_config(args,config)

	if(args.subparser_name=='full' or args.subparser_name=='fetch-data'):
		process = subprocess.Popen(generate_fetch_seq_command(args,config), shell=True)
		process.wait()

	if(args.subparser_name=='full' or args.subparser_name=='fetch-references'):
		process = subprocess.Popen(generate_fetch_ref_command(args,config), shell=True)
		process.wait()

	if(args.subparser_name=='full' or args.subparser_name=='build-index'):
		process = subprocess.Popen("bin/changeTID.sh " + config['ref_path'], shell=True)
		process.wait()
		
		tmp_dir=os.path.join(config['ref_path'],"Index")
		
		if (not os.path.isdir(tmp_dir)):
			os.makedirs(tmp_dir)
		
		if(  os.listdir(tmp_dir)!=[] and args.subparser_name=='full'):	
			pass
		else:
			if(config['mapper']=='bowtie2'):
				process = subprocess.Popen(generate_bt2_build_command(args,config), shell=True)
				process.wait()
			elif(config['mapper']=='gem'):
				process = subprocess.Popen("extra/concatenateFasta.sh " + os.path.join(config['ref_path'],'Fasta'), shell=True)
				process.wait()

				process = subprocess.Popen(generate_gem_build_command(args,config), shell=True)
				process.wait()

	if(args.subparser_name=='full' or args.subparser_name=='make'):
		generate_jobscript(config)

	if(args.subparser_name=='full' or args.subparser_name=='submit'):
		process = subprocess.Popen(generate_sbatch_command(config), shell=True)
		process.wait()

	if(args.subparser_name=='collect'):
		process = subprocess.Popen(generate_collect_command(config), shell=True)
		process.wait()
    # prov = provider.get_provider(args['provider_name'], args['project_dir'])
    # p = project.Project(args, prov)
    # p.verify_directory_structure()
    # p.read_samples()
    # p.validate_de_input()
    # generate_scripts(p, config)
    # print_instructions(p)


if __name__ == "__main__":
    main()
