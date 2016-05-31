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

from jinja2 import Environment, FileSystemLoader

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = os.path.join(os.path.dirname(CODE_DIR), "templates")
#CONFIG_DIR = os.path.join(os.path.dirname(CODE_DIR), "config")
DEFAULT_CONFIG = os.path.join(CODE_DIR, "project.conf")


def read_config(args):
    """Read project config."""
    config = configparser.ConfigParser()
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

    parser_m = subparsers.add_parser('make', help='generate an sbatch jobscript')
    parser_m.add_argument("-e", action='store_true')

    parser_fe = subparsers.add_parser('fetch', help='download references from NCBI')
    parser_fe.add_argument("-s", dest=fe_srch_file, help='A file with search strings to download.')

    parser_b = subparsers.add_parser('build-index', help='build a bowtie2 index')
    parser_b.add_argument("--opt3", action='store_true')
    parser_b.add_argument("--opt4", action='store_true')

    parser_s = subparsers.add_parser('submit', help='submit job to slurm')
    parser_s.add_argument("-opt9", action='store_true')
    parser_s.add_argument("--opt10", action='store_true')

    parser_f = subparsers.add_parser('full', help='perform all the above steps')
    parser_f.add_argument("-opt5", action='store_true')
    parser_f.add_argument("--opt6", action='store_true')

    parser_c = subparsers.add_parser('collect', help='collect output cluster data')
    parser_c.add_argument("--opt7", action='store_true')
    parser_c.add_argument("--opt8", action='store_true')

    return parser

    # parser.add_argument('-m', '--mode', dest='mode', type=str,
    #                     choices=['single', 'paired'], default='paired',
    #                     help="Sequencing read mode")

def parse_args(args):
    """Parse command line arguments"""
    args = get_parser().parse_args(args)
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
#         )
#         f.write(output)

def generate_jobscript(config):
    """Generate script for differential expression and report generation."""
    env = Environment(loader=FileSystemLoader("templates"))
    template = env.get_template("jobscript")

    with open("jobscript", "w") as f:
        output = template.render(
            project=config['Project']['ProjectID'],
            cluster=config['Project']['Cluster'],
            job_name=config['Project']['JobName'],
            job_nodes=config['Project']['JobNodes'],
            cpu_cores=config['Project']['CpuCores'],
            estimated_time=config['Project']['EstimatedTime'],
            
            node_path=config['Directories']['Node'],
            ref_path=config['Directories']['References'],
            data_path=config['Directories']['Data'],
            output_path=config['Directories']['Output'],
            doric_path=config['Directories']['DoriC'],

            ref_name=config['Other']['RefName'],
            nr_samples=config['Other']['NrSamples'],
            samples_per_node=config['Other']['SamplesPerNode'],
            email=config['Other']['Email']
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

def generate_bt2_build_command(config):
    """Generate command for building a bowtie2 index."""
    referenceList="";
    for file in os.listdir(config['References']):
        if file.endswith(".fasta"):
            referenceList=referenceList+","+file

    cmd = "bowtie2-build --large-index " + referenceList + os.join(config['References'],"Index",config['RefName'])
    return cmd

def generate_fetch_command(config):
    """Generate command for downloading reference fastas and headers."""
    cmd = "bin/fetchSeq.py -e {email} -t True -d {data_path}"
    return cmd.format(**config)

def generate_sbatch_command(config):
    """Generate command for scheduling all sample runs."""
    cmd = "sbatch --array={0}%6 run_sample.sh"
    job_range = "1-{0}".format(len(project.samples))

    return cmd.format(job_range)

def generate_collect_command(config):
    cmd="./PTRMatrix.py {data_path} {ref_path} {doric_path}"
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

    if(args.subparser_name=='full' or args.subparser_name=='fetch'):
        process = subprocess.Popen(generate_fetch_command(config), shell=True, stdout=subprocess.PIPE)
        process.wait()

    if(args.subparser_name=='full' or args.subparser_name=='build-index'):
        process = subprocess.Popen(generate_bt2_build_command(config), shell=True, stdout=subprocess.PIPE)
        process.wait()
        process = subprocess.Popen(["bin/changeTID.sh",config['Directories']['References']], shell=True, stdout=subprocess.PIPE)
        process.wait()

    if(args.subparser_name=='full' or args.subparser_name=='make'):
        generate_jobscript(config)

    if(args.subparser_name=='full' or args.subparser_name=='submit'):
        process = subprocess.Popen(generate_sbatch_command(config), shell=True, stdout=subprocess.PIPE)

    if(args.subparser_name=='collect'):
        process = subprocess.Popen(generate_collect_command(config), shell=True, stdout=subprocess.PIPE)

    # prov = provider.get_provider(args['provider_name'], args['project_dir'])
    # p = project.Project(args, prov)
    # p.verify_directory_structure()
    # p.read_samples()
    # p.validate_de_input()
    # generate_scripts(p, config)
    # print_instructions(p)


if __name__ == "__main__":
    main()
