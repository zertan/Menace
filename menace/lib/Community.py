#!/usr/bin/env python
from collections import namedtuple
from os.path import join, exists
from os import makedirs#, chdir,getcwd,fchmod
from glob import glob
import shutil
import numpy as np
#import pandas as pd
import re
#import scipy as sc
#from scipy.stats import norm
#from scipy.stats import expon
import scipy.interpolate as interpolate
import scipy.integrate
#from scipy.fftpack import fft, ifft
#from lmfit.models import Model
#from lmfit import conf_interval
#import matplotlib
#import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
#from Bio.Seq import MutableSeq
#from Bio.Alphabet import IUPAC
import docker
import configparser

# vars
#cli = docker.Client(base_url='unix://var/run/docker.sock')

class Community:
    "A community defined by genome references (Biopython SeqRecords) and corresponding growth parameters."
    def __init__(self,name,acc,growth_param,td,mapper,image,env,email):
        self.name = name
        self.conf = local_conf(join(td,name),mapper,email)
        self.d_conf = local_conf(join('/mnt/vol',name),mapper,email)
        self.args0=['ptr_pipeline.py','-c',join(self.d_conf['node_path'],'project.conf')]
        self.env=env
        self.image=image
        
        self.create_dirs(True)
        save_config(self.conf,self.d_conf)
        self.fetch_ref(acc)
        
        self.pop = self.init_population(acc,growth_param)
        self.distribution = self.community_distribution()
        self.samples = []
        
        def init_population(self,acc,growth_param):
            records=open_records(glob(join(self.conf['ref_path'],'Fasta','*.fasta')))
            population = namedtuple("population", "B C D l seq cells")

            #refId=[x.id for x in records]
            #acc=keys(comm)
            #refInd=[acc.index(x) for x in refId]

            pop = {}

            #for i,rec in enumerate(records):
            # add [:-2] to a if . removed 
            for i,a in enumerate(acc): 
                pop[a]=population(  B = growth_param[i][0],C = growth_param[i][1],
                                    D = growth_param[i][2], l = len(records[a]),
                                    seq = records[a], cells=growth_param[i][3])
            return pop
    
        def ptr(self):
            for i,a in enumerate(self.pop.keys()):
                print a+": "+str(2**growth_param[i][1])
    
        def Asnok(R,C,D,l):
            return R/(Gekv(C,D) * l)


    def community_distribution(self):
        d=np.array([Gekv(p.B,p.C)*p.l*p.cells for p in self.pop.values()])
        return d/d.sum()

    #def ab_(self):   
    #    return 
    
    def sample(self,nr_samples):
        nr_samples=np.array(nr_samples*self.distribution)
        nr_samples=nr_samples.astype(np.int)
        
        samp_run=[]
        for i,p in enumerate(self.pop.keys()):
            samp_run.append(inverse_transform_sampling(self.pop[p].C,nr_samples[i],self.pop[p].l))
        
        self.samples.append(samp_run)
    
    def write_reads(self):
        for i,samp in enumerate(self.samples):
            write_reads(samp,self,self.conf['data_path'],self.name+str(i))
    
    def compare_fit(self):
        if not self.samples:
            print "The community is not sampled, please run community.sample(nr_samples)"
            return;
        
        err_hfit=[]
        err_pfit=[]
        res_fit=[]
        
        for i,samp in enumerate(self.samples):
            for acc in self.pop.keys():
                try:
                    depth_file=join(self.conf['output_path'],self.name+str(i),'npy',acc+'.depth.npy')
                    best_file=join(self.conf['output_path'],self.name+str(i),'npy',acc+'.depth.best.npy')

                    signal=2**(np.load(depth_file))
                    signal=signal/(signal.sum()/len(signal))#*self.pop[acc].l)
                    from_ptr=np.load(best_file)

                    res=fit_signal(signal,self.pop[acc].l)
                    res_fit.append(res)
                    
                    err_hfit.append((self.pop[acc].C-res.best_values['C'])/self.pop[acc].C)

                    err_pfit.append((self.pop[acc].C-from_ptr[2]+from_ptr[3])/self.pop[acc].C)
                    print "Simulated value: "+str(self.pop[acc].C)
                    print "Error from this fit: "+str(err_hfit[-1])+ ' value: ' + str(res.best_values['C'])
                    print "Error from initial PTR fit "+str(err_pfit[-1])+' value: ' + str(from_ptr[2]-from_ptr[3])
                except Exception as Ex:
                    print "Ex"
                    pass
        return [res_fit,err_hfit,err_pfit]
    
    
    def fetch_ref(self,acc=''):
        acc_path = join(self.conf['node_path'],"acc")
        f = open(acc_path, "w")
        
        if not acc:
            f.write("\n".join(self.pop.keys()))
        else:
            f.write("\n".join(acc))
        f.close()
        
        create_mount_run(self.image,td,self.args0+['fetch-references','-s',join(self.d_conf['node_path'],'acc')],self.env)

    def build_index(self):
        create_mount_run(self.image,td,self.args0+['build-index'],self.env)

    def run_pipeline(self):
        create_mount_run(self.image,td," ".join(['/bin/bash -c "']+self.args0+['make']+[';']+self.args0+['run"']),self.env)
        
    def collect(self):
        create_mount_run(self.image,td,self.args0+['collect'],self.env)
        
    def create_dirs(self,from_init=False):
        if exists(self.conf['node_path']) and from_init:
            shutil.rmtree(self.conf['node_path'])
        
        for d in [self.conf['node_path'],self.conf['data_path'],self.conf['output_path'],self.conf['ref_path']]:
            if not exists(d):
                makedirs(d)

def create_mount_run(cli,image,mount_dir,cmd,envs):
    if envs:    
        container = cli.create_container(
            image=image, command=cmd, volumes=['/mnt/vol'],
            host_config=cli.create_host_config(binds={
                mount_dir: {
                    'bind': '/mnt/vol',
                    'mode': 'rw',
                }
            }),
            environment=envs
        )
    else:
        container = cli.create_container(
            image=image, command=cmd, volumes=['/mnt/vol'],
            host_config=cli.create_host_config(binds={
                mount_dir: {
                    'bind': '/mnt/vol',
                    'mode': 'rw',
                }
            })
        )
    ctr=container.get('Id')
    cli.start(ctr)
    cli.wait(ctr,60*60*24*10)
    return cli.logs(ctr)

def local_conf(td,mapper,email,cores):
      return {
            'project': 'ptr_simulation',
            'cluster': '',
            'job_name': 'ptr_simulation',
            'job_nodes': '1',
            'cpu_cores': cores,
            'estimated_time': '',

            'node_path': td,
            'ref_path': join(td,'References'),
            'data_path': join(td,'Data'),
            'output_path': join(td,'Out'),
            'doric_path': join(td,'DoriC'),

            'mapper': mapper,
            'ref_name': 'sim',
            'nr_samples': '1',
            'samples_per_node': '1',
            'email': email,
            'data_prefix': '',
            'start_ind': '1',
            'job_range': '1',

            'ftp_url': '',
            'data_url': ''
        }

def save_config(lconf,conf):
    Config = configparser.ConfigParser()
    Config.optionxform = str
    
    cfgfile = open(join(lconf['node_path'],'project.conf'),'w')

    Config.add_section('Project')
    Config.set('Project','ProjectID',conf['project'])
    Config.set('Project','Cluster',conf['cluster'])
    Config.set('Project','JobName',conf['job_name'])
    Config.set('Project','JobNodes',conf['job_nodes'])
    Config.set('Project','CpuCores',conf['cpu_cores'])
    Config.set('Project','EstimatedTime',conf['estimated_time'])
    
    Config.add_section('Directories')
    Config.set('Directories','Node',conf['node_path'])
    Config.set('Directories','References',conf['ref_path'])
    Config.set('Directories','Data',conf['data_path'])
    Config.set('Directories','Output',conf['output_path'])
    Config.set('Directories','DoriC',conf['doric_path'])
    
    Config.add_section('Other')
    Config.set('Other','Mapper',conf['mapper'])
    Config.set('Other','RefName',conf['ref_name'])
    Config.set('Other','NrSamples',conf['nr_samples'])
    Config.set('Other','SamplesPerNode',conf['samples_per_node'])
    Config.set('Other','Email',conf['email'])
    Config.set('Other','DataPrefix',conf['data_prefix'])
    Config.set('Other','StartInd',conf['start_ind'])
    Config.set('Other','JobRange',conf['job_range'])
    Config.set('Other','FtpURL',conf['ftp_url'])
    Config.set('Other','DataURL',conf['data_url'])
    
    Config.write(cfgfile)
    cfgfile.close()

def write_reads(samples,comm,directory,name):
    f1 = open(join(directory,name+"_1.fastq"), "w")
    f2 = open(join(directory,name+"_2.fastq"), "w")
    
    for i,p in enumerate(comm.pop.keys()):
        for j,pos in enumerate( samples[i].tolist()):
            r1,r2 = read_pair(name+"_"+p,str(j+1),comm.pop[p].seq,pos,comm.pop[p].l)

            SeqIO.write(r1, f1, "fastq")
            SeqIO.write(r2, f2, "fastq")
    
    f1.close()
    f2.close()

def open_records(fasta):
    records={};
    for f in fasta:
        handle = open(f, "rU")
        tmp=list(SeqIO.parse(handle, "fasta"))
        m=re.match('.*(N[CTZ]_([A-Z]{2})*[0-9]{6}).*',tmp[0].id)
        records[m.group(0)]=tmp[0].seq
        handle.close()
    
    return records

def read_pair(header,ind,seq,pos,l):
    def circular_yield(x,pos,sub,l):
        if (pos>l):
            return x[pos-l:pos-l+sub]
        elif (pos+sub>l):
            r=pos+sub-l
            x2=x[pos:pos+sub-r]
            return x2+x[0:r]
        else:
            return x[pos:pos+sub]

    #base_error_rate = .02
    #mutation_rate = .001
    #nuc=["A","C","G","T"]
    
    # generate a normally distributued insert size of mean 300 and sd 40 
    #insert_size = int(norm.rvs(size=1,loc=300,scale=30)[0])
    insert_size = int(500)
    
    r1 = circular_yield(seq,pos,100,l)
    r2 = circular_yield(seq,pos+insert_size,100,l)
    
    # flip reads according to seqs error and mutation rate
    #ra=[]
    #for r in [r1,r2]:
    #    rr=np.random.random(100)
    #    m_ind=[i for i,j in enumerate(rr) if j < base_error_rate]
    #    char=np.random.choice(nuc,len(m_ind))
    #    a=list(r)
    #    for i,ind in enumerate(m_ind):
    #        a[ind]=char[i] 
    #    ra.append("".join(a))
        #r_tmp=r_tmp[:m_ind[0]]
        #for i,index in enumerate(m_ind[1:]):
        #    r_tmp = r_tmp[:index] + char[i] + r[index + 1:]
    #[r1,r2]=ra

    #nrs=[str(np.random.randint(low=1000,high=2000)),str(np.random.randint(low=10**4,high=2*10**4)),str(np.random.randint(low=10**5,high=2*10**5))]
    rec1=SeqRecord(r1,id=header+"_"+str(pos)+"."+str(ind),description="RANDXXLAB"+str(pos)+"_"+str(pos+insert_size+100)+":0:0"+"/1")
    rec2=SeqRecord(r2,id=header+"_"+str(pos)+"."+str(ind),description="RANDXXLAB"+str(pos)+"_"+str(pos+insert_size+100)+":0:0"+"/2")
    #rec1=SeqRecord(r1,id=header+"."+str(ind),description="FCB00YLABXX:6:"+nrs[0]+":"+nrs[1]+":"+nrs[2]+"/1")
    #rec2=SeqRecord(r2,id=header+"."+str(ind),description="FCB00YLABXX:6:"+nrs[0]+":"+nrs[1]+":"+nrs[2]+"/2")
    
    #rec1=SeqRecord(r1,id=header+"_"+str(pos)+"_"+str(pos+insert_size)+"/1",description="")
    #rec2=SeqRecord(r2,id=header+"_"+str(pos)+"_"+str(pos+insert_size)+"/2",description="")
    rec2=rec2.reverse_complement(id=header+"_"+str(pos)+"."+str(ind),description="RANDXXLAB"+str(pos)+"_"+str(pos+insert_size+100)+":0:0"+"/2")
    #rec2=SeqRecord(r2,id=header+"."+str(ind)+"/2",description="")
    rec1.letter_annotations['phred_quality']=[17]*100
    rec2.letter_annotations['phred_quality']=[17]*100
    
    #rec1.description=
    
    return(rec1,rec2)

def inverse_transform_sampling(C,n_samples,l=1):
    x = np.linspace(0, float(1)/1.7, 100)
    y = np.linspace(pxn(.5,C,1),pxn(0,C,1),100)

    pdf=interpolate.interp1d(x,pxn(x,C,1))
    cdf = [scipy.integrate.quad(lambda x: pdf(x),0,i)[0] for i in np.linspace(0,.5,100)]
    inv_cdf = interpolate.interp1d(cdf, np.linspace(0,.5,100))
    
    r = np.random.rand(n_samples)
    v=l*np.round(inv_cdf(r[:int(n_samples/2)]),4)
    v2=l*(1-np.round(inv_cdf(r[int(n_samples/2):]),4))
    v=np.concatenate([v,v2])
    
    return v.astype(np.int)

def fit_signal(signal,l):
    x1=np.linspace(0,1,len(signal))
    
    piecewiseModel=Model(piecewise_prob)
    piecewiseModel.set_param_hint('l',     value=1,vary=False)
    piecewiseModel.set_param_hint('C',   vary=True,  value=.1,min=0,max=1)
    piecewiseModel.make_params()
    
    res=piecewiseModel.fit(signal,x=x1,weights=np.sin(1.5*x1)+1.5)
    return res

def piecewise_prob(x,C,l):
    conds=[(0 < x) & (x <= float(l)/2), x > float(l)/2]
    funcs=[lambda x: pxn(x,C,l)/float(2),lambda x: pxn((l-x),C,l)/float(2)]
    return np.piecewise(x,conds,funcs)

def Gekv(C,D):
    return 2**(C + D)/(C*np.log(2))*(1 - 2**(-C))

def pxn(x, C, l):
    return (2**(1 + C - (2*C*x)/float(l))*C*np.log(2))/(float(l)*(-1 + 2**C))

def mv_filt(L,omega):
    return (1/float(L))*(1-np.exp(-omega*1j*L))/(1-np.exp(-omega*1j))

def Asnok(R,C,D,l):
    return R/(Gekv(C,D) * l)