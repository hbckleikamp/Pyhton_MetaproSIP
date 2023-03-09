# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 11:39:26 2023

@author: ZR48SA
"""


#%% set base path
from pathlib import Path
import os
from inspect import getsourcefile
# change directory to script directory (should work on windows and mac)
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
script_dir=os.getcwd()
print(os.getcwd())
basedir=os.getcwd()

Temporary_directory  =basedir  #Writing directory for temporary indices of MSFragger and Diamond, make sure there is enough space here
Output_directory     =basedir  #Directory where all generated outputs are written to
base_command="cd "+basedir+" && "

#%% Modules
import subprocess
import pandas as pd
import subprocess
import psutil
import shutil
import pandas as pd
import numpy as np
import re
import math

#%% Functions

############ Classes #############

class subprocess_output: #monitoring outputs of functions that run a subprocess
    def __init__(self,command,stdout,stderr):
        self.command=command
        self.stdout=stdout
        self.stderr=stderr


def MSFragger_annotation(input_files,   #full_path, .mzML
                         database_path, #full_path, .fasta
                         output_folder="Spectrum_annotation"
                         ):

    MSFragger_jar_path=str(Path(basedir,"MSFragger-3.5.jar"))  
    pep_split_path=str(Path(basedir,"msfragger_pep_split.py"))
    
    #rewrite closed_fragger.params according to database path
    params_path=str(Path(basedir,"closed_fragger.params"))
    with open(params_path,"r+") as f:
        lines=f.readlines()
        lines=["database_name = "+database_path+" #database name here\n" if line.startswith("database_name =") else line for line in lines]
        f.seek(0)
        f.writelines(lines)
            
    
    os.chdir(Temporary_directory) #change to tempdir for writing indices
    stderr=""
    retry=0
    while True: #rerun with different settings untill settings are found that fit RAM
        
        #remove old peptide split indices
        if os.path.exists(str(Path(Temporary_directory,"split_peptide_index_tempdir"))): shutil.rmtree(str(Path(Temporary_directory,"split_peptide_index_tempdir"))) 
        JavaMem=int(psutil.virtual_memory().available/10**9*0.8) #check available RAM
        splits=math.ceil(os.path.getsize(database_path)/10**9*150/JavaMem) #starting number of splits for db splitting, database size*150 is an estimated scaling factor, which will be effected by the number of allowed modifications and missed cleavages
        batch_size=int(psutil.disk_usage(Path(pep_split_path).anchor).free/10**9/50) #number of files annotated per time, 50 GB is the estimated temporary writing memory allocated per file
        no_batches=int(len(input_files)/math.ceil(len(input_files)/batch_size))
        batches=[input_files[i:i + no_batches] for i in range(0, len(input_files), no_batches)]
        print("running MSFragger, total number of splits: "+str(splits)+" , total number of batches: "+str(len(batches))+", Java Heap: "+str(JavaMem))
        
        retry+=1
        if retry>1:
            print("retrying")
        
        if retry>5:
            print("retried 5 times, unknown error")
            "i"+1 #hard exit
            
        command="cd" +' "'+basedir+'" && '
        command+=" && ".join(["".join([' python ',
                                  ' "'+pep_split_path+'" ',
                                  str(splits),
                                  ' "'+"java -jar -Xmx"+str(JavaMem)+"G"+'" ',
                                  ' "'+MSFragger_jar_path+'" ',
                                  ' "'+params_path+'" ',
                                
                                  #" MSFragger-3.5.jar closed_fragger.params ", #can also be written as full paths without the "cd" 
                                  '"'+'" "'.join(batch)+'"']) for batch in batches])
        
        print(command)
        stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        
        
        if "No space left on device" in str(stderr): #retry with smaller file batch
            print("out of disk space, retrying with smaller batch size")
            batch_size-=1
            if batch_size==0:
                print("unsolvable space error, exiting")
                "i"+1 #hard exit
            continue
                
        if "java.lang.OutOfMemoryError" in str(stderr): #retry with more splits
            print("out of Java memory, retrying with more splits")
            splits+=50
            continue
        
        if "returned non-zero exit status" not in str(stderr):
            break

    os.chdir(basedir) #change back to basedir
    
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    pepxml_files=[]
    for mzml in input_files:
        pepxml=str(Path(Path(mzml).parents[0],Path(mzml).stem))+".pepXML" #files are written to the location of the mzml files, and are then moved to output folder
        shutil.move(pepxml, str(Path(Output_directory,output_folder,Path(mzml).stem+".pepXML")))
        pepxml_files.append(str(Path(Output_directory,output_folder,Path(mzml).stem+".pepXML")))
    
    com=subprocess_output(command,stdout,stderr)
    com.output_files=pepxml_files
    return com

def raw2mzML(raw_file): #

    output_folder=str(Path(Output_directory,"HB_mzML"))
    if not os.path.exists(output_folder): os.mkdir(output_folder)

    command="cd" +' "'+output_folder+'" && msconvert '
    command+='"'+raw_file+'"' 
    command+=' --mzML --filter "peakPicking vendor" --filter "zeroSamples removeExtra" --filter "titleMaker Run: <RunId>, Index: <Index>, Scan: <ScanNumber>"'
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    
    output_file=str(Path(output_folder,Path(raw_file).stem+".mzML"))
    return output_file

def MSGFPlusAdapter(file,database_path,msgf_path,output_folder="Spectrum_annotation",args=""):
    output_file=str(Path(Output_directory,output_folder,Path(file).stem+".idXML"))
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    command=base_command+"MSGFPlusAdapter -in "+'"'+file+'"'+" -out "+'"'+output_file+'"'+" -executable "+'"'+msgf_path+'"'+" -database "+'"'+database_path+'"'
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return output_file

def DecoyDatabase(database,args=""):
    output_file=str(Path(basedir,Path(database).stem+"_decoy.fa"))
    command=base_command+"DecoyDatabase -in "+'"'+database+'"'+" -out "+'"'+output_file+'"'
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return output_file    
    
def FeatureFinderMultiplex(file,output_folder="Features",args=""):
    output_file=str(Path(Output_directory,output_folder,Path(file).stem+".featureXML"))
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    command=base_command+"FeatureFinderMultiplex -in "+'"'+file+'"'+" -out "+'"'+output_file+'"'
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return output_file

def FeatureFinderMultiplex(file,output_folder="Features",args=""):
    output_file=str(Path(Output_directory,output_folder,Path(file).stem+".featureXML"))
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    command=base_command+"FeatureFinderMultiplex -in "+'"'+file+'"'+" -out "+'"'+output_file+'"'
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return output_file

def HighResPrecursorMassCorrector(file,feature_file,output_folder="Corrected_mzML",args=""):
    output_file=str(Path(Output_directory,output_folder,Path(file).name))
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    command=base_command+"HighResPrecursorMassCorrector -in "+'"'+file+'"'+" -feature:in "+'"'+feature_file+'"'+" -out "+'"'+output_file+'"'
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return output_file

def idconvert(file,output_folder="Spectrum_annotation",args=""): #proteowizard
    output_file=str(Path(Output_directory,output_folder,Path(file).stem+".mzID"))
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    command=base_command+"idconvert "+'"'+file+'"'+" -o "+'"'+output_folder+'"'
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return output_file

def IDFileConverter(file,output_folder="Spectrum_annotation",args=""): 
    output_file=str(Path(Output_directory,output_folder,Path(file).stem+".idXML"))
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    command=base_command+"IDFileConverter -in "+'"'+file+'"'+" -out "+'"'+output_file+'"'
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return output_file

def PeptideIndexer(file,database,output_folder="Spectrum_annotation",args=""):
    output_file=str(Path(Output_directory,output_folder,Path(file).name))
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    command=base_command+"PeptideIndexer -in "+'"'+file+'"'+" -fasta "+'"'+database+'"'+" -out "+'"'+output_file+'"'
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return output_file

def Proteininference(file,output_folder="Spectrum_annotation",args=""):
    output_file=str(Path(Output_directory,output_folder,Path(file).name))
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    command=base_command+"Proteininference -in "+'"'+file+'"'+" -out "+'"'+output_file+'"'
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return output_file

def FalseDiscoveryRate(file,output_folder="Spectrum_annotation",args=""):
    output_file=str(Path(Output_directory,output_folder,Path(file).name))
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    command=base_command+"FalseDiscoveryRate -in "+'"'+file+'"'+" -out "+'"'+output_file+'"'
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return output_file

def IDFilter(file,output_folder="Spectrum_annotation",args=""):
    output_file=str(Path(Output_directory,output_folder,Path(file).name))
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    command=base_command+"IDFilter -in "+'"'+file+'"'+" -out "+'"'+output_file+'"'
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return output_file

def InternalCalibration(mzML_file,idXML_file,output_folder="Corrected_mzML",args=""):
    output_file=str(Path(Output_directory,output_folder,Path(mzML_file).name))
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    command=base_command+"InternalCalibration -in "+'"'+mzML_file+'"'+" -cal:id_in  "+'"'+idXML_file+'"'+" -out "+'"'+output_file+'"'
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return output_file

def MapAlignerIdentification(files,output_folder="Alignment",args=""):
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    in_files=" ".join(['"'+str(Path(file))+'"' for file in files])
    
    output_files=  [str(Path(Output_directory,output_folder,Path(file).name))             for file in files]
    trafo_outfiles=[str(Path(Output_directory,output_folder,Path(file).stem+".trafoXML")) for file in files]
    
    # output_files=" ".join(['"'+str(Path(Output_directory,output_folder,Path(file).name))+'"' for file in files])
    # trafo_outfiles=" ".join(['"'+str(Path(Output_directory,output_folder,Path(file).stem))+".trafoXML"+'"' for file in files])
    command=base_command+"MapAlignerIdentification -in "+in_files+" -trafo_out "+'"'+'" "'.join(trafo_outfiles)+'"'+" -out "+'"'+'" "'.join(output_files)+'"'
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    #return output_files.split('" "'),trafo_outfiles.split('" "')
    return output_files,trafo_outfiles

def MapRTTransformer(file,trafofile,output_folder="Alignment",args=""):
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    output_file=str(Path(Output_directory,output_folder,Path(file).name))
    command=base_command+"MapRTTransformer -in "+'"'+file+'"'+" -trafo_in  "+'"'+trafofile+'"'+" -out "+'"'+output_file+'"'
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return output_file
    
def IDMapper(aligned_idXML_file,aligned_feature,output_folder="Alignment",args=""):
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    output_file=str(Path(Output_directory,output_folder,Path(aligned_feature).name))
    command=base_command+"IDMapper -in "+'"'+aligned_feature+'"'+" -id  "+'"'+aligned_idXML_file+'"'+" -out "+'"'+output_file+'"'
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return output_file

def FeatureLinkerUnlabeledQT(files,output_folder="Alignment",args=""):
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    in_files=" ".join(['"'+str(Path(file))+'"' for file in files])
    output_file=str(Path(Output_directory,output_folder,"features.consensusXML")) 
    command=base_command+"FeatureLinkerUnlabeledQT -in "+in_files+"  -out "+'"'+output_file+'"'
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return output_file

def IDConflictResolver(file,output_folder="Alignment",args=""):
    output_file=str(Path(Output_directory,output_folder,Path(file).name))
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    command=base_command+"IDConflictResolver -in "+file+" -out "+output_file
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return output_file

def FileConverter(file,output_folder="Alignment",args=""):
    output_file=str(Path(Output_directory,output_folder,"Pooled.featureXML"))
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    command=base_command+"FileConverter -in "+'"'+file+'"'+" -out "+'"'+output_file+'"'
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return output_file
        
    
def MetaProSIP(file,
               database,
               pooled_features,
               output_folder="MetaproSIP_results",args=""):
    
    
    out_csv=str(Path(Output_directory,output_folder,Path(file).stem+".csv"))
    out_peptide_centric_csv=str(Path(Output_directory,output_folder,Path(file).stem+"peptide_centric.csv"))           
                
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    command=base_command+"".join([" MetaProSIP ",
                                  " -in_mzML "+'"'+file+'"',
                                  " -in_fasta "+'"'+database+'"',
                                  " -in_featureXML "+'"'+pooled_features+'"',
                                  " -out_csv "+'"'+out_csv+'"',
                                  " -out_peptide_centric_csv "+'"'+out_peptide_centric_csv+'"'])
    
    
    command+=args #additional arguments
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return out_csv, out_peptide_centric_csv
   

#%%

#input variables
mzml_files=[
#"C:/Sippy/HB_mzML/Run4_Ecoli268_R1a_0_13C_400ng.mzML", this file had a parsing error
"C:/Sippy/HB_mzML/Run4_Ecoli268_R1a_0-01_13C_400ng.mzML",
"C:/Sippy/HB_mzML/Run4_Ecoli268_R1a_0-1_13C_400ng.mzML",
"C:/Sippy/HB_mzML/Run4_Ecoli268_R1a_0-025_13C_400ng.mzML",
"C:/Sippy/HB_mzML/Run4_Ecoli268_R1a_0-25_13C_400ng.mzML",
"C:/Sippy/HB_mzML/Run4_Ecoli268_R1a_1_13C_400ng.mzML",
"C:/Sippy/HB_mzML/Run4_Ecoli268_R1a_5_13C_400ng.mzML",
"C:/Sippy/HB_mzML/Run4_Ecoli268_R1a_10_13C_400ng.mzML"
]
mzml_files.sort() #does the order of files matter?

database="C:/Sippy/Datasets/PXD023693 (Ecoli C1-6)/Database/target_decoy.fa"

#Also requires R to be set up and MSFragger to be installed

#add decoy to database
#database=DecoyDatabase(database)


#Loop 1: correct Feature mz to detected features or MS1 trace peaks
#input: raw files
feature_files  =[FeatureFinderMultiplex(file)                          for    file in            mzml_files]
corrected_mzMLs=[HighResPrecursorMassCorrector(file,feature_files[ix]) for ix,file in enumerate( mzml_files)]    

#Loop 2: Peptide annotation
#MSFragger
# MSFragger_files=MSFragger_annotation(corrected_mzMLs, # MSGFPlusAdapter (replace with MSFragger)
#                                      database).output_files
# mzID_files=[idconvert(file) for file in MSFragger_files]
# idXML_files=[IDFileConverter(file) for file in mzID_files]

#MSGF
msgf_path="C:/MSGF/MSGFPlus_v20230112/MSGFPlus.jar"
idXML_files=[MSGFPlusAdapter(file,database_path=database,msgf_path=msgf_path,args=" -java_memory 20000 ") for file in corrected_mzMLs]    

Annotated_spectra=[IDFilter(FalseDiscoveryRate(
    Proteininference(
        PeptideIndexer(file,database)),args=" -FDR:PSM 0.05 "),args=" -score:pep 0.01 ") for file in idXML_files]

#filter



#Loop3: Internal Calibration
#input: loop 1 output, loop2 output, R
corrected_mzMLs=[InternalCalibration(file,Annotated_spectra[ix]) for ix,file in enumerate( mzml_files)] 

# Loop 4: Alignment and consensus features
#Align IDs, input: Loop 2 output
aligned_idXMLs,trafo_files=MapAlignerIdentification(Annotated_spectra)
#Align Features, input: features, trafoXML
aligned_features=[MapRTTransformer(file,trafo_files[ix]) for ix,file in enumerate(feature_files)] 
#Align mzMLs, input: mzML, trafoXML
aligned_mzMLs=[MapRTTransformer(file,trafo_files[ix]) for ix,file in enumerate(corrected_mzMLs)] 
#annotate aligned features with aligned idXML
annotated_features=[IDMapper(file,aligned_features[ix]) for ix,file in enumerate(aligned_idXMLs)]

#get consensus features
consensus_file=IDConflictResolver(FeatureLinkerUnlabeledQT(annotated_features))
pooled_features=FileConverter(consensus_file)

#Loop 5. MetaproSIP
#annotate aligned features with aligned idXML
result=[MetaProSIP(file,database,pooled_features) for file in aligned_mzMLs]

#%%
