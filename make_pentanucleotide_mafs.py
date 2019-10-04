import os
import subprocess

geneMafDir = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/simulatedMafs/geneMutMafs'
outFileDir = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/simulatedMafs/allPossibleGeneMutMafWithPentanucleotideContext'
for geneMaf in os.listdir(geneMafDir):
  mafPath = None
  gene = None
  if 'v2.maf' in geneMaf:
    mafPath = os.path.join(geneMafDir, geneMaf)
    gene = geneMaf.split('_')[0]
  if mafPath != None:
    cmdStr = 'python /home/friedman/friedman/mutation-signatures/make_n-nuc_maf.py'
    outfilePath = os.path.join(outFileDir, gene + '_penta.maf')
    #THE command to pass to my adjusted trinuc function
    cmd = ' '.join([cmdStr, mafPath, outfilePath])
    print 'executing command: ', cmd
    subprocess.Popen(cmd, shell=True).wait()


