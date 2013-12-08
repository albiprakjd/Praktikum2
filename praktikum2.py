#!/usr/bin/python


import time
import requests
import TreeConstruction
from Bio import AlignIO
import copy
from Bio.Phylo import BaseTree
from Bio.Align import MultipleSeqAlignment
import sys
#"home/sileadim/workspace/Praktikum2/praktikum2.fasta"


def makeMSA(fastaname):


    url = "http://www.ebi.ac.uk/Tools/services/rest/muscle"
    runUrl = url + "/run/"
#handle = open("praktikum2.fasta", "rU")


    params = {

  'email': "ehmann.christopher@gmail.com",
  'tree': None,
  'order': 'aligned',
  'format':'phyi'
  
   
    }
    files = {'sequence': open(fastaname,"r")}

    print "Requesting alignment from webserver..."
    req = requests.post(runUrl, data = params, files = files) 
    jobId = req.text
    print "JobId: "+jobId

    checkUrl = url + "/status/" + jobId

    while requests.get(checkUrl).text != "FINISHED":
        print "Waiting..."
        time.sleep(5)
    
 
 
        resultUrl = url + "/result/" + jobId + "/aln-phylip_interleaved"
        resp = requests.get(resultUrl)
        filename = jobId + 'alignment.phy'
        alnFile = open(filename , "w+")
        alnFile.write(resp.text)
        alnFile.close()
        return filename

def makeDistanceMatrix(filename):

    aln = AlignIO.read(filename, 'phylip')
    calculator = TreeConstruction.DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    return dm

def cladeHeight(clade):
    height = 0
    if clade.is_terminal():
        height = clade.branch_length
    else:
        height = height + max(cladeHeight(subclade) for subclade in clade.clades)
    return height

def upgma(distance_matrix):

    # make a copy of the distance matrix to be used
    dm = copy.deepcopy(distance_matrix)
        # init terminal clades
    clds = [BaseTree.Clade(None, name) for name in dm.names]
        # init minimum index
    minI = 0
    minJ = 0
    count = 0
    while len(dm) > 1:
        minDist = dm[1, 0]
            # find minimum index
        for i in range(1, len(dm)):
            for j in range(0, i):
                if minDist >= dm[i, j]:
                    minDist = dm[i, j]
                    minI = i
                    minJ = j

            # create clade
        clade1 = clds[minI]
        clade2 = clds[minJ]
        count += 1
        inClade = BaseTree.Clade(None, "Inner" + str(count))
        inClade.clades.append(clade1)
        inClade.clades.append(clade2)
                    #assign branch length
        if clade1.is_terminal():
            clade1.branch_length = minDist/2
        else:
            clade1.branch_length = minDist/2 - cladeHeight(clade1)
            #other clade
        if clade2.is_terminal():
            clade2.branch_length = minDist/2
        else:
            clade2.branch_length = minDist/2 - cladeHeight(clade2)

            # update node list
        clds[minJ] = inClade
        del clds[minI]

            # rebuild distance matrix,
            # set the distances of new node at the index of min_j
        for k in range(0, len(dm)):
            if k != minI and k != minJ:
                dm[minJ, k] = (dm[minI, k] + dm[minJ, k]) * 1.0 / 2

        dm.names[minJ] = "Inner" + str(count)

        del dm[minI]
    inClade.branch_length = 0
    return BaseTree.Tree(inClade)


if __name__ == "__main__":
    filename = sys.argv[1]   
    alnName = makeMSA(filename)
    dm = makeDistanceMatrix(alnName)
    print upgma(dm)
    
