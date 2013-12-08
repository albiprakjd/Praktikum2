#!/usr/bin/python


import time
import requests
import TreeConstruction
from Bio import AlignIO
from Bio import Cluster
#"home/sileadim/workspace/Praktikum2/praktikum2.fasta"





url = "http://www.ebi.ac.uk/Tools/services/rest/muscle/run/"
#handle = open("praktikum2.fasta", "rU")


params = {

  'email': "ehmann.christopher@gmail.com",
  'tree': None,
  'order': 'aligned',
  'format':'phyi'
  
   
}
files = {'sequence': open("praktikum2.fasta.txt","r")}

print "Requesting alignment from webserver..."
req = requests.post(url, data = params, files = files) 
jobId = req.text
print "jobId"
print jobId

checkUrl = "http://www.ebi.ac.uk/Tools/services/rest/muscle/status/" + jobId

while requests.get(checkUrl).text != "FINISHED":
    print "Waiting..."
    time.sleep(5)
    
 
 
resultUrl = "http://www.ebi.ac.uk/Tools/services/rest/muscle/result/" + jobId + "/aln-phylip_interleaved"
resp = requests.get(resultUrl)
print resp.text

alnFile = open('alignment.phy', "w+")
alnFile.write(resp.text)
alnFile.close()

aln = AlignIO.read("alignment.phy", 'phylip')
calculator = TreeConstruction.DistanceCalculator('identity')
dm = calculator.get_distance(aln)
print dm
