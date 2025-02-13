# txt to FASTA file
import string
import sys

wfile_name = (sys.argv[1])
rfile_name = (sys.argv[2])

wfile = open(wfile_name, "w")
rfile = open(rfile_name, "r")

count=1
while True:
  content = rfile.readline()
  if not content:
    break
  contentList=content.split()
  countstr=str(count)
  wfile.write(">kmer" + countstr + ":" + contentList[1]+ "\n" + contentList[0]+ "\n")
  count=count+1

rfile.close()

'''
file = open("python.txt", "r")
while True:
	content=file.readline()
	if not content:
		break
	print(content)
file.close()
'''
