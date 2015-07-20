#Usage notes: 
#First argument is name of the file containing all the regions, from which random ones can be selected
#Second argument is the name of the file containing known ragged regions. If you only want random ones, set this to NONE
#Third argument is a comma-separated list of bedgraph files containing methylation data for replicates
#Fourth argument is the number of random regions that should be selected. Set to 0 to take all of the regions in the argument 1 file
#Fifth argumnet is the minimum number of CpG sites allowed in each region
#Sixth argument is the minimum region length in bp

'''
allRegionsFilename = sys.argv[1]
knownRegionsFilename = sys.argv[2]
bedgraphFilenames = sys.argv[3].split(',')
numRandRegions = int(sys.argv[4])
minCpGs = int(sys.argv[5])
minRegionLength = int(sys.argv[6])
'''

import sys, os, linecache, subprocess, re, time, json, logging
import numpy as np
import numpy.random as rand
#import matplotlib.pyplot as plt
from math import *

def lineSplit(lineString):
	return re.split(re.compile("\s+"), lineString)


class StreamToLogger(object):
   """
   Fake file-like stream object that redirects writes to a logger instance.
   """
   def __init__(self, logger, log_level=logging.INFO):
      self.logger = logger
      self.log_level = log_level
      self.linebuf = ''
 
   def write(self, buf):
      for line in buf.rstrip().splitlines():
         self.logger.log(self.log_level, line.rstrip())

def bisectionSearch(fileNameToSearch, firstLine, lastLine, regionStart, regionEnd):

	def lineSplit(lineString):
		return re.split(re.compile("\s+"), lineString)

	#If firstLine or lastLine are non-positive, return None
	if firstLine <= 0 or lastLine <= 0:
		print("One of the boundary line numbers has a non-positive value.")
		return None

	#If lastLine is greater than or equal to firstLine, return None
	if lastLine <= firstLine:
		print("The last boundary line number is less than or equal to the first boundary line number")
		return None
	
	firstLineList = lineSplit(linecache.getline(fileNameToSearch, firstLine).strip())
	lastLineList = lineSplit(linecache.getline(fileNameToSearch, lastLine).strip())

	#An empty file to search will cause a None return. This is also caused is you try and get a line that does not exist in the file.
	if firstLineList == [''] or lastLineList == ['']:
		print("The bedgraph file %s appears to be empty. This may be caused by boundary line numbers that are greater than the number of lines present in the file." % fileNameToSearch)
		return None

	if regionStart <= 0 or regionEnd <= 0:
		print("Either the region start or end are non-positive.")
		return None

	if regionEnd <= regionStart:
		print("The region end is less than or equal to the region start.")
		return None

	#If the end of the region is before the first boundary line or the start of the region is after the last boundary line, return None
	if ((regionEnd < int(firstLineList[1])) or (regionStart > int(lastLineList[1]))):
		print("The selected region has no overlap with the given methylation file. Exiting.")
		return None

	firstCpGLine = 0
	lastCpGLine = 0
	if regionStart <= int(firstLineList[1]):
		firstCpGLine = firstLine
		#Region start is before the first CpG site.
		#Found first CpG in region
	else:
		alreadyFound = False
		topLine = firstLine
		bottomLine = lastLine
		midLine = int(floor((topLine+bottomLine)/2))

		while ((bottomLine != topLine+1) & (alreadyFound == False)):
			line = lineSplit(linecache.getline(fileNameToSearch, midLine).strip())

			if int(line[1]) > regionStart:
				#Taking top half
				bottomLine = midLine
				midLine = int(floor((topLine+bottomLine)/2))
			elif int(line[1]) < regionStart:
				#Taking bottom half
				topLine = midLine
				midLine = int(floor((topLine+bottomLine)/2))
			elif int(line[1]) == regionStart:
				#Exact hit
				firstCpGLine = midLine
				alreadyFound = True

		if not alreadyFound:
			firstCpGLine = bottomLine

	if regionEnd >= int(lastLineList[1]):
		lastCpGLine = lastLine
		#Region end is after the last CpG site
		#Found last CpG in region
	else:
		alreadyFound = False
		topLine = firstCpGLine
		bottomLine = lastLine
		midLine = int(floor((topLine+bottomLine)/2))
		while ((bottomLine != topLine+1) & (alreadyFound == False)):
			line = lineSplit(linecache.getline(fileNameToSearch, midLine).strip())

			if int(line[1]) > regionEnd:
				#Taking top half
				bottomLine = midLine
				midLine = int(floor((topLine+bottomLine)/2))
			elif int(line[1]) < regionEnd:
				#Taking bottom half
				topLine = midLine
				midLine = int(floor((topLine+bottomLine)/2))
			elif int(line[1]) == regionEnd:
				#Exact hit
				#Found last CPG in region
				lastCpGLine = midLine
				alreadyFound = True

		if not alreadyFound:
			lastCpGLine = topLine

	return [firstCpGLine, lastCpGLine]

def getCpGs(regionChr, regionStart, regionEnd, bedgraphFilename):
	global chrDictByRep
	if regionChr in chrDictByRep[bedgraphFilename].keys():
		print("Chromosome locations already in dictionary")
		firstLine = chrDictByRep[bedgraphFilename][regionChr][0]
		lastLine = chrDictByRep[bedgraphFilename][regionChr][1]
	else:
		print("Searching for first and last chromosome occurrences")
		nullFds = [os.open(os.devnull, os.O_RDWR) for x in xrange(2)]
		save = os.dup(1), os.dup(2)
		getChrFirstLineCmd = "cat -n %s | grep -P \'%s\s+\' | head -n 1 | cut -f 1" % (bedgraphFilename, regionChr)
		getChrLastLineCmd = "cat -n %s | grep -P \'%s\s+\' | tail -1 | cut -f 1" % (bedgraphFilename, regionChr)
		os.dup2(nullFds[0],1)
		os.dup2(nullFds[1],2)
		firstLine = int(subprocess.check_output(getChrFirstLineCmd, shell = True))
		lastLine = int(subprocess.check_output(getChrLastLineCmd, shell = True))
		print("Done")
		chrDictByRep[bedgraphFilename][regionChr] = [firstLine, lastLine]
		os.dup2(save[0],1)
		os.dup2(save[1],2)
		os.close(nullFds[0])
		os.close(nullFds[1])
	print("Searching for CpGs")
	searchResults = bisectionSearch(bedgraphFilename, firstLine, lastLine, regionStart, regionEnd)
	if searchResults != None:
		firstCpGLine, lastCpGLine = searchResults
	else:
		return None

	outputRegionLine = ''
	outputMethList = []

	n = firstCpGLine
	while n <= lastCpGLine:
		line = ' '.join(lineSplit(linecache.getline(bedgraphFilename, n).strip()))
		outputMethList.append(line)
		n+=1
	#Found all CpGs
	
	return outputMethList

class Region:

	def __init__(self, ID, chromosome, start, end, replicateList, knownStatus, minCpGs): #replicateList should be a list of filenames, knownStatus should be 1 if we know and 0 if we don't
		self.ID = ID
		#chromosome should be a string of the form /^chr[0-9XY]+$/
		self.chr = chromosome
		self.start = int(start)
		self.end = int(end)
		print("Getting methylation data for region: %s %s %d %d" % (self.ID, self.chr, self.start, self.end))
		#Each key in this dictionary is a biological replicate, and each value is a list of CpG sites with associated methylation data
		self.methylationDict = self.getMethylationDict(replicateList, minCpGs)

		#Each key is a biological replicate, and each value is the variance of the smoothed methylation levels in that replicate
		#self.withinRepVarSmoothDict = self.getWithinRepVar()
		#self.withinRepMeanVar = np.mean(self.withinRepVarSmoothDict.values())
		#Each key is a CpG site, and each value is the variance of the methylation at that CpG site across replicates
		#self.betweenRepVarDict = self.getBetweenRepVar()
		#self.betweenRepMeanVar = np.mean(self.betweenRepVarDict.values())
		if self.methylationDict != None:
			regionNameDict = {"REGION_INFO":"%s %s %d %d" % (self.ID, self.chr, self.start, self.end), "STATUS":knownStatus}
			self.methylationDict.update(regionNameDict)
			
				
	def getMethylationDict(self, replicateList, minCpGs):
		#Function to get the CpG sites for a given region from a list of bedgraph files containing methylation data
		methylationDict = dict()
		for rep in replicateList:
			print("Getting data from replicate %s" % rep)
			methList = getCpGs(self.chr, self.start, self.end, rep)
			if methList != None:
				methylationDict[rep] = methList
			else:
				print("There was a problem obtaining methylation data for this replicate. This region has not been included.")
				return None
		return methylationDict
		
	'''
	def getWithinRepVar(self):
		if self.methylationDict == None:
			return None
		varByRepDict = dict()
		for rep in self.methylationDict.keys():
			varByRepDict[rep] = np.std([float(l.split(' ')[3]) for l in self.methylationDict[rep]])**2
		return varByRepDict

	def getBetweenRepVar(self):
		if self.methylationDict == None:
			return None
		cpgPosValDict = dict()
		for rep in self.methylationDict.keys():
			cpgPosValDict[rep] = {int(l.split(' ')[1]):float(l.split(' ')[3]) for l in self.methylationDict[rep]}
		
		getCommCpGsCmd = 'list('+'&'.join(['set(cpgPosValDict["%s"].keys())' % rep for rep in cpgPosValDict.keys()])+')'
		commonCpGs = eval(getCommCpGsCmd)
		varByCpGDict = dict()

		for pos in commonCpGs:
			varByCpGDict[pos] = np.std([cpgPosValDict[rep][pos] for rep in cpgPosValDict.keys()])**2

		return varByCpGDict
	'''


def getLineCheckLength(regionFileName, lineNumber, minRegionLength):
	regionPropList = lineSplit(linecache.getline(regionFileName, lineNumber).strip())
	if int(regionPropList[2])-int(regionPropList[1]) < minRegionLength:
		#Region too short
		return None
	else:
		return regionPropList

def checkCpGNum(methylationDict, minCpGs):
	for rep in methylationDict.keys():
		if len(methylationDict[rep]) < minCpGs:
			print("Not enough CpGs in replicate %s" % rep)
			return None
		else:
			return 1

def getRegionDict(replicateList, numRandRegions, minCpGs, minRegionLength, regionsFileName, knownStatus):

	for p in [minRegionLength, numRandRegions, minCpGs]:
		if ((not isinstance(p, int)) | (p < 0)):
			print('Please check that minRegionLength, numRandRegions, and minCpGs are non-negative integers')
			return None
	try:
		f = open(regionsFileName, 'r')
	except IOError:
		print("Can't open %s. Terminating program." % regionsFileName)
		sys.exit()
	
	outRegionDict = dict()

	whitespacePattern = re.compile('\s+')

	totalNumRegions = int(subprocess.check_output("wc -l %s | awk \'{print $1}\'" % regionsFileName, shell=True))

	if numRandRegions == 0:
		for n in range(1, totalNumRegions + 1):
			lineList = getLineCheckLength(regionsFileName, n, minRegionLength)
			if (lineList == None):
				continue
			else:
				outRegionDict[lineList[3]] = Region(lineList[3], lineList[0], lineList[1], lineList[2], replicateList, knownStatus, minCpGs)

				if outRegionDict[lineList[3]].methylationDict == None:
					del outRegionDict[lineList[3]]
					continue

				cpgCheck = checkCpGNum(outRegionDict[lineList[3]].methylationDict, minCpGs)
				if cpgCheck == None:
					del outRegionDict[lineList[3]]
					continue
				else:
					sys.stdout = stdoutHolder
					print json.dumps(outRegionDict[lineList[3]].methylationDict)
					sys.stdout = slStdout

	
	else:
		availableLines = range(1, totalNumRegions+1)
		while(len(outRegionDict.keys()) < numRandRegions):
			randLineNumber = rand.choice(availableLines)
			availableLines.remove(randLineNumber)
			print("There are %d regions left to choose from." % len(availableLines))

			if len(availableLines) == 0:
				print("There are no more allowed regions in the file. Stopping search.")
				break

			lineList = getLineCheckLength(regionsFileName, randLineNumber, minRegionLength)
			if lineList == None:
				continue
			else:
				outRegionDict[lineList[3]] = Region(lineList[3], lineList[0], lineList[1], lineList[2], replicateList, knownStatus, minCpGs)

				if outRegionDict[lineList[3]].methylationDict == None:
					del outRegionDict[lineList[3]]
					continue

				cpgCheck = checkCpGNum(outRegionDict[lineList[3]].methylationDict, minCpGs)
				if cpgCheck == None:
					del outRegionDict[lineList[3]]
					continue
				else:
					sys.stdout = stdoutHolder
					print json.dumps(outRegionDict[lineList[3]].methylationDict)
					sys.stdout = slStdout

	return outRegionDict

logging.basicConfig(
   		level=logging.DEBUG,
   		format='%(asctime)s:%(levelname)s:%(name)s:%(message)s',
   		filename="out.log",
   		filemode='wa')

stdoutHolder = sys.stdout
stderrHolder = sys.stderr

stdout_logger = logging.getLogger('STDOUT')
slStdout = StreamToLogger(stdout_logger, logging.INFO)
sys.stdout = slStdout
 
stderr_logger = logging.getLogger('STDERR')
slStderr = StreamToLogger(stderr_logger, logging.ERROR)
sys.stderr = slStderr
	
print("Starting program")
allRegionsFilename = sys.argv[1]
knownRegionsFilename = sys.argv[2]
bedgraphFilenames = sys.argv[3].split(',')
numRandRegions = int(sys.argv[4])
minCpGs = int(sys.argv[5])
minRegionLength = int(sys.argv[6])

chrDictByRep = {rep:{} for rep in bedgraphFilenames}

print("Getting random regions")
randRegionDict = getRegionDict(bedgraphFilenames, numRandRegions, minCpGs, minRegionLength, allRegionsFilename, -1)

if knownRegionsFilename != "NONE":
	print("Getting known regions")
	knownRegionDict = getRegionDict(bedgraphFilenames, 0, minCpGs, minRegionLength, knownRegionsFilename, 1)

sys.stdout = stdoutHolder
sys.stderr = stderrHolder


'''
keepLog = sys.argv[7]

if keepLog in ['F', 'False', 'FALSE']:
	logFile = open('/dev/null/', 'w')
else:
	logFile = open(keepLog, 'wa')
'''
'''
x = []
y = []

for r in randRegionDict.keys():
	x.append(randRegionDict[r].betweenRepMeanVar)
	y.append(randRegionDict[r].withinRepMeanVar)
plt.plot(x,y,'ro')

if knownRegionsFilename != "NONE":
	x = []
	y = []
	for r in knownRegionDict.keys():
		x.append(knownRegionDict[r].betweenRepMeanVar)
		y.append(knownRegionDict[r].withinRepMeanVar)

	plt.plot(x,y,'bx')

plt.xlabel("Mean variance between replicates")
plt.ylabel("Mean variance within replicates")
plt.show()

print("HERE HERE HERE")

'''


























		

