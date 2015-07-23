import linecache, re, os, sys, subprocess
from math import *

def lineSplit(lineString):
		return re.split(re.compile("\s+"), lineString)

def bisectionSearch(fileNameToSearch, firstLine, lastLine, regionStart, regionEnd):

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


chrDictByRep = {"empty.txt":{'chr2':[4,19]}}
out = getCpGs("chr2", 2003, 2023, "empty.txt")
print(out)
