"""
    A python script for finding collapsed regions in genome assemblies.
    Copyright (C) 2020-2021 Jana Sperschneider  
    This program is free software; you can redistribute it and/or modify  
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation; either version 3 of the License, or     
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    Contact: jana.sperschneider@anu.edu.au
"""
import sys
import os
#--------------------------------------
HIGH_COVERAGE_THRESHOLD = 200.0
#--------------------------------------
def merge_intervals(intervals):

  sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])

  merged = []

  for higher in sorted_by_lower_bound:
      if not merged:
          merged.append(higher)
      else:
          lower = merged[-1]
          # test for intersection between lower and higher:
          # we know via sorting that lower[0] <= higher[0]
          if higher[0] <= lower[1] + 1:            
              upper_bound = max(lower[1], higher[1])
              merged[-1] = (lower[0], upper_bound)  # replace by merged interval
              
          else:
              merged.append(higher)
  return merged
#--------------------------------------
def read_in_coverage(COVERAGE_FILE):
	''' This takes the bincoverage output file of bbmap's pileup.sh as the input 
		#Mean   78.662
		#STDev  358.479
		#RefName        Cov     Pos     RunningPos
		tig00000006     6.54    1000    0
		tig00000006     9.01    2000    1000
		tig00000006     9.76    3000    2000
		tig00000006     12.36   4000    3000
		tig00000006     12.77   5000    4000
		tig00000006     10.70   6000    5000
'''
	
	COVERAGE = {}

	# Determine the bin size used in the pileup file first
	with open(COVERAGE_FILE ,'r') as f:

	    for line in f:
	    	if line.startswith('#'):
	    		pass
	    	else:
		      BINSIZE = int(line.split('\t')[2])	
		      break

	with open(COVERAGE_FILE ,'r') as f:

	    for line in f:
	    	if line.startswith('#'):
	    		pass
	    	else:
		      contig = line.split('\t')[0]
		      coverage = float(line.split('\t')[1])
		      end_position = int(line.split('\t')[2])	

		      if contig in COVERAGE:
		      	COVERAGE[contig] = COVERAGE[contig] + [(coverage, end_position)]
		      else:
		      	COVERAGE[contig] = [(coverage, end_position)]

	return COVERAGE, BINSIZE
#--------------------------------------
#--------------------------------------
try:
	COVERAGE_FILE = sys.argv[1]
except:
	print('Please provide the pileup.sh bincoverage output file.')
	print('\n-------------')
	print('Usage: python CollapsedRegions.py <bincoverage.txt> <threshold> <ouput bedfile>')
	print('-------------')
	sys.exit()

try:
	COVERAGE_THRESHOLD = int(sys.argv[2])
except:
	print('Please provide the midpoint as the collapsed region threshold from inspecting your coverage distribution plot.')
	print('\n-------------')
	print('Usage: python CollapsedRegions.py <bincoverage.txt> <threshold> <ouput bedfile>')
	print('-------------')
	sys.exit()
try:
	OUTPUT = sys.argv[3]
except:
	print('Please provide a output file name for the collapsed region BED file.')
	print('\n-------------')
	print('Usage: python CollapsedRegions.py <bincoverage.txt> <threshold> <ouput bedfile>')
	print('-------------')
	sys.exit()
#--------------------------------------
print('Read in coverage from a bincoverage file from pileup.sh')
COVERAGE, BINSIZE = read_in_coverage(COVERAGE_FILE)
print('Binsize that was used in pileup.sh:', BINSIZE)
#--------------------------------------
collapsed_regions = {}
collapsed_regions_coverage = []
#--------------------------------------
for contig, values in COVERAGE.items():

	collapsed_regions_on_contig = []

	for index, (coverage, end_position) in enumerate(values):
		start_position = index*BINSIZE + 1
				
		# To avoid counting mitochondrial and other high-coverage regions as collapsed, set a high coverage limit
		if coverage > COVERAGE_THRESHOLD and coverage < HIGH_COVERAGE_THRESHOLD: 
			# This region appears to be collapsed
			collapsed_regions_on_contig.append((start_position, end_position))
			collapsed_regions_coverage.append((contig, start_position, end_position, coverage))

	if collapsed_regions_on_contig:
		# Merge adjacent intervals to get the overall collapsed regions
		if contig in collapsed_regions:
			collapsed_regions[contig] = collapsed_regions[contig] + merge_intervals(collapsed_regions_on_contig)
		else:
			collapsed_regions[contig] = merge_intervals(collapsed_regions_on_contig)

# We now have the regions on each contig that appear to be collapsed
# Calculate the average read coverage for the merged intervals
with open(OUTPUT, 'w') as f:
	collapsed_bases = []
	for contig, regions in collapsed_regions.items():
		
		for (start, end) in regions:
			list_of_coverages = []
			# Calculate the average coverage for that region
			for (contig_ident, start_position, end_position, coverage) in collapsed_regions_coverage:
				if contig == contig_ident:
					if start_position >= start and start_position <= end and end_position >= start and end_position <= end:
						list_of_coverages.append(coverage)

			average_coverage = round(sum(list_of_coverages)/len(list_of_coverages),2)
			#collapsed_bases.append(end-start+1) 
			
			# This is more accurate, if a region is e.g. triple the coverage then add collapsed regions*3
			times_coverage = average_coverage//COVERAGE_THRESHOLD
			collapsed_bases.append(times_coverage *(end-start+1))
			f.writelines(contig + '\t' + str(start-1) + '\t' + str(end) + '\t' + 'Read-coverage:' + str(average_coverage) + '\n')

#--------------------------------------
#--------------------------------------
print
print(sum(collapsed_bases)/1000000.0, 'MB are collapsed ==', sum(collapsed_bases), 'bases.')
print('Read coverage cutoff for collapsed regions that was used:', COVERAGE_THRESHOLD)
print('Collapsed regions were saved to bedfile', OUTPUT)
#--------------------------------------
#--------------------------------------
