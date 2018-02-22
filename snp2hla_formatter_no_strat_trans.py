import sys
import os
import glob
import pandas as pd
import numpy as np
import argparse

# Future improvements
# 1. Marker QC on marker info instead of dosage file.
# 2. merge subsets using dataframes rather than PLINK - this seems slow.

VERSION = 'V1.2'

def get_options():

        parser = argparse.ArgumentParser()

        parser.add_argument('--dosage', action='store', dest='dosage_file', required=True,
                    help='Input dosage file')
        parser.add_argument('--out', action='store', dest='out_file', required=True,
                    help='Output file')
        parser.add_argument('--exclude', action='store', dest='exclude_file',
                    help='File of SNPs to exclude')
        parser.add_argument('--remove', action='store', dest='remove_file',
                    help='File of samples to exclude. Two fields; FID and IID')
        parser.add_argument('--keep', action='store', dest='keep_file',
                    help='File of samples to keep. Two fields; FID and IID')
        parser.add_argument('--maf', action='store', dest='maf_threshold', type=float,
                    help='MAF threshold for excluding SNPs')
	parser.add_argument('--r2', action='store', dest='r2_threshold', type=float,
                    help='Beagle r2 threshold for excluding SNPs')

        args =  parser.parse_args()
        return args

def make_fam(input, output):
        """
        Create an appended .fam file for dosage summary analysis. Write file to disk.
	Also create header files for individual dosage file - required when using list file
        """
        try:		
		fam_file = output + '.fam'
		df = pd.DataFrame()

		fam_files = sorted(glob.glob(input + '*.fam'))
		for (counter, file) in enumerate(fam_files, start=1):
			header_file = output + ".header%s" % counter
			fam_info = pd.read_table(file, header=None, sep=' ')
			fam_info.to_csv(header_file, header=False, index=False, sep=' ')
			df = df.append(fam_info)

		df.to_csv(fam_file, header=False, index=False, sep=' ', float_format='%.f')

        except IOError:
                sys.exit('Cannot open {input}.'.format(**locals()))
  
def write_dosage(input, output):
	"""
	Merge dosage subsets. Output written to disk.
	"""
	try:

		dosage_list = output + '.dosage.list'
		fam_file = output + '.fam'

		# create dosage file list
		f = open(dosage_list, 'w')
		dosage_files = sorted(glob.glob(input + '*.dosage'))

		for (counter, file) in enumerate(dosage_files, start=1):
			header_file = output + '.header' + str(counter)
			f.write("1 %s %s\n" % (file, header_file))
	
		f.close()

		# merge using PLINK
		plink_cmd = "/home/mlibydt3/galaxy/tools/bin/plink2 --noweb --allow-no-sex --fam %s --dosage %s list format=1 sepheader --write-dosage --out %s --silent" % (fam_file, dosage_list, output)
		os.system(plink_cmd)
				
        except IOError:
                sys.exit('Cannot open {input}.'.format(**locals()))

def make_map(input, output):
        """
        Create a map for dosage summary analysis. Write file to disk.
        """
        try:
                map_file = output + '.map'

                bim_file = glob.glob(input + '*1.bim')[0]
                df = pd.read_table(bim_file, header=None)
                df.to_csv(map_file, header=False, index=False, sep=' ', cols=[0, 1, 2, 3])

        except IOError:
                sys.exit('Cannot open {input}.'.format(**locals()))

def marker_summary(input, output):
	"""
	Generate dosage summary stats. Output written to disk.
	"""
	try:

		fam = output + '.fam'
		map = output + '.map' 
		dosage = output + '.out.dosage'

		plink_cmd = "plink --noweb --dosage %s format=1 --fam %s --map %s --out %s --silent" % (dosage, fam, map, output)
		os.system(plink_cmd)

        except IOError:
                sys.exit('Cannot open {path}.'.format(**locals()))

def average_r2(input_dir, output_dir):
	"""
	Calcualte the average Beagle r2 score across the subsets to be used in marker QC.
	"""
	try:
		info_file = output_dir + '.r2'
		list = []

                info_files = sorted(glob.glob(input_dir + '*.bgl.r2'))
                for file in info_files:
                        info = pd.read_table(file, header=None, names=['SNP','R2'])
                        list.append(info)

		ave_r2 = pd.concat(list).groupby(['SNP']).mean()
		ave_r2['SNP'] = ave_r2.index.values

		return ave_r2

	except IOError:
                sys.exit('Cannot open {input_dir}.'.format(**locals()))

def format_markers(output_dir, marker_r2):
        """
        Format marker data. Return dataframe of marker information.
        """
        try:
		marker_file = output_dir + '.assoc.dosage'

		df = pd.read_table(marker_file, sep=r'\s*')
		df.drop(['OR', 'SE', 'P'], axis=1, inplace=True)

		# create new fields
		df['MAF'] = np.where(df['FRQ'] > 0.5, 1 - df['FRQ'], df['FRQ'])
		df['RESIDUE'] = np.where(df['SNP'].str.startswith("AA_"), df['SNP'].str.split('_').str[-1], 'NA')
		df['MARKER'] = np.where(df['SNP'].str.startswith("AA_"), df['SNP'].str.split('_').str[0:3].str.join('_'), df['SNP'])		
		df['GENE'] = np.where(df['SNP'].str.contains("_"), df['SNP'].str.split('_').str[1], 'SNP')

		# include Beagle r2
		df = pd.merge(df, marker_r2, on='SNP')
		
		return df

        except IOError:
                sys.exit('Cannot open {output}.'.format(**locals()))

def format_dosage(input):
        """
        Transpose dosage data and exclude low quality variants. Output written to disk.
        """
        try:
		dosage_file = input + '.out.dosage'

		# deal with multiple entries in header.
		with open(dosage_file, 'r') as f:
			header = f.readline()

		header = header.split()
		
		header_concat = header[:3] + [i+'-'+j for i,j in zip(header[3::2], header[4::2])]
	
		# transpose dosage - use unique_header for column names.
		df = pd.read_table(dosage_file, sep=r'\s*', skiprows=1, header=None, names=header_concat, index_col=0)
		df.drop(['A1', 'A2'], axis=1, inplace=True)
		df = df.T
		df.insert(0, 'FID', header[3::2])
		df.insert(1, 'IID', header[4::2])

		df.to_csv(dosage_file, header=True, index=False, sep=' ', float_format='%.f')
		return df

        except IOError:
                sys.exit('Cannot open {input}.'.format(**locals()))

def qc_dosage(dosage_df, marker_df, maf_threshold, r2_threshold, snp_exclude, sample_remove, sample_keep):
	"""
	QC the transposed dosage file.
	"""
	try:
		snp_exclusions = []

                # get list of composite markers
		amino_acids = [x for x in list(dosage_df.columns.values) if ('AA' in x) & (len(x.split('_')) > 4)]
		composite_markers = [x for x in amino_acids if len(x.split('_')[4]) > 1]
		print ' --excluded %s composite markers' % len([x for x in composite_markers if x not in snp_exclusions])
		snp_exclusions.extend(composite_markers)

		# maf threshold
		maf_exclusions = list(marker_df['SNP'][marker_df['MAF'] < maf_threshold])
		print ' --excluded %s low MAF markers' % len([x for x in maf_exclusions if x not in snp_exclusions])
		snp_exclusions.extend(maf_exclusions)

		# r2 threshold
                r2_exclusions = list(marker_df['SNP'][marker_df['R2'] < r2_threshold])
                print ' --excluded %s low r2 markers' % len([x for x in r2_exclusions if x not in snp_exclusions])
                snp_exclusions.extend(r2_exclusions)		
		
		# append from supplied list
		if snp_exclude:
			snp_supplied = [line.strip() for line in open(snp_exclude, 'r')]
			print ' --excluded %s markers based on supplied list' % len([x for x in snp_supplied if x not in snp_exclusions])
			snp_exclusions.extend(snp_supplied)

		# exclude failed SNPs
                dosage_df.drop(snp_exclusions, axis=1, inplace=True)

		# list of samples to exclude
		if sample_remove:
			sample_exclusions = [line.replace(" ", "-").strip() for line in open(sample_remove, 'r')]
			#sample_exclusions = ["-".join(line.split()) for line in open(sample_remove, 'r')]
			dosage_df.drop(sample_exclusions, inplace=True)
			print ' --excluded %s samples based on supplied list' % len(sample_exclusions)

		# list of samples to keep
		if sample_keep:
			samples_to_keep = [line.replace(" ", "-").strip() for line in open(sample_keep, 'r')]
			dosage_df = dosage_df[dosage_df.index.isin(samples_to_keep)]
			print ' --kept only %s samples based on supplied list' % len(samples_to_keep)

		return (dosage_df, snp_exclusions)

	except IOError:
                sys.exit('Cannot open {path}.'.format(**locals()))

def merge_fam(fam, df):
	"""
	Merge fam info with dosage file. Retrun modified dataframe.
	"""
	try:
		fam_info = pd.read_table(fam, header=None, sep=' ')
		fam_info.drop([1, 2, 3], axis=1, inplace=True)
		fam_info.columns = ['ID', 'sex', 'phenotype']
		
		df = fam_info.merge(df)

		return df

        except IOError:
                sys.exit('Cannot open {path}.'.format(**locals()))

def main():

	print '\nSNP2HLA_FORMATTOR: VERSION: ' + VERSION + '\n'

	# define input/output
	options = get_options()

	# merge individual dosage files.
	print 'Merging dosage subsets and transposing file'
        make_fam(options.dosage_file, options.out_file)
	write_dosage(options.dosage_file, options.out_file)
	dosage = format_dosage(options.out_file)
	print '\nInitial dosage file contains %s markers and %s samples\n' % (len(dosage.columns) - 1, len(dosage.index))
	
	# Prepare marker information
	#print 'Preparing marker information'
	#make_map(options.dosage_file, options.out_file)
	#marker_summary(options.dosage_file, options.out_file)
	#marker_r2 = average_r2(options.dosage_file, options.out_file)
	#marker_info = format_markers(options.out_file, marker_r2)

	# QC transposed dosage file
	#print 'Excluding markers and samples'
	#dosage_qc, marker_qc = qc_dosage(dosage, marker_info, options.maf_threshold, options.r2_threshold, options.exclude_file, options.remove_file, options.keep_file)
	#dosage_outfile = options.out_file + '.qc.dosage'
	#dosage_qc.to_csv(dosage_outfile, index=False, header=True, sep=' ', float_format='%.4g')

	# write QC'd marker information	
	#marker_info.set_index('SNP', inplace=True, drop=False, append=False)
	#marker_info.drop(marker_qc, inplace=True)
	#marker_info['SNP'] = marker_info['SNP'].str.replace('1kg','X1kg')
	#marker_info['MARKER'] = marker_info['MARKER'].str.replace('1kg','X1kg')
	#marker_outfile = options.out_file + '.marker_info'
  #     marker_info.to_csv(marker_outfile, sep = ' ', index=False, float_format='%.4g')

	#print '\nFinal dosage file contains %s markers and %s samples\n' % (len(dosage_qc.columns) - 2, len(dosage_qc.index))

if __name__=='__main__':
	main()
