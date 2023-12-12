"""
expression.py
"""

import logging
import pandas as pd

from coords import read_coordinate_data
from stats import calculate_population_metrics
from stats import filter_outliers_by_sample


def queue_samples(sample_list: list, sample_folder: str, gene_data_path: str, chr_path: str, outpath: str):
    """
    Initiates sample queue for expression analysis.

    Args:
        sample_list (list): List of sample names to be processed
            EX: ['MB102', 'MB103']
        sample_folder (str): Path to folder containing gene coordinate files
        gene_data_path (str): file path for the normalized RNAseq Data
        chr_path (str): File path to chromosome len file
        outpath (str): Path to output folder. Samples will be named "{OUTPATH}_{SAMPLE_NAME}_genes.csv" 
    
    Returns:
        CSV files for each sample with elevated gene expression.
    """

    for sample in sample_list:
        try:
            coord_list = read_coordinate_data(f"{sample_folder}/{sample}_HiSV_inter_SV_result.txt")
        except FileNotFoundError:
            print(f"Coordinate file for sample {sample} not found. Skipping this sample.")
            continue

        try:
            exp = ExpressionAnalysis(gene_data_path, chr_path)
            gene_subset = exp.subset_by_range(coord_list, padding_left=2000000)
            genes_with_mean = calculate_population_metrics(gene_subset, "MB")
            outliers = filter_outliers_by_sample(genes_with_mean, sample_name="MB174")
            outliers.to_csv(f"{outpath}/{sample}_genes.csv", index=False)
        except Exception as e:
            print(f"An error occurred while processing sample {sample}: {e}")
            continue



class ExpressionAnalysis:
    """
    Takes normalized RNAseq data merged with gene coordinates and identifies
    abnormal expression flanking and inside of a predicted structural variant.

    """

    def __init__(self, gene_data_fp: str, chr_size_fp: str, logger=None) -> None:
        """ Class initiation
        Args:
            gene_data_fp (str): Path to normalized RNAseq data in csv format, should include
            both gene expression data per sample, as well as genetic coordinates.

            Ex: 
            gene	MB095	MB106	MB170	MB226	seqname	start	end	strand	gene_id
            TSPAN6	1604.586519	1360.387931	3453.889243	1668.121774	chrX	100627108	100639991	-	ENSG00000000003.15
            TNMD	5.938514135	0	6.130127966	0	chrX	100584936	100599885	+	ENSG00000000005.6
            DPM1	1340.916492	2214.562145	1980.031333	2363.08503	chr20	50934867	50959140	-	ENSG00000000419.14
        """
        if logger is not None:
            self.log = logger
        else:
            self.log = logging.getLogger(__name__)

        try:
            self.gene_data = pd.read_csv(gene_data_fp)
        except Exception as err:
            raise ValueError(self.log.error(f'Failed to read file {gene_data_fp}, Error: {err}')) from err
        try:
            self.chr_sizes = self.chromosome_sizes(chr_size_fp)
        except Exception as err:
            raise ValueError(self.log.error(f'Failed to read file {chr_size_fp}, Error: {err}')) from err

    def chromosome_sizes(self, chromosome_size_path: str) -> dict:
        """Load chromosome sizes from a TSV file into a dictionary."""
        chromosome_sizes = {}
        try:
            with open(chromosome_size_path, 'r', encoding='utf8') as file:
                for line in file:
                    parts = line.strip().split('\t')
                    chromosome_sizes[parts[0]] = int(parts[1])
            return chromosome_sizes
        except Exception as err:
            raise ValueError(self.log.error(f'Failed to read the file {chromosome_size_path}. Error: {err}')) from err

    def _truncate_coordinates(self, chromosome: str, start: int, end: int) -> int:
        """Truncate coordinates to lie within chromosome bounds."""
        max_size = self.chr_sizes.get(chromosome)
        if max_size is None:
            raise ValueError(self.log.error(f"Chromosome {chromosome} is not in the chromosome size file."))
        start = max(0, min(start, max_size))
        end = max(0, min(end, max_size))
        return start, end

    def subset_by_coordinates(self, chromosome: str, start: int, end: int):
        """Return a subset of data within the specified coordinates.

        Args:
            chromosome (str): The chromosome to filter by (e.g., 'chrX').
            start (int): The starting position of the range (inclusive).
            end (int): The ending position of the range (inclusive).

        Returns:
            pd.DataFrame: A DataFrame containing the subset within the given coordinates.
        """
        start, end = self._truncate_coordinates(chromosome, start, end)
        return self.gene_data[(self.gene_data['seqname'] == chromosome) & (self.gene_data['start'] >= start) &
                              (self.gene_data['end'] <= end)]

    def subset_by_range(self, coordinate_ranges: list, padding_left=0, padding_right=0):
        """
        Obtain a single DataFrame of data to the left and right of specified coordinates,
        removing any duplicates.

        Args:
            coordinate_ranges (list of tuple): List of tuples containing coordinate ranges.
            padding_left (int): The padding to subtract from the start of the left range.
            padding_right (int): The padding to add to the end of the left range.

        Returns:
            pandas.DataFrame: A single DataFrame containing the unique subsets from the left and right ranges.
        """
        if not isinstance(coordinate_ranges, list) or not all(isinstance(item, (list, tuple)) for item in coordinate_ranges):
            raise ValueError(self.log.error("coordinate_ranges must be a list of tuples or lists."))

        all_subsets = pd.DataFrame()
        for left_range, right_range in coordinate_ranges:
            left_start, left_end = self._truncate_coordinates(left_range[0], left_range[1] - padding_left,
                                                              left_range[2] + padding_right)
            right_start, right_end = self._truncate_coordinates(right_range[0], right_range[1], right_range[2])
            left_subset = self.subset_by_coordinates(left_range[0], left_start, left_end)
            right_subset = self.subset_by_coordinates(right_range[0], right_start, right_end)
            all_subsets = pd.concat([all_subsets, left_subset, right_subset])

        # Remove duplicates and reset the index
        unique_subsets = all_subsets.drop_duplicates().reset_index(drop=True)
        return unique_subsets


"""Early Testing"""

if __name__ == "__main__":
    sample_list = [
        'MB102',
        'MB106',
        'MB164',
        'MB174',
        'MB199',
        'MB227',
        'MB234',
        'MB244',
        'MB248',
        'MB264',
        'MB268',
        'MB274',
        'MB277'
    ]

    queue_samples(
        sample_list,
        "/Users/jkirkland/2023_chavez_rotation/analysis/expression_analysis/data/inter_sv_data",
        "/Users/jkirkland/2023_chavez_rotation/analysis/expression_analysis/data/medullo_rnaseq_annotated.csv",
        "/Users/jkirkland/2023_chavez_rotation/analysis/expression_analysis/data/hg38_no_M.len",
        "/Users/jkirkland/2023_chavez_rotation/analysis/expression_analysis/data/gene_data"
    )