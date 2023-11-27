"""analyzer.py"""

import pandas as pd


class GeneExpressionAnalyzer:
    """Class to analyze gene expression data based on structural variant coordinates called by HiSV."""

    def __init__(self, filepath):
        """
        Initialize the analyzer with a file path to the gene expression data

        Args:
            filepath (str): The file path to the CSV data.
        """
        try:
            self.gene_data = pd.read_csv(filepath)
        except Exception as e:
            raise ValueError(f"Failed to read the file {filepath}. Error: {e}")

    def subset_by_coordinates(self, chromosome, start, end):
        """Return a subset of data within the specified coordinates.

        Args:
            chromosome (str): The chromosome to filter by (e.g., 'chrX').
            start (int): The starting position of the range (inclusive).
            end (int): The ending position of the range (inclusive).

        Returns:
            pd.DataFrame: A DataFrame containing the subset within the given coordinates.
        """
        # Validation checks for chromosome, start, and end can be added here.
        subset = self.gene_data[(self.gene_data['seqname'] == chromosome) & (self.gene_data['start'] >= start) &
                                (self.gene_data['end'] <= end)]
        return subset

    def subset_by_range(self, coordinate_ranges, padding_left=0, padding_right=0):
        """
        Obtain subsets of data to the left and right of specified coordinates.

        Args:
            coordinate_ranges (list of tuple): List of tuples containing coordinate ranges.
            padding_left (int): The padding to subtract from the start of the left range.
            padding_right (int): The padding to add to the end of the left range.

        Returns:
            list of tuple: A list of DataFrame tuples containing the left and right subsets.
        """
        if not isinstance(coordinate_ranges, list) or not all(isinstance(item, (list, tuple)) for item in coordinate_ranges):
            raise ValueError("coordinate_ranges must be a list of tuples or lists.")

        subsets = []
        for left_range, right_range in coordinate_ranges:
            if padding_left < 0 or padding_right < 0:
                raise ValueError("Padding values must be non-negative.")

            left_subset = self.subset_by_coordinates(chromosome=left_range[0],
                                                     start=left_range[1] - padding_left,
                                                     end=left_range[2] + padding_right)
            right_subset = self.subset_by_coordinates(chromosome=right_range[0], start=right_range[1], end=right_range[2])
            subsets.append((left_subset, right_subset))
        return subsets

    def calculate_population_mean(self):
        """
        Calculate the mean expression across all samples and add it to the data.

        This method assumes that gene expression columns are before the 'seqname' column.
        Adds a new column 'population_mean' to the data with the calculated means.
        """
        try:
            # Assumes gene expression columns are the numeric columns before 'seqname'
            expression_data = self.gene_data.select_dtypes(include=[float, int])
            self.gene_data['population_mean'] = expression_data.mean(axis=1)
        except KeyError:
            raise ValueError("Expected 'seqname' column not found in data.")

    def compare_sample_to_population(self, sample_column):
        """
        Compare a sample's expression to the population mean.

        Args:
            sample_column (str): The name of the column for the sample to compare.

        Returns:
            pd.DataFrame: A DataFrame with the original and comparison data.
        """
        if 'population_mean' not in self.gene_data.columns:
            self.calculate_population_mean()

        if sample_column not in self.gene_data.columns:
            raise ValueError(f"Sample column '{sample_column}' not found in data.")

        comparison_column = f'{sample_column}_vs_population'
        self.gene_data[comparison_column] = self.gene_data[sample_column] - self.gene_data['population_mean']

        return self.gene_data[[sample_column, comparison_column]]
