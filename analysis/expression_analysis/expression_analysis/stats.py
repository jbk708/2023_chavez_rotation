"""stats.py"""

import pandas as pd


def calculate_population_metrics(data: pd.DataFrame, sample_prefix="MB"):
    # TODO: Fix calculate_population_metrics to better recognize sample columns
    """
    Calculate the mean and standard deviation of expression across all prefixed sample columns and add it to the data.

    This method filters columns that start with the prefix and calculates their mean and standard deviation,
    then adds new columns 'expression_mean' and 'expression_std' to the data with the calculated values.
    
    Args:
        data (pd.DataFrame): DataFrame containing normalized RNAseq data with sample names containing unique prefix.
        sample_prefix (str): String of the sample prefix to identify the relevant columns.

    Returns:
        pd.DataFrame: The DataFrame with added 'expression_mean' and 'expression_std' columns.
    """
    # Filter columns that start with the specified prefix
    columns = [col for col in data.columns if col.startswith(sample_prefix)]

    # Check if prefixed columns are present
    if not columns:
        raise ValueError(f"No columns starting with {sample_prefix} found in data.")

    # Calculate the mean for prefixed columns only
    data['expression_mean'] = data[columns].mean(axis=1)
    data['expression_std'] = data[columns].std(axis=1)

    return data


def filter_outliers_by_sample(data: pd.DataFrame, sample_name: str):
    """
    Filter rows based on a sample name and return a subset where the sample's expression
    is outside the standard deviation range.

    Args:
        data (pd.DataFrame): The DataFrame containing normalized RNAseq data with mean and std columns.
        sample_name (str): The name of the sample to filter by.

    Returns:
        pd.DataFrame: A subset of the DataFrame where the sample's expression is outside one standard deviation.
    """
    if sample_name not in data.columns:
        raise ValueError(f"Sample name {sample_name} not found in data.")

    if 'expression_mean' not in data.columns or 'expression_std' not in data.columns:
        raise ValueError("Dataframe must contain 'expression_mean' and 'expression_std' columns.")

    # Define the bounds for being within one standard deviation
    lower_bound = data['expression_mean'] - data['expression_std']
    upper_bound = data['expression_mean'] + data['expression_std']

    # Filter rows where the sample's expression is outside the standard deviation bounds
    outliers = data[(data[sample_name] < lower_bound) | (data[sample_name] > upper_bound)]

    return outliers
