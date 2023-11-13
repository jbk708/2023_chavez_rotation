"""coords.py"""


def validate_and_complete_data(line: str) -> list:
    """
    Looks for missing chr column found in intra data created by HiSV and adds it.

    Check if all six columns are present in the line. If the second chromosome
    column is missing, copy the first one and add it to the data.
    
    Args:
        line (str): A line from the coordinate data file.
        
    Returns:
        list: A list representing the validated and possibly completed data columns.
    """
    parts = line.strip().split('\t')
    # If there are exactly 5 columns, assume the second chromosome column is missing
    if len(parts) == 5:
        parts.insert(3, parts[0])  # Duplicate the first chromosome column
    elif len(parts) != 6:
        raise ValueError("Each line must contain 5 or 6 columns.")
    return parts


def read_coordinate_data(filepath: str) -> list:
    """
    Reads a file with coordinate data and returns a list of tuples.
    
    This function calls `validate_and_complete_data` to ensure data integrity.
    
    Args:
        filepath (str): The file path to the coordinate data.
        
    Returns:
        list of tuples: Each tuple contains two tuples, one for each coordinate range.
    """
    coordinate_ranges = []
    try:
        with open(filepath, 'r', encoding="utf8") as file:
            for line in file:
                parts = validate_and_complete_data(line)
                left_range = (parts[0], int(parts[1]), int(parts[2]))
                right_range = (parts[3], int(parts[4]), int(parts[5]))
                coordinate_ranges.append((left_range, right_range))
        return coordinate_ranges
    except Exception as err:
        raise ValueError(f"SV file {filepath} not found!")
