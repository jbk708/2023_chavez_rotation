# coordinate_reader.py


def read_coordinates(filepath):
    """
    Read coordinate ranges from a file and return them as a list of tuples.
    
    Args:
        filepath (str): The path to the file containing coordinates.
        
    Returns:
        list of tuples: A list where each tuple contains two tuples representing the left and right coordinate ranges.
    """
    coordinate_ranges = []
    with open(filepath, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            # Assuming each line has exactly 8 values: chr1 start1 end1 chr2 start2 end2 ...
            left_range = (parts[0], int(parts[1]), int(parts[2]))
            right_range = (parts[3], int(parts[4]), int(parts[5]))
            coordinate_ranges.append((left_range, right_range))
    return coordinate_ranges
