import os
from pathlib import Path

def save_to_downloads(filename):
    """
    Saves a text file to the user's Downloads folder dynamically.

    Parameters:
    filename (str): The name of the file (e.g., 'example.txt').
    content (str): The content to write to the file.

    Returns:
    str: The full path to the saved file.
    """
    # Get the user's home directory
    home = Path.home()

    # Construct the path to the Downloads folder
    downloads_folder = home / 'Downloads'

    # Ensure the Downloads folder exists (it should on most systems)
    if not downloads_folder.exists():
        raise FileNotFoundError("The Downloads folder could not be located.")

    # Define the full path for the file
    file_path = downloads_folder / filename

    # Return the file path
    return str(file_path)

if __name__ == '__main__':

    filename = 'untitled.txt'
    file_path = save_to_downloads(filename)
    print(file_path)