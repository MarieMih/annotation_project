"""
For downloading proteins from Uniprot, but actually
from elsewhere.
First argument of cmd is a name to record, 
second argument is a formed api query from Uniprot.
"""

import sys
import requests


def download_uniprot_data(api_url, output_file):
    try:
        response = requests.get(api_url)
        response.raise_for_status()  # Check if the request was successful
        with open(output_file, 'wb') as file:
            file.write(response.content)
        print(f'Data successfully downloaded and saved to {output_file}')
    except requests.exceptions.HTTPError as http_err:
        print(f'HTTP error occurred: {http_err}')
    except Exception as err:
        print(f'Other error occurred: {err}')


if __name__ == "__main__":
    output_file = sys.argv[1]
    api_url = sys.argv[2]
    download_uniprot_data(api_url, output_file)
