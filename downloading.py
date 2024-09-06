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
    api_url = 'https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28taxonomy_name%3A%22Escherichia+phage%22%29%29+AND+%28reviewed%3Afalse%29'
    output_file = 'uniprotkb_escherichia_phage.txt'
    download_uniprot_data(api_url, output_file)

