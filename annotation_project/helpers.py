import sys
import subprocess
import os
from datetime import datetime
import telegram_send


def check_file_exists(file):
    if os.path.exists(file) and os.path.isfile(file):
        # print(f'The file {file} exists.')
        return 0
    if os.path.isfile(file):
        # print(f'The file {file} does not exist.')
        return 1
    else:
        # print(f'The file {file} is not regular file.')
        return 2

def check_dir_exists(directory):
    if os.path.isdir(directory):
        print(f'The directory {directory} exists.')
        return 0
    print(f'The directory {directory} does not exist.')
    exit()

def check_all_prefix_unique():
    """
    Verify every genome from pool is unique (latest 24 sym)
    """
    pass

def return_str_with_date_and_time():
    now = datetime.now()
    return now.strftime("%Y-%m-%d-%H-%M-%S")

def union_files(input_list, output_file):
    with open(output_file, 'w') as outfile:
        for fname in input_list:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

def create_directory(target_or):
    target = os.path.abspath(target_or)  # отдебажить!!!
    if not os.path.exists(target):
        os.makedirs(target)
    else:
        print(f"Directory {target} exists already.")

def create_directory_with_soft_links(tsvs, target_or):
    target = os.path.abspath(target_or)  # отдебажить!!!
    create_directory(target)
    for i in tsvs:
        new_link = os.path.split(i)[1]
        os.symlink(os.path.abspath(i), os.path.join(target, new_link))  # отдебажить!!!


def create_acronym(phrase: str):
    if (phrase == "") or (phrase == "nan"):
        return "HP"
    trimmed = phrase.split('(')[0].strip()
    words = trimmed.split()
    acronym = "".join(word[0].upper() for word in words if word)
    return acronym

async def send_smth(cor_image, pan_image):
    with open(cor_image, "rb") as f:
        await telegram_send.send(images=[f])
    with open(pan_image, "rb") as f:
        await telegram_send.send(images=[f])