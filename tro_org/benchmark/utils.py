import datetime
import fnmatch
import os
import time


def find_file(path, name):
    result = []
    for root, dirs, files in os.walk(path):
        if name in files:
            path = os.path.join(root, name)
            tool = get_tool(path)
            result.append([tool, path])

    return result

def get_tool(path):
    path_split = path.split("/")
    tool = path_split[-3].split("_")[-1]
    return tool

def find_file_w_pattern(path, pattern):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
               result.append(os.path.join(root, name))

    return result

def unit_conversion_time(unix_time):
    return datetime.datetime.fromtimestamp(int(unix_time)).strftime('%Y-%m-%d %H:%M:%S')
