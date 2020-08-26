import json
def get_json_data(file_name, key):
    with open(file_name) as jsonFile:
        data = json.load(jsonFile)
        # print(data[key])
        return data[key]